#!/usr/bin/env python3
"""Measure how long the main software operations take (writes manuscript/analysis/software_benchmark.json).

Genomes used for timing are synthetic random sequences (seed 42, 50% GC): they measure speed, not biology.
"""
from __future__ import annotations

import json
import platform
import random
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from agripam import biobricks as bb, parts_io, syncom_io  # noqa: E402
from agripam.core import analyze_ttc_genome, read_pam_scores  # noqa: E402
from agripam.editor_targeting import scan_editor_targets  # noqa: E402


def best_of(fn, repeat=3):
    times = []
    result = None
    for _ in range(repeat):
        t0 = time.perf_counter()
        result = fn()
        times.append(time.perf_counter() - t0)
    return round(min(times), 3), result


def random_genome(n: int) -> str:
    rng = random.Random(42)
    return "".join(rng.choices("ACGT", k=n))


def main():
    out = {"machine": f"{platform.system()} {platform.machine()}, Python {platform.python_version()}", "operations": []}

    def add(name, seconds, note=""):
        out["operations"].append({"operation": name, "seconds": seconds, "note": note})
        print(f"{name}: {seconds} s {note}")

    # 1. community analysis on the whole bank
    workbook = ROOT / "manuscript/analysis/software_input_bank.xlsx"
    if workbook.exists():
        sheets = syncom_io.read_workbook(workbook)
        t, (bank, _) = best_of(lambda: syncom_io.build_bank(sheets))
        n = len(bank["traits"])
        add(f"Read and validate the bank ({n} strains)", t)
        t, _ = best_of(lambda: syncom_io.analyze(bank, fungi=["FI20"]))
        add("Assemble a community, rank chassis, gaps, 2,000-weighting robustness", t)
    # 2. genome scans on synthetic genomes
    scores = read_pam_scores(ROOT / "machine_learning/type_ic_all_64_pam_scores.tsv")
    for mb in (1, 5):
        genome = [("synthetic", random_genome(mb * 1_000_000))]
        t, (summary, cands) = best_of(lambda: analyze_ttc_genome(genome, 35, scores, False), repeat=1)
        add(f"5'-TTC target scan and scoring, {mb} Mb genome", t, f"{len(cands):,} candidates")
    genome1 = [("synthetic", random_genome(1_000_000))]
    t, targets = best_of(lambda: scan_editor_targets(genome1, "SpCas9", "NGG", "3prime", 20), repeat=1)
    add("SpCas9 (NGG) target scan, 1 Mb genome", t, f"{len(targets):,} sites")
    # 3. construct design
    library = {p["part_id"]: p for p in parts_io.read_parts_tsv(ROOT / "data/biobricks/starter_library.tsv")}
    genome5 = [("synthetic", random_genome(5_000_000))]

    def build_construct():
        cassette = bb.assemble([library[i] for i in ("BBa_J23100", "BBa_B0034", "BBa_E0040", "BBa_B0015")])
        arms = bb.homology_arms(genome5, "synthetic", 2_500_000, 750)
        construct = bb.integration_construct(cassette, arms, "+")
        return bb.check_sequence("construct", construct["sequence"], list(bb.STANDARDS), {"motif": "GCANNNNNNNTGC"}, "B. subtilis", "", False)
    t, findings = best_of(build_construct)
    add("Assemble a reporter construct (750-bp arms) and run all screens", t, f"{len(findings)} findings")
    # 4. reproducibility checks
    t0 = time.perf_counter()
    subprocess.run([sys.executable, str(ROOT / "scripts/heldout_pam_validation.py")], cwd=ROOT, capture_output=True, check=True)
    add("Leave-one-group-out PAM stability check (both schemes)", round(time.perf_counter() - t0, 3))
    t0 = time.perf_counter()
    r = subprocess.run([sys.executable, "-m", "unittest", "discover", "-s", "tests"], cwd=ROOT, capture_output=True, text=True)
    tail = (r.stderr or r.stdout).strip().splitlines()
    add("Full test suite", round(time.perf_counter() - t0, 3), next((l for l in reversed(tail) if l.startswith("Ran ")), ""))
    # 5. size of the code base
    def lines(pattern_dirs):
        total = 0
        for d in pattern_dirs:
            for p in (ROOT / d).rglob("*.py"):
                if "backup" in p.name or "__pycache__" in str(p):
                    continue
                total += sum(1 for _ in open(p, errors="ignore"))
        return total
    out["code_lines"] = {"library (agripam/)": lines(["agripam"]), "application (app/)": lines(["app"]), "scripts": lines(["scripts"]), "tests": lines(["tests"])}
    (ROOT / "manuscript/analysis/software_benchmark.json").write_text(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
