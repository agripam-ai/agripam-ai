#!/usr/bin/env python3
"""Run DefenseFinder reproducibly on the three 24-genome comparison panels."""

from __future__ import annotations

import concurrent.futures
import gzip
import os
from pathlib import Path
import shutil
import subprocess


ROOT = Path(__file__).resolve().parents[1]
DF_BIN = Path("/private/tmp/agripam-defense-env-v2/bin/defense-finder")
DF_PATH = str(DF_BIN.parent) + os.pathsep + os.environ.get("PATH", "")
MODELS = ROOT / "defense_atlas/models"
OUT = ROOT / "defense_atlas/diverse_24"
INPUTS = ROOT / "defense_atlas/inputs"


def accession(path: Path) -> str:
    for part in path.parts:
        if part.startswith("GCF_") and "." in part:
            return part.split("_")[0] + "_" + part.split("_")[1]
    raise ValueError(f"No accession in {path}")


def read_paths(path: Path) -> list[Path]:
    return [Path(line.strip()) for line in path.read_text().splitlines() if line.strip()]


def panels() -> list[tuple[str, str, Path]]:
    jobs: list[tuple[str, str, Path]] = []
    pp = ROOT / "data/ncbi_packages/paenibacillus_polymyxa_primary_24/ncbi_dataset/data"
    for src in sorted(pp.glob("GCF_*/*.fna.gz")):
        acc = accession(src)
        dst = INPUTS / "paenibacillus_polymyxa" / f"{acc}.fna"
        dst.parent.mkdir(parents=True, exist_ok=True)
        if not dst.exists():
            with gzip.open(src, "rb") as inp, dst.open("wb") as out:
                shutil.copyfileobj(inp, out)
        jobs.append(("paenibacillus_polymyxa", acc, dst))

    for species in ("bacillus_subtilis", "bacillus_velezensis"):
        paths = read_paths(ROOT / f"taxonomy/{species}/ani/mash_diverse_24_genomes.txt")
        jobs.extend((species, accession(src), src) for src in paths)
    return jobs


def run_one(job: tuple[str, str, Path]) -> tuple[str, str, int]:
    species, acc, src = job
    out = OUT / species / acc
    out.mkdir(parents=True, exist_ok=True)
    expected = list(out.glob("*_defense_finder_systems.tsv"))
    if expected:
        return species, acc, 0
    log = out / "defense_finder.log"
    cmd = [str(DF_BIN), "run", str(src), "-o", str(out), "--models-dir", str(MODELS),
           "-w", "1", "--skip-model-version-check"]
    env = dict(os.environ, PATH=DF_PATH)
    with log.open("w") as handle:
        proc = subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT, env=env)
    return species, acc, proc.returncode


def main() -> None:
    jobs = panels()
    print(f"Defense atlas jobs: {len(jobs)}", flush=True)
    failures = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        for species, acc, code in pool.map(run_one, jobs):
            print(f"{species}\t{acc}\t{code}", flush=True)
            if code:
                failures.append((species, acc, code))
    if failures:
        raise SystemExit(f"Failed jobs: {failures}")


if __name__ == "__main__":
    main()
