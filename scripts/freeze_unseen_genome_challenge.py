#!/usr/bin/env python3
"""Freeze an outcome-independent genome challenge for a live competition demo."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from agripam.core import analyze_ttc_genome, genome_qc, parse_fasta, read_pam_scores


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--accession", required=True)
    parser.add_argument("--organism", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--protospacer-length", type=int, default=35)
    parser.add_argument("--max-exported-candidates", type=int, default=1000)
    args = parser.parse_args()

    records = parse_fasta(args.fasta.read_text())
    scores = read_pam_scores(ROOT / "machine_learning" / "type_ic_all_64_pam_scores.tsv")
    summary, candidates = analyze_ttc_genome(records, args.protospacer_length, scores, False)
    args.output.mkdir(parents=True, exist_ok=True)
    candidate_path = args.output / "frozen_reference_ttc_candidates.tsv"
    exported = candidates[:args.max_exported_candidates]
    fields = list(exported[0]) if exported else ["candidate_id"]
    with candidate_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader(); writer.writerows(exported)

    manifest = {
        "challenge_id": "AGRIPAM_UNSEEN_GENOME_001",
        "frozen_at_utc": datetime.now(timezone.utc).isoformat(),
        "accession": args.accession,
        "organism": args.organism,
        "input_fasta_sha256": digest(args.fasta),
        "candidate_table_sha256": digest(candidate_path),
        "outcomes_available_to_analysis": False,
        "used_in_reference_model_training": False,
        "reference_model_scope": "P. polymyxa Type I-C TTC context ranking",
        "transfer_decision": "ABSTAIN from strain-specific native-system claim",
        "reason": "The genome is outside the P. polymyxa reference chassis and has no challenge-specific functional outcome.",
        "summary": summary,
        "exported_top_candidates": len(exported),
        "qc": genome_qc(records),
        "uncertainty": [
            "TTC compatibility is a reference scan, not a demonstrated PAM for this organism.",
            "Exact-copy screening is not a complete approximate off-target search.",
            "No editing efficiency or functional activity is inferred.",
        ],
        "next_experiment": [
            "Confirm the intended chassis and annotate defense, repair and delivery constraints.",
            "Compare introduced-editor and native-system options before selecting an effector.",
            "Freeze top, middle and low-ranked candidates and test them with no-guide and non-targeting controls.",
            "Sequence the target and nominated off-target loci before interpreting phenotype.",
        ],
    }
    (args.output / "challenge_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(json.dumps({"manifest": str(args.output / "challenge_manifest.json"),
                      "candidates": len(candidates), "exported": len(exported),
                      "decision": manifest["transfer_decision"]}, indent=2))


if __name__ == "__main__":
    main()
