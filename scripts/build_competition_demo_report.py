#!/usr/bin/env python3
"""Build deterministic TSV and JSON summaries used in the competition demo."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "competition_demo"


def count_rows(path: Path) -> int:
    with path.open() as handle:
        return max(sum(1 for _ in handle) - 1, 0)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    context_file = ROOT / "results" / "type_ic_context_level_evidence.tsv"
    pam_file = ROOT / "results" / "type_ic_exact_pam_evidence.tsv"
    model_file = ROOT / "machine_learning" / "model_cross_validation_summary.tsv"

    summary = {
        "project": "AgriPAM-AI",
        "organism": "Paenibacillus polymyxa",
        "system": "Type I-C",
        "inferred_pam": "5'-TTC",
        "functional_validation": "pending",
        "unique_contexts": count_rows(context_file),
        "accession_coordinate_observations": count_rows(pam_file),
        "evidence": {
            "ttc_observations": "23/23",
            "ttc_unique_contexts": "16/16",
            "ttc_spacer_groups": "11/11",
            "matched_background": "35/2100",
            "bh_adjusted_q": 6.755e-34,
        },
        "limitations": [
            "PAM is computationally inferred",
            "Model probabilities are ranking scores, not measured editing efficiencies",
            "Functional validation is pending",
        ],
    }

    json_path = OUTPUT / "agripam_competition_demo.json"
    json_path.write_text(json.dumps(summary, indent=2) + "\n")

    tsv_path = OUTPUT / "agripam_competition_demo.tsv"
    with tsv_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["field", "value"])
        writer.writerow(["organism", summary["organism"]])
        writer.writerow(["system", summary["system"]])
        writer.writerow(["inferred_pam", summary["inferred_pam"]])
        writer.writerow(["functional_validation", summary["functional_validation"]])
        writer.writerow(["unique_contexts", summary["unique_contexts"]])
        writer.writerow(["observations", summary["accession_coordinate_observations"]])

    print(f"JSON report: {json_path}")
    print(f"TSV report: {tsv_path}")


if __name__ == "__main__":
    main()
