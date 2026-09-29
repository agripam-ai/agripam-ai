#!/usr/bin/env python3
"""Score SynCom chassis candidates from data/syncom/ and write results/syncom_selector/."""
import csv
import json
import sys
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT))
from agripam import syncom  # noqa: E402

OUT = PROJECT / "results/syncom_selector"
FUNGI = ["FI20"]
SCENARIOS = {
    "laboratory_design_FI20_B1_B39_B56": ["B1", "B39", "B56"],
}


def write(name, rows):
    OUT.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    keys = list(rows[0])
    with open(OUT / name, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    bank = syncom.load_bank(PROJECT / "data/syncom")
    SCENARIOS["safety_screened_assembly"] = syncom.assemble(bank, FUNGI, size=4)
    summary = {}
    for name, members in SCENARIOS.items():
        ranking = syncom.rank_chassis(members, FUNGI, bank)
        gaps = syncom.gap_analysis(members, FUNGI, bank)
        write(f"{name}_chassis_ranking.tsv", ranking)
        write(f"{name}_gap_analysis.tsv", gaps)
        cov = syncom.coverage(members, bank)
        summary[name] = {
            "bacterial_members": members, "fungal_members": FUNGI,
            "functions_delivered": sorted(syncom.delivered(cov, bank)),
            "functions_missing": [syncom.FUNCTIONS[f] for f in syncom.FUNCTIONS if cov[f] < syncom.FUNCTION_THRESHOLD],
            "top_chassis": next((r["strain"] for r in ranking if r["rank"] == 1), None),
            "chassis_excluded_for_safety": [r["strain"] for r in ranking if not r["safety_eligible"]],
            "top_choice_frequency_under_random_weights": syncom.weight_sensitivity(members, FUNGI, bank),
        }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
