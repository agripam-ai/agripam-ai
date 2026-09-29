#!/usr/bin/env python3
"""Worked example of the SynCom selector: write results/syncom_selector/ from data/syncom/laboratory_bank_input.xlsx.

Uses the current model: fungi are scored community members (not only compatibility partners), so a fungal anchor
can also appear in the chassis ranking. The top bacterium is reported separately because editing a fungus needs
its own biosafety and precedent review.
"""
import csv
import json
import sys
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT))
from agripam import syncom_io  # noqa: E402

BANK = PROJECT / "data/syncom/laboratory_bank_input.xlsx"
OUT = PROJECT / "results/syncom_selector"
SCENARIOS = {  # name: (bacterial members or None to assemble, fungal anchors, size)
    "laboratory_design_FI20_B1_B39_B56": (["B1", "B39", "B56"], ["FI20"], 4),
    "safety_screened_assembly": (None, ["FI20"], 4),
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
    sheets = syncom_io.read_workbook(BANK)
    bank, issues = syncom_io.build_bank(sheets)
    if bank is None:
        raise SystemExit("Bank did not load: " + "; ".join(i["message"] for i in issues))
    summary = {}
    for name, (members, fungi, size) in SCENARIOS.items():
        res = syncom_io.analyze(bank, members, fungi, size)
        ranking = res["ranking"]
        write(f"{name}_chassis_ranking.tsv", ranking)
        write(f"{name}_gap_analysis.tsv", res["gaps"])
        bacteria = [m for m in res["members"] if res["kinds"][m] == "bacterium"]
        ranked = [r for r in ranking if r.get("rank")]
        top_bacterium = next((r["strain"] for r in ranked if res["kinds"][r["strain"]] == "bacterium"), None)
        summary[name] = {
            "bacterial_members": bacteria,
            "fungal_members": [m for m in res["members"] if res["kinds"][m] == "fungus"],
            "functions_delivered": res["delivered"],
            "functions_missing": res["missing"],
            "top_chassis": next((r["strain"] for r in ranked if r["rank"] == 1), None),
            "top_bacterial_chassis": top_bacterium,
            "chassis_excluded_for_safety": [r["strain"] for r in ranking if not r["safety_eligible"]],
            "top_choice_frequency_under_random_weights": res["sensitivity"],
            "warnings": res["warnings"],
        }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
