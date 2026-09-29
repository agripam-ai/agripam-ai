#!/usr/bin/env python3
"""Convert the laboratory workbook into tidy SynCom bank tables (data/syncom/).

Only tests present in the workbook are transcribed. Blank cells stay blank
(= not evaluated); they are never converted to a negative result.
Requires openpyxl (not needed by the app). Usage: python scripts/build_syncom_bank.py [workbook.xlsx]
"""
import csv
import sys
from pathlib import Path

import openpyxl

PROJECT = Path(__file__).resolve().parents[1]
OUT = PROJECT / "data/syncom"

# Column groups on the 'Pruebas bioq bacterias' sheet (7-day reading = last column of each block).
TRAITS = {
    "catalase": ["D"], "oxidase": ["E"],
    "P_solubilization": ["H", "K"], "K_solubilization": ["N"],
    "N_free_growth": ["Q", "T", "W"],
    "glucanase": ["Z", "AC"], "chitinase": ["AF", "AI"],
}
PATHOGEN_COLS = {"F7": "N", "F8": "O", "F9": "P", "F10": "Q", "F11": "R", "F14": "S",
                 "F17": "U", "F48": "V", "F53": "W", "R2": "X", "R5": "Y", "R6": "Z"}
FUNGI_COLS = {"FI1": "H", "FI2": "I", "FI20": "J", "FI62": "K"}
BASE_COLS = {"B3": "E", "B33": "F", "B34": "G", "CC25": "H", "54": "I"}


def score(value):
    """'-'=0 ... '+++++'=5 ('´- +' partial =1); None if blank/unreadable."""
    if value is None:
        return None
    s = str(value).replace("´", "").replace(" ", "").strip()
    if s == "":
        return None
    if s in ("-+", "+-"):
        return 1
    if set(s) == {"+"}:
        return len(s)
    if s in ("-", "0"):
        return 0
    return None


def write(name, header, rows):
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / name, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)
    print(f"{name}: {len(rows)} rows")


def main(path):
    wb = openpyxl.load_workbook(path, data_only=True)
    ws = wb.worksheets[6]
    rows, ids = [], {}
    for r in range(7, ws.max_row + 1):
        acr, ident = ws[f"C{r}"].value, ws[f"B{r}"].value
        if acr is None or ident is None:
            continue
        acr = str(acr).strip()
        vals = []
        for trait, cols in TRAITS.items():
            got = [score(ws[f"{c}{r}"].value) for c in cols]
            got = [g for g in got if g is not None]
            vals.append(max(got) if got else "")
        ids[acr] = str(ident).strip()
        rows.append([acr, str(ident).strip()] + vals)
    write("strain_traits.tsv", ["strain", "identification"] + list(TRAITS), rows)

    ws = wb.worksheets[1]
    rows = []
    for r in range(6, 40):
        label = ws[f"M{r}"].value
        if label is None:
            continue
        for path_id, col in PATHOGEN_COLS.items():
            s = score(ws[f"{col}{r}"].value)
            if s is not None:
                rows.append([str(label).strip(), path_id, s, "week2"])
    write("bacteria_pathogen_antagonism.tsv", ["strain", "pathogen", "score_0_4", "reading"], rows)

    ws = wb.worksheets[0]
    rows = []
    for r in range(13, 45):
        label = ws[f"G{r}"].value
        if label is None:
            continue
        for fungus, col in FUNGI_COLS.items():
            raw = ws[f"{col}{r}"].value
            if raw is not None and str(raw).strip():
                rows.append([str(label).strip(), fungus, str(raw).replace("´", "").strip(), "week2"])
    write("bacteria_fungi_compatibility.tsv", ["strain", "fungus", "symbol", "reading"], rows)

    ws = wb.worksheets[3]
    rows = []
    for r in range(4, 40):
        label = ws[f"C{r}"].value
        if label is None:
            continue
        for base, col in BASE_COLS.items():
            rows.append([str(label).strip(), base, "mixes" if ws[f"{col}{r}"].value == "X" else "not_recorded",
                         str(ws[f"K{r}"].value or "")])
    write("bacteria_pairwise_compatibility.tsv", ["strain", "base_strain", "outcome", "note"], rows)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit("Usage: python scripts/build_syncom_bank.py path/to/your_workbook.xlsx")
    main(Path(sys.argv[1]))
