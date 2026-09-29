#!/usr/bin/env python3
"""Create the SynCom input template and a small synthetic example workbook.

The example strains (S1-S6) and fungus F1 are invented for demonstration only.
Usage: python scripts/make_syncom_template.py
"""
from pathlib import Path

from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.worksheet.datavalidation import DataValidation

OUT = Path(__file__).resolve().parents[1] / "data/syncom"
HEAD = PatternFill("solid", fgColor="2E9E6B")
NOTE = PatternFill("solid", fgColor="FCF0DA")
TRAITS = ["catalase", "P_solubilization", "K_solubilization", "N_free_growth", "glucanase", "chitinase", "biocontrol_breadth"]
STRAIN_COLS = ["strain", "kind", "identification", "biosafety_hold", "pathogen_screen_score", "editing_precedent"] + TRAITS

INSTRUCTIONS = [
    ("How to use this workbook", ""),
    ("1. Strains", "One row per strain. 'strain' is your code (for example B12). 'kind' is bacterium or fungus: fungi with function scores are scored community members. 'identification' is optional and is never used in the analysis."),
    ("", "Score each function 0-5 (0 = none, 5 = strongest). Leave a cell blank if it was not tested: blank is 'not tested', never 'negative'."),
    ("", "You can rename, delete or add function columns: every extra column is treated as a function scored 0-5."),
    ("", "biosafety_hold: yes or no. Strains marked yes, or with pathogen_screen_score of 30 or more, are never ranked as editing chassis."),
    ("", "editing_precedent: species, genus or none (published genome editing in the same species, the same genus, or neither)."),
    ("2. Fungi_compatibility", "Growth of each bacterium with each beneficial fungus: +++ or ++ compatible, + touch only, -+ partial inhibition, - inhibited."),
    ("3. Bacteria_compatibility", "Pairs tested together: outcome 'mixes' or 'inhibits'. Untested pairs stay blank."),
    ("", "Fungus_compatibility: two fungi against each other (same symbols as Fungi_compatibility). Two fungi that inhibit each other trigger a warning."),
    ("4. Community (optional)", "List the members of a planned community (type bacterium or fungus). Leave empty to let the software assemble one."),
    ("Rows starting with EXAMPLE", "are ignored. Delete them and enter your own data."),
    ("Interpretation", "Results are decision support: not measured editing efficiency, safety clearance or field performance."),
]


def sheet(wb, title, header, rows, widths=None):
    ws = wb.create_sheet(title)
    ws.append(header)
    for cell in ws[1]:
        cell.fill, cell.font = HEAD, Font(bold=True, color="FFFFFF")
        cell.alignment = Alignment(wrap_text=True, vertical="center")
    for row in rows:
        ws.append(row)
    for i, col in enumerate(header, 1):
        ws.column_dimensions[ws.cell(1, i).column_letter].width = (widths or {}).get(col, max(14, len(col) + 3))
    ws.freeze_panes = "B2"
    return ws


def build(path, strains, fungi, pairs, community, example_note, fungus_pairs=None):
    wb = Workbook()
    ws = wb.active
    ws.title = "Instructions"
    for row in INSTRUCTIONS:
        ws.append(row)
    if example_note:
        ws.append(("Example data", example_note))
    ws.column_dimensions["A"].width = 28
    ws.column_dimensions["B"].width = 120
    ws["A1"].font = Font(bold=True, size=14)
    for r in ws.iter_rows(min_row=2):
        r[0].font = Font(bold=True)
        r[1].alignment = Alignment(wrap_text=True, vertical="top")
    st = sheet(wb, "Strains", STRAIN_COLS, strains, {"strain": 16, "identification": 24, "editing_precedent": 18})
    dv_yn = DataValidation(type="list", formula1='"yes,no"', allow_blank=True)
    dv_pre = DataValidation(type="list", formula1='"species,genus,none"', allow_blank=True)
    dv_score = DataValidation(type="whole", operator="between", formula1="0", formula2="5", allow_blank=True)
    for dv in (dv_yn, dv_pre, dv_score):
        st.add_data_validation(dv)
    dv_kind = DataValidation(type="list", formula1='"bacterium,fungus"', allow_blank=True)
    st.add_data_validation(dv_kind)
    dv_kind.add("B2:B500")
    dv_yn.add("D2:D500")
    dv_pre.add("F2:F500")
    dv_score.add("G2:Z500")
    fc = sheet(wb, "Fungi_compatibility", ["strain", "fungus", "symbol"], fungi, {"strain": 16, "fungus": 16, "symbol": 10})
    dv_sym = DataValidation(type="list", formula1='"+++,++,+,-+,-"', allow_blank=True)
    fc.add_data_validation(dv_sym)
    dv_sym.add("C2:C500")
    bc = sheet(wb, "Bacteria_compatibility", ["strain_a", "strain_b", "outcome"], pairs, {"strain_a": 16, "strain_b": 16, "outcome": 12})
    dv_out = DataValidation(type="list", formula1='"mixes,inhibits"', allow_blank=True)
    bc.add_data_validation(dv_out)
    dv_out.add("C2:C500")
    ffs = sheet(wb, "Fungus_compatibility", ["fungus_a", "fungus_b", "symbol"], fungus_pairs or [], {"fungus_a": 16, "fungus_b": 16, "symbol": 10})
    dv_sym2 = DataValidation(type="list", formula1='"+++,++,+,-+,-"', allow_blank=True)
    ffs.add_data_validation(dv_sym2)
    dv_sym2.add("C2:C500")
    sheet(wb, "Community", ["member", "type"], community, {"member": 16, "type": 14})
    wb.save(path)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    blank_strains = [["EXAMPLE_delete_this_row", "bacterium", "optional free text", "no", 0, "none", 3, 4, 3, 2, 0, 1, 2],
                     ["", "", "", "", "", "", "", "", "", "", "", "", ""]]
    build(OUT / "SynCom_input_template.xlsx", blank_strains,
          [["EXAMPLE_delete_this_row", "F1", "++"]], [["EXAMPLE_delete_this_row", "EXAMPLE_other", "mixes"]],
          [["EXAMPLE_delete_this_row", "bacterium"]], "", [["EXAMPLE_delete_this_row", "EXAMPLE_other", "-"]])
    b = "bacterium"
    strains = [
        ["S1", b, "synthetic", "no", 5, "genus", 3, 4, 3, 3, 0, 1, 0],
        ["S2", b, "synthetic", "no", 2, "none", 0, 4, 3, 3, 1, 3, 1],
        ["S3", b, "synthetic", "no", 0, "species", 3, 5, 3, None, 0, 1, 0],
        ["S4", b, "synthetic", "no", 4, "species", 1, 0, 1, 0, 5, 5, 4],
        ["S5", b, "synthetic", "yes", 54, "none", 3, 4, 4, 3, 0, 1, 2],
        ["S6", b, "synthetic", "no", 3, "species", 2, 1, 1, 1, 4, 2, 2],
        ["F1", "fungus", "synthetic", "no", 2, "genus", None, None, None, None, None, None, 5],
        ["F2", "fungus", "synthetic", "no", 3, "none", None, None, None, None, None, None, 4],
    ]
    fungi = [[s, "F1", sym] for s, sym in (("S1", "++"), ("S2", "++"), ("S3", "++"), ("S4", "-"), ("S5", "++"), ("S6", "++"))]
    fungi += [[s, "F2", sym] for s, sym in (("S1", "++"), ("S2", "+"), ("S3", "++"))]
    pairs = [["S1", "S2", "mixes"], ["S2", "S3", "mixes"], ["S1", "S6", "mixes"], ["S3", "S6", "inhibits"]]
    build(OUT / "SynCom_example_filled.xlsx", strains, fungi, pairs, [],
          "S1-S6 (bacteria) and F1, F2 (fungi) are invented strains for demonstration only. They are not real data.",
          [["F1", "F2", "-"]])
    print("wrote", OUT / "SynCom_input_template.xlsx", "and", OUT / "SynCom_example_filled.xlsx")


if __name__ == "__main__":
    main()
