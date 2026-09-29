#!/usr/bin/env python3
"""Build the results tables of the manuscript (Tables 1-7) as one Excel workbook.

Every number is read from the project's result files (the same sources as the manuscript draft).
Usage: python scripts/build_manuscript_tables.py [output.xlsx]
"""
import json
import sys
from pathlib import Path

import pandas as pd
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_manuscript_draft as draft  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
WINDOWS = [
    ("1 Start here", "Orient the user; show why a PAM alone is not enough; a ten-step project plan", "Choices on screen", "On-screen plan"),
    ("2 Genome evaluation", "One genome from assembly to ranked edit sites and a multi-objective design portfolio", "FASTA, reads or NCBI accession; optional mobile-element reference", "Tables, logs and version record in one ZIP"),
    ("3 SynCom candidate bank", "Community design, gap analysis, chassis ranking; native-bank evidence and PAM reference panels", "Bank workbook (strains, traits, compatibility)", "Community, gaps, chassis ranking, robustness, results workbook"),
    ("4 Agricultural editing knowledgebase", "What is known about editing agricultural microbes; protocol planner", "Search terms; optional reviewed records", "Filtered evidence, protocol"),
    ("5 Parts & constructs", "Assemble and screen a reporter or CRISPRi construct; follow the experiment", "Parts (starter set, Excel, iGEM Registry); neutral region or guide", "FASTA, GenBank, PCR sizes, tracker workbook"),
    ("6 External validation", "Check predictions against other laboratories' data; evaluate a frozen held-out table", "Predictions with outcomes joined afterwards", "Correlations, controls, report"),
    ("7 Software & reproducibility", "Editing-system atlas; reproducibility contract", "None", "Versions, seeds, checksums, commands"),
]


def main(out: Path):
    f = draft.facts()
    n = f["names"]
    t1 = pd.DataFrame(WINDOWS, columns=["Window", "Purpose", "Input", "Output"])
    s2 = f["s2"].head(10).copy()
    s2.insert(1, "species", [n.get(s, "") for s in s2["strain"]])
    t2 = s2.rename(columns={"rank": "Rank", "strain": "Strain", "species": "Species", "total_score": "Total score", "mean_intensity": "Mean intensity",
                             "positive_function_groups (of 7)": "Function groups (of 7)"})[["Rank", "Strain", "Species", "Total score", "Mean intensity", "Function groups (of 7)"]]
    t3 = f["summary"].reset_index()[["scenario", "description", "members", "functions_missing", "top_chassis", "not_ranked_biosafety", "robustness", "warnings"]]
    t3.columns = ["Scenario", "Description", "Members", "Functions missing", "Top chassis", "Not ranked (biosafety flag)", "Robustness (2,000 weightings)", "Warnings"]
    rows = []
    for tag in "ABCD":
        rk = f[f"{tag}_ranking"]
        for r in rk.itertuples():
            rows.append((tag, r.strain, n.get(r.strain, ""), getattr(r, "kind", ""), r.fit_score, r.dispensability_score, r.editability_precedent_score, r.composite,
                         "" if pd.isna(r.rank) else int(r.rank), "yes" if r.safety_eligible else "no (held)"))
    t4 = pd.DataFrame(rows, columns=["Scenario", "Strain", "Species", "Kind", "Fit", "Dispensability", "Editing precedent", "Composite", "Rank", "Eligible as chassis"])
    t5 = pd.DataFrame([
        ("Target observations with 5'-TTC", f["ttc_observations"], "accession-coordinate observations"),
        ("Unique contexts", f"{f['contexts']}/{f['contexts']}", ""),
        ("Spacer groups", f"{f['spacer_groups']}/{f['spacer_groups']}", "conservative unit"),
        ("Matched random loci with 5'-TTC", f["background"], ""),
        ("Benjamini-Hochberg q", f"{f['q_value']:.2e}", "over 64 triplets"),
        ("ROC AUC (grouped cross-validation)", f["auc"], ""),
        ("Average precision", round(f["ap"], 2), "0.99% positive prevalence"),
        ("Leave-one-out folds (spacer groups, genomes)", f"{f['heldout_folds'][0]}, {f['heldout_folds'][1]}", ""),
        ("Held-out observations matched", f["heldout"][0], "in both schemes"),
    ], columns=["Measure", "Value", "Note"])
    t6 = pd.DataFrame([
        ("Yu et al. CRISPRi screen, E. coli", f["yu_n"], f["yu_rho"], f"within-gene median {f['yu_median']}; permutation {f['yu_null']}"),
        ("Transfer E. coli to B. subtilis (GFP)", f["hk_n"], f["hk_rho"], f"95% interval {f['hk_lo']} to {f['hk_hi']}"),
        ("Bacterial Cas9 released holdout (crisprHAL)", f["ch_n"], f["ch_rho"], f"GC-only baseline {f['ch_base']}"),
    ], columns=["Benchmark", "Guides", "Spearman rho", "Controls and interval"])
    sb = json.loads((ROOT / "manuscript/analysis/software_benchmark.json").read_text())
    ops = pd.DataFrame(sb["operations"]).rename(columns={"operation": "Operation", "seconds": "Seconds", "note": "Note"})
    t7 = ops[["Operation", "Seconds", "Note"]]
    sheets = [("Table1_windows", t1), ("Table2_biochem_top10", t2), ("Table3_communities", t3), ("Table4_chassis_ranking", t4),
              ("Table5_PAM_evidence", t5), ("Table6_validation", t6), ("Table7_performance", t7)]
    out.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(out, engine="openpyxl") as w:
        legends = pd.DataFrame([
            ("Table1_windows", "The seven windows of the application: purpose, input and output."),
            ("Table2_biochem_top10", "Top ten bacterial isolates by total biochemical score (Appendix S2)."),
            ("Table3_communities", "Communities assembled by the software in four scenarios (A: fungal anchor plus assembled bacteria; B: bacteria only; C: laboratory-planned community; D: two fungal anchors)."),
            ("Table4_chassis_ranking", "Chassis ranking in each scenario: fit, dispensability, editing precedent (0-100), composite = 0.40 fit + 0.35 dispensability + 0.25 precedent. Strains with a biosafety flag are not ranked."),
            ("Table5_PAM_evidence", "Evidence for the inferred native 5'-TTC PAM in the P. polymyxa reference panel (24 genomes, 10 with CRISPR-Cas)."),
            ("Table6_validation", "Frozen-prediction benchmarks of the guide-ranking components on deposited data."),
            ("Table7_performance", f"Timings on {sb['machine']}."),
        ], columns=["Sheet", "Legend"])
        legends.to_excel(w, sheet_name="README_legends", index=False)
        for name, df in sheets:
            df.to_excel(w, sheet_name=name, index=False)
        for ws in w.book.worksheets:
            for c in ws[1]:
                c.font = Font(bold=True, color="FFFFFF"); c.fill = PatternFill("solid", fgColor="0E3B36")
            for i, col in enumerate(ws.columns, 1):
                width = min(60, max(len(str(c.value)) if c.value is not None else 0 for c in col) + 2)
                ws.column_dimensions[get_column_letter(i)].width = max(10, width)
                for c in col:
                    c.alignment = Alignment(wrap_text=True, vertical="top")
            ws.freeze_panes = "A2"
    print("wrote", out, [n for n, _ in sheets])


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "manuscript/tables/Manuscript_results_tables.xlsx")
