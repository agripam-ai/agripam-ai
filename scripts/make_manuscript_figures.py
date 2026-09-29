#!/usr/bin/env python3
"""Draw the manuscript figures from the project's result files (manuscript/figures/*.png, 300 dpi).

Needs matplotlib and pandas. Figures 1 and 6 are schematics; 2, 3 and 5 plot result tables; 4 reuses the PAM evidence figure.
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from final_check import source_facts  # noqa: E402
from agripam import biobricks as bb, parts_io  # noqa: E402

OUT = ROOT / "manuscript/figures"
DEEP, LEAF, AMBER, BLUE, GREY, PALE = "#0E3B36", "#2E9E6B", "#E8A33D", "#2B7BB9", "#9FB3AD", "#E6F3EC"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False, "axes.spines.right": False})


def box(ax, x, y, w, h, text, fc="white", ec=DEEP, fs=8.5, bold=False, color=DEEP):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.04", fc=fc, ec=ec, lw=1.1))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color=color, fontweight="bold" if bold else "normal", wrap=True)


def arrow(ax, x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=12, color=DEEP, lw=1.1))


def fig1():
    fig, ax = plt.subplots(figsize=(11, 6.4))
    ax.set_xlim(0, 11); ax.set_ylim(0, 6.4); ax.axis("off")
    ax.text(0.1, 6.15, "Inputs", fontweight="bold", color=LEAF); ax.text(3.85, 6.15, "Application windows", fontweight="bold", color=LEAF)
    ax.text(8.3, 6.15, "Outputs", fontweight="bold", color=LEAF)
    for i, t in enumerate(["Bank workbook\n(strains, traits, compatibility)", "Genome\n(FASTA, accession, reads)", "Parts\n(starter set, Excel, iGEM Registry)", "Held-out table\n(predictions + outcomes)"]):
        box(ax, 0.1, 5.05 - i * 1.25, 2.6, 0.85, t, PALE)
    wins = ["1  Start here", "2  Genome evaluation", "3  SynCom candidate bank", "4  Editing knowledgebase", "5  Parts & constructs", "6  External validation", "7  Software & reproducibility"]
    for i, t in enumerate(wins):
        box(ax, 3.85, 5.5 - i * 0.62, 3.6, 0.46, t, "white", fs=8.5, bold=True)
    for i, t in enumerate(["Community, gaps and\nchassis ranking (tables)", "Targets, gene context and\nportfolio (ZIP + manifest)", "Construct maps: FASTA,\nGenBank, tracker (xlsx)", "Validation reports and\nprovenance (json, tsv)"]):
        box(ax, 8.3, 5.05 - i * 1.25, 2.6, 0.85, t, PALE)
    for i in range(4):
        arrow(ax, 2.72, 5.47 - i * 1.25, 3.83, [4.49, 5.11, 3.25, 2.63][i]); arrow(ax, 7.47, [4.49, 5.11, 3.25, 2.63][i], 8.28, 5.47 - i * 1.25)
    ax.add_patch(Rectangle((0.1, 0.05), 10.8, 1.0, fc="#F7FAF8", ec=DEEP, lw=1.1))
    ax.text(0.25, 0.85, "Python library  agripam/", fontweight="bold", color=LEAF, fontsize=9)
    mods = ["syncom", "syncom_io", "core", "workflow", "editor_targeting", "pam_discovery", "multiobjective", "validation", "biobricks", "parts_io"]
    for i, m in enumerate(mods):
        box(ax, 0.25 + i * 1.06, 0.15, 1.0, 0.38, m, "white", ec=GREY, fs=7)
    ax.text(10.8, 0.85, "60 automated tests · fixed seed 42 · SHA-256 checksums · versioned outputs", ha="right", fontsize=8.5, color=DEEP, style="italic")
    fig.savefig(OUT / "Fig1_architecture.png", dpi=300, bbox_inches="tight"); plt.close(fig)


def fig2():
    s3 = pd.read_csv(ROOT / "manuscript/appendices/tsv/S3_bact_vs_pathogens.tsv", sep="\t")
    s3 = s3[s3.score.notna()].copy()
    order = {"week 1": 1, "week 2": 2, "week 3": 3}
    s3["r"] = s3.reading.map(order)
    last = s3.sort_values("r").groupby(["strain", "pathogen"]).tail(1)
    mat = last.pivot(index="strain", columns="pathogen", values="score")
    mat = mat.loc[mat.sum(axis=1, skipna=True).sort_values(ascending=False).index]
    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(mat.values, cmap="YlGn", vmin=0, vmax=4, aspect="auto")
    ax.set_xticks(range(mat.shape[1])); ax.set_xticklabels(mat.columns, rotation=45, ha="right")
    ax.set_yticks(range(mat.shape[0])); ax.set_yticklabels(mat.index, fontsize=7.5)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat.values[i, j]
            if np.isnan(v):
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fc="#EEEEEE", ec="white"))
            else:
                ax.text(j, i, int(v), ha="center", va="center", fontsize=7, color="white" if v >= 3 else DEEP)
    ax.set_title("Antagonism of fungal pathogens by bacterial isolates (last reading, 0–4; grey = not tested)", fontsize=9.5)
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="score")
    for sp in ax.spines.values():
        sp.set_visible(False)
    fig.savefig(OUT / "Fig2_antagonism_heatmap.png", dpi=300, bbox_inches="tight"); plt.close(fig)


def fig3():
    bank = pd.read_excel(ROOT / "manuscript/analysis/software_input_bank.xlsx", sheet_name="Strains").set_index("strain")
    summary = pd.read_excel(ROOT / "manuscript/analysis/analysis_results.xlsx", sheet_name="summary").set_index("scenario")
    funcs = ["P_solubilization", "K_solubilization", "N_free_growth", "catalase", "glucanase", "chitinase", "biocontrol_breadth", "biocontrol_peak"]
    labels = ["P", "K", "N-free", "catalase", "glucanase", "chitinase", "biocontrol\nbreadth", "biocontrol\npeak"]
    fig, axes = plt.subplots(2, 2, figsize=(11, 6.6), gridspec_kw={"width_ratios": [2.3, 1]})
    for row, tag in enumerate(["A", "C"]):
        members = [m.strip() for m in summary.loc[tag, "members"].split(",")]
        mat = bank.loc[members, funcs].apply(pd.to_numeric, errors="coerce")
        ax = axes[row, 0]
        ax.imshow(mat.values.astype(float), cmap="YlGn", vmin=0, vmax=5, aspect="auto")
        ax.set_xticks(range(len(funcs))); ax.set_xticklabels(labels, fontsize=7.5)
        ax.set_yticks(range(len(members))); ax.set_yticklabels([f"{m} ({bank.loc[m, 'kind'][0]})" for m in members])
        for i in range(len(members)):
            for j in range(len(funcs)):
                v = mat.values[i, j]
                ax.text(j, i, "–" if np.isnan(v) else int(v), ha="center", va="center", fontsize=8, color="white" if (not np.isnan(v) and v >= 4) else DEEP)
        ax.set_title(f"Scenario {tag}: {summary.loc[tag, 'description']}\nmissing: {summary.loc[tag, 'functions_missing']}", fontsize=8.5, loc="left")
        for sp in ax.spines.values():
            sp.set_visible(False)
        rk = pd.read_excel(ROOT / "manuscript/analysis/analysis_results.xlsx", sheet_name=f"{tag}_ranking").set_index("strain")
        axb = axes[row, 1]
        comp = rk["composite"].reindex(members)
        held = ~rk["safety_eligible"].reindex(members).astype(bool)
        axb.barh(range(len(members)), comp.values, color=[GREY if h else LEAF for h in held.values])
        axb.set_yticks(range(len(members))); axb.set_yticklabels([]); axb.invert_yaxis(); axb.set_xlim(0, 105)
        axb.set_xlabel("composite score" if row == 1 else "")
        for i, (v, h) in enumerate(zip(comp.values, held.values)):
            axb.text(v + 1.5, i, "not ranked (biosafety)" if h else f"{v:.0f}", va="center", fontsize=7.5)
    fig.suptitle("Community design and chassis ranking (bars: green = eligible; grey = held by the biosafety gate)", fontsize=10, y=0.995)
    fig.tight_layout()
    fig.savefig(OUT / "Fig3_community_design.png", dpi=300, bbox_inches="tight"); plt.close(fig)


def fig5():
    f = source_facts()
    names = ["GSE196911\n(742 guides)", "Cross-chassis\nGFP (33)", "crisprHAL\nholdout (5,049)", "GC-only\nbaseline", "Permutation\ncontrol"]
    vals = [f["yu_rho"], f["hk_rho"], f["ch_rho"], f["ch_base"], f["yu_null"]]
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    ax.bar(names, vals, color=[LEAF, LEAF, LEAF, GREY, GREY])
    for i, v in enumerate(vals):
        ax.text(i, v + (0.02 if v >= 0 else -0.05), f"{v:.3f}", ha="center", fontsize=8.5)
    ax.axhline(0, color="black", lw=0.8); ax.set_ylabel("Spearman ρ (frozen predictions)"); ax.set_ylim(-0.12, 0.75)
    ax.set_title("Independent, frozen-prediction benchmarks", fontsize=10)
    fig.savefig(OUT / "Fig5_validation.png", dpi=300, bbox_inches="tight"); plt.close(fig)


def fig6():
    lib = {p["part_id"]: p for p in parts_io.read_parts_tsv(ROOT / "data/biobricks/starter_library.tsv")}
    cassette = bb.assemble([lib[i] for i in ("BBa_J23100", "BBa_B0034", "BBa_E0040", "BBa_B0015")])
    fig, ax = plt.subplots(figsize=(10, 2.6))
    colors = {"promoter": LEAF, "rbs": AMBER, "reporter": "#5DBA8A", "cds": "#5DBA8A", "terminator": BLUE, "other": GREY}
    total = len(cassette["sequence"])
    for n in cassette["annotations"]:
        w = (n["end"] - n["start"] + 1)
        ax.add_patch(Rectangle((n["start"], 0.4), w, 0.5, fc=colors.get(n["role"], GREY), ec="white"))
        if w > 60:
            ax.text(n["start"] + w / 2, 0.65, n["label"].split(" ")[0] if n["role"] != "other" else "", ha="center", va="center", fontsize=8, color="white")
    labels = [("BBa_J23100\npromoter", 18), ("B0034\nRBS", 62), ("BBa_E0040\nGFP reporter (720 bp)", 450), ("BBa_B0015\nterminator (129 bp)", 855)]
    for t, x in labels:
        y = 0.15 if not t.startswith("B0034") else -0.3
        ax.text(x, y, t, ha="center", va="top", fontsize=7.5, color=DEEP)
        if t.startswith("B0034"):
            ax.plot([x, x], [0.4, y], color=GREY, lw=0.8)
    ax.text(total / 2, 1.15, f"Reporter cassette assembled from four iGEM Registry parts: {total} bp including three 8-bp scars", ha="center", fontsize=9)
    ax.set_xlim(-10, total + 10); ax.set_ylim(-0.9, 1.4); ax.axis("off")
    fig.savefig(OUT / "Fig6_cassette.png", dpi=300, bbox_inches="tight"); plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for fn in (fig1, fig2, fig3, fig5, fig6):
        fn()
        print("drew", fn.__name__)
    shutil.copy2(ROOT / "deliverables/assets/type_ic_pam_summary_v2.png", OUT / "Fig4_pam_evidence.png")


if __name__ == "__main__":
    main()
