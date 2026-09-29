#!/usr/bin/env python3
"""Create a compact competition figure from frozen validation outputs."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
GSE = ROOT / "data/external_validation/processed/YU2024_GSE196911_PURINE/validation_bundle"
HAL = ROOT / "data/external_validation/processed/CRISPRHAL_TEV_RELEASED_HOLDOUT"
OUT = ROOT / "figures"


def main() -> None:
    gse_metrics = json.loads((GSE / "validation_metrics.json").read_text())
    hal_metrics = json.loads((HAL / "validation_metrics.json").read_text())
    hal = pd.read_csv(HAL / "frozen_predictions_and_outcomes.tsv", sep="\t")
    per_gene = pd.DataFrame(gse_metrics["analysis"]["per_gene_spearman"]).sort_values("spearman_rho")

    navy, green, orange, red, gray = "#173B57", "#00A67D", "#E66A00", "#C73E3A", "#647481"
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.25), gridspec_kw={"width_ratios": [1.05, 1.35, 0.95]})

    colors = [red if value < 0 else green for value in per_gene.spearman_rho]
    axes[0].barh(per_gene.gene_id, per_gene.spearman_rho, color=colors)
    axes[0].axvline(0, color=gray, lw=0.8)
    axes[0].axvline(gse_metrics["random_effects_meta_analysis"]["pooled_spearman"], color=navy, lw=2, ls="--")
    axes[0].set_xlim(-0.4, 0.8); axes[0].set_xlabel("Within-gene Spearman ρ")
    axes[0].set_title("A  Heterogeneity is visible", loc="left", color=navy, fontweight="bold")
    axes[0].text(0.02, 0.02, "Dashed: random-effects ρ=0.484\nI²=75.2%; purE retained", transform=axes[0].transAxes, fontsize=8, color=gray)

    hb = axes[1].hexbin(hal.prediction, hal.outcome, gridsize=45, mincnt=1, cmap="viridis", bins="log")
    axes[1].set_xlabel("Frozen sequence-model score"); axes[1].set_ylabel("Held-out activity")
    axes[1].set_title("B  5,049 outcome-separated guides", loc="left", color=navy, fontweight="bold")
    axes[1].text(0.03, 0.95, "ρ=0.633  |  95% CI 0.613–0.650\npermutation p=1/2,001", transform=axes[1].transAxes,
                 va="top", fontsize=9, color=navy, bbox={"facecolor": "white", "alpha": .86, "edgecolor": "none"})
    cb = fig.colorbar(hb, ax=axes[1], fraction=.05, pad=.03); cb.set_label("log₁₀ bin count", fontsize=8)

    similarity = hal_metrics["sequence_similarity_audit"]
    labels = ["GC-only\nbaseline", "All test\nguides", "Distant\nsubset"]
    values = [hal_metrics["gc_baseline"]["spearman_rho"], hal_metrics["primary"]["spearman_rho"],
              similarity["distant_subset_metrics"]["spearman_rho"]]
    bars = axes[2].bar(labels, values, color=[gray, orange, green], width=.66)
    axes[2].set_ylim(0, .72); axes[2].set_ylabel("Spearman ρ")
    axes[2].set_title("C  Leakage stress test", loc="left", color=navy, fontweight="bold")
    for bar, value in zip(bars, values):
        axes[2].text(bar.get_x()+bar.get_width()/2, value+.018, f"{value:.3f}", ha="center", fontweight="bold", color=navy)
    axes[2].text(.5, .06, "Distant subset: 5,031 guides\nnearest training identity ≤17/20", transform=axes[2].transAxes,
                 ha="center", fontsize=8, color=gray)

    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False); ax.grid(axis="y", alpha=.18); ax.tick_params(labelsize=8)
    fig.suptitle("AgriPAM-AI validation: performance, heterogeneity and leakage controls", x=.055, ha="left",
                 fontsize=15, fontweight="bold", color=navy)
    fig.text(.055, .005, "Frozen predictions; deposited outcomes joined afterward. No panel represents B26 wet-lab validation.", fontsize=8, color=gray)
    fig.tight_layout(rect=[0.03, .045, 1, .92])
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / "external_validation_robustness.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "external_validation_robustness.svg", bbox_inches="tight")


if __name__ == "__main__":
    main()
