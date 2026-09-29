#!/usr/bin/env python3
"""Redraw the 3-panel 5'-TTC evidence figure from the counts in results/ (clean labels)."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "deliverables/assets/type_ic_pam_summary_v2.png"
GREEN, BLUE, GREY, RED = "#2E9E6B", "#2B7BB9", "#B8C2BF", "#D9534F"
fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.6), gridspec_kw={"width_ratios": [0.9, 1.25, 1.25]})
ax = axes[0]
ax.axis("off")
ax.set_title("A   Inferred Type I-C PAM", loc="left", fontsize=13, fontweight="bold")
for i, (letter, colr) in enumerate((("T", RED), ("T", RED), ("C", BLUE))):
    ax.add_patch(plt.Rectangle((0.22 + i * 0.2, 0.5), 0.17, 0.22, color=colr, transform=ax.transAxes))
    ax.text(0.305 + i * 0.2, 0.61, letter, color="white", fontsize=24, fontweight="bold", ha="center", va="center", transform=ax.transAxes)
ax.text(0.05, 0.61, "5′", fontsize=14, va="center", transform=ax.transAxes)
ax.text(0.85, 0.61, "→", fontsize=16, va="center", transform=ax.transAxes)
ax.text(0.5, 0.3, "16/16 unique contexts\n11/11 spacer groups", ha="center", fontsize=12, transform=ax.transAxes)
ax.text(0.5, 0.12, "orientation normalised", ha="center", fontsize=10, color="#5E7A74", transform=ax.transAxes)


def bars(ax, title, labels, vals, texts, colr, ylabel):
    ax.bar(range(4), vals, color=[colr] * 3 + [GREY], width=0.7)
    for i, t in enumerate(texts):
        ax.text(i, vals[i] + 2, t, ha="center", fontsize=10.5)
    ax.set_xticks(range(4))
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylim(0, 118)
    ax.set_ylabel(ylabel, fontsize=10.5)
    ax.set_title(title, loc="left", fontsize=13, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)


cats = ["Accession\nobservations", "Unique\ncontexts", "Unique\nspacer groups", "Random\nbackground"]
bars(axes[1], "B   TTC frequency", cats, [100, 100, 100, 100 * 35 / 2100], ["23/23", "16/16", "11/11", "35/2100"], BLUE, "Sequences with TTC (%)")
bars(axes[2], "C   Mobile-region association", cats[:3] + ["Random\nloci"], [100 * 20 / 21, 100 * 14 / 15, 90, 100 * 502 / 2100],
     ["20/21", "14/15", "9/10", "502/2100"], GREEN, "Near mobile-element features (%)")
fig.text(0.5, 0.005, "P. polymyxa reference panel · BH-adjusted TTC enrichment q = 6.755e-34 · matched random controls, seed 42",
         ha="center", fontsize=10, color="#5E7A74")
plt.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig(OUT, dpi=200)
print(OUT)
