#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
FIGURES = PROJECT / "figures"

ROC_INPUT = FIGURES / "prism_ml_roc_curves.tsv"
PR_INPUT = FIGURES / "prism_ml_precision_recall_curves.tsv"

ROC_OUTPUT = FIGURES / "prism_ml_roc_xy.tsv"
PR_OUTPUT = FIGURES / "prism_ml_precision_recall_xy.tsv"

MODELS = [
    "pam_only",
    "mobile_only",
    "pam_plus_mobile",
]

DISPLAY_NAMES = {
    "pam_only": "PAM_only",
    "mobile_only": "Mobile_only",
    "pam_plus_mobile": "PAM_plus_mobile",
}

roc = pd.read_csv(ROC_INPUT, sep="\t")
pr = pd.read_csv(PR_INPUT, sep="\t")


def padded_series(values, length):
    values = list(values)
    return pd.Series(
        values + [np.nan] * (length - len(values))
    )


roc_lengths = [
    len(roc[roc["model"] == model])
    for model in MODELS
]
roc_length = max(max(roc_lengths), 2)

roc_wide = pd.DataFrame(index=range(roc_length))

for model in MODELS:
    subset = (
        roc[roc["model"] == model]
        .sort_values("point")
    )
    name = DISPLAY_NAMES[model]

    roc_wide[f"{name}_FPR"] = padded_series(
        subset["false_positive_rate"],
        roc_length,
    )
    roc_wide[f"{name}_TPR"] = padded_series(
        subset["true_positive_rate"],
        roc_length,
    )

roc_wide["Chance_FPR"] = padded_series(
    [0.0, 1.0],
    roc_length,
)
roc_wide["Chance_TPR"] = padded_series(
    [0.0, 1.0],
    roc_length,
)

pr_lengths = [
    len(pr[pr["model"] == model])
    for model in MODELS
]
pr_length = max(max(pr_lengths), 2)

pr_wide = pd.DataFrame(index=range(pr_length))

for model in MODELS:
    subset = (
        pr[pr["model"] == model]
        .sort_values("recall")
    )
    name = DISPLAY_NAMES[model]

    pr_wide[f"{name}_Recall"] = padded_series(
        subset["recall"],
        pr_length,
    )
    pr_wide[f"{name}_Precision"] = padded_series(
        subset["precision"],
        pr_length,
    )

baseline = 15 / 1515

pr_wide["Prevalence_Recall"] = padded_series(
    [0.0, 1.0],
    pr_length,
)
pr_wide["Prevalence_Precision"] = padded_series(
    [baseline, baseline],
    pr_length,
)

roc_wide.to_csv(
    ROC_OUTPUT,
    sep="\t",
    index=False,
    na_rep="",
)

pr_wide.to_csv(
    PR_OUTPUT,
    sep="\t",
    index=False,
    na_rep="",
)

print("ROC XY rows:", len(roc_wide))
print("ROC XY columns:", len(roc_wide.columns))
print("PR XY rows:", len(pr_wide))
print("PR XY columns:", len(pr_wide.columns))
print("ROC Prism table:", ROC_OUTPUT)
print("PR Prism table:", PR_OUTPUT)
