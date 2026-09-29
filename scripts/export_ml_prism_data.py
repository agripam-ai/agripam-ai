#!/usr/bin/env python3

from pathlib import Path

import pandas as pd
from sklearn.metrics import (
    auc,
    precision_recall_curve,
    roc_curve,
)


PROJECT = Path(__file__).resolve().parents[1]

PREDICTIONS = (
    PROJECT
    / "machine_learning/model_oof_predictions.tsv"
)

SUMMARY = (
    PROJECT
    / "machine_learning/model_cross_validation_summary.tsv"
)

OUTPUT_DIR = PROJECT / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

ROC_OUTPUT = OUTPUT_DIR / "prism_ml_roc_curves.tsv"
PR_OUTPUT = OUTPUT_DIR / "prism_ml_precision_recall_curves.tsv"
SUMMARY_OUTPUT = OUTPUT_DIR / "prism_ml_model_summary.tsv"

predictions = pd.read_csv(PREDICTIONS, sep="\t")
summary = pd.read_csv(SUMMARY, sep="\t")

model_order = [
    "pam_only",
    "mobile_only",
    "pam_plus_mobile",
]

roc_rows = []
pr_rows = []

for model_name in model_order:
    subset = predictions[
        predictions["model"] == model_name
    ].copy()

    y_true = subset["label"].astype(int).to_numpy()
    probability = subset["oof_probability"].to_numpy()

    false_positive_rate, true_positive_rate, roc_threshold = (
        roc_curve(y_true, probability)
    )

    for index, (fpr, tpr, threshold) in enumerate(
        zip(
            false_positive_rate,
            true_positive_rate,
            roc_threshold,
        ),
        start=1,
    ):
        roc_rows.append(
            {
                "model": model_name,
                "point": index,
                "false_positive_rate": fpr,
                "true_positive_rate": tpr,
                "threshold": threshold,
            }
        )

    precision, recall, pr_threshold = precision_recall_curve(
        y_true,
        probability,
    )

    thresholds = list(pr_threshold) + [float("nan")]

    for index, (rec, prec, threshold) in enumerate(
        zip(recall, precision, thresholds),
        start=1,
    ):
        pr_rows.append(
            {
                "model": model_name,
                "point": index,
                "recall": rec,
                "precision": prec,
                "threshold": threshold,
            }
        )

roc_table = pd.DataFrame(roc_rows)
pr_table = pd.DataFrame(pr_rows)

roc_table.to_csv(ROC_OUTPUT, sep="\t", index=False)
pr_table.to_csv(PR_OUTPUT, sep="\t", index=False)

summary_table = summary[
    [
        "model",
        "roc_auc",
        "average_precision",
        "balanced_accuracy",
        "sensitivity_at_0.5",
        "specificity_at_0.5",
        "true_positive",
        "false_positive",
        "true_negative",
        "false_negative",
    ]
].copy()

summary_table["positive_prevalence"] = 15 / 1515
summary_table["random_roc_baseline"] = 0.5
summary_table["random_pr_baseline"] = 15 / 1515

summary_table.to_csv(
    SUMMARY_OUTPUT,
    sep="\t",
    index=False,
)

print("Models exported:", len(model_order))
print("ROC points:", len(roc_table))
print("Precision-recall points:", len(pr_table))
print("Positive-prevalence PR baseline:", f"{15/1515:.6f}")
print("ROC data:", ROC_OUTPUT)
print("Precision-recall data:", PR_OUTPUT)
print("Summary data:", SUMMARY_OUTPUT)
