#!/usr/bin/env python3

from pathlib import Path
import sys

import numpy as np
import pandas as pd

from sklearn.compose import ColumnTransformer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    confusion_matrix,
    roc_auc_score,
)
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder


SEED = 42
N_SPLITS = 5

PROJECT = Path(__file__).resolve().parents[1]
INPUT = PROJECT / "machine_learning/type_ic_context_model_dataset.tsv"
OUTPUT_DIR = PROJECT / "machine_learning"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

FOLD_METRICS = OUTPUT_DIR / "model_cross_validation_metrics.tsv"
SUMMARY_METRICS = OUTPUT_DIR / "model_cross_validation_summary.tsv"
OOF_FILE = OUTPUT_DIR / "model_oof_predictions.tsv"
COEFFICIENT_FILE = OUTPUT_DIR / "model_coefficients.tsv"

MODEL_SPECS = {
    "pam_only": {
        "categorical": ["pam_minus3", "pam_minus2", "pam_minus1"],
        "numeric": [],
    },
    "mobile_only": {
        "categorical": [],
        "numeric": ["mobile_nearby"],
    },
    "pam_plus_mobile": {
        "categorical": ["pam_minus3", "pam_minus2", "pam_minus1"],
        "numeric": ["mobile_nearby"],
    },
}


def build_model(categorical, numeric):
    transformers = []

    if categorical:
        transformers.append(
            (
                "categorical",
                OneHotEncoder(
                    categories=[["A", "C", "G", "T"]] * len(categorical),
                    drop="first",
                    handle_unknown="ignore",
                    sparse_output=False,
                ),
                categorical,
            )
        )

    if numeric:
        transformers.append(
            (
                "numeric",
                "passthrough",
                numeric,
            )
        )

    preprocessing = ColumnTransformer(
        transformers=transformers,
        remainder="drop",
        verbose_feature_names_out=False,
    )

    classifier = LogisticRegression(
        class_weight="balanced",
        solver="liblinear",
        max_iter=5000,
        random_state=SEED,
    )

    return Pipeline(
        [
            ("preprocessing", preprocessing),
            ("classifier", classifier),
        ]
    )


def calculate_metrics(y_true, probabilities):
    predictions = (probabilities >= 0.5).astype(int)
    tn, fp, fn, tp = confusion_matrix(
        y_true,
        predictions,
        labels=[0, 1],
    ).ravel()

    sensitivity = tp / (tp + fn) if (tp + fn) else np.nan
    specificity = tn / (tn + fp) if (tn + fp) else np.nan

    return {
        "roc_auc": roc_auc_score(y_true, probabilities),
        "average_precision": average_precision_score(
            y_true, probabilities
        ),
        "balanced_accuracy": balanced_accuracy_score(
            y_true, predictions
        ),
        "sensitivity_at_0.5": sensitivity,
        "specificity_at_0.5": specificity,
        "true_positive": tp,
        "false_positive": fp,
        "true_negative": tn,
        "false_negative": fn,
    }


data = pd.read_csv(INPUT, sep="\t")

required = {
    "example_id",
    "label",
    "query_group",
    "context_id",
    "pam_minus3",
    "pam_minus2",
    "pam_minus1",
    "mobile_nearby",
}

missing = required - set(data.columns)
if missing:
    raise SystemExit(
        "Missing required columns: " + ", ".join(sorted(missing))
    )

if data["example_id"].duplicated().any():
    raise SystemExit("Duplicate example_id values detected.")

data["label"] = data["label"].astype(int)
data["mobile_nearby"] = data["mobile_nearby"].astype(int)

if sorted(data["label"].unique()) != [0, 1]:
    raise SystemExit("The label column must contain both 0 and 1.")

if data["query_group"].nunique() < N_SPLITS:
    raise SystemExit("Insufficient query groups for five-fold CV.")

cv = GroupKFold(n_splits=N_SPLITS)

all_fold_metrics = []
all_predictions = []
all_coefficients = []

for model_name, specification in MODEL_SPECS.items():
    categorical = specification["categorical"]
    numeric = specification["numeric"]
    feature_columns = categorical + numeric

    oof_probability = np.full(len(data), np.nan)
    fold_assignment = np.full(len(data), -1, dtype=int)

    for fold, (train_index, test_index) in enumerate(
        cv.split(
            data[feature_columns],
            data["label"],
            groups=data["query_group"],
        ),
        start=1,
    ):
        train_groups = set(
            data.iloc[train_index]["query_group"]
        )
        test_groups = set(
            data.iloc[test_index]["query_group"]
        )

        if train_groups & test_groups:
            raise SystemExit(
                f"Group leakage detected in {model_name}, fold {fold}."
            )

        model = build_model(categorical, numeric)
        model.fit(
            data.iloc[train_index][feature_columns],
            data.iloc[train_index]["label"],
        )

        probabilities = model.predict_proba(
            data.iloc[test_index][feature_columns]
        )[:, 1]

        oof_probability[test_index] = probabilities
        fold_assignment[test_index] = fold

        fold_result = calculate_metrics(
            data.iloc[test_index]["label"].to_numpy(),
            probabilities,
        )

        fold_result.update(
            {
                "model": model_name,
                "fold": fold,
                "training_rows": len(train_index),
                "test_rows": len(test_index),
                "training_groups": len(train_groups),
                "test_groups": len(test_groups),
                "test_positive": int(
                    data.iloc[test_index]["label"].sum()
                ),
                "test_negative": int(
                    len(test_index)
                    - data.iloc[test_index]["label"].sum()
                ),
            }
        )
        all_fold_metrics.append(fold_result)

    if np.isnan(oof_probability).any():
        raise SystemExit(
            f"Missing OOF predictions for model {model_name}."
        )

    aggregate = calculate_metrics(
        data["label"].to_numpy(),
        oof_probability,
    )
    aggregate.update(
        {
            "model": model_name,
            "fold": "all_oof",
            "training_rows": "NA",
            "test_rows": len(data),
            "training_groups": "NA",
            "test_groups": data["query_group"].nunique(),
            "test_positive": int(data["label"].sum()),
            "test_negative": int(
                len(data) - data["label"].sum()
            ),
        }
    )
    all_fold_metrics.append(aggregate)

    prediction_table = data[
        [
            "example_id",
            "label",
            "query_group",
            "context_id",
            "example_type",
            "pam3",
            "mobile_nearby",
        ]
    ].copy()

    prediction_table["model"] = model_name
    prediction_table["cv_fold"] = fold_assignment
    prediction_table["oof_probability"] = oof_probability
    prediction_table["predicted_label_0.5"] = (
        oof_probability >= 0.5
    ).astype(int)

    all_predictions.append(prediction_table)

    final_model = build_model(categorical, numeric)
    final_model.fit(data[feature_columns], data["label"])

    feature_names = (
        final_model.named_steps["preprocessing"]
        .get_feature_names_out()
    )
    coefficients = (
        final_model.named_steps["classifier"].coef_[0]
    )

    for feature, coefficient in zip(
        feature_names, coefficients
    ):
        all_coefficients.append(
            {
                "model": model_name,
                "feature": feature,
                "coefficient": coefficient,
                "odds_ratio": np.exp(coefficient),
            }
        )

    intercept = final_model.named_steps[
        "classifier"
    ].intercept_[0]

    all_coefficients.append(
        {
            "model": model_name,
            "feature": "intercept",
            "coefficient": intercept,
            "odds_ratio": np.exp(intercept),
        }
    )

metrics = pd.DataFrame(all_fold_metrics)
metrics.to_csv(FOLD_METRICS, sep="\t", index=False)

summary = metrics.loc[
    metrics["fold"].astype(str) == "all_oof"
].copy()

summary[
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
].to_csv(SUMMARY_METRICS, sep="\t", index=False)

pd.concat(
    all_predictions,
    ignore_index=True,
).to_csv(OOF_FILE, sep="\t", index=False)

pd.DataFrame(all_coefficients).to_csv(
    COEFFICIENT_FILE,
    sep="\t",
    index=False,
)

print("Input rows:", len(data))
print("Positive contexts:", int(data["label"].sum()))
print("Query groups:", data["query_group"].nunique())
print("Cross-validation:", f"{N_SPLITS}-fold GroupKFold")
print()
print(
    summary[
        [
            "model",
            "roc_auc",
            "average_precision",
            "balanced_accuracy",
            "sensitivity_at_0.5",
            "specificity_at_0.5",
        ]
    ].to_string(index=False)
)
print()
print("Metrics:", FOLD_METRICS)
print("Summary:", SUMMARY_METRICS)
print("OOF predictions:", OOF_FILE)
print("Coefficients:", COEFFICIENT_FILE)
