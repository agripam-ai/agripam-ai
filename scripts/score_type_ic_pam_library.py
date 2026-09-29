#!/usr/bin/env python3

from itertools import product
from pathlib import Path

import pandas as pd

from sklearn.compose import ColumnTransformer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder


SEED = 42

PROJECT = Path(__file__).resolve().parents[1]
INPUT = PROJECT / "machine_learning/type_ic_context_model_dataset.tsv"
OUTPUT = PROJECT / "machine_learning/type_ic_all_64_pam_scores.tsv"

FEATURES = [
    "pam_minus3",
    "pam_minus2",
    "pam_minus1",
    "mobile_nearby",
]

data = pd.read_csv(INPUT, sep="\t")
data["label"] = data["label"].astype(int)
data["mobile_nearby"] = data["mobile_nearby"].astype(int)

preprocessing = ColumnTransformer(
    transformers=[
        (
            "pam",
            OneHotEncoder(
                categories=[["A", "C", "G", "T"]] * 3,
                drop="first",
                handle_unknown="ignore",
                sparse_output=False,
            ),
            ["pam_minus3", "pam_minus2", "pam_minus1"],
        ),
        (
            "mobile",
            "passthrough",
            ["mobile_nearby"],
        ),
    ],
    remainder="drop",
)

model = Pipeline(
    [
        ("preprocessing", preprocessing),
        (
            "classifier",
            LogisticRegression(
                class_weight="balanced",
                solver="liblinear",
                max_iter=5000,
                random_state=SEED,
            ),
        ),
    ]
)

model.fit(data[FEATURES], data["label"])

rows = []

for bases in product("ACGT", repeat=3):
    pam = "".join(bases)

    for mobile in (0, 1):
        rows.append(
            {
                "pam3": pam,
                "pam_minus3": bases[0],
                "pam_minus2": bases[1],
                "pam_minus1": bases[2],
                "mobile_nearby": mobile,
            }
        )

scores = pd.DataFrame(rows)
scores["predicted_probability"] = model.predict_proba(
    scores[FEATURES]
)[:, 1]

scores["predicted_percent"] = (
    100 * scores["predicted_probability"]
)

scores["rank_within_context"] = (
    scores.groupby("mobile_nearby")["predicted_probability"]
    .rank(method="min", ascending=False)
    .astype(int)
)

scores = scores.sort_values(
    ["mobile_nearby", "predicted_probability"],
    ascending=[True, False],
)

scores.to_csv(OUTPUT, sep="\t", index=False)

print("PAM-context combinations:", len(scores))
print("Unique PAMs:", scores["pam3"].nunique())
print()

for mobile in (0, 1):
    label = "mobile-associated" if mobile else "non-mobile"
    subset = scores[scores["mobile_nearby"] == mobile].head(10)

    print(f"TOP 10: {label}")
    print(
        subset[
            [
                "rank_within_context",
                "pam3",
                "predicted_probability",
                "predicted_percent",
            ]
        ].to_string(index=False)
    )
    print()

print("Output:", OUTPUT)
