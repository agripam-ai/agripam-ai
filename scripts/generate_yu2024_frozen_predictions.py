#!/usr/bin/env python3
"""Generate outcome-blind predictions with the published Yu et al. forest.

The released figure table supplies precomputed sequence/thermodynamic inputs.
Experimental columns in that table are never selected into the model matrix.
Guide IDs come from the already blinded GSE196911 prediction-input file.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import pickle
from pathlib import Path

import pandas as pd
import numpy as np

from agripam.sklearn_compat import load_legacy_sklearn_forest


BASE_FEATURES = [
    "guide_GC_content", "distance_start_codon", "distance_start_codon_perc", "homopolymers",
    "MFE_hybrid_full", "MFE_hybrid_seed", "MFE_homodimer_guide", "MFE_monomer_guide",
    "if_promoter",
]
FORBIDDEN_TOKENS = ("logfc", "pvalue", "edger", "outcome", "score")
FEATURE_GROUPS = {
    "sequence": lambda columns: [column for column in columns if column.startswith("sequence_")],
    "thermodynamics": lambda columns: [column for column in columns if column.startswith("MFE_")],
    "position": lambda columns: [column for column in columns if column.startswith("distance_start_codon")],
    "composition": lambda columns: [column for column in ("guide_GC_content", "homopolymers") if column in columns],
    "promoter": lambda columns: ["if_promoter"] if "if_promoter" in columns else [],
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_model_matrix(blinded: pd.DataFrame, released_features: pd.DataFrame,
                       headers: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    if blinded["sequence"].duplicated().any():
        raise ValueError("Blinded guide sequences must be unique.")
    safe = released_features[["sequence", "sequence_30nt", *BASE_FEATURES]].copy()
    if safe["sequence"].duplicated().any():
        raise ValueError("Released feature-table sequences must be unique.")
    joined = blinded[["guides", "sequence"]].merge(safe, on="sequence", how="left", validate="one_to_one")
    if joined["sequence_30nt"].isna().any():
        raise ValueError("Not every blinded guide has released predictor features.")
    encoded = {
        f"sequence_{position + 1}_{nt}": (joined["sequence_30nt"].str[position] == nt).astype(int)
        for position in range(30) for nt in "ATCG"
    }
    matrix = pd.concat([joined[BASE_FEATURES], pd.DataFrame(encoded, index=joined.index)], axis=1)
    matrix = matrix[headers].astype(float)
    forbidden = [column for column in matrix if any(token in column.lower() for token in FORBIDDEN_TOKENS)]
    if forbidden:
        raise ValueError(f"Outcome-like columns entered the model matrix: {forbidden}")
    return joined, matrix


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("blinded_inputs", type=Path)
    parser.add_argument("released_features", type=Path)
    parser.add_argument("model", type=Path)
    parser.add_argument("headers", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    blinded = pd.read_csv(args.blinded_inputs, sep="\t")
    released = pd.read_csv(args.released_features, sep="\t")
    with args.headers.open("rb") as handle:
        headers = pickle.load(handle)
    joined, matrix = build_model_matrix(blinded, released, headers)
    model = load_legacy_sklearn_forest(args.model)
    predicted_log2fc = model.predict(matrix)

    result = pd.DataFrame({
        "guide_id": joined["guides"],
        "prediction": -predicted_log2fc,  # higher = stronger predicted depletion
        "predicted_log2fc": predicted_log2fc,
        "prediction_status": "frozen_before_outcome_join",
    })
    rng = np.random.default_rng(42)
    sensitivity_groups = {}
    for group, selector in FEATURE_GROUPS.items():
        columns = selector(list(matrix.columns))
        if not columns:
            continue
        permuted = matrix.copy()
        order = rng.permutation(len(permuted))
        permuted.loc[:, columns] = matrix.loc[order, columns].to_numpy()
        result[f"prediction_sensitivity_permuted_{group}"] = -model.predict(permuted)
        sensitivity_groups[group] = columns
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, sep="\t", index=False)
    manifest = {
        "n_predictions": len(result),
        "prediction_orientation": "higher prediction = stronger depletion; equals -predicted_log2fc",
        "outcome_columns_in_model_matrix": [],
        "feature_count": len(headers),
        "sensitivity_analysis": {
            "method": "fixed-seed joint row permutation within each feature group; model not retrained",
            "interpretation": "model-reliance sensitivity, not a retrained feature ablation",
            "seed": 42, "groups": sensitivity_groups,
        },
        "model_conversion": "scikit-learn 0.24 tree nodes loaded in 1.3; missing_go_to_left added as zero; no retraining",
        "inputs": {name: {"path": str(path), "sha256": sha256(path)} for name, path in {
            "blinded_inputs": args.blinded_inputs, "released_features": args.released_features,
            "model": args.model, "headers": args.headers,
        }.items()},
        "output_sha256": sha256(args.output),
    }
    args.output.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
