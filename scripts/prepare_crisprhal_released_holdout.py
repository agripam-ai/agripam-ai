#!/usr/bin/env python3
"""Outcome-separated benchmark on the crisprHAL TevSpCas9 released holdout.

The released training table is the only source used to fit models. The testing
table is first separated into sequence-only inputs and outcomes; predictions
are written and checksummed before the outcome table is joined.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression, RidgeCV
from sklearn.metrics import mean_absolute_error, mean_squared_error

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from agripam.validation import ValidationConfig, evaluate_predictions, top_k_permutation_test  # noqa: E402

DNA = "ACGT"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_release(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, header=None, names=["context_37nt", "activity"])
    frame["context_37nt"] = frame["context_37nt"].astype(str).str.upper()
    frame["activity"] = pd.to_numeric(frame["activity"], errors="raise")
    bad = frame[~frame["context_37nt"].str.fullmatch("[ACGT]{37}")]
    if not bad.empty:
        raise ValueError(f"{path} contains {len(bad)} invalid 37-nt contexts")
    return frame


def sequence_features(sequences: pd.Series) -> tuple[np.ndarray, list[str]]:
    names = [f"pos_{i+1}_{base}" for i in range(37) for base in DNA]
    names += ["gc_fraction", "guide_gc_fraction"]
    names += [f"dinucleotide_{a}{b}" for a in DNA for b in DNA]
    rows = []
    for sequence in sequences:
        one_hot = [float(sequence[i] == base) for i in range(37) for base in DNA]
        gc = (sequence.count("G") + sequence.count("C")) / 37
        guide = sequence[3:23]
        guide_gc = (guide.count("G") + guide.count("C")) / 20
        dinucs = [sum(sequence[i:i+2] == a+b for i in range(36)) / 36 for a in DNA for b in DNA]
        rows.append(one_hot + [gc, guide_gc] + dinucs)
    return np.asarray(rows, dtype=float), names


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    return float(pd.Series(x).corr(pd.Series(y), method="spearman"))


def bootstrap_and_permutation(prediction: np.ndarray, outcome: np.ndarray, seed: int, resamples: int) -> dict:
    rng = np.random.default_rng(seed); n = len(outcome)
    boot = []
    for _ in range(resamples):
        idx = rng.integers(0, n, n); boot.append(spearman(prediction[idx], outcome[idx]))
    observed = spearman(prediction, outcome)
    null = [spearman(prediction, rng.permutation(outcome)) for _ in range(resamples)]
    return {
        "unit": "guide", "resamples": resamples, "observed_spearman": observed,
        "bootstrap_ci95_low": float(np.nanquantile(boot, 0.025)),
        "bootstrap_ci95_high": float(np.nanquantile(boot, 0.975)),
        "permutation_null_mean": float(np.nanmean(null)),
        "permutation_two_sided_p": float((1 + np.sum(np.abs(null) >= abs(observed))) / (resamples + 1)),
    }


def nearest_sequence_identity(x_train: np.ndarray, x_test: np.ndarray, batch_size: int = 128) -> tuple[np.ndarray, np.ndarray]:
    """Return maximum positional identity to training for 37-nt context and 20-nt guide."""
    train_context = x_train[:, :148].astype(np.uint8)
    test_context = x_test[:, :148].astype(np.uint8)
    train_guide = train_context[:, 12:92]
    test_guide = test_context[:, 12:92]
    best_context, best_guide = [], []
    for start in range(0, len(test_context), batch_size):
        stop = min(len(test_context), start + batch_size)
        best_context.extend(np.max(test_context[start:stop] @ train_context.T, axis=1).tolist())
        best_guide.extend(np.max(test_guide[start:stop] @ train_guide.T, axis=1).tolist())
    return np.asarray(best_context, dtype=int), np.asarray(best_guide, dtype=int)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", type=Path, default=ROOT / "data/external_validation/raw/CRISPRHAL_TEV_RELEASED_HOLDOUT")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "data/external_validation/processed/CRISPRHAL_TEV_RELEASED_HOLDOUT")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--resamples", type=int, default=2000)
    args = parser.parse_args(); args.output_dir.mkdir(parents=True, exist_ok=True)
    train_path = args.raw_dir / "TevSpCas9_training_data.csv"
    test_path = args.raw_dir / "TevSpCas9_testing_data.csv"
    train, test = read_release(train_path), read_release(test_path)
    overlap = set(train.context_37nt) & set(test.context_37nt)
    blinded = pd.DataFrame({"guide_id": [f"CRISPRHAL_TEV_{i:05d}" for i in range(len(test))],
                            "context_37nt": test.context_37nt})
    outcomes = pd.DataFrame({"guide_id": blinded.guide_id, "outcome": test.activity})
    blinded_path = args.output_dir / "prediction_inputs_blinded.tsv"
    outcome_path = args.output_dir / "heldout_outcomes.tsv"
    blinded.to_csv(blinded_path, sep="\t", index=False); outcomes.to_csv(outcome_path, sep="\t", index=False)

    x_train, features = sequence_features(train.context_37nt)
    x_test, _ = sequence_features(blinded.context_37nt)
    nearest_context_identity, nearest_guide_identity = nearest_sequence_identity(x_train, x_test)
    ridge = RidgeCV(alphas=np.logspace(-3, 3, 25), cv=5).fit(x_train, train.activity)
    gc_model = LinearRegression().fit(x_train[:, -18:-16], train.activity)
    predictions = pd.DataFrame({"guide_id": blinded.guide_id,
                                "prediction": ridge.predict(x_test),
                                "prediction_gc_baseline": gc_model.predict(x_test[:, -18:-16])})
    prediction_path = args.output_dir / "frozen_predictions.tsv"
    predictions.to_csv(prediction_path, sep="\t", index=False)
    frozen_at = datetime.now(timezone.utc).isoformat()

    joined = blinded.merge(predictions, on="guide_id", validate="one_to_one").merge(outcomes, on="guide_id", validate="one_to_one")
    joined.insert(0, "dataset_id", "CRISPRHAL_TEV_RELEASED_HOLDOUT")
    joined.insert(1, "organism", "Citrobacter rodentium")
    joined.insert(2, "genome_id", "C_rodentium_TevSpCas9_release")
    joined["nearest_training_context_identity_37"] = nearest_context_identity
    joined["nearest_training_guide_identity_20"] = nearest_guide_identity
    joined_path = args.output_dir / "frozen_predictions_and_outcomes.tsv"
    joined.to_csv(joined_path, sep="\t", index=False)
    config = ValidationConfig(top_fractions=(0.01, 0.05, 0.10, 0.20), seed=args.seed, resamples=args.resamples)
    full = evaluate_predictions(joined, config)
    baseline_frame = joined.copy(); baseline_frame["prediction"] = baseline_frame["prediction_gc_baseline"]
    baseline = evaluate_predictions(baseline_frame, config)
    y = joined.outcome.to_numpy(); pred = joined.prediction.to_numpy()
    distant = joined[joined["nearest_training_guide_identity_20"] <= 17].copy()
    distant_metrics = evaluate_predictions(distant, config) if len(distant) >= 2 else {"status": "not_available"}
    metrics = {
        "benchmark_role": "predefined released holdout; same study and chassis, not external-study or held-out-genome validation",
        "separation": "models fitted on released training table only; testing outcomes stored separately until predictions were frozen",
        "n_training": len(train), "n_testing": len(test), "exact_sequence_overlap_train_test": len(overlap),
        "ridge_alpha_selected_by_training_only_5fold_cv": float(ridge.alpha_),
        "primary": full, "gc_baseline": baseline,
        "sequence_similarity_audit": {
            "method": "maximum position-wise identity to any released training sequence",
            "context_identity_37_percentiles": {str(q): float(np.quantile(nearest_context_identity, q)) for q in (0, .25, .5, .75, .95, 1)},
            "guide_identity_20_percentiles": {str(q): float(np.quantile(nearest_guide_identity, q)) for q in (0, .25, .5, .75, .95, 1)},
            "test_guides_with_context_identity_at_least_34_of_37": int(np.sum(nearest_context_identity >= 34)),
            "test_guides_with_guide_identity_at_least_18_of_20": int(np.sum(nearest_guide_identity >= 18)),
            "distant_subset_rule": "nearest training guide identity <=17/20",
            "distant_subset_n": int(len(distant)),
            "distant_subset_metrics": distant_metrics,
        },
        "delta_spearman_over_gc_baseline": float(full["spearman_rho"] - baseline["spearman_rho"]),
        "rmse": float(mean_squared_error(y, pred) ** 0.5), "mae": float(mean_absolute_error(y, pred)),
        "uncertainty_and_null": bootstrap_and_permutation(pred, y, args.seed, args.resamples),
        "top_k_permutation": top_k_permutation_test(joined, "gene_id", config),
        "limitations": [
            "The release split is guide-level and comes from one study and chassis.",
            "Guide bootstrap treats guides as independent because gene or target-cluster identifiers are absent.",
            "This transparent ridge model is an AgriPAM benchmark model, not the authors' neural-network reproduction.",
            "No activity threshold was pre-specified, so ROC, PR and calibration are not reported.",
        ],
        "source": {"repository": "https://github.com/tbrowne5/crisprHAL",
                   "paper_doi": "10.1038/s41467-023-41143-7",
                   "training_sha256": sha(train_path), "testing_sha256": sha(test_path)},
        "artifacts": {"blinded_input_sha256": sha(blinded_path), "heldout_outcome_sha256": sha(outcome_path),
                      "frozen_prediction_sha256": sha(prediction_path), "joined_sha256": sha(joined_path),
                      "frozen_at_utc": frozen_at},
        "feature_count": len(features), "seed": args.seed,
    }
    (args.output_dir / "validation_metrics.json").write_text(json.dumps(metrics, indent=2, allow_nan=True))
    print(json.dumps({"n_training": len(train), "n_testing": len(test), "overlap": len(overlap),
                      "rho": full["spearman_rho"], "gc_baseline_rho": baseline["spearman_rho"],
                      "alpha": ridge.alpha_}, indent=2))


if __name__ == "__main__":
    main()
