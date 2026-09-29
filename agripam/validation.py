"""Leakage-resistant retrospective validation for deposited CRISPR datasets.

The evaluator is deliberately model-agnostic: predictions must already be frozen in
``prediction`` before the held-out ``outcome`` column is supplied.  This prevents the
dashboard from silently training on the outcomes it claims to validate against.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd


REQUIRED_COLUMNS = {"dataset_id", "organism", "genome_id", "guide_id", "prediction", "outcome"}


@dataclass(frozen=True)
class ValidationConfig:
    positive_threshold: float | None = None
    top_fractions: tuple[float, ...] = (0.01, 0.05, 0.10)
    calibration_bins: int = 10
    seed: int = 42
    resamples: int = 2000


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_table(frame: pd.DataFrame) -> pd.DataFrame:
    missing = REQUIRED_COLUMNS - set(frame.columns)
    if missing:
        raise ValueError("Validation table is missing: " + ", ".join(sorted(missing)))
    clean = frame.copy()
    clean["prediction"] = pd.to_numeric(clean["prediction"], errors="coerce")
    clean["outcome"] = pd.to_numeric(clean["outcome"], errors="coerce")
    clean = clean.dropna(subset=["prediction", "outcome"])
    if clean.empty:
        raise ValueError("Validation table contains no finite prediction/outcome pairs.")
    if clean["guide_id"].duplicated().any():
        raise ValueError("guide_id must be unique within a validation file.")
    return clean


def _rank(values: pd.Series) -> np.ndarray:
    return values.rank(method="average").to_numpy(dtype=float)


def _spearman_values(predictions, outcomes) -> float:
    predictions = pd.Series(predictions)
    outcomes = pd.Series(outcomes)
    if len(predictions) < 2 or predictions.nunique() < 2 or outcomes.nunique() < 2:
        return float("nan")
    return float(np.corrcoef(_rank(predictions), _rank(outcomes))[0, 1])


def _within_group_median(frame: pd.DataFrame, group_column: str) -> float:
    values = []
    for _, group in frame.groupby(group_column):
        rho = _spearman_values(group["prediction"], group["outcome"])
        if np.isfinite(rho):
            values.append(rho)
    return float(np.median(values)) if values else float("nan")


def _roc_auc(labels: np.ndarray, scores: np.ndarray) -> float:
    positives, negatives = labels.sum(), len(labels) - labels.sum()
    if positives == 0 or negatives == 0:
        return float("nan")
    ranks = pd.Series(scores).rank(method="average").to_numpy()
    return float((ranks[labels == 1].sum() - positives * (positives + 1) / 2) / (positives * negatives))


def _average_precision(labels: np.ndarray, scores: np.ndarray) -> float:
    positives = int(labels.sum())
    if positives == 0:
        return float("nan")
    order = np.argsort(-scores, kind="stable")
    ordered = labels[order]
    precision = np.cumsum(ordered) / np.arange(1, len(ordered) + 1)
    return float((precision * ordered).sum() / positives)


def evaluate_predictions(frame: pd.DataFrame, config: ValidationConfig = ValidationConfig()) -> dict:
    data = validate_table(frame)
    predictions = data["prediction"].to_numpy(dtype=float)
    outcomes = data["outcome"].to_numpy(dtype=float)
    spearman = _spearman_values(data["prediction"], data["outcome"])
    result: dict = {"n_guides": len(data), "n_genomes": int(data["genome_id"].nunique()), "spearman_rho": spearman}
    continuous_enrichment = []
    overall_mean = float(np.mean(outcomes))
    for fraction in config.top_fractions:
        k = max(1, int(np.ceil(len(data) * fraction)))
        selected = outcomes[np.argsort(-predictions, kind="stable")[:k]]
        continuous_enrichment.append({
            "top_fraction": fraction, "k": k, "mean_outcome": float(np.mean(selected)),
            "median_outcome": float(np.median(selected)),
            "mean_lift_over_all": float(np.mean(selected) - overall_mean),
        })
    result["continuous_top_k"] = continuous_enrichment
    if "gene_id" in data.columns:
        grouped = []
        for gene_id, group in data.groupby("gene_id"):
            if len(group) < 3 or group["prediction"].nunique() < 2 or group["outcome"].nunique() < 2:
                continue
            rho = float(np.corrcoef(_rank(group["prediction"]), _rank(group["outcome"]))[0, 1])
            grouped.append({"gene_id": str(gene_id), "n_guides": len(group), "spearman_rho": rho})
        if grouped:
            result["per_gene_spearman"] = grouped
            result["median_within_gene_spearman"] = float(np.median([row["spearman_rho"] for row in grouped]))

    if config.positive_threshold is not None:
        labels = (outcomes >= config.positive_threshold).astype(int)
        prevalence = float(labels.mean())
        result.update({"positive_threshold": config.positive_threshold, "positive_prevalence": prevalence,
                       "roc_auc": _roc_auc(labels, predictions), "average_precision": _average_precision(labels, predictions)})
        enrichment = []
        for fraction in config.top_fractions:
            k = max(1, int(np.ceil(len(data) * fraction)))
            chosen = labels[np.argsort(-predictions, kind="stable")[:k]].mean()
            enrichment.append({"top_fraction": fraction, "k": k, "positive_rate": float(chosen),
                               "fold_enrichment": float(chosen / prevalence) if prevalence else float("nan")})
        result["top_k_enrichment"] = enrichment

        # Calibration is meaningful only when prediction values are probabilities.
        if np.all((predictions >= 0) & (predictions <= 1)):
            bins = pd.cut(predictions, np.linspace(0, 1, config.calibration_bins + 1), include_lowest=True)
            calibration = (pd.DataFrame({"bin": bins, "prediction": predictions, "label": labels})
                           .groupby("bin", observed=True).agg(n=("label", "size"), mean_prediction=("prediction", "mean"), observed_rate=("label", "mean"))
                           .reset_index().astype({"bin": str}).to_dict("records"))
            result["brier_score"] = float(np.mean((predictions - labels) ** 2))
            result["calibration"] = calibration
        else:
            result["calibration_note"] = "Skipped: predictions are ranks/scores, not calibrated probabilities."
    else:
        result["classification_note"] = "ROC/PR, enrichment and calibration require a pre-specified biological threshold."
    return result


def permuted_outcome_control(frame: pd.DataFrame, config: ValidationConfig = ValidationConfig()) -> dict:
    clean = validate_table(frame)
    rng = np.random.default_rng(config.seed)
    shuffled = clean.copy()
    shuffled["outcome"] = rng.permutation(shuffled["outcome"].to_numpy())
    return evaluate_predictions(shuffled, config)


def cluster_bootstrap_intervals(frame: pd.DataFrame, cluster_column: str = "gene_id",
                                config: ValidationConfig = ValidationConfig()) -> dict:
    """Bootstrap whole biological groups rather than treating guides as independent."""
    data = validate_table(frame)
    if cluster_column not in data or data[cluster_column].nunique() < 2:
        return {"status": "not_available", "reason": f"At least two {cluster_column} groups are required."}
    groups = {str(key): group for key, group in data.groupby(cluster_column)}
    keys = np.array(list(groups), dtype=object)
    rng = np.random.default_rng(config.seed)
    pooled, within = [], []
    for _ in range(config.resamples):
        sampled_keys = rng.choice(keys, size=len(keys), replace=True)
        sample = pd.concat([groups[str(key)] for key in sampled_keys], ignore_index=True)
        pooled.append(_spearman_values(sample["prediction"], sample["outcome"]))
        # Resampling group-level correlations directly keeps the gene as the unit.
        per_group = [_spearman_values(groups[str(key)]["prediction"], groups[str(key)]["outcome"])
                     for key in sampled_keys]
        within.append(float(np.nanmedian(per_group)))
    def interval(values):
        finite = np.asarray(values, dtype=float)
        finite = finite[np.isfinite(finite)]
        return {"median": float(np.median(finite)), "ci95_low": float(np.quantile(finite, 0.025)),
                "ci95_high": float(np.quantile(finite, 0.975)), "valid_resamples": int(len(finite))}
    return {"status": "complete", "cluster_unit": cluster_column, "resamples": config.resamples,
            "pooled_spearman": interval(pooled), "median_within_cluster_spearman": interval(within)}


def permutation_null_test(frame: pd.DataFrame, group_column: str = "gene_id",
                          config: ValidationConfig = ValidationConfig()) -> dict:
    """Permute outcomes within genes, preserving gene-level outcome distributions."""
    data = validate_table(frame)
    if group_column not in data:
        return {"status": "not_available", "reason": f"Missing {group_column}."}
    observed_pooled = _spearman_values(data["prediction"], data["outcome"])
    observed_within = _within_group_median(data, group_column)
    rng = np.random.default_rng(config.seed)
    pooled_null, within_null = [], []
    group_indices = [group.index.to_numpy() for _, group in data.groupby(group_column)]
    for _ in range(config.resamples):
        shuffled = data.copy()
        for indices in group_indices:
            shuffled.loc[indices, "outcome"] = rng.permutation(data.loc[indices, "outcome"].to_numpy())
        pooled_null.append(_spearman_values(shuffled["prediction"], shuffled["outcome"]))
        within_null.append(_within_group_median(shuffled, group_column))
    def summarize(observed, null):
        values = np.asarray(null, dtype=float)
        values = values[np.isfinite(values)]
        p = (1 + int(np.sum(np.abs(values) >= abs(observed)))) / (len(values) + 1)
        return {"observed": observed, "null_mean": float(np.mean(values)),
                "null_sd": float(np.std(values, ddof=1)), "empirical_two_sided_p": float(p)}
    return {"status": "complete", "permutation_unit": f"within {group_column}",
            "resamples": config.resamples, "pooled_spearman": summarize(observed_pooled, pooled_null),
            "median_within_cluster_spearman": summarize(observed_within, within_null)}


def leave_one_group_out_influence(frame: pd.DataFrame, group_column: str = "gene_id") -> dict:
    """Quantify whether one biological group drives the reported association."""
    data = validate_table(frame)
    if group_column not in data or data[group_column].nunique() < 3:
        return {"status": "not_available", "reason": f"At least three {group_column} groups are required."}
    full = _spearman_values(data["prediction"], data["outcome"])
    rows = []
    for key, _ in data.groupby(group_column):
        retained = data[data[group_column] != key]
        rho = _spearman_values(retained["prediction"], retained["outcome"])
        rows.append({"omitted_group": str(key), "n_retained": int(len(retained)),
                     "pooled_spearman": rho, "delta_from_full": float(rho - full)})
    values = np.asarray([row["pooled_spearman"] for row in rows], dtype=float)
    return {"status": "complete", "group_unit": group_column, "full_pooled_spearman": full,
            "minimum_leave_one_out_spearman": float(np.nanmin(values)),
            "maximum_leave_one_out_spearman": float(np.nanmax(values)), "rows": rows}


def random_effects_spearman_meta(frame: pd.DataFrame, group_column: str = "gene_id") -> dict:
    """Exploratory heterogeneity-aware synthesis of within-group correlations."""
    data = validate_table(frame)
    if group_column not in data:
        return {"status": "not_available", "reason": f"Missing {group_column}."}
    rows = []
    for key, group in data.groupby(group_column):
        n = len(group)
        rho = _spearman_values(group["prediction"], group["outcome"])
        if n < 4 or not np.isfinite(rho):
            continue
        clipped = float(np.clip(rho, -0.999999, 0.999999))
        rows.append({"group": str(key), "n": int(n), "spearman_rho": float(rho),
                     "z": float(np.arctanh(clipped)), "variance": float(1 / (n - 3))})
    if len(rows) < 2:
        return {"status": "not_available", "reason": "At least two estimable groups are required."}
    z = np.asarray([row["z"] for row in rows]); variances = np.asarray([row["variance"] for row in rows])
    fixed_w = 1 / variances; fixed_z = float(np.sum(fixed_w * z) / np.sum(fixed_w))
    q = float(np.sum(fixed_w * (z - fixed_z) ** 2)); df = len(rows) - 1
    c = float(np.sum(fixed_w) - np.sum(fixed_w ** 2) / np.sum(fixed_w))
    tau2 = max(0.0, (q - df) / c) if c > 0 else 0.0
    random_w = 1 / (variances + tau2)
    pooled_z = float(np.sum(random_w * z) / np.sum(random_w)); se = float(np.sqrt(1 / np.sum(random_w)))
    return {"status": "complete", "method": "DerSimonian-Laird random effects on Fisher-z transformed Spearman correlations",
            "groups": len(rows), "pooled_spearman": float(np.tanh(pooled_z)),
            "ci95_low": float(np.tanh(pooled_z - 1.96 * se)),
            "ci95_high": float(np.tanh(pooled_z + 1.96 * se)), "tau_squared": float(tau2),
            "cochran_q": q, "q_df": df,
            "i_squared_percent": float(max(0.0, (q - df) / q) * 100) if q > 0 else 0.0,
            "note": "Exploratory robustness analysis; Fisher-z variances are approximate for Spearman correlations.",
            "per_group": [{k: v for k, v in row.items() if k not in {"z", "variance"}} for row in rows]}


def top_k_permutation_test(frame: pd.DataFrame, group_column: str = "gene_id",
                           config: ValidationConfig = ValidationConfig()) -> dict:
    """Test top-k mean lift against permutations preserving group distributions."""
    data = validate_table(frame)
    predictions = data["prediction"].to_numpy(dtype=float)
    outcomes = data["outcome"].to_numpy(dtype=float)
    order = np.argsort(-predictions, kind="stable")
    overall = float(np.mean(outcomes))
    observed = []
    for fraction in config.top_fractions:
        k = max(1, int(np.ceil(len(data) * fraction)))
        observed.append((fraction, k, float(np.mean(outcomes[order[:k]]) - overall)))
    rng = np.random.default_rng(config.seed)
    group_indices = ([group.index.to_numpy() for _, group in data.groupby(group_column)]
                     if group_column in data else [data.index.to_numpy()])
    null = [[] for _ in observed]
    base = data["outcome"].copy()
    for _ in range(config.resamples):
        shuffled = base.copy()
        for indices in group_indices:
            shuffled.loc[indices] = rng.permutation(base.loc[indices].to_numpy())
        shuffled_values = shuffled.to_numpy(dtype=float)
        shuffled_overall = float(np.mean(shuffled_values))
        for index, (_, k, _) in enumerate(observed):
            null[index].append(float(np.mean(shuffled_values[order[:k]]) - shuffled_overall))
    rows = []
    for (fraction, k, lift), values in zip(observed, null):
        values = np.asarray(values, dtype=float)
        p = (1 + int(np.sum(values >= lift))) / (len(values) + 1)
        rows.append({"top_fraction": fraction, "k": k, "observed_mean_lift": lift,
                     "null_mean_lift": float(np.mean(values)), "empirical_one_sided_p": float(p)})
    return {"status": "complete", "permutation_unit": f"within {group_column}" if group_column in data else "all rows",
            "resamples": config.resamples, "rows": rows}


def grouped_holdout_manifest(frame: pd.DataFrame, held_out_genomes: list[str]) -> dict:
    clean = validate_table(frame)
    held = set(map(str, held_out_genomes))
    test = clean[clean["genome_id"].astype(str).isin(held)]
    train = clean[~clean["genome_id"].astype(str).isin(held)]
    if test.empty or train.empty:
        raise ValueError("Held-out genomes must produce non-empty, disjoint train and test partitions.")
    return {"group_key": "genome_id", "train_genomes": sorted(train["genome_id"].astype(str).unique()),
            "test_genomes": sorted(test["genome_id"].astype(str).unique()), "train_n": len(train), "test_n": len(test),
            "outcomes_hidden_until_predictions_frozen": True}


def evaluate_ablations(frame: pd.DataFrame, prediction_columns: list[str],
                       config: ValidationConfig = ValidationConfig()) -> list[dict]:
    """Evaluate precomputed ablation predictions on one unchanged held-out set."""
    rows = []
    for column in prediction_columns:
        if column not in frame:
            raise ValueError(f"Missing ablation prediction column: {column}")
        candidate = frame.copy()
        candidate["prediction"] = candidate[column]
        metrics = evaluate_predictions(candidate, config)
        rows.append({"model": column, **{key: value for key, value in metrics.items()
                                        if not isinstance(value, (list, dict))}})
    if rows:
        baseline_rho = rows[0].get("spearman_rho", float("nan"))
        baseline_within = rows[0].get("median_within_gene_spearman", float("nan"))
        for row in rows:
            row["delta_spearman_vs_full"] = float(row.get("spearman_rho", float("nan")) - baseline_rho)
            row["delta_median_within_gene_vs_full"] = float(
                row.get("median_within_gene_spearman", float("nan")) - baseline_within
            )
    return rows


def write_validation_bundle(frame: pd.DataFrame, output: Path, config: ValidationConfig = ValidationConfig()) -> dict:
    output.mkdir(parents=True, exist_ok=True)
    clean = validate_table(frame)
    clean.to_csv(output / "frozen_predictions_and_outcomes.tsv", sep="\t", index=False)
    sensitivity_columns = [column for column in clean if column.startswith("prediction_sensitivity_")]
    payload = {"analysis": evaluate_predictions(clean, config), "negative_control": permuted_outcome_control(clean, config),
               "cluster_bootstrap": cluster_bootstrap_intervals(clean, "gene_id", config),
               "permutation_null": permutation_null_test(clean, "gene_id", config),
               "leave_one_group_out": leave_one_group_out_influence(clean, "gene_id"),
               "random_effects_meta_analysis": random_effects_spearman_meta(clean, "gene_id"),
               "top_k_permutation": top_k_permutation_test(clean, "gene_id", config),
               "feature_group_sensitivity": evaluate_ablations(clean, ["prediction", *sensitivity_columns], config)
               if sensitivity_columns else [],
               "config": asdict(config), "input_sha256": sha256_file(output / "frozen_predictions_and_outcomes.tsv")}
    (output / "validation_metrics.json").write_text(json.dumps(payload, indent=2, allow_nan=True))
    return payload
