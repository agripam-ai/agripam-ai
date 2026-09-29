#!/usr/bin/env python3
"""Build the Hawkins 2020 E. coli -> B. subtilis transfer benchmark.

E. coli measurements are frozen as a sequence-matched transfer score before
the B. subtilis measurements are joined. This is an experimental cross-chassis
baseline, not a trained AgriPAM model and not evidence for B26.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
from pathlib import Path

FILES = ("eco_biorep1.csv", "eco_biorep2.csv", "bsu_biorep1.csv", "bsu_biorep2.csv")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_original_map(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        return {row["variant"]: row["original"] for row in rows if row["nmm"] in {"1", "2"}}


def parent_by_original(path: Path, mapping: dict[str, str]) -> dict[str, float]:
    values: dict[str, list[float]] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            original = mapping.get(row[""])
            if not original:
                continue
            try:
                activity = -float(row["parent"])
            except (TypeError, ValueError):
                continue
            values.setdefault(original, []).append(activity)
    return {guide: sum(group) / len(group) for guide, group in values.items()}


def ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    result = [0.0] * len(values)
    cursor = 0
    while cursor < len(order):
        end = cursor + 1
        while end < len(order) and values[order[end]] == values[order[cursor]]:
            end += 1
        average = (cursor + 1 + end) / 2
        for index in order[cursor:end]:
            result[index] = average
        cursor = end
    return result


def pearson(x: list[float], y: list[float]) -> float:
    mx, my = sum(x) / len(x), sum(y) / len(y)
    numerator = sum((a - mx) * (b - my) for a, b in zip(x, y))
    denominator = math.sqrt(sum((a - mx) ** 2 for a in x) * sum((b - my) ** 2 for b in y))
    return numerator / denominator if denominator else float("nan")


def spearman(x: list[float], y: list[float]) -> float:
    return pearson(ranks(x), ranks(y))


def interval(values: list[float]) -> dict[str, float]:
    ordered = sorted(value for value in values if math.isfinite(value))

    def quantile(q: float) -> float:
        position = (len(ordered) - 1) * q
        lower = int(position)
        upper = min(lower + 1, len(ordered) - 1)
        weight = position - lower
        return ordered[lower] * (1 - weight) + ordered[upper] * weight

    return {"ci95_low": quantile(0.025), "ci95_high": quantile(0.975)}


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build(raw: Path, output: Path, map_path: Path, resamples: int = 2000, seed: int = 42) -> dict:
    mapping = read_original_map(map_path)
    eco1 = parent_by_original(raw / "eco_biorep1.csv", mapping)
    eco2 = parent_by_original(raw / "eco_biorep2.csv", mapping)
    bsu1 = parent_by_original(raw / "bsu_biorep1.csv", mapping)
    bsu2 = parent_by_original(raw / "bsu_biorep2.csv", mapping)
    guides = sorted(set(eco1) & set(eco2) & set(bsu1) & set(bsu2))
    if len(guides) != 33:
        raise ValueError(f"Expected 33 fully matched guides; reconstructed {len(guides)}")

    prediction_rows, outcome_rows, joined = [], [], []
    for index, guide in enumerate(guides, 1):
        prediction = (eco1[guide] + eco2[guide]) / 2
        outcome = (bsu1[guide] + bsu2[guide]) / 2
        guide_id = f"HAWKINS_GFP_{index:02d}"
        prediction_rows.append({"dataset_id": "HAWKINS2020_GFP_TRANSFER", "guide_id": guide_id,
                                "guide_sequence": guide, "prediction": prediction,
                                "source_chassis": "Escherichia coli"})
        outcome_rows.append({"dataset_id": "HAWKINS2020_GFP_TRANSFER", "guide_id": guide_id,
                             "outcome": outcome, "heldout_chassis": "Bacillus subtilis"})
        joined.append({"dataset_id": "HAWKINS2020_GFP_TRANSFER", "organism": "Bacillus subtilis",
                       "genome_id": "B_subtilis_168_reporter", "guide_id": guide_id,
                       "guide_sequence": guide, "prediction": prediction, "outcome": outcome,
                       "ecoli_activity_rep1": eco1[guide], "ecoli_activity_rep2": eco2[guide],
                       "bsub_activity_rep1": bsu1[guide], "bsub_activity_rep2": bsu2[guide]})

    output.mkdir(parents=True, exist_ok=True)
    write_tsv(output / "frozen_ecoli_transfer_scores.tsv", prediction_rows, list(prediction_rows[0]))
    write_tsv(output / "heldout_bsubtilis_outcomes.tsv", outcome_rows, list(outcome_rows[0]))
    write_tsv(output / "frozen_scores_and_heldout_outcomes.tsv", joined, list(joined[0]))

    predictions = [row["prediction"] for row in joined]
    outcomes = [row["outcome"] for row in joined]
    observed = spearman(predictions, outcomes)
    rng = random.Random(seed)
    bootstrap, null = [], []
    for _ in range(resamples):
        indices = [rng.randrange(len(guides)) for _ in guides]
        bootstrap.append(spearman([predictions[i] for i in indices], [outcomes[i] for i in indices]))
        shuffled = outcomes[:]
        rng.shuffle(shuffled)
        null.append(spearman(predictions, shuffled))
    p_value = (1 + sum(abs(value) >= abs(observed) for value in null)) / (resamples + 1)

    top_k = []
    overall = sum(outcomes) / len(outcomes)
    for fraction in (0.10, 0.20, 0.25):
        k = max(1, math.ceil(len(guides) * fraction))
        selected = [outcomes[i] for i in sorted(range(len(guides)), key=lambda i: predictions[i], reverse=True)[:k]]
        top_k.append({"top_fraction": fraction, "k": k, "mean_outcome": sum(selected) / k,
                      "mean_lift_over_all": sum(selected) / k - overall})

    metrics = {
        "benchmark_role": "sequence-matched experimental organism-transfer baseline",
        "prediction_definition": "mean E. coli fully matched-guide repression activity across two biological replicates",
        "heldout_outcome_definition": "mean B. subtilis fully matched-guide repression activity across two biological replicates",
        "n_guides": len(guides), "spearman_rho": observed,
        "bootstrap": {"unit": "guide", "resamples": resamples, **interval(bootstrap)},
        "permutation": {"resamples": resamples, "empirical_two_sided_p": p_value,
                        "null_mean": sum(null) / len(null)},
        "replicate_reliability": {"ecoli_spearman": spearman([eco1[g] for g in guides], [eco2[g] for g in guides]),
                                  "bsubtilis_spearman": spearman([bsu1[g] for g in guides], [bsu2[g] for g in guides])},
        "continuous_top_k": top_k,
        "classification_note": "ROC/PR and calibration are not reported: no binary threshold or probability output was pre-specified.",
        "limitations": ["E. coli activity is a transfer score, not an AgriPAM-trained prediction.",
                        "The 33 guides target one GFP reporter, so this is not genome-wide validation.",
                        "The target is shared; this tests chassis transfer, not different genomic contexts.",
                        "Guide bootstrap does not remove dependence from overlapping GFP targets."],
        "raw_sha256": {name: sha256(raw / name) for name in FILES}, "map_sha256": sha256(map_path),
        "frozen_prediction_sha256": sha256(output / "frozen_ecoli_transfer_scores.tsv"),
        "heldout_outcome_sha256": sha256(output / "heldout_bsubtilis_outcomes.tsv"),
        "joined_sha256": sha256(output / "frozen_scores_and_heldout_outcomes.tsv"), "seed": seed,
    }
    (output / "validation_metrics.json").write_text(json.dumps(metrics, indent=2, allow_nan=False) + "\n")
    return metrics


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw", type=Path, default=root / "data/external_validation/raw/HAWKINS2020_GFP_CROSS_SPECIES")
    parser.add_argument("--map", type=Path, default=root / "data/external_validation/raw/HAWKINS2020_GFP_CROSS_SPECIES/gfp.origmap.tsv")
    parser.add_argument("--output", type=Path, default=root / "data/external_validation/processed/HAWKINS2020_GFP_TRANSFER")
    parser.add_argument("--resamples", type=int, default=2000)
    args = parser.parse_args()
    print(json.dumps(build(args.raw, args.output, args.map, args.resamples), indent=2))


if __name__ == "__main__":
    main()
