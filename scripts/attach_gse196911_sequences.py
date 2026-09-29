#!/usr/bin/env python3
"""Attach supplementary guide sequences to normalized GSE196911 outcomes.

The GEO file lacks sequences. Table S13 of the article supplement contains the
same three time-point measurements plus sequence context. Matching uses gene and
all three measurements, never model predictions.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


def _key(frame: pd.DataFrame, gene_column: str) -> pd.Series:
    return frame.apply(lambda row: "|".join([
        str(row[gene_column]),
        *[f"{float(row[column]):.12f}" for column in ("OD02.logFC", "OD06.logFC", "OD1.logFC")],
    ]), axis=1)


def attach(outcomes: pd.DataFrame, supplement: pd.DataFrame) -> pd.DataFrame:
    geo = outcomes.copy()
    geo["gene_name"] = geo["guides"].astype(str).str.rsplit("_", n=1).str[0]
    geo["match_key"] = _key(geo, "gene_name")
    sup = supplement.copy()
    sup["match_key"] = _key(sup, "gene_name")
    if geo["match_key"].duplicated().any():
        raise ValueError("GEO outcome keys are not unique.")
    usable = sup[sup["match_key"].isin(set(geo["match_key"]))]
    if usable["match_key"].duplicated().any():
        raise ValueError("Matched supplementary keys are not unique.")
    fields = ["match_key", "sequence", "sequence_30nt", "seq_60nt", "PAM", "genome_pos",
              "guide_strand", "guide_GC_content", "distance_start_codon", "distance_start_codon_perc"]
    joined = geo.merge(usable[fields], on="match_key", how="left", validate="one_to_one")
    if joined["sequence"].isna().any():
        raise ValueError(f"Only {joined['sequence'].notna().sum()} of {len(joined)} guides matched Table S13.")
    return joined.drop(columns=["gene_name", "match_key"])


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("geo", type=Path)
    parser.add_argument("supplement", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    outcomes = pd.read_csv(args.geo, sep="\t")
    supplement = pd.read_excel(args.supplement, sheet_name="TableS13_purine_screen_guides")
    joined = attach(outcomes, supplement)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    joined.to_csv(args.output, sep="\t", index=False)
    outcome_columns = [column for column in joined if column.endswith(".logFC") or column.endswith(".pvalue")]
    features_path = args.output.with_name("prediction_inputs_blinded.tsv")
    joined.drop(columns=outcome_columns).to_csv(features_path, sep="\t", index=False)
    manifest = {
        "geo_sha256": sha256(args.geo), "supplement_sha256": sha256(args.supplement),
        "output_sha256": sha256(args.output), "matched_guides": len(joined),
        "blinded_features_sha256": sha256(features_path),
        "unique_sequences": int(joined["sequence"].nunique()),
        "matching_fields": ["gene", "OD02.logFC", "OD06.logFC", "OD1.logFC"],
        "prediction_columns_present": [column for column in joined if column.startswith("prediction")],
        "blinded_file_outcome_columns": [column for column in pd.read_csv(features_path, sep="\t", nrows=1)
                                          if column.endswith(".logFC") or column.endswith(".pvalue")],
    }
    args.output.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
