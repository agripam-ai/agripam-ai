#!/usr/bin/env python3
"""Normalize the deposited GSE196911 purine screen without generating predictions."""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


EXPECTED = {"guides", "type", "OD02.logFC", "OD02.pvalue", "OD06.logFC", "OD06.pvalue", "OD1.logFC", "OD1.pvalue"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalize(source: Path) -> pd.DataFrame:
    raw = pd.read_csv(source, sep="\t")
    missing = EXPECTED - set(raw.columns)
    if missing:
        raise ValueError("Unexpected GSE196911 schema; missing: " + ", ".join(sorted(missing)))
    if raw["guides"].duplicated().any():
        raise ValueError("Deposited guide identifiers are not unique.")
    if set(raw["type"]) != {"targeting"}:
        raise ValueError("This release is expected to contain targeting guides only.")
    normalized = pd.DataFrame({
        "dataset_id": "YU2024_GSE196911_PURINE",
        "deposit": "GSE196911",
        "organism": "Escherichia coli K-12 MG1655",
        "genome_id": "NC_000913.3",
        "study_partition": "external_purine_screen",
        "guide_id": raw["guides"].astype(str),
        "gene_id": raw["guides"].astype(str).str.rsplit("_", n=1).str[0],
        "guide_type": raw["type"].astype(str),
        # More positive means stronger depletion/activity; the sign change is declared here.
        "outcome": -pd.to_numeric(raw["OD1.logFC"]),
        "primary_raw_log2fc_od1": pd.to_numeric(raw["OD1.logFC"]),
        "primary_pvalue_od1": pd.to_numeric(raw["OD1.pvalue"]),
        "secondary_raw_log2fc_od02": pd.to_numeric(raw["OD02.logFC"]),
        "secondary_pvalue_od02": pd.to_numeric(raw["OD02.pvalue"]),
        "secondary_raw_log2fc_od06": pd.to_numeric(raw["OD06.logFC"]),
        "secondary_pvalue_od06": pd.to_numeric(raw["OD06.pvalue"]),
        "prediction_status": "withheld — generate without outcome access",
    })
    if normalized["outcome"].isna().any():
        raise ValueError("Primary outcome contains missing or non-numeric values.")
    return normalized.sort_values(["gene_id", "guide_id"]).reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    table = normalize(args.source)
    table_path = args.output_dir / "heldout_outcomes.tsv"
    table.to_csv(table_path, sep="\t", index=False)
    manifest = {
        "dataset_id": "YU2024_GSE196911_PURINE",
        "deposit_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196911",
        "source_filename": args.source.name,
        "source_sha256": sha256(args.source),
        "normalized_sha256": sha256(table_path),
        "normalized_utc": datetime.now(timezone.utc).isoformat(),
        "n_guides": len(table),
        "n_genes": int(table["gene_id"].nunique()),
        "genes": sorted(table["gene_id"].unique()),
        "primary_outcome": "negative OD1 log2 fold-change (higher = stronger depletion)",
        "partition": "fully external; no rows permitted in model fitting or threshold selection",
        "prediction_file": None,
        "validation_status": "OUTCOMES NORMALIZED; PREDICTIONS NOT YET GENERATED",
        "interpretation": "Deposited CRISPRi depletion is an indirect activity proxy, not genome editing efficiency.",
    }
    (args.output_dir / "provenance_manifest.json").write_text(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
