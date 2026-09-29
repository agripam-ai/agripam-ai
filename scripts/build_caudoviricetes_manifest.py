#!/usr/bin/env python3

from pathlib import Path
import csv
import json

PROJECT = Path(__file__).resolve().parents[1]

INPUT = (
    PROJECT / "metadata" / "ncbi_raw" /
    "caudoviricetes_refseq_complete.jsonl"
)
MANIFEST = (
    PROJECT / "metadata" /
    "caudoviricetes_refseq_complete_manifest.tsv"
)
ACCESSIONS = (
    PROJECT / "config" /
    "caudoviricetes_refseq_complete_accessions.txt"
)

records = []

with INPUT.open() as handle:
    for line in handle:
        if line.strip():
            records.append(json.loads(line))

rows = []

for record in records:
    virus = record.get("virus") or {}
    host = record.get("host") or {}

    row = {
        "accession": record["accession"],
        "virus_name": virus.get("organism_name", "NA"),
        "virus_taxid": virus.get("tax_id", "NA"),
        "host_name": host.get("organism_name", "NA"),
        "host_taxid": host.get("tax_id", "NA"),
        "completeness": record.get("completeness", "NA"),
        "length": record.get("length", "NA"),
        "source_database": record.get(
            "source_database", "NA"
        ),
        "release_date": record.get("release_date", "NA"),
        "query_taxon": "Caudoviricetes",
        "query_complete_only": "true",
        "query_refseq_only": "true",
        "access_date": "2026-09-10",
    }
    rows.append(row)

rows.sort(key=lambda row: row["accession"])

if len(rows) != 5497:
    raise SystemExit(
        f"Expected 5497 records, found {len(rows)}"
    )

if len({row["accession"] for row in rows}) != len(rows):
    raise SystemExit("Duplicate accessions detected")

fields = list(rows[0])

with MANIFEST.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)

with ACCESSIONS.open("w") as handle:
    for row in rows:
        handle.write(row["accession"] + "\n")

print(f"Manifest rows: {len(rows)}")
print(f"Unique accessions: {len(rows)}")
print(f"Manifest: {MANIFEST}")
print(f"Accession list: {ACCESSIONS}")
