#!/usr/bin/env python3

from pathlib import Path
import csv
import json

PROJECT = Path(__file__).resolve().parents[1]
RAW = PROJECT / "metadata" / "ncbi_raw"

TITLE_JSON = RAW / "paenibacillus_phage_title_accessions.json"
LARVAE_JSONL = RAW / "paenibacillus_larvae_host_complete.jsonl"
GLUCANO_JSONL = RAW / "paenibacillus_glucanolyticus_host_complete.jsonl"

MANIFEST = PROJECT / "metadata" / "paenibacillus_phage_query_manifest.tsv"
ACCESSIONS = PROJECT / "config" / "paenibacillus_phage_accessions.txt"


def read_jsonl(path):
    records = []

    with path.open() as handle:
        for line in handle:
            if line.strip():
                records.append(json.loads(line))

    return records


with TITLE_JSON.open() as handle:
    title_data = json.load(handle)

title_accessions = set(
    title_data["esearchresult"]["idlist"]
)

larvae_records = read_jsonl(LARVAE_JSONL)
glucano_records = read_jsonl(GLUCANO_JSONL)

larvae_by_accession = {
    record["accession"]: record
    for record in larvae_records
}

glucano_by_accession = {
    record["accession"]: record
    for record in glucano_records
}

larvae_accessions = set(larvae_by_accession)
glucano_accessions = set(glucano_by_accession)

all_accessions = sorted(
    title_accessions |
    larvae_accessions |
    glucano_accessions
)

rows = []

for accession in all_accessions:
    record = (
        larvae_by_accession.get(accession)
        or glucano_by_accession.get(accession)
        or {}
    )

    virus = record.get("virus") or {}
    host = record.get("host") or {}

    rows.append({
        "accession": accession,
        "in_title_query": int(accession in title_accessions),
        "in_larvae_host_query": int(accession in larvae_accessions),
        "in_glucanolyticus_host_query": int(
            accession in glucano_accessions
        ),
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
    })

fields = [
    "accession",
    "in_title_query",
    "in_larvae_host_query",
    "in_glucanolyticus_host_query",
    "virus_name",
    "virus_taxid",
    "host_name",
    "host_taxid",
    "completeness",
    "length",
    "source_database",
    "release_date",
]

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
    for accession in all_accessions:
        handle.write(accession + "\n")

print(f"Title-query accessions: {len(title_accessions)}")
print(f"P. larvae host accessions: {len(larvae_accessions)}")
print(
    "P. glucanolyticus host accessions:",
    len(glucano_accessions),
)
print(
    "Title/larvae overlap:",
    len(title_accessions & larvae_accessions),
)
print(
    "Title-only accessions:",
    len(title_accessions - larvae_accessions - glucano_accessions),
)
print(
    "Host-only accessions:",
    len(
        (larvae_accessions | glucano_accessions)
        - title_accessions
    ),
)
print(f"Union accessions: {len(all_accessions)}")
print(f"Manifest: {MANIFEST}")
print(f"Accession list: {ACCESSIONS}")
