#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict
import csv

PROJECT = Path(__file__).resolve().parents[1]
SPACER_ROOT = PROJECT / "spacers"

INPUT_FILE = SPACER_ROOT / "spacers.tsv"
UNIQUE_TABLE = SPACER_ROOT / "unique_spacers.tsv"
MEMBERSHIP_TABLE = SPACER_ROOT / "spacer_membership.tsv"
UNIQUE_FASTA = SPACER_ROOT / "unique_spacers.fasta"


def joined_unique(rows, field):
    values = {
        row[field]
        for row in rows
        if row.get(field) and row[field] != "NA"
    }
    return ";".join(sorted(values)) if values else "NA"


with INPUT_FILE.open(newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

groups = defaultdict(list)

for row in rows:
    groups[row["spacer_sequence"]].append(row)

sequence_to_query = {
    sequence: f"ASP{index:06d}"
    for index, sequence in enumerate(sorted(groups), start=1)
}

unique_rows = []
membership_rows = []

for sequence in sorted(groups):
    query_id = sequence_to_query[sequence]
    members = groups[sequence]

    unique_rows.append({
        "query_id": query_id,
        "query_sequence": sequence,
        "spacer_length": len(sequence),
        "n_spacer_records": len(members),
        "n_accessions": len({
            row["accession"] for row in members
        }),
        "accessions": joined_unique(members, "accession"),
        "strains": joined_unique(members, "strain"),
        "crispr_arrays": joined_unique(members, "crispr"),
        "extraction_sources": joined_unique(
            members, "extraction_source"
        ),
        "repeat_subtypes": joined_unique(
            members, "repeat_subtype"
        ),
        "n_near_cas_records": sum(
            int(row["is_near_cas"]) for row in members
        ),
        "n_orphan_records": sum(
            int(row["is_orphan"]) for row in members
        ),
    })

    for member in members:
        membership_rows.append({
            "query_id": query_id,
            "spacer_id": member["spacer_id"],
            "accession": member["accession"],
            "strain": member["strain"],
            "contig": member["contig"],
            "crispr": member["crispr"],
            "spacer_index": member["spacer_index"],
            "extraction_source": member["extraction_source"],
            "repeat_subtype": member["repeat_subtype"],
            "is_near_cas": member["is_near_cas"],
            "is_orphan": member["is_orphan"],
        })

unique_fields = [
    "query_id",
    "query_sequence",
    "spacer_length",
    "n_spacer_records",
    "n_accessions",
    "accessions",
    "strains",
    "crispr_arrays",
    "extraction_sources",
    "repeat_subtypes",
    "n_near_cas_records",
    "n_orphan_records",
]

membership_fields = [
    "query_id",
    "spacer_id",
    "accession",
    "strain",
    "contig",
    "crispr",
    "spacer_index",
    "extraction_source",
    "repeat_subtype",
    "is_near_cas",
    "is_orphan",
]

with UNIQUE_TABLE.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=unique_fields,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(unique_rows)

with MEMBERSHIP_TABLE.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=membership_fields,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(membership_rows)

with UNIQUE_FASTA.open("w") as handle:
    for row in unique_rows:
        handle.write(
            f'>{row["query_id"]}'
            f' length={row["spacer_length"]}'
            f'|records={row["n_spacer_records"]}\n'
        )
        handle.write(row["query_sequence"] + "\n")

if len(rows) != 378:
    raise SystemExit(f"Expected 378 memberships, found {len(rows)}")

if len(unique_rows) != 249:
    raise SystemExit(
        f"Expected 249 exact unique spacers, found "
        f"{len(unique_rows)}"
    )

if len(membership_rows) != len(rows):
    raise SystemExit("Membership-row count mismatch")

print(f"Input spacer records: {len(rows)}")
print(f"Unique query sequences: {len(unique_rows)}")
print(f"Membership rows: {len(membership_rows)}")
print(f"Unique table: {UNIQUE_TABLE}")
print(f"Membership table: {MEMBERSHIP_TABLE}")
print(f"Query FASTA: {UNIQUE_FASTA}")
