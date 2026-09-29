#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict
import csv
import hashlib
import json

PROJECT = Path(__file__).resolve().parents[1]

PHAGE_ROOT = PROJECT / "data" / "phages"
PACKAGE_DATA = (
    PHAGE_ROOT / "ncbi_paenibacillus_union_129" /
    "ncbi_dataset" / "data"
)
INPUT_FASTA = PACKAGE_DATA / "genomic.fna"
METADATA_JSONL = PACKAGE_DATA / "data_report.jsonl"
QUERY_MANIFEST = (
    PROJECT / "metadata" /
    "paenibacillus_phage_query_manifest.tsv"
)

PROCESSED_ROOT = PHAGE_ROOT / "processed"
OUTPUT_FASTA = (
    PROCESSED_ROOT /
    "paenibacillus_phages_exact_deduplicated.fna"
)
SEQUENCE_MANIFEST = (
    PROJECT / "metadata" /
    "paenibacillus_phage_sequence_manifest.tsv"
)
GROUP_TABLE = (
    PROJECT / "metadata" /
    "paenibacillus_phage_sequence_groups.tsv"
)


def read_fasta(path):
    records = {}
    name = None
    description = None
    parts = []

    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")

            if line.startswith(">"):
                if name is not None:
                    records[name] = {
                        "description": description,
                        "sequence": "".join(parts).upper(),
                    }

                description = line[1:]
                name = description.split()[0]
                parts = []
            else:
                parts.append(line.strip())

    if name is not None:
        records[name] = {
            "description": description,
            "sequence": "".join(parts).upper(),
        }

    return records


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fieldnames):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


PROCESSED_ROOT.mkdir(parents=True, exist_ok=True)

fasta_records = read_fasta(INPUT_FASTA)

metadata_by_accession = {}

with METADATA_JSONL.open() as handle:
    for line in handle:
        if not line.strip():
            continue

        record = json.loads(line)

        record["source_database"] = record.get(
            "source_database",
            record.get("sourceDatabase", "NA"),
        )
        record["release_date"] = record.get(
            "release_date",
            record.get("releaseDate", "NA"),
        )

        virus = record.get("virus") or {}
        virus["organism_name"] = virus.get(
            "organism_name",
            virus.get("organismName", "NA"),
        )
        virus["tax_id"] = virus.get(
            "tax_id",
            virus.get("taxId", "NA"),
        )
        record["virus"] = virus

        host = record.get("host") or {}
        host["organism_name"] = host.get(
            "organism_name",
            host.get("organismName", "NA"),
        )
        host["tax_id"] = host.get(
            "tax_id",
            host.get("taxId", "NA"),
        )
        record["host"] = host

        metadata_by_accession[record["accession"]] = record

query_by_accession = {
    row["accession"]: row
    for row in read_tsv(QUERY_MANIFEST)
}

if set(fasta_records) != set(metadata_by_accession):
    raise SystemExit(
        "FASTA and metadata accession sets do not match"
    )

groups_by_hash = defaultdict(list)

for accession, fasta_record in fasta_records.items():
    digest = hashlib.sha256(
        fasta_record["sequence"].encode()
    ).hexdigest()
    groups_by_hash[digest].append(accession)

group_rows = []
manifest_rows = []
representative_records = []

sorted_groups = sorted(
    groups_by_hash.items(),
    key=lambda item: sorted(item[1])[0],
)

for index, (digest, accessions) in enumerate(
    sorted_groups,
    start=1,
):
    group_id = f"PHG{index:05d}"

    accessions = sorted(
        accessions,
        key=lambda accession: (
            0 if metadata_by_accession[accession].get(
                "source_database"
            ) == "RefSeq" else 1,
            accession,
        ),
    )

    representative = accessions[0]
    representative_metadata = metadata_by_accession[
        representative
    ]
    representative_sequence = fasta_records[
        representative
    ]["sequence"]

    representative_records.append({
        "accession": representative,
        "group_id": group_id,
        "sequence": representative_sequence,
        "member_count": len(accessions),
    })

    group_rows.append({
        "sequence_group": group_id,
        "representative_accession": representative,
        "n_accessions": len(accessions),
        "all_accessions": ";".join(sorted(accessions)),
        "sequence_length": len(representative_sequence),
        "sha256": digest,
        "representative_source_database":
            representative_metadata.get(
                "source_database", "NA"
            ),
    })

    for accession in accessions:
        metadata = metadata_by_accession[accession]
        query = query_by_accession.get(accession, {})
        virus = metadata.get("virus") or {}
        host = metadata.get("host") or {}

        manifest_rows.append({
            "accession": accession,
            "sequence_group": group_id,
            "representative_accession": representative,
            "is_representative": int(
                accession == representative
            ),
            "is_exact_duplicate": int(
                len(accessions) > 1
            ),
            "n_accessions_in_group": len(accessions),
            "sequence_length": len(
                fasta_records[accession]["sequence"]
            ),
            "sha256": digest,
            "source_database": metadata.get(
                "source_database", "NA"
            ),
            "virus_name": virus.get(
                "organism_name", "NA"
            ),
            "virus_taxid": virus.get("tax_id", "NA"),
            "host_name": host.get(
                "organism_name", "NA"
            ),
            "host_taxid": host.get("tax_id", "NA"),
            "completeness": metadata.get(
                "completeness", "NA"
            ),
            "release_date": metadata.get(
                "release_date", "NA"
            ),
            "in_title_query": query.get(
                "in_title_query", "0"
            ),
            "in_larvae_host_query": query.get(
                "in_larvae_host_query", "0"
            ),
            "in_glucanolyticus_host_query": query.get(
                "in_glucanolyticus_host_query", "0"
            ),
        })

with OUTPUT_FASTA.open("w") as handle:
    for record in representative_records:
        handle.write(
            f'>{record["accession"]}'
            f' sequence_group={record["group_id"]}'
            f' accessions={record["member_count"]}\n'
        )

        sequence = record["sequence"]
        for start in range(0, len(sequence), 80):
            handle.write(sequence[start:start + 80] + "\n")

write_tsv(
    GROUP_TABLE,
    group_rows,
    [
        "sequence_group",
        "representative_accession",
        "n_accessions",
        "all_accessions",
        "sequence_length",
        "sha256",
        "representative_source_database",
    ],
)

write_tsv(
    SEQUENCE_MANIFEST,
    sorted(
        manifest_rows,
        key=lambda row: (
            row["sequence_group"],
            -int(row["is_representative"]),
            row["accession"],
        ),
    ),
    [
        "accession",
        "sequence_group",
        "representative_accession",
        "is_representative",
        "is_exact_duplicate",
        "n_accessions_in_group",
        "sequence_length",
        "sha256",
        "source_database",
        "virus_name",
        "virus_taxid",
        "host_name",
        "host_taxid",
        "completeness",
        "release_date",
        "in_title_query",
        "in_larvae_host_query",
        "in_glucanolyticus_host_query",
    ],
)

if len(fasta_records) != 129:
    raise SystemExit(
        f"Expected 129 input genomes, found "
        f"{len(fasta_records)}"
    )

if len(representative_records) != 98:
    raise SystemExit(
        f"Expected 98 unique sequences, found "
        f"{len(representative_records)}"
    )

if len(manifest_rows) != 129:
    raise SystemExit(
        f"Expected 129 manifest rows, found "
        f"{len(manifest_rows)}"
    )

print(f"Input accessions: {len(fasta_records)}")
print(f"Exact sequence groups: {len(group_rows)}")
print(
    "Duplicate accession records removed from FASTA:",
    len(fasta_records) - len(representative_records),
)
print(
    "RefSeq representatives:",
    sum(
        row["representative_source_database"] == "RefSeq"
        for row in group_rows
    ),
)
print(f"Deduplicated FASTA: {OUTPUT_FASTA}")
print(f"Sequence manifest: {SEQUENCE_MANIFEST}")
print(f"Sequence groups: {GROUP_TABLE}")
