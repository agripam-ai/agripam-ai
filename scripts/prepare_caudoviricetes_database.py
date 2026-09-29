#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict
import csv
import hashlib
import zipfile

PROJECT = Path(__file__).resolve().parents[1]

PACKAGE_ROOT = (
    PROJECT / "data" / "ncbi_packages" /
    "caudoviricetes_refseq_chunks"
)
ACCESSION_FILE = (
    PROJECT / "config" /
    "caudoviricetes_refseq_complete_accessions.txt"
)
METADATA_FILE = (
    PROJECT / "metadata" /
    "caudoviricetes_refseq_complete_manifest.tsv"
)

PROCESSED_ROOT = PROJECT / "data" / "phages" / "processed"
OUTPUT_FASTA = (
    PROCESSED_ROOT /
    "caudoviricetes_refseq_exact_deduplicated.fna"
)
SEQUENCE_MANIFEST = (
    PROJECT / "metadata" /
    "caudoviricetes_refseq_sequence_manifest.tsv"
)
GROUP_TABLE = (
    PROJECT / "metadata" /
    "caudoviricetes_refseq_sequence_groups.tsv"
)

MAIN_CHUNKS = [
    f"chunk_{index:02d}"
    for index in range(11)
    if index != 4
]

PACKAGE_FILES = [
    PACKAGE_ROOT / f"{chunk}.zip"
    for chunk in MAIN_CHUNKS
]

PACKAGE_FILES.extend(
    PACKAGE_ROOT / "chunk_04_parts" / f"part_{index}.zip"
    for index in range(5)
)


def parse_fasta(text):
    records = {}
    accession = None
    description = None
    parts = []

    for line in text.splitlines():
        if line.startswith(">"):
            if accession is not None:
                records[accession] = {
                    "description": description,
                    "sequence": "".join(parts).upper(),
                }

            description = line[1:]
            accession = description.split()[0]
            parts = []
        elif line.strip():
            parts.append(line.strip())

    if accession is not None:
        records[accession] = {
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

expected_accessions = {
    line.strip()
    for line in ACCESSION_FILE.open()
    if line.strip()
}

metadata_by_accession = {
    row["accession"]: row
    for row in read_tsv(METADATA_FILE)
}

all_records = {}

for package_file in PACKAGE_FILES:
    if not package_file.is_file():
        raise SystemExit(f"Missing package: {package_file}")

    with zipfile.ZipFile(package_file) as archive:
        bad_member = archive.testzip()

        if bad_member is not None:
            raise SystemExit(
                f"Corrupt ZIP member in {package_file}: "
                f"{bad_member}"
            )

        fasta_members = [
            name for name in archive.namelist()
            if name.endswith("/genomic.fna")
        ]

        if len(fasta_members) != 1:
            raise SystemExit(
                f"Expected one genomic.fna in {package_file}, "
                f"found {len(fasta_members)}"
            )

        text = archive.read(
            fasta_members[0]
        ).decode("utf-8")

    package_records = parse_fasta(text)

    for accession, record in package_records.items():
        if accession in all_records:
            raise SystemExit(
                f"Duplicate accession across packages: "
                f"{accession}"
            )

        record["source_package"] = package_file.name
        all_records[accession] = record

observed_accessions = set(all_records)

if observed_accessions != expected_accessions:
    missing = sorted(expected_accessions - observed_accessions)
    unexpected = sorted(observed_accessions - expected_accessions)

    raise SystemExit(
        f"Accession mismatch. Missing={len(missing)}, "
        f"unexpected={len(unexpected)}"
    )

groups_by_hash = defaultdict(list)

for accession, record in all_records.items():
    digest = hashlib.sha256(
        record["sequence"].encode()
    ).hexdigest()
    groups_by_hash[digest].append(accession)

group_rows = []
manifest_rows = []
representatives = []

for index, (digest, accessions) in enumerate(
    sorted(
        groups_by_hash.items(),
        key=lambda item: sorted(item[1])[0],
    ),
    start=1,
):
    group_id = f"CRV{index:06d}"
    accessions = sorted(accessions)
    representative = accessions[0]
    representative_sequence = all_records[
        representative
    ]["sequence"]

    representatives.append({
        "accession": representative,
        "group_id": group_id,
        "sequence": representative_sequence,
        "member_count": len(accessions),
    })

    group_rows.append({
        "sequence_group": group_id,
        "representative_accession": representative,
        "n_accessions": len(accessions),
        "all_accessions": ";".join(accessions),
        "sequence_length": len(representative_sequence),
        "sha256": digest,
    })

    for accession in accessions:
        metadata = metadata_by_accession[accession]
        record = all_records[accession]

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
            "sequence_length": len(record["sequence"]),
            "sha256": digest,
            "source_package": record["source_package"],
            "virus_name": metadata.get(
                "virus_name", "NA"
            ),
            "virus_taxid": metadata.get(
                "virus_taxid", "NA"
            ),
            "host_name": metadata.get(
                "host_name", "NA"
            ),
            "host_taxid": metadata.get(
                "host_taxid", "NA"
            ),
            "release_date": metadata.get(
                "release_date", "NA"
            ),
        })

with OUTPUT_FASTA.open("w") as handle:
    for record in representatives:
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
        "source_package",
        "virus_name",
        "virus_taxid",
        "host_name",
        "host_taxid",
        "release_date",
    ],
)

if len(PACKAGE_FILES) != 15:
    raise SystemExit(
        f"Expected 15 ZIP packages, found "
        f"{len(PACKAGE_FILES)}"
    )

if len(all_records) != 5497:
    raise SystemExit(
        f"Expected 5497 accessions, found "
        f"{len(all_records)}"
    )

if len(manifest_rows) != 5497:
    raise SystemExit(
        "Sequence manifest row-count mismatch"
    )

print(f"ZIP packages validated: {len(PACKAGE_FILES)}")
print(f"Input accessions validated: {len(all_records)}")
print(
    "Input total bases:",
    sum(len(record["sequence"]) for record in all_records.values()),
)
print(f"Exact sequence groups: {len(group_rows)}")
print(
    "Exact duplicate accession records removed:",
    len(all_records) - len(representatives),
)
print(f"Deduplicated FASTA: {OUTPUT_FASTA}")
print(f"Sequence manifest: {SEQUENCE_MANIFEST}")
print(f"Sequence groups: {GROUP_TABLE}")
