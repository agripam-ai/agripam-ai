#!/usr/bin/env python3

from pathlib import Path
import csv
import re

PROJECT = Path(__file__).resolve().parents[1]
CCTYPER_ROOT = PROJECT / "crispr_detection" / "cctyper_1.8.0"
RESULT_ROOT = PROJECT / "results" / "cctyper_1.8.0"
SPACER_ROOT = PROJECT / "spacers"

ARRAY_FILE = RESULT_ROOT / "crisprs_all.tsv"
NEAR_FILE = RESULT_ROOT / "crisprs_near_cas.tsv"
ORPHAN_FILE = RESULT_ROOT / "crisprs_orphan.tsv"
SUMMARY_FILE = RESULT_ROOT / "genome_summary_enriched.tsv"

SPACER_TABLE = SPACER_ROOT / "spacers.tsv"
SPACER_FASTA = SPACER_ROOT / "spacers.fasta"
QC_FILE = SPACER_ROOT / "spacer_extraction_qc.tsv"

SEQUENCE_RE = re.compile(r"^Sequence '([^']+)'")
CRISPR_RE = re.compile(
    r"^CRISPR\s+(\d+)\s+Range:\s+(\d+)\s+-\s+(\d+)"
)
SPACER_RE = re.compile(
    r"^\s*(\d+)\s+([A-Za-z]+)\s+([A-Za-z]+)\s+\["
)


def read_tsv(path):
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return []
        return list(reader)


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


SPACER_ROOT.mkdir(parents=True, exist_ok=True)

array_rows = read_tsv(ARRAY_FILE)
near_rows = read_tsv(NEAR_FILE)
orphan_rows = read_tsv(ORPHAN_FILE)
summary_rows = read_tsv(SUMMARY_FILE)

strain_by_accession = {
    row["accession"]: row["strain"]
    for row in summary_rows
}

array_by_key = {
    (row["accession"], row["CRISPR"]): row
    for row in array_rows
}

near_keys = {
    (row["accession"], row["CRISPR"])
    for row in near_rows
}

orphan_keys = {
    (row["accession"], row["CRISPR"])
    for row in orphan_rows
}

spacer_rows = []
observed_by_array = {}

for genome_dir in sorted(CCTYPER_ROOT.glob("GCF_*")):
    if not genome_dir.is_dir():
        continue

    accession = genome_dir.name
    minced_file = genome_dir / "minced.out"

    if not minced_file.is_file():
        continue

    contig = None
    array_number = None
    array_start = None
    array_end = None
    spacer_index = 0

    with minced_file.open(errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")

            sequence_match = SEQUENCE_RE.match(line)
            if sequence_match:
                contig = sequence_match.group(1)
                array_number = None
                continue

            crispr_match = CRISPR_RE.match(line)
            if crispr_match:
                if contig is None:
                    raise SystemExit(
                        f"CRISPR encountered before sequence in {minced_file}"
                    )

                array_number = int(crispr_match.group(1))
                array_start = int(crispr_match.group(2))
                array_end = int(crispr_match.group(3))
                spacer_index = 0
                continue

            spacer_match = SPACER_RE.match(line)
            if spacer_match and contig is not None and array_number is not None:
                repeat_start = int(spacer_match.group(1))
                repeat_sequence = spacer_match.group(2).upper()
                spacer_sequence = spacer_match.group(3).upper()

                if not set(repeat_sequence) <= set("ACGTN"):
                    raise SystemExit(
                        f"Unexpected repeat symbols in {accession}: "
                        f"{repeat_sequence}"
                    )

                if not set(spacer_sequence) <= set("ACGTN"):
                    raise SystemExit(
                        f"Unexpected spacer symbols in {accession}: "
                        f"{spacer_sequence}"
                    )

                spacer_index += 1
                crispr_id = f"{contig}_{array_number}"
                key = (accession, crispr_id)

                if key not in array_by_key:
                    raise SystemExit(
                        f"MinCED array absent from crisprs_all.tsv: "
                        f"{accession} {crispr_id}"
                    )

                annotation = array_by_key[key]
                spacer_start = repeat_start + len(repeat_sequence)
                spacer_end = spacer_start + len(spacer_sequence) - 1
                spacer_id = (
                    f"{accession}|{crispr_id}|S{spacer_index:03d}"
                )

                spacer_rows.append({
                    "spacer_id": spacer_id,
                    "accession": accession,
                    "strain": strain_by_accession.get(accession, "NA"),
                    "contig": contig,
                    "crispr": crispr_id,
                    "array_number": array_number,
                    "array_start": array_start,
                    "array_end": array_end,
                    "spacer_index": spacer_index,
                    "repeat_start": repeat_start,
                    "repeat_sequence": repeat_sequence,
                    "spacer_start": spacer_start,
                    "spacer_end": spacer_end,
                    "spacer_length": len(spacer_sequence),
                    "spacer_sequence": spacer_sequence,
                    "trusted": annotation.get("Trusted", "NA"),
                    "repeat_prediction": annotation.get("Prediction", "NA"),
                    "repeat_subtype": annotation.get("Subtype", "NA"),
                    "subtype_probability": annotation.get(
                        "Subtype_probability", "NA"
                    ),
                    "is_near_cas": int(key in near_keys),
                    "is_orphan": int(key in orphan_keys),
                })

                observed_by_array[key] = (
                    observed_by_array.get(key, 0) + 1
                )

qc_rows = []

for key, annotation in sorted(array_by_key.items()):
    accession, crispr_id = key
    expected = int(annotation["N_repeats"]) - 1
    observed = observed_by_array.get(key, 0)

    qc_rows.append({
        "accession": accession,
        "crispr": crispr_id,
        "n_repeats": annotation["N_repeats"],
        "expected_spacers": expected,
        "observed_spacers": observed,
        "status": "PASS" if expected == observed else "FAIL",
    })

failed_arrays = [
    row for row in qc_rows
    if row["status"] != "PASS"
]

expected_total = sum(
    int(row["N_repeats"]) - 1
    for row in array_rows
)

if failed_arrays:
    details = ", ".join(
        f'{row["accession"]}:{row["crispr"]}'
        for row in failed_arrays
    )
    raise SystemExit(
        f"Spacer-count validation failed for: {details}"
    )

if len(spacer_rows) != expected_total:
    raise SystemExit(
        f"Expected {expected_total} total spacers, "
        f"extracted {len(spacer_rows)}"
    )

if len({row["spacer_id"] for row in spacer_rows}) != len(spacer_rows):
    raise SystemExit("Duplicate spacer identifiers detected")

spacer_fields = [
    "spacer_id",
    "accession",
    "strain",
    "contig",
    "crispr",
    "array_number",
    "array_start",
    "array_end",
    "spacer_index",
    "repeat_start",
    "repeat_sequence",
    "spacer_start",
    "spacer_end",
    "spacer_length",
    "spacer_sequence",
    "trusted",
    "repeat_prediction",
    "repeat_subtype",
    "subtype_probability",
    "is_near_cas",
    "is_orphan",
]

write_tsv(SPACER_TABLE, spacer_rows, spacer_fields)

write_tsv(
    QC_FILE,
    qc_rows,
    [
        "accession",
        "crispr",
        "n_repeats",
        "expected_spacers",
        "observed_spacers",
        "status",
    ],
)

with SPACER_FASTA.open("w") as handle:
    for row in spacer_rows:
        strain_safe = re.sub(
            r"[^A-Za-z0-9_.-]+", "_", row["strain"]
        )
        header = (
            f'>{row["spacer_id"]}'
            f'|strain={strain_safe}'
            f'|location={row["contig"]}:'
            f'{row["spacer_start"]}-{row["spacer_end"]}'
        )
        handle.write(header + "\n")
        handle.write(row["spacer_sequence"] + "\n")

print(f"Arrays validated: {len(qc_rows)}")
print(f"Arrays failing validation: {len(failed_arrays)}")
print(f"Spacers extracted: {len(spacer_rows)}")
print(f"Spacer table: {SPACER_TABLE}")
print(f"Spacer FASTA: {SPACER_FASTA}")
print(f"QC table: {QC_FILE}")
