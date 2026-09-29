#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict
import csv
import re

PROJECT = Path(__file__).resolve().parents[1]
CCTYPER_ROOT = PROJECT / "crispr_detection" / "cctyper_1.8.0"
GENOME_ROOT = PROJECT / "data" / "genomes_raw"
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

COMPLEMENT = str.maketrans(
    "ACGTRYKMSWBDHVNacgtrykmswbdhvn",
    "TGCAYRMKSWVHDBNtgcayrmkswvhdbn",
)


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


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


def read_fasta(path):
    sequences = {}
    name = None
    parts = []

    with path.open() as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(parts).upper()

                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)

    if name is not None:
        sequences[name] = "".join(parts).upper()

    return sequences


def parse_attributes(text):
    attributes = {}

    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attributes[key] = value

    return attributes


def native_minced_array_ids(path):
    native_ids = set()
    contig = None

    if not path.is_file():
        return native_ids

    with path.open(errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")

            sequence_match = SEQUENCE_RE.match(line)
            if sequence_match:
                contig = sequence_match.group(1)
                continue

            crispr_match = CRISPR_RE.match(line)
            if crispr_match and contig is not None:
                native_ids.add(
                    f"{contig}_{int(crispr_match.group(1))}"
                )

    return native_ids


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
qc_rows = []
total_repeat_features = 0

for genome_dir in sorted(CCTYPER_ROOT.glob("GCF_*")):
    if not genome_dir.is_dir():
        continue

    accession = genome_dir.name
    gff_file = genome_dir / "crisprs.gff"
    genome_file = GENOME_ROOT / f"{accession}.fna"

    if not gff_file.is_file():
        continue

    if not genome_file.is_file():
        raise SystemExit(f"Missing genome FASTA: {genome_file}")

    sequences = read_fasta(genome_file)
    native_ids = native_minced_array_ids(
        genome_dir / "minced.out"
    )

    regions = {}
    repeats_by_parent = defaultdict(list)

    with gff_file.open() as handle:
        for raw_line in handle:
            if raw_line.startswith("#") or not raw_line.strip():
                continue

            fields = raw_line.rstrip("\n").split("\t")

            if len(fields) != 9:
                raise SystemExit(
                    f"Malformed GFF line in {gff_file}: {raw_line}"
                )

            contig, _, feature_type, start, end, _, strand, _, attrs = fields
            attributes = parse_attributes(attrs)

            feature = {
                "contig": contig,
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "attributes": attributes,
            }

            if feature_type == "repeat_region":
                region_id = attributes.get("ID")
                if not region_id:
                    raise SystemExit(
                        f"repeat_region without ID in {gff_file}"
                    )
                regions[region_id] = feature

            elif feature_type == "direct_repeat":
                parent = attributes.get("Parent")
                if not parent:
                    raise SystemExit(
                        f"direct_repeat without Parent in {gff_file}"
                    )
                repeats_by_parent[parent].append(feature)
                total_repeat_features += 1

    for crispr_id, region in sorted(
        regions.items(),
        key=lambda item: (
            item[1]["contig"],
            item[1]["start"],
            item[0],
        ),
    ):
        key = (accession, crispr_id)

        if key not in array_by_key:
            raise SystemExit(
                f"GFF array absent from crisprs_all.tsv: "
                f"{accession} {crispr_id}"
            )

        annotation = array_by_key[key]
        repeats = repeats_by_parent.get(crispr_id, [])
        expected_repeats = int(annotation["N_repeats"])
        strand = region["strand"]

        if strand not in {"+", "-"}:
            raise SystemExit(
                f"Unsupported strand for {accession} {crispr_id}: "
                f"{strand}"
            )

        if len(repeats) != expected_repeats:
            raise SystemExit(
                f"Repeat-count mismatch for {accession} {crispr_id}: "
                f"expected {expected_repeats}, observed {len(repeats)}"
            )

        repeats.sort(
            key=lambda repeat: repeat["start"],
            reverse=(strand == "-"),
        )

        extraction_source = (
            "minced"
            if crispr_id in native_ids
            else "blast_recovered"
        )

        observed_spacers = 0
        spacer_lengths = []

        for index in range(len(repeats) - 1):
            current_repeat = repeats[index]
            next_repeat = repeats[index + 1]

            if strand == "+":
                spacer_start = current_repeat["end"] + 1
                spacer_end = next_repeat["start"] - 1
            else:
                spacer_start = next_repeat["end"] + 1
                spacer_end = current_repeat["start"] - 1

            if spacer_end < spacer_start:
                raise SystemExit(
                    f"Invalid spacer interval for "
                    f"{accession} {crispr_id} spacer {index + 1}"
                )

            contig = region["contig"]

            if contig not in sequences:
                raise SystemExit(
                    f"Contig {contig} missing from {genome_file}"
                )

            genomic_sequence = sequences[contig][
                spacer_start - 1:spacer_end
            ]

            if strand == "-":
                spacer_sequence = reverse_complement(
                    genomic_sequence
                )
            else:
                spacer_sequence = genomic_sequence

            if not spacer_sequence:
                raise SystemExit(
                    f"Empty spacer for "
                    f"{accession} {crispr_id} spacer {index + 1}"
                )

            if not set(spacer_sequence) <= set(
                "ACGTRYKMSWBDHVN"
            ):
                raise SystemExit(
                    f"Unexpected spacer symbols for "
                    f"{accession} {crispr_id} spacer {index + 1}"
                )

            spacer_index = index + 1
            spacer_id = (
                f"{accession}|{crispr_id}|S{spacer_index:03d}"
            )

            spacer_rows.append({
                "spacer_id": spacer_id,
                "accession": accession,
                "strain": strain_by_accession.get(accession, "NA"),
                "contig": contig,
                "crispr": crispr_id,
                "array_start": region["start"],
                "array_end": region["end"],
                "gff_strand": strand,
                "spacer_index": spacer_index,
                "left_repeat_start": min(
                    current_repeat["start"],
                    next_repeat["start"],
                ),
                "right_repeat_end": max(
                    current_repeat["end"],
                    next_repeat["end"],
                ),
                "spacer_start": spacer_start,
                "spacer_end": spacer_end,
                "spacer_length": len(spacer_sequence),
                "spacer_sequence": spacer_sequence,
                "extraction_source": extraction_source,
                "trusted": annotation.get("Trusted", "NA"),
                "repeat_prediction": annotation.get(
                    "Prediction", "NA"
                ),
                "repeat_subtype": annotation.get("Subtype", "NA"),
                "subtype_probability": annotation.get(
                    "Subtype_probability", "NA"
                ),
                "is_near_cas": int(key in near_keys),
                "is_orphan": int(key in orphan_keys),
            })

            observed_spacers += 1
            spacer_lengths.append(len(spacer_sequence))

        expected_spacers = expected_repeats - 1

        qc_rows.append({
            "accession": accession,
            "crispr": crispr_id,
            "contig": region["contig"],
            "gff_strand": strand,
            "extraction_source": extraction_source,
            "expected_repeats": expected_repeats,
            "observed_repeats": len(repeats),
            "expected_spacers": expected_spacers,
            "observed_spacers": observed_spacers,
            "minimum_spacer_length": min(spacer_lengths),
            "maximum_spacer_length": max(spacer_lengths),
            "status": (
                "PASS"
                if expected_spacers == observed_spacers
                else "FAIL"
            ),
        })

expected_total = sum(
    int(row["N_repeats"]) - 1
    for row in array_rows
)

if len(qc_rows) != 41:
    raise SystemExit(
        f"Expected 41 arrays, processed {len(qc_rows)}"
    )

if total_repeat_features != 419:
    raise SystemExit(
        f"Expected 419 repeat features, found "
        f"{total_repeat_features}"
    )

failed_arrays = [
    row for row in qc_rows
    if row["status"] != "PASS"
]

if failed_arrays:
    raise SystemExit(
        f"{len(failed_arrays)} arrays failed spacer validation"
    )

if len(spacer_rows) != expected_total:
    raise SystemExit(
        f"Expected {expected_total} spacers, "
        f"extracted {len(spacer_rows)}"
    )

if len({row["spacer_id"] for row in spacer_rows}) != len(
    spacer_rows
):
    raise SystemExit("Duplicate spacer identifiers detected")

spacer_fields = [
    "spacer_id",
    "accession",
    "strain",
    "contig",
    "crispr",
    "array_start",
    "array_end",
    "gff_strand",
    "spacer_index",
    "left_repeat_start",
    "right_repeat_end",
    "spacer_start",
    "spacer_end",
    "spacer_length",
    "spacer_sequence",
    "extraction_source",
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
        "contig",
        "gff_strand",
        "extraction_source",
        "expected_repeats",
        "observed_repeats",
        "expected_spacers",
        "observed_spacers",
        "minimum_spacer_length",
        "maximum_spacer_length",
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
            f'|source={row["extraction_source"]}'
            f'|location={row["contig"]}:'
            f'{row["spacer_start"]}-{row["spacer_end"]}'
        )
        handle.write(header + "\n")
        handle.write(row["spacer_sequence"] + "\n")

native_count = sum(
    row["extraction_source"] == "minced"
    for row in spacer_rows
)
recovered_count = sum(
    row["extraction_source"] == "blast_recovered"
    for row in spacer_rows
)

print(f"Arrays validated: {len(qc_rows)}")
print(f"Direct-repeat features validated: {total_repeat_features}")
print(f"Spacers extracted: {len(spacer_rows)}")
print(f"MinCED spacers: {native_count}")
print(f"BLAST-recovered spacers: {recovered_count}")
print(f"Arrays failing validation: {len(failed_arrays)}")
print(f"Spacer table: {SPACER_TABLE}")
print(f"Spacer FASTA: {SPACER_FASTA}")
print(f"QC table: {QC_FILE}")
