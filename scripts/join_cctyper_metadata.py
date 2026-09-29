#!/usr/bin/env python3

from pathlib import Path
import csv

PROJECT = Path(__file__).resolve().parents[1]

SUMMARY_FILE = (
    PROJECT / "results" / "cctyper_1.8.0" / "genome_summary.tsv"
)
METADATA_FILE = (
    PROJECT / "taxonomy" /
    "paenibacillus_polymyxa_primary_quality_pass.tsv"
)
OUTPUT_FILE = (
    PROJECT / "results" / "cctyper_1.8.0" /
    "genome_summary_enriched.tsv"
)

METADATA_FIELDS = [
    ("organism", "Organism Name"),
    ("strain", "Organism Infraspecific Names Strain"),
    ("assembly_level", "Assembly Level"),
    ("assembly_name", "Assembly Name"),
    ("bioproject", "Assembly BioProject Accession"),
    ("biosample", "Assembly BioSample Accession"),
    ("isolation_source", "Assembly BioSample Isolation source"),
    ("geographic_location", "Assembly BioSample Geographic location"),
    ("release_date", "Assembly Release Date"),
    ("genome_length_bp", "Assembly Stats Total Sequence Length"),
    ("gc_percent", "Assembly Stats GC Percent"),
    ("number_of_contigs", "Assembly Stats Number of Contigs"),
    ("contig_n50", "Assembly Stats Contig N50"),
    ("checkm_completeness", "CheckM completeness"),
    ("checkm_contamination", "CheckM contamination"),
    ("ani_best_match", "ANI Best ANI match ANI"),
    ("ani_reference_assembly", "ANI Best ANI match Assembly"),
    ("ani_status", "ANI Best match status"),
    ("ani_check_status", "ANI Check status"),
]


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


summary_rows = read_tsv(SUMMARY_FILE)
metadata_rows = read_tsv(METADATA_FILE)

metadata_by_accession = {}

for row in metadata_rows:
    accession = row["Assembly Accession"].strip()
    if accession in metadata_by_accession:
        raise SystemExit(f"Duplicate metadata accession: {accession}")
    metadata_by_accession[accession] = row

missing = [
    row["accession"]
    for row in summary_rows
    if row["accession"] not in metadata_by_accession
]

if missing:
    raise SystemExit(
        "Missing metadata for: " + ", ".join(missing)
    )

enriched_rows = []

for summary in summary_rows:
    accession = summary["accession"]
    metadata = metadata_by_accession[accession]

    enriched = {"accession": accession}

    for output_name, metadata_name in METADATA_FIELDS:
        value = metadata.get(metadata_name, "").strip()
        enriched[output_name] = value if value else "NA"

    for field, value in summary.items():
        if field != "accession":
            enriched[field] = value if value else "NA"

    enriched_rows.append(enriched)

if len(enriched_rows) != 24:
    raise SystemExit(
        f"Expected 24 enriched rows, found {len(enriched_rows)}"
    )

with OUTPUT_FILE.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=list(enriched_rows[0]),
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(enriched_rows)

print(f"Metadata rows available: {len(metadata_rows)}")
print(f"Summary rows joined: {len(enriched_rows)}")
print("Missing metadata rows: 0")
print(f"Output: {OUTPUT_FILE}")
