#!/usr/bin/env python3

from pathlib import Path
import csv

PROJECT = Path(__file__).resolve().parents[1]
INPUT_ROOT = PROJECT / "crispr_detection" / "cctyper_1.8.0"
RESULT_ROOT = PROJECT / "results" / "cctyper_1.8.0"
STATUS_FILE = PROJECT / "results" / "cctyper_batch_status.tsv"

TABLES = {
    "cas_operons": "cas_operons.tab",
    "cas_operons_putative": "cas_operons_putative.tab",
    "crisprs_all": "crisprs_all.tab",
    "crisprs_near_cas": "crisprs_near_cas.tab",
    "crisprs_orphan": "crisprs_orphan.tab",
    "crispr_cas_links": "CRISPR_Cas.tab",
}


def read_tsv(path):
    if not path.is_file() or path.stat().st_size == 0:
        return []

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return []
        return [dict(row) for row in reader]


def unique_values(rows, field):
    values = {
        row.get(field, "").strip()
        for row in rows
        if row.get(field, "").strip()
    }
    return ";".join(sorted(values)) if values else "NA"


def write_rows(path, rows):
    first_fields = ["accession"]
    remaining_fields = []

    for row in rows:
        for field in row:
            if field not in first_fields and field not in remaining_fields:
                remaining_fields.append(field)

    fieldnames = first_fields + remaining_fields

    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


RESULT_ROOT.mkdir(parents=True, exist_ok=True)

status_rows = read_tsv(STATUS_FILE)
status_by_accession = {
    row["accession"]: row.get("status", "NA")
    for row in status_rows
}

combined = {name: [] for name in TABLES}
summary_rows = []

genome_dirs = sorted(
    path for path in INPUT_ROOT.glob("GCF_*")
    if path.is_dir()
)

for genome_dir in genome_dirs:
    accession = genome_dir.name
    genome_tables = {}

    for table_name, filename in TABLES.items():
        rows = read_tsv(genome_dir / filename)
        genome_tables[table_name] = rows

        for row in rows:
            combined[table_name].append({
                "accession": accession,
                **row,
            })

    cas_rows = genome_tables["cas_operons"]
    putative_rows = genome_tables["cas_operons_putative"]
    crispr_rows = genome_tables["crisprs_all"]
    near_rows = genome_tables["crisprs_near_cas"]
    orphan_rows = genome_tables["crisprs_orphan"]
    link_rows = genome_tables["crispr_cas_links"]

    summary_rows.append({
        "accession": accession,
        "batch_status": status_by_accession.get(accession, "NA"),
        "n_cas_operons": len(cas_rows),
        "n_putative_operons": len(putative_rows),
        "n_crispr_arrays": len(crispr_rows),
        "n_crisprs_near_cas": len(near_rows),
        "n_orphan_crisprs": len(orphan_rows),
        "n_crispr_cas_links": len(link_rows),
        "has_cas_operon": int(bool(cas_rows)),
        "has_crispr_array": int(bool(crispr_rows)),
        "has_linked_system": int(bool(link_rows)),
        "cas_predictions": unique_values(cas_rows, "Prediction"),
        "crispr_repeat_subtypes": unique_values(crispr_rows, "Subtype"),
    })

if len(summary_rows) != 24:
    raise SystemExit(
        f"Expected 24 genome directories, found {len(summary_rows)}"
    )

if len({row["accession"] for row in summary_rows}) != 24:
    raise SystemExit("Duplicate genome accessions detected")

for table_name, rows in combined.items():
    write_rows(RESULT_ROOT / f"{table_name}.tsv", rows)

write_rows(RESULT_ROOT / "genome_summary.tsv", summary_rows)

print(f"Genome summaries: {len(summary_rows)}")
for table_name, rows in combined.items():
    print(f"{table_name}: {len(rows)} rows")
print(f"Results directory: {RESULT_ROOT}")
