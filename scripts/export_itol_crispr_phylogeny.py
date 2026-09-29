#!/usr/bin/env python3

from pathlib import Path

import pandas as pd
from Bio import Phylo


PROJECT = Path(__file__).resolve().parents[1]

TREE_INPUT = (
    PROJECT
    / "taxonomy/core_genome/iqtree/"
      "paenibacillus_polymyxa_core_snps_midpoint_rooted.treefile"
)

METADATA_INPUT = (
    PROJECT
    / "taxonomy/core_genome/tree_tip_crispr_metadata.tsv"
)

OUTPUT_DIR = PROJECT / "figures/itol_crispr_phylogeny"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

TREE_OUTPUT = OUTPUT_DIR / "core_snp_midpoint_rooted.treefile"
LABELS_OUTPUT = OUTPUT_DIR / "itol_tip_labels.txt"
BINARY_OUTPUT = OUTPUT_DIR / "itol_crispr_subtypes_binary.txt"
ARRAYS_OUTPUT = OUTPUT_DIR / "itol_crispr_array_counts.txt"
STRIP_OUTPUT = OUTPUT_DIR / "itol_crispr_presence_strip.txt"

metadata = pd.read_csv(METADATA_INPUT, sep="\t")

tree = Phylo.read(TREE_INPUT, "newick")
tree_tips = {
    terminal.name
    for terminal in tree.get_terminals()
}

metadata_tips = set(metadata["accession"])

if tree_tips != metadata_tips:
    raise SystemExit(
        "Tree and metadata accessions do not match."
    )

Phylo.write(tree, TREE_OUTPUT, "newick")

metadata = metadata.sort_values("accession")


with LABELS_OUTPUT.open("w") as handle:
    handle.write("LABELS\n")
    handle.write("SEPARATOR TAB\n")
    handle.write("DATA\n")

    for row in metadata.itertuples(index=False):
        display = f"{row.strain} | {row.accession}"
        handle.write(
            f"{row.accession}\t{display}\n"
        )


with BINARY_OUTPUT.open("w") as handle:
    handle.write("DATASET_BINARY\n")
    handle.write("SEPARATOR TAB\n")
    handle.write("DATASET_LABEL\tCRISPR-Cas subtypes\n")
    handle.write("COLOR\t#333333\n")
    handle.write("FIELD_SHAPES\t1\t1\n")
    handle.write("FIELD_LABELS\tType I-C\tType III-B\n")
    handle.write("FIELD_COLORS\t#0072B2\t#D55E00\n")
    handle.write("SHOW_LABELS\t1\n")
    handle.write("DATA\n")

    for row in metadata.itertuples(index=False):
        handle.write(
            f"{row.accession}\t"
            f"{int(row.has_type_ic)}\t"
            f"{int(row.has_type_iiib)}\n"
        )


maximum_arrays = int(
    metadata["n_crispr_arrays"].max()
)

with ARRAYS_OUTPUT.open("w") as handle:
    handle.write("DATASET_SIMPLEBAR\n")
    handle.write("SEPARATOR TAB\n")
    handle.write("DATASET_LABEL\tCRISPR arrays\n")
    handle.write("COLOR\t#009E73\n")
    handle.write("WIDTH\t80\n")
    handle.write("MARGIN\t5\n")
    handle.write("HEIGHT_FACTOR\t0.8\n")
    handle.write(f"DATASET_SCALE\t0-{maximum_arrays}\n")
    handle.write("DATA\n")

    for row in metadata.itertuples(index=False):
        handle.write(
            f"{row.accession}\t"
            f"{int(row.n_crispr_arrays)}\n"
        )


with STRIP_OUTPUT.open("w") as handle:
    handle.write("DATASET_COLORSTRIP\n")
    handle.write("SEPARATOR TAB\n")
    handle.write("DATASET_LABEL\tCRISPR-Cas presence\n")
    handle.write("COLOR\t#000000\n")
    handle.write("STRIP_WIDTH\t30\n")
    handle.write("MARGIN\t5\n")
    handle.write("BORDER_WIDTH\t0\n")
    handle.write("SHOW_INTERNAL\t0\n")
    handle.write("LEGEND_TITLE\tCRISPR-Cas detected\n")
    handle.write("LEGEND_SHAPES\t1\t1\n")
    handle.write("LEGEND_COLORS\t#6A3D9A\t#D9D9D9\n")
    handle.write("LEGEND_LABELS\tPresent\tNot detected\n")
    handle.write("DATA\n")

    for row in metadata.itertuples(index=False):
        color = (
            "#6A3D9A"
            if int(row.has_cas) == 1
            else "#D9D9D9"
        )

        label = (
            "Present"
            if int(row.has_cas) == 1
            else "Not detected"
        )

        handle.write(
            f"{row.accession}\t{color}\t{label}\n"
        )


print("Tree tips:", len(tree_tips))
print("Maximum CRISPR-array count:", maximum_arrays)
print("Tree:", TREE_OUTPUT)
print("Tip labels:", LABELS_OUTPUT)
print("Subtype symbols:", BINARY_OUTPUT)
print("Array-count bars:", ARRAYS_OUTPUT)
print("Presence strip:", STRIP_OUTPUT)
