#!/usr/bin/env python3

import re
from pathlib import Path

import pandas as pd
from Bio import Phylo


PROJECT = Path(__file__).resolve().parents[1]

TREE_INPUT = (
    PROJECT
    / "taxonomy/core_genome/iqtree/"
      "paenibacillus_polymyxa_core_snps.treefile"
)

METADATA_INPUT = (
    PROJECT
    / "results/cctyper_1.8.0/"
      "genome_summary_enriched.tsv"
)

ROOTED_OUTPUT = (
    PROJECT
    / "taxonomy/core_genome/iqtree/"
      "paenibacillus_polymyxa_core_snps_midpoint_rooted.treefile"
)

LABELLED_OUTPUT = (
    PROJECT
    / "taxonomy/core_genome/iqtree/"
      "paenibacillus_polymyxa_core_snps_midpoint_rooted_labelled.treefile"
)

ANNOTATION_OUTPUT = (
    PROJECT
    / "taxonomy/core_genome/"
      "tree_tip_crispr_metadata.tsv"
)

SUPPORT_OUTPUT = (
    PROJECT
    / "taxonomy/core_genome/"
      "tree_support_summary.tsv"
)


def safe_label(value):
    value = str(value).strip()
    value = re.sub(r"\s+", "_", value)
    value = re.sub(r"[^A-Za-z0-9_.|+-]", "_", value)
    return value


metadata = pd.read_csv(METADATA_INPUT, sep="\t")

required = {
    "accession",
    "strain",
    "n_cas_operons",
    "n_crispr_arrays",
    "n_crisprs_near_cas",
    "n_orphan_crisprs",
    "cas_predictions",
    "crispr_repeat_subtypes",
}

missing = required - set(metadata.columns)
if missing:
    raise SystemExit(
        "Missing metadata columns: "
        + ", ".join(sorted(missing))
    )

for column in [
    "n_cas_operons",
    "n_crispr_arrays",
    "n_crisprs_near_cas",
    "n_orphan_crisprs",
]:
    metadata[column] = pd.to_numeric(
        metadata[column],
        errors="raise",
    ).astype(int)

metadata["strain"] = metadata["strain"].fillna("unknown")
metadata["tip_label"] = metadata.apply(
    lambda row: (
        f"{safe_label(row['strain'])}|{row['accession']}"
    ),
    axis=1,
)

metadata["has_crispr"] = (
    metadata["n_crispr_arrays"] > 0
).astype(int)

metadata["has_cas"] = (
    metadata["n_cas_operons"] > 0
).astype(int)

metadata["has_type_ic"] = (
    metadata["cas_predictions"]
    .fillna("")
    .str.split(";")
    .map(lambda values: int("I-C" in values))
)

metadata["has_type_iiib"] = (
    metadata["cas_predictions"]
    .fillna("")
    .str.split(";")
    .map(lambda values: int("III-B" in values))
)

annotation_columns = [
    "accession",
    "strain",
    "tip_label",
    "n_cas_operons",
    "n_crispr_arrays",
    "n_crisprs_near_cas",
    "n_orphan_crisprs",
    "has_crispr",
    "has_cas",
    "has_type_ic",
    "has_type_iiib",
    "cas_predictions",
    "crispr_repeat_subtypes",
]

annotations = metadata[annotation_columns].sort_values(
    "accession"
)

tree = Phylo.read(TREE_INPUT, "newick")
tree_accessions = {
    terminal.name
    for terminal in tree.get_terminals()
}

metadata_accessions = set(annotations["accession"])

if tree_accessions != metadata_accessions:
    raise SystemExit(
        "Tree/metadata mismatch. "
        f"Missing metadata={sorted(tree_accessions - metadata_accessions)}; "
        f"missing tree tips={sorted(metadata_accessions - tree_accessions)}"
    )

tree.root_at_midpoint()
Phylo.write(tree, ROOTED_OUTPUT, "newick")

label_lookup = dict(
    zip(
        annotations["accession"],
        annotations["tip_label"],
    )
)

for terminal in tree.get_terminals():
    terminal.name = label_lookup[terminal.name]

Phylo.write(tree, LABELLED_OUTPUT, "newick")

annotations.to_csv(
    ANNOTATION_OUTPUT,
    sep="\t",
    index=False,
)

tree_text = TREE_INPUT.read_text()

support_pairs = [
    (float(sh), float(ufboot))
    for sh, ufboot in re.findall(
        r"\)([0-9.]+)/([0-9.]+):",
        tree_text,
    )
]

if not support_pairs:
    raise SystemExit(
        "No SH-aLRT/UFBoot support pairs found."
    )

support = pd.DataFrame(
    support_pairs,
    columns=["sh_alrt", "ultrafast_bootstrap"],
)

summary_rows = [
    {
        "metric": "supported_internal_branches",
        "value": len(support),
    },
    {
        "metric": "sh_alrt_ge_80",
        "value": int((support["sh_alrt"] >= 80).sum()),
    },
    {
        "metric": "ufboot_ge_95",
        "value": int(
            (support["ultrafast_bootstrap"] >= 95).sum()
        ),
    },
    {
        "metric": "both_sh80_ufboot95",
        "value": int(
            (
                (support["sh_alrt"] >= 80)
                & (support["ultrafast_bootstrap"] >= 95)
            ).sum()
        ),
    },
    {
        "metric": "ufboot_below_70",
        "value": int(
            (support["ultrafast_bootstrap"] < 70).sum()
        ),
    },
]

pd.DataFrame(summary_rows).to_csv(
    SUPPORT_OUTPUT,
    sep="\t",
    index=False,
)

print("Tree tips:", len(tree_accessions))
print("CRISPR-positive genomes:", int(annotations["has_crispr"].sum()))
print("Cas-positive genomes:", int(annotations["has_cas"].sum()))
print("Type I-C-positive genomes:", int(annotations["has_type_ic"].sum()))
print("Type III-B-positive genomes:", int(annotations["has_type_iiib"].sum()))
print("Internal branches with dual support:", len(support))
print(
    "Branches with SH-aLRT >=80 and UFBoot >=95:",
    int(
        (
            (support["sh_alrt"] >= 80)
            & (support["ultrafast_bootstrap"] >= 95)
        ).sum()
    ),
)
print(
    "Branches with UFBoot <70:",
    int((support["ultrafast_bootstrap"] < 70).sum()),
)
print("Rooted tree:", ROOTED_OUTPUT)
print("Labelled tree:", LABELLED_OUTPUT)
print("Tip annotations:", ANNOTATION_OUTPUT)
print("Support summary:", SUPPORT_OUTPUT)
