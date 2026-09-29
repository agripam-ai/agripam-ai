#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]

SPACER_INPUT = PROJECT / "results/spacer_sharing_by_spacer.tsv"
PAIR_INPUT = PROJECT / "results/spacer_sharing_positive_genome_pairs.tsv"
METADATA_INPUT = (
    PROJECT
    / "taxonomy/core_genome/tree_tip_crispr_metadata.tsv"
)

SUMMARY_OUTPUT = (
    PROJECT / "figures/prism_spacer_private_shared.tsv"
)
ANI_OUTPUT = (
    PROJECT / "figures/prism_spacer_ani_scatter.tsv"
)
TREE_OUTPUT = (
    PROJECT / "figures/prism_spacer_tree_distance_scatter.tsv"
)


spacers = pd.read_csv(SPACER_INPUT, sep="\t")
pairs = pd.read_csv(PAIR_INPUT, sep="\t")
metadata = pd.read_csv(METADATA_INPUT, sep="\t")

strain_lookup = dict(
    zip(metadata["accession"], metadata["strain"])
)

private_count = int(
    (spacers["genomes_present"] == 1).sum()
)

shared_count = int(
    (spacers["genomes_present"] > 1).sum()
)

total = private_count + shared_count

summary = pd.DataFrame(
    [
        {
            "category": "Genome-private",
            "spacer_sequences": private_count,
            "percent": 100 * private_count / total,
            "total_unique_spacers": total,
        },
        {
            "category": "Shared_between_genomes",
            "spacer_sequences": shared_count,
            "percent": 100 * shared_count / total,
            "total_unique_spacers": total,
        },
    ]
)

summary.to_csv(
    SUMMARY_OUTPUT,
    sep="\t",
    index=False,
)

pairs["strain_a"] = pairs["accession_a"].map(
    strain_lookup
)

pairs["strain_b"] = pairs["accession_b"].map(
    strain_lookup
)

pairs["pair_label"] = (
    pairs["strain_a"].astype(str)
    + " vs "
    + pairs["strain_b"].astype(str)
)

ani_table = pairs[
    [
        "pair_label",
        "symmetric_ani",
        "jaccard_similarity",
        "shared_spacers",
        "union_spacers",
    ]
].copy()

ani_table = ani_table.rename(
    columns={
        "symmetric_ani": "X_ANI_percent",
        "jaccard_similarity": "Y_spacer_Jaccard",
    }
)

ani_table = ani_table.sort_values(
    ["X_ANI_percent", "Y_spacer_Jaccard"]
)

ani_table.to_csv(
    ANI_OUTPUT,
    sep="\t",
    index=False,
)

tree_table = pairs[
    [
        "pair_label",
        "tree_patristic_distance",
        "jaccard_similarity",
        "shared_spacers",
        "union_spacers",
    ]
].copy()

tree_table = tree_table.rename(
    columns={
        "tree_patristic_distance": "X_tree_distance",
        "jaccard_similarity": "Y_spacer_Jaccard",
    }
)

tree_table = tree_table.sort_values(
    ["X_tree_distance", "Y_spacer_Jaccard"]
)

tree_table.to_csv(
    TREE_OUTPUT,
    sep="\t",
    index=False,
)

if len(ani_table) != 45:
    raise SystemExit(
        f"Expected 45 positive-genome pairs, found {len(ani_table)}."
    )

print("Unique spacer sequences:", total)
print(
    "Genome-private:",
    f"{private_count} ({100 * private_count / total:.2f}%)",
)
print(
    "Shared between genomes:",
    f"{shared_count} ({100 * shared_count / total:.2f}%)",
)
print("ANI scatter points:", len(ani_table))
print("Tree-distance scatter points:", len(tree_table))
print("Summary data:", SUMMARY_OUTPUT)
print("ANI scatter data:", ANI_OUTPUT)
print("Tree-distance data:", TREE_OUTPUT)
