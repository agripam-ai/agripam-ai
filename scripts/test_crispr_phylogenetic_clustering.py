#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo


SEED = 42
PERMUTATIONS = 10000

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

SUMMARY_OUTPUT = (
    PROJECT
    / "results/crispr_phylogenetic_clustering.tsv"
)

NULL_OUTPUT = (
    PROJECT
    / "controls/crispr_phylogenetic_clustering_seed42.tsv"
)

TRAITS = {
    "any_crispr": "has_crispr",
    "any_cas": "has_cas",
    "type_I_C": "has_type_ic",
    "type_III_B": "has_type_iiib",
}


def fitch_score(tree, states):
    score = 0

    def visit(clade):
        nonlocal score

        if clade.is_terminal():
            return {int(states[clade.name])}

        child_sets = [
            visit(child)
            for child in clade.clades
        ]

        current = child_sets[0]

        for child_set in child_sets[1:]:
            intersection = current & child_set

            if intersection:
                current = intersection
            else:
                current = current | child_set
                score += 1

        return current

    visit(tree.root)
    return score


tree = Phylo.read(TREE_INPUT, "newick")
metadata = pd.read_csv(METADATA_INPUT, sep="\t")

tips = [
    terminal.name
    for terminal in tree.get_terminals()
]

if set(tips) != set(metadata["accession"]):
    raise SystemExit("Tree tips and metadata accessions differ.")

metadata = metadata.set_index("accession").loc[tips]

rng = np.random.default_rng(SEED)

summary_rows = []
null_rows = []

for trait_name, column in TRAITS.items():
    observed_values = (
        metadata[column].astype(int).to_numpy()
    )

    states = dict(zip(tips, observed_values))
    observed_score = fitch_score(tree, states)

    null_scores = np.empty(
        PERMUTATIONS,
        dtype=int,
    )

    for permutation in range(PERMUTATIONS):
        shuffled = rng.permutation(observed_values)

        null_scores[permutation] = fitch_score(
            tree,
            dict(zip(tips, shuffled)),
        )

    lower_tail_p = (
        1
        + np.sum(null_scores <= observed_score)
    ) / (PERMUTATIONS + 1)

    upper_tail_p = (
        1
        + np.sum(null_scores >= observed_score)
    ) / (PERMUTATIONS + 1)

    summary_rows.append(
        {
            "trait": trait_name,
            "positive_tips": int(
                observed_values.sum()
            ),
            "total_tips": len(observed_values),
            "observed_fitch_changes": observed_score,
            "null_mean_changes": null_scores.mean(),
            "null_sd_changes": null_scores.std(
                ddof=1
            ),
            "null_min_changes": null_scores.min(),
            "null_max_changes": null_scores.max(),
            "clustering_lower_tail_p": lower_tail_p,
            "dispersion_upper_tail_p": upper_tail_p,
            "permutations": PERMUTATIONS,
            "random_seed": SEED,
        }
    )

    for index, value in enumerate(
        null_scores,
        start=1,
    ):
        null_rows.append(
            {
                "trait": trait_name,
                "permutation": index,
                "fitch_changes": int(value),
                "random_seed": SEED,
            }
        )

summary = pd.DataFrame(summary_rows)
null = pd.DataFrame(null_rows)

summary.to_csv(
    SUMMARY_OUTPUT,
    sep="\t",
    index=False,
)

null.to_csv(
    NULL_OUTPUT,
    sep="\t",
    index=False,
)

print(
    summary[
        [
            "trait",
            "positive_tips",
            "observed_fitch_changes",
            "null_mean_changes",
            "clustering_lower_tail_p",
        ]
    ].to_string(index=False)
)

print()
print("Permutations per trait:", PERMUTATIONS)
print("Summary:", SUMMARY_OUTPUT)
print("Permutation distributions:", NULL_OUTPUT)
