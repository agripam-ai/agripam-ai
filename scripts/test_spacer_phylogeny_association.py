#!/usr/bin/env python3

from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo
from scipy.stats import spearmanr


SEED = 42
PERMUTATIONS = 10000

PROJECT = Path(__file__).resolve().parents[1]

PRESENCE_INPUT = (
    PROJECT / "spacers/genome_by_spacer_presence.tsv"
)

ANI_INPUT = (
    PROJECT / "taxonomy/ani/fastani_symmetric_matrix.tsv"
)

TREE_INPUT = (
    PROJECT
    / "taxonomy/core_genome/iqtree/"
      "paenibacillus_polymyxa_core_snps_midpoint_rooted.treefile"
)

SUMMARY_OUTPUT = (
    PROJECT
    / "results/spacer_sharing_phylogeny_association.tsv"
)

PAIR_OUTPUT = (
    PROJECT
    / "results/spacer_sharing_positive_genome_pairs.tsv"
)

NULL_OUTPUT = (
    PROJECT
    / "controls/spacer_phylogeny_permutations_seed42.tsv"
)


presence = pd.read_csv(
    PRESENCE_INPUT,
    sep="\t",
    index_col=0,
)

ani = pd.read_csv(
    ANI_INPUT,
    sep="\t",
    index_col=0,
)

tree = Phylo.read(TREE_INPUT, "newick")

positive_genomes = sorted(
    presence.index[
        presence.sum(axis=1) > 0
    ]
)

if len(positive_genomes) != 10:
    raise SystemExit(
        f"Expected 10 spacer-positive genomes, "
        f"found {len(positive_genomes)}."
    )

profiles = {
    accession: set(
        presence.columns[
            presence.loc[accession].to_numpy() == 1
        ]
    )
    for accession in positive_genomes
}

pairs = list(combinations(positive_genomes, 2))

pair_rows = []

for accession_a, accession_b in pairs:
    set_a = profiles[accession_a]
    set_b = profiles[accession_b]

    intersection = len(set_a & set_b)
    union = len(set_a | set_b)

    pair_rows.append(
        {
            "accession_a": accession_a,
            "accession_b": accession_b,
            "symmetric_ani": ani.loc[
                accession_a,
                accession_b,
            ],
            "ani_distance": (
                100
                - ani.loc[accession_a, accession_b]
            ),
            "tree_patristic_distance": tree.distance(
                accession_a,
                accession_b,
            ),
            "shared_spacers": intersection,
            "union_spacers": union,
            "jaccard_similarity": intersection / union,
        }
    )

pair_table = pd.DataFrame(pair_rows)
pair_table.to_csv(
    PAIR_OUTPUT,
    sep="\t",
    index=False,
)

observed_ani_rho = spearmanr(
    pair_table["symmetric_ani"],
    pair_table["jaccard_similarity"],
).statistic

observed_tree_rho = spearmanr(
    pair_table["tree_patristic_distance"],
    pair_table["jaccard_similarity"],
).statistic

rng = np.random.default_rng(SEED)
null_rows = []

genome_array = np.array(positive_genomes)

for permutation in range(1, PERMUTATIONS + 1):
    permuted_sources = rng.permutation(
        genome_array
    )

    permuted_profile = {
        target: profiles[source]
        for target, source in zip(
            positive_genomes,
            permuted_sources,
        )
    }

    jaccard_values = []

    for accession_a, accession_b in pairs:
        set_a = permuted_profile[accession_a]
        set_b = permuted_profile[accession_b]

        union = len(set_a | set_b)

        jaccard_values.append(
            len(set_a & set_b) / union
        )

    ani_rho = spearmanr(
        pair_table["symmetric_ani"],
        jaccard_values,
    ).statistic

    tree_rho = spearmanr(
        pair_table["tree_patristic_distance"],
        jaccard_values,
    ).statistic

    null_rows.append(
        {
            "permutation": permutation,
            "ani_similarity_rho": ani_rho,
            "tree_distance_rho": tree_rho,
            "random_seed": SEED,
        }
    )

null = pd.DataFrame(null_rows)

ani_p = (
    1
    + int(
        (
            null["ani_similarity_rho"]
            >= observed_ani_rho
        ).sum()
    )
) / (PERMUTATIONS + 1)

tree_p = (
    1
    + int(
        (
            null["tree_distance_rho"]
            <= observed_tree_rho
        ).sum()
    )
) / (PERMUTATIONS + 1)

summary = pd.DataFrame(
    [
        {
            "comparison": "ANI_similarity_vs_spacer_Jaccard",
            "genomes": len(positive_genomes),
            "genome_pairs": len(pairs),
            "spearman_rho": observed_ani_rho,
            "permutation_p_one_sided": ani_p,
            "expected_direction": "positive",
            "permutations": PERMUTATIONS,
            "random_seed": SEED,
        },
        {
            "comparison": "tree_distance_vs_spacer_Jaccard",
            "genomes": len(positive_genomes),
            "genome_pairs": len(pairs),
            "spearman_rho": observed_tree_rho,
            "permutation_p_one_sided": tree_p,
            "expected_direction": "negative",
            "permutations": PERMUTATIONS,
            "random_seed": SEED,
        },
    ]
)

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

print("Spacer-positive genomes:", len(positive_genomes))
print("Genome pairs tested:", len(pairs))
print(
    "ANI similarity vs spacer Jaccard:",
    f"rho={observed_ani_rho:.4f}, "
    f"permutation p={ani_p:.6f}",
)
print(
    "Tree distance vs spacer Jaccard:",
    f"rho={observed_tree_rho:.4f}, "
    f"permutation p={tree_p:.6f}",
)
print("Summary:", SUMMARY_OUTPUT)
print("Pairwise data:", PAIR_OUTPUT)
print("Permutation data:", NULL_OUTPUT)
