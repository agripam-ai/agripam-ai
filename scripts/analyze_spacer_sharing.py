#!/usr/bin/env python3

from itertools import combinations
from pathlib import Path

import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]

MEMBERSHIP_INPUT = (
    PROJECT / "spacers/spacer_membership.tsv"
)

METADATA_INPUT = (
    PROJECT
    / "taxonomy/core_genome/tree_tip_crispr_metadata.tsv"
)

MATRIX_OUTPUT = (
    PROJECT
    / "spacers/genome_by_spacer_presence.tsv"
)

PAIRWISE_OUTPUT = (
    PROJECT
    / "results/spacer_sharing_pairwise.tsv"
)

SPACER_SUMMARY_OUTPUT = (
    PROJECT
    / "results/spacer_sharing_by_spacer.tsv"
)

GENOME_SUMMARY_OUTPUT = (
    PROJECT
    / "results/spacer_sharing_by_genome.tsv"
)


membership = pd.read_csv(
    MEMBERSHIP_INPUT,
    sep="\t",
)

metadata = pd.read_csv(
    METADATA_INPUT,
    sep="\t",
)

required_membership = {
    "query_id",
    "accession",
    "strain",
    "crispr",
    "spacer_id",
}

missing = required_membership - set(membership.columns)

if missing:
    raise SystemExit(
        "Missing membership columns: "
        + ", ".join(sorted(missing))
    )

accessions = sorted(metadata["accession"].unique())
query_ids = sorted(membership["query_id"].unique())

if len(accessions) != 24:
    raise SystemExit(
        f"Expected 24 genomes, found {len(accessions)}."
    )

if len(query_ids) != 249:
    raise SystemExit(
        f"Expected 249 unique spacer sequences, "
        f"found {len(query_ids)}."
    )

presence = (
    membership[
        ["accession", "query_id"]
    ]
    .drop_duplicates()
    .assign(present=1)
    .pivot(
        index="accession",
        columns="query_id",
        values="present",
    )
    .reindex(
        index=accessions,
        columns=query_ids,
        fill_value=0,
    )
    .fillna(0)
    .astype(int)
)

presence.index.name = "accession"
presence.to_csv(MATRIX_OUTPUT, sep="\t")


spacer_rows = []

for query_id in query_ids:
    rows = membership[
        membership["query_id"] == query_id
    ]

    genome_values = sorted(
        rows["accession"].unique()
    )

    strain_values = sorted(
        rows["strain"].astype(str).unique()
    )

    spacer_rows.append(
        {
            "query_id": query_id,
            "spacer_records": len(rows),
            "genomes_present": len(genome_values),
            "arrays_present": rows["crispr"].nunique(),
            "accessions": ";".join(genome_values),
            "strains": ";".join(strain_values),
        }
    )

spacer_summary = pd.DataFrame(spacer_rows).sort_values(
    [
        "genomes_present",
        "spacer_records",
        "query_id",
    ],
    ascending=[False, False, True],
)

spacer_summary.to_csv(
    SPACER_SUMMARY_OUTPUT,
    sep="\t",
    index=False,
)


strain_lookup = dict(
    zip(metadata["accession"], metadata["strain"])
)

pair_rows = []

for accession_a, accession_b in combinations(
    accessions,
    2,
):
    set_a = set(
        presence.columns[
            presence.loc[accession_a].to_numpy() == 1
        ]
    )

    set_b = set(
        presence.columns[
            presence.loc[accession_b].to_numpy() == 1
        ]
    )

    intersection = len(set_a & set_b)
    union = len(set_a | set_b)

    jaccard = (
        intersection / union
        if union > 0
        else float("nan")
    )

    pair_rows.append(
        {
            "accession_a": accession_a,
            "strain_a": strain_lookup[accession_a],
            "accession_b": accession_b,
            "strain_b": strain_lookup[accession_b],
            "spacers_a": len(set_a),
            "spacers_b": len(set_b),
            "shared_spacers": intersection,
            "union_spacers": union,
            "jaccard_similarity": jaccard,
            "jaccard_distance": (
                1 - jaccard
                if union > 0
                else float("nan")
            ),
        }
    )

pairwise = pd.DataFrame(pair_rows).sort_values(
    [
        "shared_spacers",
        "jaccard_similarity",
    ],
    ascending=False,
)

pairwise.to_csv(
    PAIRWISE_OUTPUT,
    sep="\t",
    index=False,
)


genome_rows = []

for accession in accessions:
    genome_membership = membership[
        membership["accession"] == accession
    ]

    present_queries = set(
        genome_membership["query_id"]
    )

    shared_queries = set(
        spacer_summary.loc[
            spacer_summary["genomes_present"] > 1,
            "query_id",
        ]
    )

    genome_rows.append(
        {
            "accession": accession,
            "strain": strain_lookup[accession],
            "spacer_records": len(genome_membership),
            "unique_spacer_sequences": len(
                present_queries
            ),
            "spacers_shared_with_other_genomes": len(
                present_queries & shared_queries
            ),
            "genome_private_spacers": len(
                present_queries - shared_queries
            ),
        }
    )

genome_summary = pd.DataFrame(genome_rows).sort_values(
    [
        "unique_spacer_sequences",
        "accession",
    ],
    ascending=[False, True],
)

genome_summary.to_csv(
    GENOME_SUMMARY_OUTPUT,
    sep="\t",
    index=False,
)


shared_spacers = int(
    (spacer_summary["genomes_present"] > 1).sum()
)

private_spacers = int(
    (spacer_summary["genomes_present"] == 1).sum()
)

pairs_with_shared = int(
    (pairwise["shared_spacers"] > 0).sum()
)

print("Genomes:", len(accessions))
print("Exact unique spacer sequences:", len(query_ids))
print("Genome-private spacer sequences:", private_spacers)
print(
    "Spacer sequences found in multiple genomes:",
    shared_spacers,
)
print("Genome pairs:", len(pairwise))
print(
    "Genome pairs sharing at least one spacer:",
    pairs_with_shared,
)
print()
print("TOP FIVE GENOME PAIRS BY SHARED SPACERS")
print(
    pairwise.head(5)[
        [
            "strain_a",
            "strain_b",
            "shared_spacers",
            "jaccard_similarity",
        ]
    ].to_string(index=False)
)
print()
print("Presence matrix:", MATRIX_OUTPUT)
print("Pairwise sharing:", PAIRWISE_OUTPUT)
print("Spacer summary:", SPACER_SUMMARY_OUTPUT)
print("Genome summary:", GENOME_SUMMARY_OUTPUT)
