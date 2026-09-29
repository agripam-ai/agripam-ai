#!/usr/bin/env python3

from pathlib import Path

import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
INPUT = PROJECT / "taxonomy/ani/fastani_all_vs_all.tsv"
OUTPUT_DIR = PROJECT / "taxonomy/ani"

PAIRWISE_OUTPUT = OUTPUT_DIR / "fastani_unordered_pairs.tsv"
MATRIX_OUTPUT = OUTPUT_DIR / "fastani_symmetric_matrix.tsv"

columns = [
    "query_path",
    "reference_path",
    "ani",
    "matched_fragments",
    "total_query_fragments",
]

data = pd.read_csv(
    INPUT,
    sep="\t",
    header=None,
    names=columns,
)

for column in [
    "ani",
    "matched_fragments",
    "total_query_fragments",
]:
    data[column] = pd.to_numeric(data[column])

data["query"] = (
    data["query_path"]
    .map(lambda value: Path(value).stem)
)

data["reference"] = (
    data["reference_path"]
    .map(lambda value: Path(value).stem)
)

data["query_coverage"] = (
    data["matched_fragments"]
    / data["total_query_fragments"]
)

accessions = sorted(
    set(data["query"]) | set(data["reference"])
)

if len(accessions) != 24:
    raise SystemExit(
        f"Expected 24 accessions, found {len(accessions)}."
    )

expected_rows = len(accessions) ** 2
if len(data) != expected_rows:
    raise SystemExit(
        f"Expected {expected_rows} directed comparisons, "
        f"found {len(data)}."
    )

lookup = {
    (row.query, row.reference): row
    for row in data.itertuples(index=False)
}

pair_rows = []

for index, accession_a in enumerate(accessions):
    for accession_b in accessions[index + 1:]:
        forward = lookup[(accession_a, accession_b)]
        reverse = lookup[(accession_b, accession_a)]

        pair_rows.append(
            {
                "accession_a": accession_a,
                "accession_b": accession_b,
                "ani_a_to_b": forward.ani,
                "ani_b_to_a": reverse.ani,
                "symmetric_ani": (
                    forward.ani + reverse.ani
                ) / 2,
                "ani_directional_difference": abs(
                    forward.ani - reverse.ani
                ),
                "coverage_a_to_b": forward.query_coverage,
                "coverage_b_to_a": reverse.query_coverage,
                "minimum_reciprocal_coverage": min(
                    forward.query_coverage,
                    reverse.query_coverage,
                ),
            }
        )

pairs = pd.DataFrame(pair_rows).sort_values(
    ["symmetric_ani", "accession_a", "accession_b"],
    ascending=[False, True, True],
)

pairs.to_csv(
    PAIRWISE_OUTPUT,
    sep="\t",
    index=False,
)

matrix = pd.DataFrame(
    index=accessions,
    columns=accessions,
    dtype=float,
)

for accession in accessions:
    matrix.loc[accession, accession] = 100.0

for row in pairs.itertuples(index=False):
    matrix.loc[
        row.accession_a,
        row.accession_b,
    ] = row.symmetric_ani

    matrix.loc[
        row.accession_b,
        row.accession_a,
    ] = row.symmetric_ani

matrix.index.name = "accession"
matrix.to_csv(
    MATRIX_OUTPUT,
    sep="\t",
)

print("Accessions:", len(accessions))
print("Directed comparisons:", len(data))
print("Unique unordered pairs:", len(pairs))
print(
    "Minimum off-diagonal ANI:",
    f"{pairs['symmetric_ani'].min():.4f}%",
)
print(
    "Maximum off-diagonal ANI:",
    f"{pairs['symmetric_ani'].max():.4f}%",
)
print(
    "Pairs below 95% ANI:",
    int((pairs["symmetric_ani"] < 95).sum()),
)
print(
    "Minimum reciprocal alignment coverage:",
    f"{pairs['minimum_reciprocal_coverage'].min():.4f}",
)
print(
    "Maximum directional ANI difference:",
    f"{pairs['ani_directional_difference'].max():.4f}",
)
print()
print("FIVE MOST DIVERGENT PAIRS")
print(
    pairs.sort_values("symmetric_ani")
    .head(5)[
        [
            "accession_a",
            "accession_b",
            "symmetric_ani",
            "minimum_reciprocal_coverage",
        ]
    ]
    .to_string(index=False)
)
print()
print("Pairwise table:", PAIRWISE_OUTPUT)
print("Symmetric matrix:", MATRIX_OUTPUT)
