#!/usr/bin/env python3

from collections import Counter
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]

INPUT = (
    PROJECT
    / "taxonomy/core_genome/parsnp_24/parsnp.snps.mblocks"
)

EXPECTED_LIST = (
    PROJECT
    / "taxonomy/ani/paenibacillus_polymyxa_24_genomes.txt"
)

OUTPUT = (
    PROJECT
    / "taxonomy/core_genome/core_snps_24_clean.fasta"
)

STATS_OUTPUT = (
    PROJECT
    / "taxonomy/core_genome/core_snps_24_alignment_stats.tsv"
)


def read_fasta(path):
    records = []
    identifier = None
    sequence_parts = []

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if identifier is not None:
                    records.append(
                        (identifier, "".join(sequence_parts).upper())
                    )

                identifier = line[1:].split()[0]
                sequence_parts = []
            else:
                sequence_parts.append(line)

    if identifier is not None:
        records.append(
            (identifier, "".join(sequence_parts).upper())
        )

    return records


def clean_identifier(identifier):
    for suffix in [
        ".fna.ref",
        ".fasta.ref",
        ".fa.ref",
        ".fna",
        ".fasta",
        ".fa",
        ".ref",
    ]:
        if identifier.endswith(suffix):
            identifier = identifier[
                : -len(suffix)
            ]
            break

    return identifier


records = read_fasta(INPUT)

if len(records) != 24:
    raise SystemExit(
        f"Expected 24 sequences, found {len(records)}."
    )

cleaned = [
    (clean_identifier(identifier), sequence)
    for identifier, sequence in records
]

identifiers = [identifier for identifier, _ in cleaned]

if len(set(identifiers)) != len(identifiers):
    raise SystemExit(
        "Duplicate identifiers after label cleaning."
    )

expected = {
    Path(line.strip()).stem
    for line in EXPECTED_LIST.read_text().splitlines()
    if line.strip()
}

observed = set(identifiers)

if observed != expected:
    raise SystemExit(
        "Alignment accession mismatch. "
        f"Missing={sorted(expected - observed)}; "
        f"unexpected={sorted(observed - expected)}"
    )

lengths = {
    len(sequence)
    for _, sequence in cleaned
}

if len(lengths) != 1:
    raise SystemExit(
        f"Sequences have unequal lengths: {sorted(lengths)}"
    )

alignment_length = lengths.pop()
allowed = set("ACGTN-?")

unexpected_characters = sorted(
    {
        character
        for _, sequence in cleaned
        for character in sequence
        if character not in allowed
    }
)

if unexpected_characters:
    raise SystemExit(
        "Unexpected sequence characters: "
        + ",".join(unexpected_characters)
    )

variable_sites = 0
parsimony_informative_sites = 0
constant_sites = 0
missing_cells = 0

for position in range(alignment_length):
    column = [
        sequence[position]
        for _, sequence in cleaned
    ]

    missing_cells += sum(
        base in {"N", "-", "?"}
        for base in column
    )

    observed_bases = [
        base
        for base in column
        if base in {"A", "C", "G", "T"}
    ]

    counts = Counter(observed_bases)

    if len(counts) <= 1:
        constant_sites += 1
    else:
        variable_sites += 1

        if sum(count >= 2 for count in counts.values()) >= 2:
            parsimony_informative_sites += 1

with OUTPUT.open("w") as handle:
    for identifier, sequence in sorted(cleaned):
        handle.write(f">{identifier}\n")

        for start in range(0, len(sequence), 80):
            handle.write(sequence[start:start + 80] + "\n")

total_cells = len(cleaned) * alignment_length
missing_fraction = (
    missing_cells / total_cells
    if total_cells
    else 0.0
)

stats = [
    ("sequences", len(cleaned)),
    ("alignment_length", alignment_length),
    ("constant_sites", constant_sites),
    ("variable_sites", variable_sites),
    (
        "parsimony_informative_sites",
        parsimony_informative_sites,
    ),
    ("missing_cells", missing_cells),
    ("missing_fraction", f"{missing_fraction:.8f}"),
]

with STATS_OUTPUT.open("w") as handle:
    handle.write("metric\tvalue\n")

    for metric, value in stats:
        handle.write(f"{metric}\t{value}\n")

print("Sequences:", len(cleaned))
print("Alignment length:", alignment_length)
print("Variable sites:", variable_sites)
print(
    "Parsimony-informative sites:",
    parsimony_informative_sites,
)
print("Constant sites:", constant_sites)
print("Missing fraction:", f"{missing_fraction:.6%}")
print("Clean alignment:", OUTPUT)
print("Statistics:", STATS_OUTPUT)
