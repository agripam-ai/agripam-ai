"""Small, testable analysis core shared by the app and command-line reports."""

from __future__ import annotations

import csv
import hashlib
import re
from collections import Counter
from pathlib import Path
from typing import Iterable


DNA = set("ACGT")
FASTA_DNA = set("ACGTNRYKMSWBDHV")
DNA_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def normalize_dna(sequence: str) -> str:
    sequence = "".join(sequence.upper().split())
    invalid = sorted(set(sequence) - DNA)
    if invalid:
        raise ValueError(f"Sequence contains unsupported symbols: {', '.join(invalid)}")
    return sequence


def parse_fasta(text: str) -> list[tuple[str, str]]:
    """Parse FASTA records while retaining standard IUPAC ambiguity codes."""
    records, chunks = [], []
    name = None
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(chunks).upper()))
            name, chunks = line[1:].strip() or f"sequence_{len(records)+1}", []
        else:
            if name is None:
                name = "sequence_1"
            chunks.append(re.sub(r"\s+", "", line))
    if name is not None:
        records.append((name, "".join(chunks).upper()))
    if not records:
        raise ValueError("No FASTA sequence records were found")
    for record_name, sequence in records:
        if not sequence:
            raise ValueError(f"FASTA record '{record_name}' is empty")
        invalid = sorted(set(sequence) - FASTA_DNA)
        if invalid:
            raise ValueError(f"FASTA record '{record_name}' contains unsupported symbols: {', '.join(invalid)}")
    return records


def genome_summary(records: Iterable[tuple[str, str]]) -> dict[str, float | int | str]:
    records = list(records)
    sequence = "".join(seq for _, seq in records)
    acgt = sum(sequence.count(base) for base in "ACGT")
    gc = sequence.count("G") + sequence.count("C")
    lengths = sorted((len(seq) for _, seq in records), reverse=True)
    cumulative, n50, l50 = 0, 0, 0
    for l50, length in enumerate(lengths, 1):
        cumulative += length
        if cumulative >= len(sequence) / 2:
            n50 = length
            break
    return {"contigs": len(records), "total_length_bp": len(sequence),
            "gc_percent": 100.0 * gc / acgt if acgt else 0.0,
            "ambiguous_bases": len(sequence) - acgt,
            "ambiguous_percent": 100.0 * (len(sequence) - acgt) / len(sequence) if sequence else 0.0,
            "n50_bp": n50, "l50_contigs": l50,
            "longest_contig_bp": lengths[0], "shortest_contig_bp": lengths[-1],
            "sha256": hashlib.sha256(sequence.encode()).hexdigest()}


def genome_qc(records: Iterable[tuple[str, str]]) -> list[dict[str, float | int | str]]:
    """Return taxon-agnostic assembly checks without inventing completeness claims."""
    summary = genome_summary(records)
    return [
        {"check": "sequence_content", "value": round(float(summary["ambiguous_percent"]), 4),
         "unit": "% ambiguous bases", "status": "pass" if summary["ambiguous_percent"] <= 0.1 else "review",
         "interpretation": "High ambiguity can hide PAMs and invalidate guide sequences."},
        {"check": "contiguity", "value": int(summary["contigs"]), "unit": "contigs",
         "status": "pass" if summary["contigs"] <= 20 else "review",
         "interpretation": "Fragmented assemblies can split defense loci and conceal off-target copies."},
        {"check": "N50", "value": int(summary["n50_bp"]), "unit": "bp",
         "status": "pass" if summary["n50_bp"] >= 100_000 else "review",
         "interpretation": "A descriptive contiguity statistic; it is not a completeness estimate."},
        {"check": "assembly_span", "value": int(summary["total_length_bp"]), "unit": "bp",
         "status": "reported",
         "interpretation": "Compare with a trusted close reference; no taxon-specific expected size is assumed."},
    ]


def exact_target_copy_counts(records: Iterable[tuple[str, str]], targets: Iterable[str]) -> dict[str, int]:
    """Count exact target occurrences across both DNA strands using overlapping windows."""
    records = list(records)
    targets = list(dict.fromkeys(targets))
    by_length: dict[int, set[str]] = {}
    for target in targets:
        by_length.setdefault(len(target), set()).add(target)
    result = Counter()
    for _, sequence in records:
        for length, wanted in by_length.items():
            for index in range(0, len(sequence) - length + 1):
                window = sequence[index:index + length]
                if window in wanted:
                    result[window] += 1
                reverse = window.translate(DNA_COMPLEMENT)[::-1]
                if reverse in wanted and reverse != window:
                    result[reverse] += 1
    return {target: result[target] for target in targets}


def _longest_homopolymer(sequence: str) -> int:
    return max((len(m.group(0)) for m in re.finditer(r"(A+|C+|G+|T+)", sequence)), default=0)


def reverse_complement(sequence: str) -> str:
    sequence = normalize_dna(sequence)
    return sequence.translate(DNA_COMPLEMENT)[::-1]


def read_pam_scores(path: Path | str) -> list[dict[str, str]]:
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def score_pam(
    pam: str,
    mobile_nearby: bool,
    score_rows: Iterable[dict[str, str]],
) -> dict[str, float | int | str]:
    pam = normalize_dna(pam)
    if len(pam) != 3:
        raise ValueError("PAM must contain exactly three nucleotides")
    mobile_value = "1" if mobile_nearby else "0"
    for row in score_rows:
        if row["pam3"] == pam and row["mobile_nearby"] == mobile_value:
            return {
                "pam": pam,
                "mobile_nearby": int(mobile_nearby),
                "probability": float(row["predicted_probability"]),
                "percent": float(row["predicted_percent"]),
                "rank": int(row["rank_within_context"]),
            }
    raise KeyError(f"No score available for PAM={pam}, mobile={mobile_value}")


def _scan_oriented(
    oriented: str,
    strand: str,
    original_length: int,
    protospacer_length: int,
) -> list[dict[str, int | str]]:
    hits: list[dict[str, int | str]] = []
    for index in range(0, len(oriented) - protospacer_length - 2):
        if oriented[index:index + 3] != "TTC":
            continue
        proto_start = index + 3
        proto_end = proto_start + protospacer_length
        protospacer = oriented[proto_start:proto_end]
        if len(protospacer) != protospacer_length:
            continue
        if strand == "+":
            start = proto_start + 1
            end = proto_end
        else:
            start = original_length - proto_end + 1
            end = original_length - proto_start
        hits.append({
            "strand": strand,
            "pam": "TTC",
            "protospacer": protospacer,
            "start_1based": start,
            "end_1based": end,
        })
    return hits


def scan_ttc_targets(sequence: str, protospacer_length: int = 35) -> list[dict[str, int | str]]:
    """Find 5'-TTC-protospacer candidates on both sequence orientations."""
    sequence = normalize_dna(sequence)
    if protospacer_length < 20 or protospacer_length > 60:
        raise ValueError("Protospacer length must be between 20 and 60 nt")
    hits = _scan_oriented(sequence, "+", len(sequence), protospacer_length)
    hits.extend(
        _scan_oriented(reverse_complement(sequence), "-", len(sequence), protospacer_length)
    )
    return sorted(hits, key=lambda row: (int(row["start_1based"]), str(row["strand"])))


def scan_ttc_targets_with_ambiguity(sequence: str, protospacer_length: int = 35) -> list[dict[str, int | str]]:
    sequence = sequence.upper()
    invalid = sorted(set(sequence) - FASTA_DNA)
    if invalid:
        raise ValueError(f"Sequence contains unsupported symbols: {', '.join(invalid)}")
    if protospacer_length < 20 or protospacer_length > 60:
        raise ValueError("Protospacer length must be between 20 and 60 nt")
    hits = []
    for start in range(len(sequence) - protospacer_length - 2):
        window = sequence[start:start + 3 + protospacer_length]
        if set(window) <= DNA and window.startswith("TTC"):
            hits.append({"strand": "+", "pam": "TTC", "protospacer": window[3:],
                         "start_1based": start + 4, "end_1based": start + 3 + protospacer_length})
    reverse = sequence.translate(str.maketrans("ACGTNRYKMSWBDHV", "TGCANYRMKSWVHDB"))[::-1]
    for start in range(len(reverse) - protospacer_length - 2):
        window = reverse[start:start + 3 + protospacer_length]
        if set(window) <= DNA and window.startswith("TTC"):
            proto_start, proto_end = start + 3, start + 3 + protospacer_length
            hits.append({"strand": "-", "pam": "TTC", "protospacer": window[3:],
                         "start_1based": len(sequence) - proto_end + 1,
                         "end_1based": len(sequence) - proto_start})
    return sorted(hits, key=lambda row: (int(row["start_1based"]), str(row["strand"])))


def scan_pam_targets_with_ambiguity(sequence: str, pams: Iterable[str],
                                    protospacer_length: int = 35) -> list[dict[str, int | str]]:
    """Find targets for an explicit set of 3-nt PAMs on both orientations."""
    motifs = {normalize_dna(pam) for pam in pams}
    if not motifs or any(len(pam) != 3 for pam in motifs):
        raise ValueError("Provide at least one three-nucleotide PAM")
    sequence = sequence.upper()
    invalid = sorted(set(sequence) - FASTA_DNA)
    if invalid:
        raise ValueError(f"Sequence contains unsupported symbols: {', '.join(invalid)}")
    if protospacer_length < 20 or protospacer_length > 60:
        raise ValueError("Protospacer length must be between 20 and 60 nt")
    hits = []
    orientations = [('+', sequence), ('-', sequence.translate(
        str.maketrans("ACGTNRYKMSWBDHV", "TGCANYRMKSWVHDB"))[::-1])]
    for strand, oriented in orientations:
        for start in range(len(oriented) - protospacer_length - 2):
            window = oriented[start:start + 3 + protospacer_length]
            if set(window) <= DNA and window[:3] in motifs:
                proto_start, proto_end = start + 3, start + 3 + protospacer_length
                if strand == '+':
                    target_start, target_end = proto_start + 1, proto_end
                else:
                    target_start = len(sequence) - proto_end + 1
                    target_end = len(sequence) - proto_start
                hits.append({"strand": strand, "pam": window[:3], "protospacer": window[3:],
                             "start_1based": target_start, "end_1based": target_end})
    return sorted(hits, key=lambda row: (int(row["start_1based"]), str(row["strand"]), str(row["pam"])))


def analyze_ttc_genome(records, protospacer_length, score_rows, mobile_nearby=False):
    """Scan a genome and rank TTC targets using transparent sequence checks."""
    records = list(records)
    summary = genome_summary(records)
    pam_score = score_pam("TTC", mobile_nearby, score_rows)
    candidates, counts = [], {}
    for contig, sequence in records:
        for hit in scan_ttc_targets_with_ambiguity(sequence, protospacer_length):
            hit["contig"] = contig
            candidates.append(hit)
            counts[hit["protospacer"]] = counts.get(hit["protospacer"], 0) + 1
    for index, row in enumerate(candidates, 1):
        proto = str(row["protospacer"])
        gc = 100.0 * (proto.count("G") + proto.count("C")) / len(proto)
        copies, homopolymer = counts[proto], _longest_homopolymer(proto)
        gc_quality = max(0.0, 1.0 - abs(gc - 50.0) / 50.0)
        sequence_quality = max(0.0, 1.0 - max(0, homopolymer - 4) * 0.12)
        design_score = 0.60 * float(pam_score["probability"]) + 0.20 * gc_quality + 0.10 / copies + 0.10 * sequence_quality
        row.update({"candidate_id": f"AGP{index:06d}", "gc_percent": round(gc, 2),
                    "exact_genome_copies": copies, "longest_homopolymer": homopolymer,
                    "pam_model_percent": round(float(pam_score["percent"]), 2),
                    "design_score": round(100 * design_score, 2),
                    "review_status": "review approximate off-targets and gene context"})
    candidates.sort(key=lambda row: (-float(row["design_score"]), str(row["contig"]), int(row["start_1based"])))
    summary.update({"ttc_candidates": len(candidates), "unique_candidate_spacers": len(counts),
                    "mobile_context_assumption": int(mobile_nearby)})
    return summary, candidates
