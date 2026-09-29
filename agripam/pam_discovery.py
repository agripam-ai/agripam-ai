"""Evidence-led, strain-specific PAM inference.

The module deliberately separates PAM inference from target-site scanning.  PAMs
are inferred only from spacer matches in an external mobile-sequence collection;
an uploaded host genome by itself is not PAM evidence.
"""
from __future__ import annotations

import csv
import hashlib
import itertools
import math
import re
from collections import Counter
from pathlib import Path

from .core import DNA, parse_fasta, reverse_complement


def extract_cctyper_spacers(output_dir: Path) -> tuple[list[dict], list[dict]]:
    """Read trusted arrays, subtypes and spacer FASTAs from CRISPRCasTyper."""
    output_dir = Path(output_dir)
    arrays: dict[str, dict] = {}
    table = output_dir / "crisprs_all.tab"
    if table.exists():
        with table.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                trusted = str(row.get("Trusted", "")).lower() == "true"
                arrays[row.get("CRISPR", "")] = {
                    "array_id": row.get("CRISPR", ""),
                    "contig": row.get("Contig", ""),
                    "trusted": trusted,
                    "array_prediction": row.get("Prediction", ""),
                    "subtype": row.get("Subtype", "") or row.get("Prediction", ""),
                    "subtype_probability": row.get("Subtype_probability", ""),
                }
    spacers = []
    for fasta in sorted((output_dir / "spacers").glob("*.fa")):
        array_id = fasta.stem
        metadata = arrays.get(array_id, {"array_id": array_id, "contig": "", "trusted": False,
                                         "array_prediction": "", "subtype": "", "subtype_probability": ""})
        for spacer_id, sequence in parse_fasta(fasta.read_text()):
            if set(sequence) <= DNA:
                spacers.append({**metadata, "spacer_id": spacer_id, "spacer": sequence})
    return list(arrays.values()), spacers


def _wilson(successes: int, total: int, z: float = 1.96) -> tuple[float, float]:
    if total == 0:
        return 0.0, 0.0
    p = successes / total
    denominator = 1 + z * z / total
    centre = (p + z * z / (2 * total)) / denominator
    margin = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denominator
    return max(0.0, centre - margin), min(1.0, centre + margin)


def _binomial_tail(k: int, n: int, p: float) -> float:
    if n == 0 or k == 0:
        return 1.0
    if p <= 0:
        return 0.0
    if p >= 1:
        return 1.0
    return min(1.0, sum(math.comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(k, n + 1)))


def _bh_adjust(rows: list[dict]) -> None:
    ordered = sorted(enumerate(rows), key=lambda item: item[1]["p_value"])
    running = 1.0
    for reverse_rank, (original, row) in enumerate(reversed(ordered), 1):
        rank = len(rows) - reverse_rank + 1
        running = min(running, row["p_value"] * len(rows) / rank)
        rows[original]["bh_q_value"] = running


def _find_exact(sequence: str, query: str):
    start = sequence.find(query)
    while start >= 0:
        yield start
        start = sequence.find(query, start + 1)


def infer_pams_from_mobile_references(
    spacers: list[dict], reference_fasta: str, pam_length: int = 3,
) -> tuple[dict, list[dict], list[dict], list[dict]]:
    """Infer upstream PAMs from exact spacer–protospacer matches on both strands.

    This first auditable implementation uses exact full-length matches.  Each
    spacer/reference/orientation/locus is an observation; duplicated observations
    remain visible and can be collapsed by users during sensitivity analyses.
    """
    if pam_length != 3:
        raise ValueError("This release ranks the complete 64-member 3-nt PAM space")
    references = parse_fasta(reference_fasta)
    observations = []
    reference_triplets = Counter()
    reference_triplet_total = 0
    for reference_id, sequence in references:
        tier_match = re.search(r"evidence_tier=(\d+)", reference_id)
        scope_match = re.search(r"scope=([^\s]+)", reference_id)
        evidence_tier = int(tier_match.group(1)) if tier_match else 5
        reference_scope = scope_match.group(1).replace("_", " ") if scope_match else "user-supplied or broad collection"
        clean = "".join(base for base in sequence if base in DNA)
        for index in range(0, len(clean) - pam_length + 1):
            triplet = clean[index:index + pam_length]
            if len(triplet) == pam_length:
                reference_triplets[triplet] += 1
                reference_triplet_total += 1
        reverse = sequence.translate(str.maketrans("ACGTNRYKMSWBDHV", "TGCANYRMKSWVHDB"))[::-1]
        oriented_references = [("+", sequence), ("-", reverse)]
        for spacer in spacers:
            query = str(spacer["spacer"])
            for strand, oriented in oriented_references:
                for start in _find_exact(oriented, query):
                    if start < pam_length:
                        continue
                    pam = oriented[start - pam_length:start]
                    if set(pam) != DNA.intersection(set(pam)) or len(pam) != pam_length:
                        continue
                    if strand == "+":
                        reference_start = start + 1
                        reference_end = start + len(query)
                    else:
                        reference_start = len(sequence) - (start + len(query)) + 1
                        reference_end = len(sequence) - start
                    observations.append({
                        "spacer_id": spacer["spacer_id"], "array_id": spacer["array_id"],
                        "crispr_subtype": spacer.get("subtype", ""), "reference_id": reference_id,
                        "reference_strand": strand, "protospacer_start_1based": reference_start,
                        "protospacer_end_1based": reference_end, "match_identity_percent": 100.0,
                        "match_length": len(query), "pam3_normalized": pam,
                        "evidence_tier": evidence_tier, "reference_scope": reference_scope,
                        "orientation_rule": "5-prime upstream in protospacer orientation",
                    })

    counts = Counter(row["pam3_normalized"] for row in observations)
    total = len(observations)
    ranking = []
    for pam_tuple in itertools.product("ACGT", repeat=3):
        pam = "".join(pam_tuple)
        count = counts[pam]
        expected = reference_triplets[pam] / reference_triplet_total if reference_triplet_total else 1 / 64
        lower, upper = _wilson(count, total)
        ranking.append({
            "pam3": pam, "observations": count, "total_observations": total,
            "observed_percent": 100 * count / total if total else 0.0,
            "wilson_95ci_lower_percent": 100 * lower,
            "wilson_95ci_upper_percent": 100 * upper,
            "reference_background_percent": 100 * expected,
            "fold_enrichment": (count / total) / expected if total and expected else 0.0,
            "p_value": _binomial_tail(count, total, expected),
            "best_evidence_tier": min((int(obs["evidence_tier"]) for obs in observations
                                       if obs["pam3_normalized"] == pam), default=""),
        })
    _bh_adjust(ranking)
    ranking.sort(key=lambda row: (-row["observations"], row["bh_q_value"], -row["fold_enrichment"], row["pam3"]))
    for rank, row in enumerate(ranking, 1):
        row["rank"] = rank
        n = int(row["observations"])
        q = float(row["bh_q_value"])
        tier = int(row["best_evidence_tier"]) if row["best_evidence_tier"] != "" else 9
        row["confidence"] = ("high computational support" if n >= 10 and q < 0.05 and tier <= 2
                             else "moderate computational support" if n >= 5 and q < 0.10 and tier <= 3
                             else "exploratory taxonomically transferred evidence" if n and tier == 4
                             else "limited evidence" if n else "not observed")
        row["required_validation"] = ("independent computational replication; locus confirmation; "
                                      "PAM-library or interference reporter assay with PAM-mutant and no-guide controls"
                                      if n else "not applicable")

    logo = []
    for position in range(3):
        position_counts = Counter(row["pam3_normalized"][position] for row in observations)
        for nucleotide in "ACGT":
            logo.append({"position": position + 1, "nucleotide": nucleotide,
                         "count": position_counts[nucleotide],
                         "frequency": position_counts[nucleotide] / total if total else 0.0})
    digest = hashlib.sha256("".join(seq for _, seq in references).encode()).hexdigest()
    summary = {
        "trusted_arrays": len({s["array_id"] for s in spacers if s.get("trusted")}),
        "spacers_extracted": len(spacers), "mobile_reference_records": len(references),
        "exact_protospacer_matches": total, "unique_spacers_matched": len({r["spacer_id"] for r in observations}),
        "reference_sha256": digest, "match_policy": "exact full-length spacer match",
        "pam_orientation": "5-prime upstream in normalized protospacer orientation",
        "interpretation": "computational PAM inference; experimental PAM-library validation pending",
    }
    return summary, observations, ranking, logo


def supported_pams(ranking: list[dict], minimum_observations: int = 3) -> list[str]:
    """Return motifs with evidence, never silently falling back to TTC."""
    return [str(row["pam3"]) for row in ranking
            if int(row["observations"]) >= minimum_observations and float(row["bh_q_value"]) < 0.10]
