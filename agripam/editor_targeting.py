"""Editor-specific PAM scanning for introduced tools (not native inference)."""
from __future__ import annotations
from collections import Counter

IUPAC = {"A":"A", "C":"C", "G":"G", "T":"T", "N":"ACGT", "R":"AG", "Y":"CT", "S":"GC", "W":"AT", "K":"GT", "M":"AC", "B":"CGT", "D":"AGT", "H":"ACT", "V":"ACG"}
EDITOR_PRESETS = {
    "SpCas9": {"pam_pattern":"NGG", "pam_side":"3prime", "protospacer_length":20, "application":"DNA cleavage"},
    "Cas12a/Cpf1": {"pam_pattern":"TTTV", "pam_side":"5prime", "protospacer_length":23, "application":"DNA cleavage"},
    "dCas9/CRISPRi": {"pam_pattern":"NGG", "pam_side":"3prime", "protospacer_length":20, "application":"transcriptional repression"},
}

def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

def pam_matches(seq: str, pattern: str) -> bool:
    return len(seq) == len(pattern) and all(base in IUPAC.get(code, code) for base, code in zip(seq.upper(), pattern.upper()))

def scan_editor_targets(records, editor_name: str, pam_pattern: str, pam_side: str = "3prime", protospacer_length: int = 20):
    """Scan both strands and return target coordinates (1-based, inclusive)."""
    out = []
    p = pam_pattern.upper().strip()
    for contig, raw in records:
        seq = raw.upper()
        for strand, oriented in (("+", seq), ("-", reverse_complement(seq))):
            for i in range(len(oriented) - len(p) - protospacer_length + 1):
                if pam_side == "3prime":
                    prot = oriented[i:i+protospacer_length]; pam = oriented[i+protospacer_length:i+protospacer_length+len(p)]
                    a, b = i, i + protospacer_length + len(p) - 1
                else:
                    pam = oriented[i:i+len(p)]; prot = oriented[i+len(p):i+len(p)+protospacer_length]
                    a, b = i, i + len(p) + protospacer_length - 1
                if not pam_matches(pam, p): continue
                if strand == "+":
                    ps, pe = a + 1, b + 1
                    pam_start = (i + protospacer_length + 1) if pam_side == "3prime" else i + 1
                else:
                    ps, pe = len(seq) - b, len(seq) - a
                    pam_start = len(seq) - (i + len(p)) + 1 if pam_side == "3prime" else len(seq) - i - len(p) + 1
                out.append({"editor": editor_name, "pam_pattern": p, "pam": pam, "protospacer": prot, "strand": strand,
                            "contig": contig, "start": ps, "end": pe, "pam_start": pam_start,
                            "pam_end": pam_start + len(p) - 1, "application": EDITOR_PRESETS.get(editor_name, {}).get("application", "custom editing hypothesis")})
    counts = Counter(x["protospacer"] for x in out)
    for x in out:
        x["guide_occurrences"] = counts[x["protospacer"]]
        x["guide_uniqueness"] = "unique" if counts[x["protospacer"]] == 1 else "repeated"
    return out


def approximate_offtarget_counts(targets: list[dict], max_mismatches: int = 2,
                                 max_queries: int = 2000,
                                 max_reference_sequences: int = 500_000) -> dict[str, dict]:
    """Screen top guide sequences against PAM-compatible sites with <=2 mismatches.

    Three exact seed partitions guarantee retrieval of every sequence with at most
    two substitutions. Indels/bulges are deliberately outside this bounded screen.
    """
    if max_mismatches != 2:
        raise ValueError("The audited seed scheme currently supports exactly two mismatches.")
    results: dict[str, dict] = {}
    by_editor: dict[str, list[dict]] = {}
    for row in targets:
        by_editor.setdefault(str(row["editor"]), []).append(row)
    for editor, rows in by_editor.items():
        site_counts = Counter(str(row["protospacer"]) for row in rows)
        queries = list(dict.fromkeys(str(row["protospacer"]) for row in rows[:max_queries]))
        if len(site_counts) > max_reference_sequences:
            for query in queries:
                results[f"{editor}\t{query}"] = {"status": "skipped_reference_limit", "count": None,
                                                    "nearest_mismatches": None}
            continue
        length = len(next(iter(site_counts), ""))
        cuts = (0, length // 3, 2 * length // 3, length)
        index: dict[tuple[int, str], set[str]] = {}
        for sequence in site_counts:
            for part in range(3):
                index.setdefault((part, sequence[cuts[part]:cuts[part + 1]]), set()).add(sequence)
        for query in queries:
            candidates: set[str] = set()
            for part in range(3):
                candidates.update(index.get((part, query[cuts[part]:cuts[part + 1]]), set()))
            mismatch_sites = Counter()
            for candidate in candidates:
                distance = sum(a != b for a, b in zip(query, candidate))
                if 0 < distance <= 2:
                    mismatch_sites[distance] += site_counts[candidate]
            nearest = min(mismatch_sites) if mismatch_sites else None
            results[f"{editor}\t{query}"] = {
                "status": "complete_substitutions_only", "count": int(sum(mismatch_sites.values())),
                "nearest_mismatches": nearest,
            }
    return results
