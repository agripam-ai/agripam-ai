"""BioBrick parts, in-silico construct assembly, compatibility checks and an experiment tracker.

The module turns a chosen edit (reporter insertion at a neutral site, or a CRISPRi guide cassette)
into an annotated DNA construct, checks it against standard-assembly rules and the chassis's
restriction motifs, and tracks what happens at the bench.

Every check is a computational screen. A compatible part is not a working part: activity in the
chosen chassis has to be measured.
"""
from __future__ import annotations

import datetime as _dt
import re
import urllib.parse
import urllib.request
from collections import Counter
from typing import Callable

ROLES = ["promoter", "rbs", "cds", "reporter", "terminator", "scaffold", "backbone", "other"]
REGISTRY_ROLE = {"regulatory": "promoter", "rbs": "rbs", "coding": "cds", "reporter": "reporter",
                 "terminator": "terminator", "plasmid_backbone": "backbone", "plasmid": "backbone"}

RFC10 = {"EcoRI": "GAATTC", "XbaI": "TCTAGA", "SpeI": "ACTAGT", "PstI": "CTGCAG"}
RFC10_EXTRA = {"NotI": "GCGGCCGC"}
RFC25_EXTRA = {"AgeI": "ACCGGT", "NgoMIV": "GCCGGC"}
GOLDEN_GATE = {"BsaI": "GGTCTC", "BsmBI": "CGTCTC"}
STANDARDS = {
    "RFC10 (BioBrick)": {**RFC10, **RFC10_EXTRA},
    "RFC25 (Freiburg)": {**RFC10, **RFC10_EXTRA, **RFC25_EXTRA},
    "Golden Gate (BsaI/BsmBI)": dict(GOLDEN_GATE),
}
SCARS = {"8-bp scar TACTAGAG (standard)": "TACTAGAG", "6-bp scar TACTAG (in-frame fusions)": "TACTAG",
         "Seamless (no scar)": ""}
PREFIX = "GAATTCGCGGCCGCTTCTAGAG"
SUFFIX = "TACTAGTAGCGGCCGCTGCAG"

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT",
         "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
STAGES = ["designed", "assembled", "sequence_confirmed", "delivered", "edit_confirmed", "phenotype_measured"]
REGISTRY_UA = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) "
               "Chrome/124.0 Safari/537.36")


# ------------------------------------------------------------------ sequences
def clean_sequence(text: str) -> str:
    """Upper-case DNA without whitespace or digits; raises ValueError on non-ACGT characters."""
    seq = re.sub(r"[\s\d]", "", text or "").upper()
    bad = sorted(set(seq) - set("ACGT"))
    if bad:
        raise ValueError(f"Sequence contains characters other than A, C, G, T: {''.join(bad)}")
    return seq


def revcomp(seq: str) -> str:
    return seq.upper().translate(COMPLEMENT)[::-1]


def gc_percent(seq: str) -> float:
    return round(100 * (seq.count("G") + seq.count("C")) / len(seq), 1) if seq else 0.0


def longest_homopolymer(seq: str) -> int:
    return max((len(m.group(0)) for m in re.finditer(r"(.)\1*", seq)), default=0)


def find_motif(seq: str, motif: str) -> list[dict]:
    """Every occurrence of an IUPAC motif on both strands (1-based start on the forward strand)."""
    motif = motif.upper().strip()
    if not motif or set(motif) - set(IUPAC):
        raise ValueError(f"Invalid IUPAC motif: '{motif}'")
    pattern = "".join(f"[{IUPAC[c]}]" if len(IUPAC[c]) > 1 else IUPAC[c] for c in motif)
    hits = [{"start": m.start() + 1, "strand": "+"} for m in re.finditer(f"(?=({pattern}))", seq)]
    rc = revcomp(motif)
    if rc != motif:
        rc_pattern = "".join(f"[{IUPAC[c]}]" if len(IUPAC[c]) > 1 else IUPAC[c] for c in rc)
        hits += [{"start": m.start() + 1, "strand": "-"} for m in re.finditer(f"(?=({rc_pattern}))", seq)]
    return hits


def repeated_kmers(seq: str, k: int = 24) -> list[tuple[int, int]]:
    """Pairs of positions (1-based) where a k-mer is repeated: a recombination risk."""
    seen: dict[str, int] = {}
    pairs = []
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        if kmer in seen and i - seen[kmer] >= k:
            pairs.append((seen[kmer] + 1, i + 1))
            if len(pairs) >= 5:
                break
        else:
            seen.setdefault(kmer, i)
    return pairs


# ------------------------------------------------------------------ parts and checks
def norm_taxon(name: str) -> str:
    """'Bacillus subtilis' and 'B. subtilis' both become 'b. subtilis'."""
    words = re.sub(r"[^A-Za-z. ]", "", name or "").replace(".", ". ").split()
    if len(words) >= 2:
        return f"{words[0][0].lower()}. {words[1].lower()}"
    return (name or "").strip().lower()


def check_sequence(name: str, seq: str, standards: list[str], chassis_motifs: dict[str, str] | None = None,
                   chassis: str = "", hosts: str = "", is_part: bool = True) -> list[dict]:
    """Screen one sequence. Returns findings: {level: ok|warning|error, check, detail}."""
    findings: list[dict] = []

    def add(level, check, detail):
        findings.append({"item": name, "level": level, "check": check, "detail": detail})

    for standard in standards:
        sites = STANDARDS[standard]
        hits = [(enzyme, h["start"]) for enzyme, motif in sites.items() for h in find_motif(seq, motif)]
        if hits:
            add("error", standard, "internal sites: " + ", ".join(f"{e} at {p}" for e, p in hits[:8]))
        else:
            add("ok", standard, "no forbidden internal restriction sites")
    for label, motif in (chassis_motifs or {}).items():
        hits = find_motif(seq, motif)
        if hits:
            add("warning", "chassis restriction motif", f"{label} ({motif}) occurs {len(hits)}x, e.g. at {hits[0]['start']}: "
                "the construct may be cut after delivery unless the sites are protected, removed or methylated")
    homopolymer = longest_homopolymer(seq)
    if homopolymer >= 8:
        add("warning", "homopolymer", f"run of {homopolymer} identical bases (synthesis and sequencing risk)")
    gc = gc_percent(seq)
    if len(seq) >= 100 and (gc < 25 or gc > 75):
        add("warning", "GC content", f"{gc}% is extreme")
    repeats = repeated_kmers(seq)
    if repeats:
        add("warning", "repeats", f"a 24-bp sequence repeats (positions {repeats[0][0]} and {repeats[0][1]}): recombination risk")
    if is_part and chassis:
        declared = [norm_taxon(h) for h in re.split(r"[;,]", hosts or "") if h.strip()]
        if not declared:
            add("warning", "host range", "no host range declared for this part; activity in the chassis is unknown")
        elif norm_taxon(chassis) not in declared:
            add("warning", "host range", f"declared hosts ({hosts}) do not include {chassis}; test activity there")
        else:
            add("ok", "host range", f"declared for {chassis}")
    return findings


def parse_registry_xml(xml: str) -> dict:
    """Extract the fields we need from an iGEM Registry part XML (plain-text extraction, no XML parser)."""
    def tag(name):
        m = re.search(rf"<{name}>(.*?)</{name}>", xml, re.S)
        return m.group(1).strip() if m else ""
    seq = re.sub(r"[^acgtACGT]", "", tag("seq_data")).upper()
    if not seq or not tag("part_name"):
        raise ValueError("The registry response has no part name or sequence.")
    rtype = tag("part_type")
    return {"part_id": tag("part_name"), "name": tag("part_short_desc") or tag("part_name"),
            "role": REGISTRY_ROLE.get(rtype.lower(), "other"), "sequence": seq, "hosts": "",
            "source": f"iGEM Registry ({rtype}); {tag('release_status')}; results: {tag('part_results') or 'n/a'}; "
                      f"author: {tag('part_author') or 'n/a'}",
            "url": f"https://parts.igem.org/Part:{tag('part_name')}"}


def fetch_registry_part(part_id: str, opener: Callable | None = None, timeout: int = 25) -> dict:
    """Download one part by ID (for example BBa_J23100). `opener` is injectable for tests."""
    part_id = part_id.strip()
    if not re.fullmatch(r"BBa_[A-Za-z0-9]{2,20}", part_id):
        raise ValueError(f"'{part_id}' is not a valid BioBrick ID (expected something like BBa_J23100).")
    url = "https://parts.igem.org/cgi/xml/part.cgi?part=" + urllib.parse.quote(part_id)
    request = urllib.request.Request(url, headers={"User-Agent": REGISTRY_UA})
    try:
        with (opener or urllib.request.urlopen)(request, timeout=timeout) as response:
            xml = response.read().decode("utf-8", "replace")
    except Exception as error:  # offline, blocked or server error
        raise RuntimeError(f"Could not reach the iGEM Registry for {part_id}: {error}") from error
    part = parse_registry_xml(xml)
    part["fetched"] = _dt.date.today().isoformat()
    return part


# ------------------------------------------------------------------ assembly
def assemble(parts: list[dict], scar: str = "TACTAGAG", add_biobrick_ends: bool = False,
             junction_scars: list[str] | None = None) -> dict:
    """Join parts in order with scars. Returns {sequence, annotations}; positions are 1-based.

    `junction_scars` (one entry per junction) overrides `scar` where a junction must be exact,
    for example between a guide spacer and its scaffold.
    """
    if not parts:
        raise ValueError("Choose at least one part.")
    if junction_scars is not None and len(junction_scars) != len(parts) - 1:
        raise ValueError("junction_scars needs one entry per junction.")
    seq, notes = "", []
    if add_biobrick_ends:
        seq = PREFIX
        notes.append({"label": "BioBrick prefix", "role": "other", "start": 1, "end": len(PREFIX)})
    for i, part in enumerate(parts):
        joint = (junction_scars[i - 1] if junction_scars is not None else scar) if i else ""
        if joint:
            notes.append({"label": f"scar {joint}", "role": "other", "start": len(seq) + 1, "end": len(seq) + len(joint)})
            seq += joint
        s = clean_sequence(part["sequence"])
        notes.append({"label": part.get("name") or part["part_id"], "part_id": part["part_id"],
                      "role": part.get("role", "other"), "start": len(seq) + 1, "end": len(seq) + len(s)})
        seq += s
    if add_biobrick_ends:
        notes.append({"label": "BioBrick suffix", "role": "other", "start": len(seq) + 1, "end": len(seq) + len(SUFFIX)})
        seq += SUFFIX
    return {"sequence": seq, "annotations": notes}


def homology_arms(records: list[tuple[str, str]], contig: str, position: int, arm_length: int) -> dict:
    """Left and right arms flanking an insertion point placed after base `position` (1-based)."""
    lookup = dict(records)
    if contig not in lookup:
        raise ValueError(f"Contig '{contig}' is not in the genome.")
    genome = lookup[contig].upper()
    if not 1 <= position < len(genome):
        raise ValueError(f"Position must be between 1 and {len(genome) - 1} for contig '{contig}'.")
    if position < arm_length or position + arm_length > len(genome):
        raise ValueError(f"The site is closer than {arm_length} bp to a contig end: use shorter arms.")
    return {"left": genome[position - arm_length:position], "right": genome[position:position + arm_length],
            "left_span": (position - arm_length + 1, position), "right_span": (position + 1, position + arm_length)}


def integration_construct(cassette: dict, arms: dict, orientation: str = "+") -> dict:
    """Left arm + cassette + right arm; the cassette may be inserted in reverse orientation."""
    seq = cassette["sequence"] if orientation == "+" else revcomp(cassette["sequence"])
    notes, offset = [], len(arms["left"])
    notes.append({"label": "left homology arm", "role": "other", "start": 1, "end": offset})
    length = len(cassette["sequence"])
    for n in cassette["annotations"]:
        start, end = (n["start"], n["end"]) if orientation == "+" else (length - n["end"] + 1, length - n["start"] + 1)
        notes.append({**n, "start": start + offset, "end": end + offset,
                      "strand": "+" if orientation == "+" else "-"})
    right_start = offset + len(seq) + 1
    notes.append({"label": "right homology arm", "role": "other", "start": right_start, "end": right_start + len(arms["right"]) - 1})
    notes.sort(key=lambda n: n["start"])
    return {"sequence": arms["left"] + seq + arms["right"], "annotations": notes}


def diagnostic_pcr(arms: dict, cassette_length: int, outside: int = 150) -> dict:
    """Expected product sizes for primers placed `outside` bp beyond each arm (wild type vs edited)."""
    wild = len(arms["left"]) + len(arms["right"]) + 2 * outside
    return {"outside_bp": outside, "wild_type_bp": wild, "edited_bp": wild + cassette_length}


def sgrna_cassette(promoter: dict, spacer: str, scaffold: dict, terminator: dict | None, scar: str = "") -> dict:
    """Promoter + 20-nt spacer + sgRNA scaffold (+ terminator).

    The spacer is the protospacer without the PAM. Promoter, spacer and scaffold are always joined
    without a scar (a scar inside the guide RNA would break it); `scar` applies only before the terminator.
    """
    spacer = clean_sequence(spacer)
    if len(spacer) != 20:
        raise ValueError(f"The spacer must be 20 nt; got {len(spacer)}.")
    spacer_part = {"part_id": "spacer", "name": "CRISPRi spacer", "role": "scaffold", "sequence": spacer}
    order = [promoter, spacer_part, scaffold] + ([terminator] if terminator else [])
    return assemble(order, junction_scars=["", ""] + ([scar] if terminator else []))


def to_fasta(name: str, seq: str, width: int = 70) -> str:
    return f">{name} {len(seq)} bp\n" + "\n".join(seq[i:i + width] for i in range(0, len(seq), width)) + "\n"


def to_genbank(name: str, seq: str, annotations: list[dict], definition: str = "designed construct") -> str:
    """Minimal GenBank flat file so the construct opens in common sequence editors."""
    key = {"promoter": "promoter", "rbs": "RBS", "cds": "CDS", "reporter": "CDS", "terminator": "terminator"}
    date = _dt.date.today().strftime("%d-%b-%Y").upper()
    lines = [f"LOCUS       {name[:16]:<16} {len(seq)} bp    DNA     linear   SYN {date}",
             f"DEFINITION  {definition}.", "FEATURES             Location/Qualifiers",
             f"     source          1..{len(seq)}", '                     /organism="synthetic construct"']
    for n in annotations:
        loc = f"{n['start']}..{n['end']}"
        if n.get("strand") == "-":
            loc = f"complement({loc})"
        lines += [f"     {key.get(n.get('role', 'other'), 'misc_feature'):<15} {loc}", f'                     /label="{n["label"]}"']
    lines.append("ORIGIN")
    lower = seq.lower()
    for i in range(0, len(lower), 60):
        chunks = " ".join(lower[j:j + 10] for j in range(i, min(i + 60, len(lower)), 10))
        lines.append(f"{i + 1:>9} {chunks}")
    lines.append("//")
    return "\n".join(lines) + "\n"


# ------------------------------------------------------------------ tracker
TRACKER_COLUMNS = (["construct_id", "group", "role_in_experiment", "edit_type", "target", "parts", "length_bp",
                    "predicted_score"] + [f"stage_{s}" for s in STAGES] +
                   ["assay", "replicates", "measured_value", "units", "normalized_effect", "date", "notes"])


def tracker_rows(constructs: list[dict], group: str) -> list[dict]:
    """Rows for the built constructs plus the controls every experiment needs."""
    rows = []
    for c in constructs:
        rows.append({**{col: "" for col in TRACKER_COLUMNS}, "construct_id": c["construct_id"], "group": group,
                     "role_in_experiment": "test", "edit_type": c["edit_type"], "target": c["target"],
                     "parts": c["parts"], "length_bp": c["length_bp"], "predicted_score": c.get("predicted_score", ""),
                     "stage_designed": _dt.date.today().isoformat()})
    controls = [("parental (no construct)", "parental"), ("empty cassette / no-guide control", "negative")]
    if any(c["edit_type"] == "CRISPRi" for c in constructs):
        controls.append(("low-ranked guide", "low-ranked"))
    for label, role in controls:
        rows.append({**{col: "" for col in TRACKER_COLUMNS}, "construct_id": f"CONTROL: {label}", "group": group,
                     "role_in_experiment": role, "edit_type": constructs[0]["edit_type"] if constructs else "",
                     "target": "", "parts": "", "length_bp": ""})
    return rows


def analyze_tracker(rows: list[dict]) -> dict:
    """Funnel counts, control completeness and predicted-vs-measured agreement."""
    tests = [r for r in rows if str(r.get("role_in_experiment", "")) == "test"]
    funnel = {s: sum(1 for r in tests if str(r.get(f"stage_{s}", "")).strip() not in ("", "nan", "None")) for s in STAGES}
    warnings = []
    groups = sorted({str(r.get("group", "")) for r in rows if r.get("group")})
    for g in groups:
        members = [r for r in rows if str(r.get("group", "")) == g]
        measured = lambda r: str(r.get("measured_value", "")).strip() not in ("", "nan", "None")
        have = {str(r.get("role_in_experiment")) for r in members if measured(r)}
        if any(measured(r) for r in members if r.get("role_in_experiment") == "test"):
            for needed in ("parental", "negative"):
                if needed not in have:
                    warnings.append(f"Group '{g}': the {needed} control has no measured value, so effects cannot be normalised.")
        for r in members:
            if measured(r):
                try:
                    if float(r.get("replicates") or 0) < 3:
                        warnings.append(f"{r['construct_id']}: fewer than 3 replicates recorded.")
                except (TypeError, ValueError):
                    warnings.append(f"{r['construct_id']}: replicate count is not a number.")
    pairs = []
    for r in tests:
        try:
            pairs.append((float(r["predicted_score"]), float(r["normalized_effect"])))
        except (TypeError, ValueError, KeyError):
            continue
    rho = None
    if len(pairs) >= 5:
        from .validation import _spearman_values  # local import keeps this module light
        import pandas as pd
        rho = float(_spearman_values(pd.Series([p[0] for p in pairs]), pd.Series([p[1] for p in pairs])))
    return {"funnel": funnel, "warnings": warnings, "n_paired": len(pairs), "spearman_rho": rho,
            "counts": Counter(str(r.get("role_in_experiment")) for r in rows)}
