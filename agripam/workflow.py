"""Reproducible, tool-aware genome annotation workflow for RhizoForge-Select.

Every result records whether it came from a specialist external tool, a transparent
screen, or was unavailable.  The workflow never converts a missing database into a
negative biological result.
"""
from __future__ import annotations

import csv
import json
import os
import re
import shutil
import subprocess
import tempfile
import time
import zipfile
from dataclasses import asdict, dataclass
from bisect import bisect_right
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

from .core import (exact_target_copy_counts, genome_qc, genome_summary, parse_fasta,
                   scan_pam_targets_with_ambiguity, scan_ttc_targets_with_ambiguity)
from .pam_discovery import extract_cctyper_spacers, infer_pams_from_mobile_references, supported_pams
from .multiobjective import experimental_results_template, rank_editing_candidates
from .system_assessment import assess_editing_systems, assess_introduced_editors
from .editor_targeting import EDITOR_PRESETS, approximate_offtarget_counts, scan_editor_targets


FUNCTIONAL_TERMS = {
    "antifungal_compounds": ("antifungal", "iturin", "fengycin", "bacillomycin", "surfactin", "mycosubtilin"),
    "lipopeptide_synthesis": ("non-ribosomal peptide", "nonribosomal peptide", "nrps", "lipopeptide", "srf", "fen", "itu"),
    "secretion": ("secretion", "secretory", "type vii", "type iv", "sec transloc", "tat ", "signal peptid"),
    "hydrolytic_enzymes": ("chitinase", "glucanase", "cellulase", "xylanase", "protease", "lipase", "hydrolase"),
    "phosphate_solubilization": ("phosphatase", "phytase", "gluconate dehydrogenase", "phosphate transport", "pqq"),
    "nitrogen_metabolism": ("nitrogenase", "nitrate reductase", "nitrite reductase", "glutamine synthetase", "nif", "nar", "nir"),
    "siderophore_production": ("siderophore", "bacillibactin", "petrobactin", "iron chelat", "dhb"),
    "root_colonization": ("biofilm", "chemotaxis", "flagellar", "swarming", "adhesin", "exopolysaccharide", "tasA"),
}
RM_TERMS = ("restriction", "methyltransferase", "dna methylase", "restriction-modification")
MOBILE_TERMS = ("transposase", "integrase", "recombinase", "prophage", "phage", "insertion sequence", "conjug")
GROWTH_TERMS = (
    "ribosomal protein", "dna polymerase", "dna gyrase", "rna polymerase",
    "cell division", "chromosome partition", "peptidoglycan", "cell wall",
    "replication initiator", "translation initiation", "elongation factor",
    "atp synthase", "tricarboxylic", "glycolysis", "amino acid biosynthesis",
)
GROWTH_GENE_PREFIXES = ("dnaA", "dnaN", "gyrA", "gyrB", "rpo", "rpl", "rps", "fts", "mur", "atp")
EDITING_SYSTEM_TERMS = {
    "recombinase_integrase": ("recombinase", "integrase", "xerC", "xerD", "resolvase"),
    "retron_reverse_transcriptase": ("retron", "reverse transcriptase", "msr-msd"),
    "crispr_associated_transposase": ("tnsA", "tnsB", "tnsC", "tnsD", "cas12k", "crispr-associated transpos"),
    "homology_directed_repair": ("recA", "recF", "recO", "recR", "ruvA", "ruvB", "ruvC"),
    "rna_guided_nuclease": ("cas9", "cas12", "cas13", "cas3", "cas10", "cpf1", "rna-guided"),
}


@dataclass
class Stage:
    name: str
    status: str
    method: str
    output: str = ""
    message: str = ""


def _which(root: Path, *names: str) -> str | None:
    candidates = []
    for name in names:
        found = shutil.which(name)
        if found:
            candidates.append(found)
        candidates.extend(str(p) for p in root.glob(f".envs/*/bin/{name}") if p.is_file())
    return candidates[0] if candidates else None


def tool_inventory(root: Path) -> list[dict[str, str]]:
    tools = {
        "Prodigal": ("prodigal",), "Bakta": ("bakta",), "Prokka": ("prokka",),
        "CRISPRCasTyper": ("cctyper",), "MinCED": ("minced",),
        "DefenseFinder": ("defense-finder", "defensefinder"),
        "MOB-suite": ("mob_recon",), "geNomad": ("genomad",),
        "SPAdes": ("spades.py",), "Unicycler": ("unicycler",), "Flye": ("flye",),
    }
    return [{"component": label, "available": "yes" if (p := _which(root, *names)) else "no", "path": p or ""}
            for label, names in tools.items()]


def _run(
    command: list[str],
    log: Path,
    timeout: int = 3600,
) -> subprocess.CompletedProcess[str]:
    """Run an external tool without allowing one module to block the workflow."""
    env = os.environ.copy()
    executable_dir = str(Path(command[0]).resolve().parent)
    env["PATH"] = executable_dir + os.pathsep + env.get("PATH", "")

    try:
        result = subprocess.run(
            command,
            text=True,
            capture_output=True,
            timeout=timeout,
            check=False,
            env=env,
        )
    except subprocess.TimeoutExpired as error:
        stdout = error.stdout or ""
        stderr = error.stderr or ""

        if isinstance(stdout, bytes):
            stdout = stdout.decode(errors="replace")
        if isinstance(stderr, bytes):
            stderr = stderr.decode(errors="replace")

        message = (
            f"External tool exceeded the {timeout}-second timeout. "
            "RhizoForge continued with the remaining modules."
        )

        log.write_text(
            "COMMAND\n"
            + " ".join(command)
            + "\n\nSTATUS\nTIMEOUT\n\nSTDOUT\n"
            + stdout
            + "\nSTDERR\n"
            + stderr
            + "\n\nRHIZOFORGE\n"
            + message
        )

        return subprocess.CompletedProcess(
            args=command,
            returncode=124,
            stdout=stdout,
            stderr=(stderr + "\n" + message).strip(),
        )
    except OSError as error:
        message = (
            f"External tool could not be executed: {error}. "
            "RhizoForge continued with the remaining modules."
        )

        log.write_text(
            "COMMAND\n"
            + " ".join(command)
            + "\n\nSTATUS\nEXECUTION_ERROR\n\nRHIZOFORGE\n"
            + message
        )

        return subprocess.CompletedProcess(
            args=command,
            returncode=127,
            stdout="",
            stderr=message,
        )

    log.write_text(
        "COMMAND\n"
        + " ".join(command)
        + "\n\nSTDOUT\n"
        + result.stdout
        + "\nSTDERR\n"
        + result.stderr
    )
    return result


def _write_fasta(text: str, path: Path) -> list[tuple[str, str]]:
    parsed_records = parse_fasta(text)
    records = [(name.split()[0], seq) for name, seq in parsed_records]
    with path.open("w") as handle:
        for name, seq in records:
            handle.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i + 80] + "\n")
    return records


def read_gff(path: Path) -> list[dict[str, object]]:
    rows = []
    if not path.exists():
        return rows
    with path.open(errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError:
                continue
            attrs = {k: v for k, _, v in (item.partition("=") for item in fields[8].split(";")) if k}
            rows.append({"contig": fields[0], "source": fields[1], "type": fields[2], "start": start,
                         "end": end, "strand": fields[6], "attributes": fields[8],
                         "id": attrs.get("ID", attrs.get("locus_tag", "")),
                         "gene": attrs.get("gene", attrs.get("Name", "")),
                         "product": attrs.get("product", attrs.get("Name", "")).replace("%20", " ")})
    return rows



def _build_cds_index(features: list[dict[str, object]]) -> dict[str, tuple[list[int], list[dict[str, object]]]]:
    """Index CDS features by contig for fast coordinate-overlap lookup."""
    by_contig: dict[str, list[dict[str, object]]] = {}
    for feature in features:
        if str(feature.get("type", "")).upper() == "CDS":
            by_contig.setdefault(str(feature["contig"]), []).append(feature)

    index = {}
    for contig, rows in by_contig.items():
        rows.sort(key=lambda row: int(row["start"]))
        index[contig] = ([int(row["start"]) for row in rows], rows)
    return index


def _find_overlapping_cds(
    cds_index: dict[str, tuple[list[int], list[dict[str, object]]]],
    contig: str,
    start: int,
    end: int,
) -> dict[str, object] | None:
    """Return the first CDS overlapping a genomic interval, or None."""
    indexed = cds_index.get(str(contig))
    if not indexed:
        return None

    starts, rows = indexed
    pos = bisect_right(starts, int(end)) - 1

    while pos >= 0:
        feature = rows[pos]
        if int(feature["end"]) >= int(start):
            return feature
        break

    return None


def _write_tsv(rows: Iterable[dict], path: Path, fields: list[str] | None = None) -> None:
    rows = list(rows)
    fields = fields or (list(rows[0]) if rows else ["status"])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        if rows:
            writer.writerows(rows)


def _keyword_hits(features: list[dict[str, object]], terms: dict[str, tuple[str, ...]]) -> list[dict]:
    hits = []
    for feature in features:
        text = (str(feature.get("product", "")) + " " + str(feature.get("attributes", ""))).lower()
        for category, words in terms.items():
            matched = []
            for word in words:
                term = word.lower().strip()
                # Short gene symbols must match complete tokens: recO must not match recombinase.
                if len(term) <= 4 and term.isalnum():
                    found = re.search(rf"(?<![a-z0-9]){re.escape(term)}(?![a-z0-9])", text)
                else:
                    found = term in text
                if found:
                    matched.append(word)
            matched = sorted(set(matched))
            if matched:
                hits.append({"category": category, "contig": feature["contig"], "start": feature["start"],
                             "end": feature["end"], "strand": feature["strand"], "feature_id": feature["id"],
                             "product": feature["product"], "matched_terms": ", ".join(matched),
                             "evidence_level": "annotation-keyword screen"})
    return hits


def _interval_hits(features: list[dict[str, object]], terms: tuple[str, ...], category: str) -> list[dict]:
    return _keyword_hits(features, {category: terms})


def _classify_feature(feature: dict[str, object]) -> tuple[str, str]:
    """Return conservative growth and agronomic screening labels."""
    text = (str(feature.get("gene", "")) + " " + str(feature.get("product", ""))).lower()
    growth_matches = [term for term in GROWTH_TERMS if term in text]
    gene = str(feature.get("gene", ""))
    if gene and any(gene.lower().startswith(prefix.lower()) for prefix in GROWTH_GENE_PREFIXES):
        growth_matches.append(gene)
    agronomic = [category for category, terms in FUNCTIONAL_TERMS.items()
                 if any(term.lower() in text for term in terms)]
    return (", ".join(sorted(set(growth_matches))) or "not flagged",
            ", ".join(sorted(set(agronomic))) or "not flagged")


def _pam_gene_context(records: list[tuple[str, str]], features: list[dict[str, object]],
                      protospacer_length: int = 35, pams: list[str] | None = None) -> list[dict]:
    """Annotate targets for explicit PAMs with overlapping or nearest CDS context."""
    by_contig: dict[str, list[dict[str, object]]] = {}
    for feature in features:
        if str(feature.get("type", "")).upper() == "CDS":
            by_contig.setdefault(str(feature["contig"]), []).append(feature)
    for contig in by_contig:
        by_contig[contig].sort(key=lambda row: int(row["start"]))

    rows = []
    index = 0
    for contig, sequence in records:
        contig_features = by_contig.get(contig, [])
        motifs = pams or ["TTC"]
        for hit in scan_pam_targets_with_ambiguity(sequence, motifs, protospacer_length):
            index += 1
            target_start, target_end = int(hit["start_1based"]), int(hit["end_1based"])
            if hit["strand"] == "+":
                pam_start, pam_end = target_start - 3, target_start - 1
            else:
                pam_start, pam_end = target_end + 1, target_end + 3
            overlaps = [f for f in contig_features
                        if int(f["start"]) <= target_end and int(f["end"]) >= target_start]
            if overlaps:
                feature = max(overlaps, key=lambda f: min(target_end, int(f["end"])) - max(target_start, int(f["start"])) + 1)
                overlap_bp = min(target_end, int(feature["end"])) - max(target_start, int(feature["start"])) + 1
                relation, distance = "within CDS", 0
                consequence = "may disrupt coding sequence; consequence depends on edit type and repair outcome"
            elif contig_features:
                feature = min(contig_features, key=lambda f: min(abs(target_start - int(f["end"])), abs(int(f["start"]) - target_end)))
                distance = max(0, max(int(feature["start"]) - target_end, target_start - int(feature["end"])))
                relation = "promoter-proximal" if distance <= 250 else "intergenic"
                overlap_bp = 0
                consequence = ("may alter nearby gene regulation; promoter mapping is required" if distance <= 250
                               else "no direct CDS disruption predicted; distal regulatory effects remain possible")
            else:
                feature = {"id": "", "gene": "", "product": "", "start": "", "end": "", "strand": ""}
                relation, distance, overlap_bp = "unannotated", "", 0
                consequence = "gene consequence unavailable because no CDS annotation was available"
            growth, agronomic = _classify_feature(feature)
            rows.append({
                "candidate_id": f"RFS{index:07d}", "contig": contig, "pam": hit["pam"],
                "pam_start_1based": pam_start, "pam_end_1based": pam_end,
                "target_start_1based": target_start, "target_end_1based": target_end,
                "target_strand": hit["strand"], "protospacer": hit["protospacer"],
                "gene_relation": relation, "overlap_bp": overlap_bp, "distance_to_gene_bp": distance,
                "feature_id": feature.get("id", ""), "gene": feature.get("gene", ""),
                "product": feature.get("product", ""), "gene_start": feature.get("start", ""),
                "gene_end": feature.get("end", ""), "gene_strand": feature.get("strand", ""),
                "growth_relevance_screen": growth, "agronomic_category_screen": agronomic,
                "potential_consequence": consequence,
                "evidence_level": "coordinate overlap + annotation screen; experimental validation required",
            })
    return rows


def _neutral_windows(records: list[tuple[str, str]], features: list[dict[str, object]], excluded: list[dict],
                     window: int = 1000, flank: int = 500) -> list[dict]:
    blocked: dict[str, list[tuple[int, int, str]]] = {}
    for row in [*features, *excluded]:
        contig = str(row.get("contig", ""))
        if not contig:
            continue
        blocked.setdefault(contig, []).append((max(1, int(row["start"]) - flank), int(row["end"]) + flank,
                                               str(row.get("type", row.get("category", "feature")))))
    candidates = []
    for contig, seq in records:
        spans = sorted(blocked.get(contig, []))
        cursor = 1
        for start, end, _ in spans + [(len(seq) + 1, len(seq) + 1, "end")]:
            if start - cursor >= window:
                region_start, region_end = cursor, start - 1
                segment = seq[region_start - 1:region_end]
                gc = 100 * (segment.count("G") + segment.count("C")) / max(1, len(segment))
                candidates.append({"contig": contig, "start": region_start, "end": region_end,
                                   "length_bp": region_end - region_start + 1, "gc_percent": round(gc, 2),
                                   "status": "candidate only", "basis": f">={flank}-bp from annotated/excluded features"})
            cursor = max(cursor, end + 1)
    return sorted(candidates, key=lambda r: -int(r["length_bp"]))


def assemble_reads(read_paths: list[Path], technology: str, root: Path, out: Path) -> tuple[Path, Stage]:
    if technology == "Illumina paired-end":
        exe = _which(root, "spades.py")
        if not exe or len(read_paths) != 2:
            raise RuntimeError("Paired-end assembly requires SPAdes and exactly two FASTQ files. Install the full pipeline environment or upload an assembled FASTA.")
        command = [exe, "-1", str(read_paths[0]), "-2", str(read_paths[1]), "-o", str(out / "spades")]
        assembly = out / "spades" / "contigs.fasta"
    elif technology == "Long reads":
        exe = _which(root, "flye")
        if not exe or len(read_paths) != 1:
            raise RuntimeError("Long-read assembly requires Flye and one FASTQ file. Install the full pipeline environment or upload an assembled FASTA.")
        command = [exe, "--nano-raw", str(read_paths[0]), "--out-dir", str(out / "flye")]
        assembly = out / "flye" / "assembly.fasta"
    else:
        exe = _which(root, "unicycler")
        if not exe or len(read_paths) < 3:
            raise RuntimeError("Hybrid assembly requires Unicycler, paired short reads and long reads.")
        command = [exe, "-1", str(read_paths[0]), "-2", str(read_paths[1]), "-l", str(read_paths[2]), "-o", str(out / "unicycler")]
        assembly = out / "unicycler" / "assembly.fasta"
    result = _run(command, out / "assembly.log", timeout=14400)
    if result.returncode or not assembly.exists():
        raise RuntimeError("Read assembly failed; download the assembly log for details.")
    return assembly, Stage("read_assembly", "complete", Path(exe).name, str(assembly.relative_to(out)))


def run_genome_workflow(fasta_text: str, source: str, root: Path, include_heavy: bool = True,
                        mobile_reference_fasta: str | None = None) -> tuple[Path, dict]:
    """Run available modules in a temporary directory and return directory + manifest."""
    work = Path(tempfile.mkdtemp(prefix="rhizoforge_"))
    genome = work / "genome.fna"
    records = _write_fasta(fasta_text, genome)
    stages: list[Stage] = [Stage("input_validation", "complete", "RhizoForge FASTA validator", "genome.fna")]
    qc_rows = genome_qc(records)
    _write_tsv(qc_rows, work / "genome_qc.tsv")
    qc_review = [str(row["check"]) for row in qc_rows if row["status"] == "review"]
    stages.append(Stage("genome_qc", "review" if qc_review else "complete",
                        "taxon-agnostic FASTA metrics", "genome_qc.tsv",
                        ("Review: " + ", ".join(qc_review)) if qc_review else "No generic QC threshold was triggered"))

    # Annotation: prefer Bakta/Prokka; always retain Prodigal as a database-free fallback.
    gff, proteins = work / "annotation.gff", work / "proteins.faa"
    bakta, prokka, prodigal = _which(root, "bakta"), _which(root, "prokka"), _which(root, "prodigal")
    bakta_db = os.environ.get("BAKTA_DB")
    annotation_method = "none"
    annotation_message = ""
    if include_heavy and bakta and bakta_db and Path(bakta_db).is_dir():
        result = _run([bakta, "--db", bakta_db, "--output", str(work / "bakta"), "--prefix", "genome", str(genome)], work / "annotation.log")
        src_gff, src_faa = work / "bakta" / "genome.gff3", work / "bakta" / "genome.faa"
        if result.returncode == 0 and src_gff.exists():
            shutil.copy2(src_gff, gff); shutil.copy2(src_faa, proteins)
            annotation_method = "Bakta"
        else:
            annotation_message = "Bakta failed; "
    elif include_heavy and prokka:
        result = _run([prokka, "--outdir", str(work / "prokka"), "--prefix", "genome", str(genome)], work / "annotation.log")
        src_gff, src_faa = work / "prokka" / "genome.gff", work / "prokka" / "genome.faa"
        if result.returncode == 0 and src_gff.exists():
            shutil.copy2(src_gff, gff); shutil.copy2(src_faa, proteins)
            annotation_method = "Prokka"
        else:
            annotation_message = "Prokka failed; "

    # A failed database-backed annotator must not turn an ordinary genome into a false zero-gene result.
    if not gff.exists() and prodigal:
        prodigal_mode = "single" if sum(len(seq) for _, seq in records) >= 20_000 else "meta"
        result = _run([prodigal, "-i", str(genome), "-a", str(proteins), "-f", "gff", "-o", str(gff), "-p", prodigal_mode], work / "prodigal_fallback.log")
        if result.returncode == 0 and gff.exists():
            annotation_method = "Prodigal fallback"
            annotation_message += "CDS coordinates and proteins recovered; named products require Bakta/Prokka"
    annotation_status = "complete" if gff.exists() else ("failed" if annotation_method != "none" or annotation_message else "unavailable")
    if annotation_status == "unavailable":
        annotation_message = "Install Bakta, Prokka or Prodigal"
    stages.append(Stage("coding_sequence_annotation", annotation_status, annotation_method,
                        "annotation.gff" if gff.exists() else "", annotation_message or "Database-supported products"))

    features = read_gff(gff)
    functional = _keyword_hits(features, FUNCTIONAL_TERMS)
    rm_hits = _interval_hits(features, RM_TERMS, "restriction_modification")
    mobile_hits = _interval_hits(features, MOBILE_TERMS, "mobile_element")
    editing_hits = _keyword_hits(features, EDITING_SYSTEM_TERMS)
    introduced_timing = {}
    introduced_start = time.perf_counter()

    introduced_targets = []
    scan_start = time.perf_counter()

    # Editors with identical PAM geometry share the same genomic target space.
    # Scan each unique geometry once, then materialize editor-specific records
    # so downstream biological interpretation remains separate.
    target_space_cache = {}

    for editor, spec in EDITOR_PRESETS.items():
        geometry = (
            str(spec["pam_pattern"]),
            str(spec["pam_side"]),
            int(spec["protospacer_length"]),
        )

        if geometry not in target_space_cache:
            scanned = scan_editor_targets(
                records,
                editor,
                spec["pam_pattern"],
                spec["pam_side"],
                spec["protospacer_length"],
            )
            target_space_cache[geometry] = scanned
            introduced_targets.extend(scanned)
        else:
            application = spec.get("application", "custom editing hypothesis")
            for source in target_space_cache[geometry]:
                cloned = dict(source)
                cloned["editor"] = editor
                cloned["application"] = application
                introduced_targets.append(cloned)

    introduced_timing["pam_scan_seconds"] = time.perf_counter() - scan_start
    introduced_timing["unique_target_geometries"] = len(target_space_cache)

    copy_start = time.perf_counter()
    exact_copies = exact_target_copy_counts(
        records,
        (str(target["protospacer"]) for target in introduced_targets),
    )
    introduced_timing["exact_copy_seconds"] = time.perf_counter() - copy_start

    mapping_start = time.perf_counter()
    cds_index = _build_cds_index(features)
    for target in introduced_targets:
        feature = _find_overlapping_cds(
            cds_index,
            str(target["contig"]),
            int(target["start"]),
            int(target["end"]),
        )
        protospacer = str(target["protospacer"])
        gc = 100 * (protospacer.count("G") + protospacer.count("C")) / len(protospacer)
        copies = exact_copies.get(protospacer, 0)
        homopolymer = max((len(match.group()) for match in re.finditer(r"(A+|C+|G+|T+)", protospacer)), default=0)
        priority = max(0.0, 100 - abs(gc - 50) * 1.5 - max(0, homopolymer - 4) * 8 - max(0, copies - 1) * 25)
        target.update({"gene_relation": "within CDS" if feature else "intergenic/uncalled", "feature_id": feature.get("id", "") if feature else "", "product": feature.get("product", "") if feature else "", "possible_consequence": "candidate gene disruption or regulation" if feature else "candidate non-coding edit",
                       "gc_percent": round(gc, 2), "longest_homopolymer": homopolymer,
                       "exact_genome_copies": copies, "exact_uniqueness": "unique" if copies == 1 else "repeated",
                       "guide_priority_score": round(priority, 2),
                       "score_interpretation": "sequence/uniqueness heuristic; not an activity probability",
                       "off_target_warning": "approximate mismatch and bulge search not run",
                       "evidence_level": "editor-specific PAM scan; computational targetability only", "validation_required": "approximate off-target search, delivery, nuclease activity, repair and sequence/phenotype confirmation"})
    introduced_timing["annotation_mapping_seconds"] = time.perf_counter() - mapping_start

    introduced_targets.sort(
        key=lambda row: (
            -float(row["guide_priority_score"]),
            int(row["exact_genome_copies"]),
            str(row["contig"]),
            int(row["start"]),
        )
    )

    offtarget_start = time.perf_counter()
    approximate = approximate_offtarget_counts(introduced_targets)
    introduced_timing["offtarget_seconds"] = time.perf_counter() - offtarget_start
    for target in introduced_targets:
        result = approximate.get(f"{target['editor']}\t{target['protospacer']}")
        if result:
            target["off_target_screen_status"] = result["status"]
            target["approx_offtargets_le2"] = result["count"] if result["count"] is not None else ""
            target["nearest_approx_offtarget_mismatches"] = result["nearest_mismatches"] or ""
            if result["status"] == "complete_substitutions_only":
                penalty = min(40, 5 * int(result["count"] or 0))
                target["offtarget_adjusted_priority_score"] = round(
                    max(0, float(target["guide_priority_score"]) - penalty), 2
                )
                target["off_target_warning"] = (
                    "substitution screen complete; DNA/RNA bulges and noncanonical PAMs not evaluated"
                )
            else:
                target["offtarget_adjusted_priority_score"] = target["guide_priority_score"]
                target["off_target_warning"] = (
                    "approximate off-target screen not performed because the "
                    "PAM-compatible reference space exceeded the configured limit; "
                    "off-target risk remains unresolved"
                )
        else:
            target.update({"off_target_screen_status": "not_prioritized",
                           "approx_offtargets_le2": "", "nearest_approx_offtarget_mismatches": "",
                           "offtarget_adjusted_priority_score": target["guide_priority_score"],
                           "off_target_warning": "not among the top 2,000 unique guide queries; approximate mismatch/bulge search not run"})
    introduced_targets.sort(key=lambda row: (-float(row["offtarget_adjusted_priority_score"]),
                                               int(row["exact_genome_copies"]), str(row["contig"]), int(row["start"])))
    _write_tsv(introduced_targets, work / "introduced_editor_targets.tsv")

    introduced_timing["total_seconds"] = time.perf_counter() - introduced_start
    introduced_timing["total_candidates"] = len(introduced_targets)
    (work / "introduced_editor_timing.json").write_text(
        json.dumps(introduced_timing, indent=2)
    )

    stages.append(Stage("introduced_editor_targetability", "complete", "editor-specific IUPAC PAM scanning + coordinate/gene mapping + bounded mismatch screen", "introduced_editor_targets.tsv", "Up to 2,000 unique guides per editor are screened against PAM-compatible sites at 1–2 substitutions when the compatible reference space contains no more than 500,000 unique sequences; larger reference spaces are explicitly reported as not screened. Bulges and noncanonical PAMs remain unresolved."))
    editing_applications = {
        "recombinase_integrase": "site-specific integration or rearrangement",
        "retron_reverse_transcriptase": "template-generating editing hypothesis",
        "crispr_associated_transposase": "RNA-guided insertion hypothesis",
        "homology_directed_repair": "template-based repair compatibility",
        "rna_guided_nuclease": "alternative RNA-guided targeting",
    }
    for hit in editing_hits:
        hit.update({"possible_application": editing_applications.get(hit["category"], "editing-system hypothesis"),
                    "confidence": "exploratory annotation screen",
                    "limitations": "keyword/homology annotation alone does not establish a complete or active system",
                    "recommended_validation": "confirm gene identity, locus architecture and activity with an appropriate system-specific assay"})
    _write_tsv(functional, work / "agricultural_function_genes.tsv")
    _write_tsv(rm_hits, work / "restriction_modification_systems.tsv")
    _write_tsv(mobile_hits, work / "mobile_elements.tsv")
    _write_tsv(editing_hits, work / "editing_systems_screen.tsv")
    annotation_has_products = any(row.get("product") for row in features)
    method_note = "annotation-keyword screen" if annotation_has_products else "not interpretable: CDS products unavailable"
    stages.extend([Stage("agricultural_function_screen", "complete" if annotation_has_products else "limited", method_note, "agricultural_function_genes.tsv"),
                   Stage("restriction_modification", "complete" if annotation_has_products else "limited", method_note, "restriction_modification_systems.tsv"),
                   Stage("mobile_elements", "complete" if annotation_has_products else "limited", method_note, "mobile_elements.tsv"),
                   Stage("alternative_editing_systems", "complete" if annotation_has_products else "limited", method_note, "editing_systems_screen.tsv", "Hits are candidates, not active editing systems")])

    pam_context = _pam_gene_context(records, features)
    _write_tsv(pam_context, work / "pam_gene_context.tsv")
    stages.append(Stage("pam_gene_context", "complete" if features else "limited",
                        "coordinate overlap and nearest-CDS annotation", "pam_gene_context.tsv",
                        "Growth and agronomic labels are screening flags, not essentiality or phenotype proof"))

    # Plasmids: specialist call if MOB-suite exists; otherwise a clearly labelled contig screen.
    mob = _which(root, "mob_recon")
    plasmids = []
    if include_heavy and mob:
        result = _run([mob, "--infile", str(genome), "--outdir", str(work / "mob_suite")], work / "plasmids.log")
        report = work / "mob_suite" / "contig_report.txt"
        if report.exists(): shutil.copy2(report, work / "plasmid_predictions.tsv")
        stages.append(Stage("plasmid_prediction", "complete" if result.returncode == 0 else "failed", "MOB-suite", "plasmid_predictions.tsv"))
    else:
        for name, seq in records:
            header_signal = bool(re.search(r"plasmid|circular", name, re.I))
            size_signal = len(seq) < 350_000
            plasmids.append({"contig": name, "length_bp": len(seq), "header_signal": int(header_signal),
                             "small_contig_signal": int(size_signal), "classification": "candidate" if header_signal else "unresolved",
                             "evidence_level": "heuristic; not a MOB-suite call"})
        _write_tsv(plasmids, work / "plasmid_predictions.tsv")
        stages.append(Stage("plasmid_prediction", "limited", "transparent contig screen", "plasmid_predictions.tsv", "Install MOB-suite for database-supported classification"))

    # CRISPR-Cas and defense systems.
    cctyper = _which(root, "cctyper")
    if include_heavy and cctyper:
        cctyper_command = [cctyper, str(genome), str(work / "cctyper"), "--prodigal", "meta"]
        cctyper_db = root / "data" / "databases" / "cctyper_1.8.0"
        if cctyper_db.is_dir():
            cctyper_command.extend(["--db", str(cctyper_db)])
        result = _run(cctyper_command, work / "crispr_cas.log", timeout=300)
        stages.append(Stage("crispr_cas_annotation", "complete" if result.returncode == 0 else ("timeout" if result.returncode == 124 else "failed"), "CRISPRCasTyper", "cctyper/", "Completed successfully" if result.returncode == 0 else ("Specialist analysis exceeded the configured time limit; biological status remains unresolved" if result.returncode == 124 else "Specialist analysis failed; inspect crispr_cas.log")))
    else:
        stages.append(Stage("crispr_cas_annotation", "unavailable", "none", message="CRISPRCasTyper not installed or specialist analyses disabled"))

    # Strain-specific inference requires both native spacers and external mobile references.
    pam_ranking: list[dict] = []
    cctyper_output = work / "cctyper"
    if mobile_reference_fasta and cctyper_output.is_dir():
        arrays, spacers = extract_cctyper_spacers(cctyper_output)
        _write_tsv(arrays, work / "crispr_arrays.tsv")
        _write_tsv(spacers, work / "crispr_spacers.tsv")
        pam_summary, pam_observations, pam_ranking, pam_logo = infer_pams_from_mobile_references(
            spacers, mobile_reference_fasta
        )
        _write_tsv(pam_observations, work / "pam_protospacer_observations.tsv")
        _write_tsv(pam_ranking, work / "pam_64_triplet_ranking.tsv")
        _write_tsv(pam_logo, work / "pam_logo_frequencies.tsv")
        (work / "pam_discovery_summary.json").write_text(json.dumps(pam_summary, indent=2))
        discovered_motifs = supported_pams(pam_ranking)
        # When discovery was requested, never silently retain the legacy TTC assumption.
        pam_context = _pam_gene_context(records, features, pams=discovered_motifs) if discovered_motifs else []
        _write_tsv(pam_context, work / "pam_gene_context.tsv")
        context_stage = next(stage for stage in stages if stage.name == "pam_gene_context")
        context_stage.status = "complete" if discovered_motifs and features else "insufficient_evidence"
        context_stage.method = "strain-specific supported PAMs + coordinate overlap"
        context_stage.message = (f"Mapped supported motifs: {', '.join(discovered_motifs)}" if discovered_motifs
                                 else "No motif passed the prespecified support threshold; TTC was not substituted")
        pam_status = "complete" if pam_observations else "insufficient_evidence"
        stages.append(Stage("strain_specific_pam_discovery", pam_status,
                            "exact spacer–protospacer matching + background enrichment",
                            "pam_64_triplet_ranking.tsv",
                            "Experimental PAM-library validation remains required"))
    elif mobile_reference_fasta:
        stages.append(Stage("strain_specific_pam_discovery", "unavailable", "none",
                            message="CRISPRCasTyper output was unavailable; no PAM was inferred"))
    else:
        stages.append(Stage("strain_specific_pam_discovery", "not_run", "none",
                            message="Provide a versioned phage/plasmid/mobile-element FASTA collection"))
    defense = _which(root, "defense-finder", "defensefinder")
    if include_heavy and defense and proteins.exists():
        result = _run([defense, "run", str(proteins), "-o", str(work / "defensefinder")], work / "defensefinder.log", timeout=300)
        stages.append(Stage("defensefinder", "complete" if result.returncode == 0 else ("timeout" if result.returncode == 124 else "failed"), "DefenseFinder", "defensefinder/", "Completed successfully" if result.returncode == 0 else ("Specialist analysis exceeded the configured time limit; defense status remains unresolved" if result.returncode == 124 else "Specialist analysis failed; inspect defensefinder.log")))
    else:
        stages.append(Stage("defensefinder", "unavailable", "none", message="Install DefenseFinder and its models; absence is not a biological negative"))

    crispr_status = next((stage.status for stage in stages if stage.name == "crispr_cas_annotation"), "unavailable")
    defense_status = next((stage.status for stage in stages if stage.name == "defensefinder"), "unavailable")
    system_assessments = assess_editing_systems(
        features, {name: len(sequence) for name, sequence in records}, len(rm_hits),
        defense_status, pam_ranking, crispr_status,
    )
    _write_tsv(system_assessments, work / "editing_system_assessments.tsv")
    stages.append(Stage(
        "editing_system_assessment", "complete" if features else "unavailable",
        "explicit component-completeness and locus-organization models",
        "editing_system_assessments.tsv",
        "Expression and activity remain not evaluated without omics or experimental evidence",
    ))

    repair_supported = any(hit.get("category") == "homology_directed_repair" for hit in editing_hits)
    editor_compatibility = assess_introduced_editors(
        introduced_targets, len(rm_hits), defense_status, repair_supported, bool(features)
    )
    _write_tsv(editor_compatibility, work / "effector_chassis_compatibility.tsv")
    stages.append(Stage("effector_chassis_compatibility", "complete", "explicit targeting/delivery/repair evidence matrix",
                        "effector_chassis_compatibility.tsv",
                        "Compatibility tiers are decision support, not measured editing efficiency"))

    neutral = _neutral_windows(records, features, [*rm_hits, *mobile_hits, *functional]) if features else []
    _write_tsv(neutral, work / "candidate_neutral_regions.tsv")
    stages.append(Stage("candidate_neutral_regions", "complete" if features else "unavailable", "exclusion-window screen", "candidate_neutral_regions.tsv", "Candidates require synteny, essentiality and experimental validation"))

    portfolio = rank_editing_candidates(pam_context, pam_ranking, len(rm_hits), defense_status, repair_supported) if pam_ranking else []
    feedback_template = experimental_results_template(portfolio)
    _write_tsv(portfolio, work / "multiobjective_design_portfolio.tsv")
    _write_tsv(feedback_template, work / "experimental_results_template.tsv")
    stages.append(Stage("multiobjective_portfolio", "complete" if portfolio else "insufficient_evidence",
                        "transparent four-objective Pareto ranking", "multiobjective_design_portfolio.tsv",
                        "Requires supported strain-specific PAM evidence; scores are not editing probabilities"))

    summary = genome_summary(records)
    summary.update({"coding_sequences": sum(1 for r in features if str(r["type"]).upper() == "CDS"),
                    "functional_gene_hits": len(functional), "rm_hits": len(rm_hits),
                    "mobile_feature_hits": len(mobile_hits), "neutral_region_candidates": len(neutral),
                    "pam_targets_annotated": len(pam_context),
                    "pam_targets_within_cds": sum(r["gene_relation"] == "within CDS" for r in pam_context),
                    "pam_targets_growth_flagged": sum(r["growth_relevance_screen"] != "not flagged" for r in pam_context),
                    "pam_targets_agronomic_flagged": sum(r["agronomic_category_screen"] != "not flagged" for r in pam_context)})
    summary.update({"editing_system_candidates": len(editing_hits), "portfolio_candidates": len(portfolio),
                    "pareto_candidates": sum(r.get("pareto_optimal") == "yes" for r in portfolio)})
    confidence_rows = []
    specialist_methods = {"Bakta", "Prokka", "CRISPRCasTyper", "DefenseFinder", "MOB-suite"}
    for stage in stages:
        confidence = ("high" if stage.status == "complete" and stage.method in specialist_methods
                      else "moderate" if stage.status == "complete" else "limited")
        confidence_rows.append({"module": stage.name, "status": stage.status, "confidence": confidence,
                                "method": stage.method,
                                "uncertainty_or_next_action": stage.message or "Review output and validate experimentally."})
    _write_tsv(confidence_rows, work / "uncertainty_summary.tsv")
    manifest = {"product": "RhizoForge-Select", "engine": "AgriPAM-AI", "source": source,
                "run_timestamp_utc": datetime.now(timezone.utc).isoformat(), "summary": summary,
                "stages": [asdict(s) for s in stages], "interpretation_boundary":
                "Candidate neutral regions, editing-system hits and multi-objective scores are hypotheses, not validated functions, efficiencies or safe insertion sites.",
                "decision_model": {"objectives": ["editability", "deliverability", "agronomic value", "preservation and safety"],
                                   "ranking": "Pareto front plus transparent balanced score",
                                   "feedback_file": "experimental_results_template.tsv"},
                "quality_control": {"review_required": bool(qc_review), "triggered_checks": qc_review,
                                    "report": "genome_qc.tsv"},
                "uncertainty_report": "uncertainty_summary.tsv"}
    (work / "run_manifest.json").write_text(json.dumps(manifest, indent=2))
    _write_tsv([asdict(s) for s in stages], work / "pipeline_status.tsv")
    return work, manifest


def zip_results(directory: Path) -> bytes:
    archive = directory / "rhizoforge_select_results.zip"
    with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as bundle:
        for path in sorted(directory.rglob("*")):
            if path.is_file() and path != archive:
                bundle.write(path, path.relative_to(directory))
    return archive.read_bytes()
