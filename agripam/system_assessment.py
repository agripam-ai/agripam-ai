"""Conservative system-level assessment of genome-engineering machinery.

The output distinguishes genome-derived evidence from measurements that cannot be
obtained from an assembly (expression, activity and editing efficiency).
"""
from __future__ import annotations

import re

from .editor_targeting import EDITOR_PRESETS


SYSTEM_MODELS = {
    "CRISPR interference machinery": {
        "groups": {
            "effector": ("cas9", "cas12", "cpf1", "cas3", "cas10", "cas13"),
            "adaptation": ("cas1",),
            "adaptation_partner": ("cas2",),
        },
        "minimum": ("effector",),
        "application": "RNA-guided cleavage, interference or regulation",
        "validation": "confirm locus sequence and subtype; determine guide architecture; test candidate and PAM-mutant targets with no-guide controls",
    },
    "Recombinase / integrase route": {
        "groups": {"catalytic_component": ("recombinase", "integrase", "resolvase", "xerc", "xerd")},
        "minimum": ("catalytic_component",),
        "application": "site-specific integration, excision or rearrangement",
        "validation": "confirm catalytic gene and attachment-site architecture; verify predicted integration junctions against negative controls",
    },
    "Retron / reverse-transcriptase route": {
        "groups": {
            "reverse_transcriptase": ("retron", "reverse transcriptase"),
            "retron_ncrna_or_accessory": ("msr", "msd", "retron accessory", "effector"),
        },
        "minimum": ("reverse_transcriptase", "retron_ncrna_or_accessory"),
        "application": "intracellular donor-template generation or retron-derived engineering",
        "validation": "confirm RT–ncRNA locus boundaries and product formation before testing a sequence-confirmed editing readout",
    },
    "CRISPR-associated transposase route": {
        "groups": {
            "rna_guided_component": ("cas12k", "cas6", "cas7", "cas8"),
            "transposase_b": ("tnsb", "transposase b"),
            "transposition_atpase": ("tnsc", "transposition protein c"),
            "target_selector": ("tnsd", "tnse", "target selector"),
        },
        "minimum": ("rna_guided_component", "transposase_b", "transposition_atpase"),
        "application": "RNA-guided DNA insertion hypothesis",
        "validation": "confirm complete co-localized CAST architecture, then map insertion junctions and orientation with guide-free controls",
    },
    "Homology-directed repair route": {
        "groups": {
            "strand_exchange": ("reca", "recombinase a"),
            "presynaptic_repair": ("recf", "reco", "recr"),
            "junction_processing": ("ruva", "ruvb", "ruvc"),
        },
        "minimum": ("strand_exchange", "presynaptic_repair"),
        "application": "template-directed repair or recombineering support",
        "validation": "confirm repair-gene integrity and compare sequence-confirmed donor incorporation with donor-free and no-editor controls",
    },
    "Alternative RNA-guided nuclease route": {
        "groups": {"nuclease": ("cas9", "cas12", "cpf1", "cas13", "cas3", "cas10")},
        "minimum": ("nuclease",),
        "application": "non-Type-I-C RNA-guided targeting",
        "validation": "confirm nuclease identity and guide/PAM requirements, followed by a matched-target reporter assay",
    },
}


def _feature_text(feature: dict) -> str:
    return " ".join(str(feature.get(key, "")) for key in ("gene", "product", "attributes")).lower()


def _group_hits(features: list[dict], terms: tuple[str, ...]) -> list[dict]:
    patterns = tuple(term.lower() for term in terms)
    return [feature for feature in features if any(re.search(rf"\b{re.escape(term)}\b", _feature_text(feature))
                                                   for term in patterns)]


def assess_editing_systems(
    features: list[dict], contig_lengths: dict[str, int], rm_hit_count: int,
    defense_status: str, pam_ranking: list[dict], crispr_status: str,
) -> list[dict]:
    """Return one auditable assessment per editing-system model."""
    supported_pams = [str(row.get("pam3")) for row in pam_ranking
                      if int(row.get("observations", 0) or 0) >= 3 and float(row.get("bh_q_value", 1) or 1) < 0.10]
    rows = []
    for system, model in SYSTEM_MODELS.items():
        grouped = {name: _group_hits(features, terms) for name, terms in model["groups"].items()}
        detected_groups = [name for name, hits in grouped.items() if hits]
        required = list(model["minimum"])
        missing = [name for name in required if not grouped.get(name)]
        all_hits = [hit for hits in grouped.values() for hit in hits]
        unique_hits = {(str(hit.get("contig", "")), int(hit.get("start", 0)), int(hit.get("end", 0))): hit
                       for hit in all_hits}
        all_hits = list(unique_hits.values())

        by_contig: dict[str, list[dict]] = {}
        for hit in all_hits:
            by_contig.setdefault(str(hit.get("contig", "")), []).append(hit)
        best_locus = max(by_contig.values(), key=len, default=[])
        if len(best_locus) >= 2:
            locus_span = max(int(hit["end"]) for hit in best_locus) - min(int(hit["start"]) for hit in best_locus) + 1
            organization = "co-localized candidate locus" if locus_span <= 50_000 else "components dispersed over >50 kb"
        elif best_locus:
            locus_span, organization = 0, "single component; organization not assessable"
        else:
            locus_span, organization = 0, "no component detected"

        integrity_flags = []
        for hit in all_hits:
            text = _feature_text(hit)
            if any(flag in text for flag in ("pseudo", "partial", "truncated", "frameshift")):
                integrity_flags.append(f"{hit.get('id') or hit.get('gene') or 'feature'} annotated as partial/pseudogene")
            contig_length = contig_lengths.get(str(hit.get("contig", "")), 0)
            if contig_length and (int(hit.get("start", 0)) <= 30 or int(hit.get("end", 0)) >= contig_length - 30):
                integrity_flags.append(f"{hit.get('id') or hit.get('gene') or 'feature'} touches a contig boundary")

        if not all_hits:
            completeness = "not detected"
        elif not missing:
            completeness = "minimum component model satisfied"
        else:
            completeness = "partial candidate"
        if all_hits and not integrity_flags:
            integrity = "no annotation-level truncation flag; protein/domain integrity not proven"
        elif integrity_flags:
            integrity = "; ".join(sorted(set(integrity_flags)))
        else:
            integrity = "not assessable"

        if rm_hit_count:
            delivery = f"potential restriction barrier: {rm_hit_count} annotation hit(s); transformation testing required"
        elif defense_status == "complete":
            delivery = "no restriction annotation hit in this screen; other defence barriers may remain"
        else:
            delivery = "unresolved because defence-system analysis was incomplete or unavailable"

        is_crispr = "CRISPR" in system or "RNA-guided" in system
        targeting = (f"supported candidate PAM(s): {', '.join(supported_pams)}" if is_crispr and supported_pams
                     else "PAM/guide requirement unresolved" if is_crispr
                     else "system-specific recognition/attachment requirements must be established")
        confidence = ("moderate computational candidate" if completeness == "minimum component model satisfied" and not integrity_flags
                      else "limited computational candidate" if all_hits else "no genomic candidate evidence")
        if system == "CRISPR interference machinery" and crispr_status != "complete":
            confidence = "limited: specialist CRISPR system call unavailable or incomplete"

        coordinates = "; ".join(
            f"{hit.get('contig')}:{hit.get('start')}-{hit.get('end')}({hit.get('strand')})"
            for hit in sorted(all_hits, key=lambda item: (str(item.get("contig")), int(item.get("start", 0))))
        ) or "none"
        products = "; ".join(sorted({str(hit.get("product") or hit.get("gene") or hit.get("id")) for hit in all_hits})) or "none"
        rows.append({
            "system_model": system, "completeness": completeness,
            "required_component_groups": ", ".join(required),
            "detected_component_groups": ", ".join(detected_groups) or "none",
            "missing_required_groups": ", ".join(missing) or "none",
            "detected_features": products, "genomic_coordinates": coordinates,
            "locus_organization": organization, "candidate_locus_span_bp": locus_span,
            "sequence_integrity_assessment": integrity,
            "expression_status": "not evaluated from genome sequence",
            "functional_activity_status": "not evaluated; requires experimental evidence",
            "pam_or_target_requirement": targeting, "delivery_compatibility": delivery,
            "possible_application": model["application"], "confidence": confidence,
            "evidence_source": "database-backed genome annotation + explicit component model",
            "limitations": "component models are screening rules; domain homology, transcription and biochemical activity are not proven",
            "recommended_validation": model["validation"],
        })
    return rows


def assess_introduced_editors(targets: list[dict], rm_hit_count: int, defense_status: str,
                              repair_supported: bool, annotation_available: bool) -> list[dict]:
    """Summarize chassis compatibility without converting evidence into efficiency."""
    rows = []
    for editor in ("SpCas9", "Cas12a/Cpf1", "dCas9/CRISPRi"):
        candidates = [row for row in targets if row.get("editor") == editor]
        unique = sum(row.get("exact_uniqueness") == "unique" for row in candidates)
        screened = [
            row for row in candidates
            if row.get("off_target_screen_status") == "complete_substitutions_only"
        ]
        skipped_reference_limit = any(
            row.get("off_target_screen_status") == "skipped_reference_limit"
            for row in candidates
        )

        if screened:
            screening_status = "screened"
            screening_interpretation = (
                "Bounded substitution screen completed for the selected top guides; "
                "bulges and noncanonical PAMs remain unresolved."
            )
            clean = sum(
                row.get("approx_offtargets_le2") == 0
                for row in screened
            )
        elif skipped_reference_limit:
            screening_status = "not_screened_reference_limit"
            screening_interpretation = (
                "Not screened: the PAM-compatible reference space exceeded the "
                "configured screening limit; off-target status remains unresolved."
            )
            clean = None
        else:
            screening_status = "not_screened"
            screening_interpretation = (
                "Off-target screening was not performed for the available candidates."
            )
            clean = None

        if rm_hit_count:
            delivery = f"elevated delivery risk: {rm_hit_count} restriction/modification annotation hit(s)"
        elif defense_status == "complete":
            delivery = "no restriction annotation hit; other detected or uncharacterized defenses may remain"
        else:
            delivery = "delivery risk unresolved because specialist defense analysis was unavailable"
        repair = ("native repair-support annotations detected" if repair_supported
                  else "repair support unresolved; CRISPRi avoids double-strand-break repair dependence"
                  if editor == "dCas9/CRISPRi" else "repair support unresolved")
        if not candidates:
            tier = "not targetable under scanned PAM rule"
        elif rm_hit_count or defense_status != "complete":
            tier = "conditional — targeting sites exist; delivery/defense unresolved"
        elif editor != "dCas9/CRISPRi" and not repair_supported:
            tier = "conditional — targeting supported; repair route unresolved"
        else:
            tier = "computationally compatible candidate"
        rows.append({
            "editor": editor, "pam_rule": EDITOR_PRESETS[editor]["pam_pattern"],
            "application": EDITOR_PRESETS[editor]["application"], "pam_compatible_sites": len(candidates),
            "exact_unique_sites": unique, "approximately_screened_top_sites": len(screened),
            "screened_sites_without_1_or_2_mismatch_hit": clean,
            "off_target_screening_status": screening_status,
            "off_target_screening_interpretation": screening_interpretation, 
            "delivery_assessment": delivery, "repair_assessment": repair,
            "gene_context_status": "available" if annotation_available else "unavailable",
            "compatibility_tier": tier,
            "uncertainty": "not an editing-efficiency prediction; expression, toxicity, delivery, mismatch/bulge tolerance and phenotype require testing",
        })
    return rows
