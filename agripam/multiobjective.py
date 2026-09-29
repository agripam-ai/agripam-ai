"""Transparent multi-objective candidate ranking for RhizoForge-Select.

Scores are decision-support heuristics, not estimates of editing efficiency.
Every component is retained in the output so a user can audit the ranking.
"""
from __future__ import annotations

from collections import Counter


SCORE_FIELDS = ("editability_score", "deliverability_score", "agronomic_value_score", "preservation_safety_score")


def _clamp(value: float) -> float:
    return round(max(0.0, min(100.0, value)), 2)


def _pareto_front(rows: list[dict]) -> set[str]:
    front: set[str] = set()
    for candidate in rows:
        dominated = False
        for other in rows:
            if other is candidate:
                continue
            at_least = all(float(other[field]) >= float(candidate[field]) for field in SCORE_FIELDS)
            better = any(float(other[field]) > float(candidate[field]) for field in SCORE_FIELDS)
            if at_least and better:
                dominated = True
                break
        if not dominated:
            front.add(str(candidate["candidate_id"]))
    return front


def rank_editing_candidates(
    candidates: list[dict], pam_ranking: list[dict], rm_hit_count: int,
    defense_status: str, repair_supported: bool,
) -> list[dict]:
    """Return auditable scores, Pareto status, portfolio roles and rejection reasons."""
    if not candidates:
        return []
    pam_support = {str(row.get("pam3", "")): float(row.get("observed_percent", 0) or 0)
                   for row in pam_ranking}
    copy_counts = Counter(str(row.get("protospacer", "")) for row in candidates)
    scored = []
    for row in candidates:
        pam = str(row.get("pam", ""))
        copies = max(1, copy_counts[str(row.get("protospacer", ""))])
        uniqueness = 100.0 if copies == 1 else max(0.0, 100.0 / copies)
        pam_percent = pam_support.get(pam, 0.0)
        repair = 100.0 if repair_supported else 45.0
        editability = _clamp(0.50 * pam_percent + 0.35 * uniqueness + 0.15 * repair)

        defense_penalty = 15.0 if defense_status in {"complete", "available"} else 5.0
        rm_penalty = min(55.0, 8.0 * rm_hit_count)
        deliverability = _clamp(100.0 - rm_penalty - defense_penalty)

        agronomic_flag = str(row.get("agronomic_category_screen", "not flagged"))
        growth_flag = str(row.get("growth_relevance_screen", "not flagged"))
        relation = str(row.get("gene_relation", ""))
        agronomic_value = 80.0 if agronomic_flag != "not flagged" else (45.0 if relation == "promoter-proximal" else 25.0)

        safety = 90.0
        reasons = []
        if growth_flag != "not flagged":
            safety -= 65.0
            reasons.append("growth/essentiality-related annotation flag")
        if agronomic_flag != "not flagged":
            safety -= 35.0
            reasons.append("may disrupt a beneficial agricultural function")
        if relation == "within CDS":
            safety -= 15.0
        if copies > 1:
            safety -= min(30.0, 5.0 * (copies - 1))
            reasons.append(f"protospacer occurs {copies} times in candidate set")
        safety = _clamp(safety)
        rejected = safety < 35.0 or editability < 25.0
        if editability < 25.0:
            reasons.append("weak or missing PAM/guide support")

        out = dict(row)
        out.update({
            "pam_evidence_percent": round(pam_percent, 2), "candidate_exact_copies": copies,
            "repair_compatibility_basis": "repair genes flagged" if repair_supported else "repair support unresolved",
            "editability_score": editability, "deliverability_score": deliverability,
            "agronomic_value_score": _clamp(agronomic_value), "preservation_safety_score": safety,
            "balanced_score": _clamp((editability * deliverability * max(1.0, agronomic_value) * max(1.0, safety)) ** 0.25),
            "decision": "reject" if rejected else "retain",
            "rejection_reasons": "; ".join(reasons) if reasons else "none",
            "possible_application": "candidate gene disruption/regulation; edit type must be specified",
            "limitations": "heuristic scores; off-target, expression, delivery and phenotype require validation",
            "recommended_validation": "sequence-confirmed reporter or phenotype assay with parental, no-guide and lower-ranked controls",
        })
        scored.append(out)

    retained = [row for row in scored if row["decision"] == "retain"]
    front = _pareto_front(retained) if retained else set()
    for row in scored:
        row["pareto_optimal"] = "yes" if row["candidate_id"] in front else "no"
        row["portfolio_role"] = "rejected" if row["decision"] == "reject" else "retained candidate"
    if retained:
        best = max(retained, key=lambda r: (r["balanced_score"], r["editability_score"]))
        safest = max(retained, key=lambda r: (r["preservation_safety_score"], r["balanced_score"]))
        easiest = max(retained, key=lambda r: (min(r["editability_score"], r["deliverability_score"]), r["balanced_score"]))
        best["portfolio_role"] = "best predicted design"
        if safest is not best:
            safest["portfolio_role"] = "biologically safer alternative"
        if easiest is not best and easiest is not safest:
            easiest["portfolio_role"] = "easier experimental alternative"
    return sorted(scored, key=lambda r: (r["decision"] == "reject", -float(r["pareto_optimal"] == "yes"), -float(r["balanced_score"])))


def experimental_results_template(portfolio: list[dict]) -> list[dict]:
    """Create a prospective ledger that can later be joined to measured outcomes."""
    selected = [r for r in portfolio if r.get("portfolio_role") not in {"retained candidate", "rejected"}][:3]
    return [{
        "candidate_id": row["candidate_id"], "portfolio_role": row["portfolio_role"],
        "experimental_status": "not tested", "construct_or_edit_id": "",
        "sequence_confirmed": "", "reporter_or_phenotype": "", "measured_value": "",
        "units": "", "replicate_count": "", "control_normalized_effect": "",
        "notes": "",
    } for row in selected]
