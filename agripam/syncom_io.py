"""Read a user-supplied SynCom workbook, validate it, and run the chassis-selection analysis.

The workbook has these sheets (see data/syncom/SynCom_input_template.xlsx):

* ``Strains``: one row per strain: ``strain`` (code), optional ``identification``, ``biosafety_hold``
  (yes/no), ``pathogen_screen_score`` (0-100), ``editing_precedent`` (species/genus/none), then any
  number of function columns scored 0-5. Blank cells mean "not tested", never "negative".
  ``kind`` is ``bacterium`` (default) or ``fungus``: fungi with function scores are scored community members.
* ``Fungi_compatibility``: ``strain``, ``fungus``, ``symbol`` (+++, ++, +, -+, -): a bacterium against a fungus.
* ``Fungus_compatibility`` (optional): ``fungus_a``, ``fungus_b``, ``symbol``: two fungi against each other.
* ``Bacteria_compatibility``: ``strain_a``, ``strain_b``, ``outcome`` (mixes / inhibits).
* ``Community`` (optional): ``member``, ``type`` (bacterium / fungus). If it is empty, a compatible
  community is assembled automatically.
"""
from __future__ import annotations

import io
from itertools import combinations
from typing import Any

import pandas as pd

from . import syncom

META_COLUMNS = {"strain", "kind", "identification", "biosafety_hold", "pathogen_screen_score",
                "editing_precedent", "notes"}
YES = {"yes", "y", "true", "1"}
NO = {"no", "n", "false", "0"}
PRECEDENT = {"species": "species-level editing precedent declared by the user",
             "genus": "genus-level precedent declared by the user",
             "none": "no editing precedent declared"}
PATHOGEN_LIMIT = 30.0


def _clean(value: Any) -> str:
    if value is None or (isinstance(value, float) and value != value):
        return ""
    return str(value).strip()


def read_workbook(source) -> dict[str, pd.DataFrame]:
    """Read every sheet as text so nothing is silently coerced (requires openpyxl)."""
    sheets = pd.read_excel(source, sheet_name=None, dtype=object, engine="openpyxl")
    return {name.strip(): frame for name, frame in sheets.items()}


def build_bank(sheets: dict[str, pd.DataFrame]) -> tuple[dict | None, list[dict]]:
    """Return (bank, issues). bank is None if the workbook cannot be analysed."""
    issues: list[dict] = []

    def problem(level: str, where: str, message: str):
        issues.append({"level": level, "where": where, "message": message})

    if "Strains" not in sheets:
        problem("error", "Strains", "The sheet 'Strains' is missing.")
        return None, issues
    frame = sheets["Strains"].copy()
    frame.columns = [_clean(c).lower().replace(" ", "_") for c in frame.columns]
    if "strain" not in frame.columns:
        problem("error", "Strains", "The column 'strain' is missing.")
        return None, issues
    frame = frame[frame["strain"].map(_clean) != ""]
    frame = frame[~frame["strain"].map(_clean).str.upper().str.startswith("EXAMPLE")]
    if frame.empty:
        problem("error", "Strains", "No strains found. Replace the EXAMPLE rows with your own.")
        return None, issues
    codes = frame["strain"].map(_clean)
    for dup in sorted(set(codes[codes.duplicated()])):
        problem("error", "Strains", f"Strain code '{dup}' appears more than once.")
    trait_cols = [c for c in frame.columns if c not in META_COLUMNS and not c.startswith("unnamed")]
    if not trait_cols:
        problem("error", "Strains", "No function columns found (add columns scored 0-5).")
        return None, issues

    traits, meta, kinds = {}, {}, {}
    for _, row in frame.iterrows():
        code = _clean(row["strain"])
        kind = _clean(row.get("kind", "")).lower() or "bacterium"
        if kind not in {"bacterium", "fungus"}:
            problem("error", "Strains", f"{code}: kind '{kind}' must be bacterium or fungus.")
            kind = "bacterium"
        kinds[code] = kind
        scores = {}
        for col in trait_cols:
            raw = _clean(row[col])
            if raw == "":
                scores[col] = None
                continue
            try:
                value = float(raw)
            except ValueError:
                problem("error", "Strains", f"{code}: '{col}' = '{raw}' is not a number.")
                scores[col] = None
                continue
            if not 0 <= value <= 5:
                problem("error", "Strains", f"{code}: '{col}' = {raw} is outside 0-5.")
                value = max(0.0, min(5.0, value))
            scores[col] = int(round(value))
        traits[code] = scores

        hold_raw = _clean(row.get("biosafety_hold", "")).lower()
        score_raw = _clean(row.get("pathogen_screen_score", ""))
        reasons = []
        if hold_raw in YES:
            reasons.append("biosafety hold declared")
        elif hold_raw not in NO:
            problem("warning", "Strains", f"{code}: biosafety_hold is blank, so this strain is treated as eligible. "
                    "Assess biosafety before any use or editing.")
        if score_raw:
            try:
                if float(score_raw) >= PATHOGEN_LIMIT:
                    reasons.append(f"pathogen-screen score {score_raw} >= {PATHOGEN_LIMIT:g}")
            except ValueError:
                problem("error", "Strains", f"{code}: pathogen_screen_score '{score_raw}' is not a number.")
        tier = _clean(row.get("editing_precedent", "")).lower() or "none"
        if tier not in PRECEDENT:
            problem("error", "Strains", f"{code}: editing_precedent '{tier}' must be species, genus or none.")
            tier = "none"
        meta[code] = {"biosafety_hold": "yes" if reasons else "no",
                      "biosafety_reason": "; ".join(reasons) or (
                          "no hold declared" if hold_raw in NO else "biosafety not assessed (treated as eligible)"),
                      "editing_precedent_tier": tier, "editing_precedent_note": PRECEDENT[tier],
                      "identification": _clean(row.get("identification", ""))}

    fungi: dict[str, dict[str, str]] = {}
    if "Fungi_compatibility" in sheets:
        f = sheets["Fungi_compatibility"].copy()
        f.columns = [_clean(c).lower() for c in f.columns]
        for _, row in f.iterrows():
            strain, fungus, symbol = _clean(row.get("strain")), _clean(row.get("fungus")), _clean(row.get("symbol"))
            if not (strain and fungus and symbol) or strain.upper().startswith("EXAMPLE"):
                continue
            if strain not in traits:
                problem("warning", "Fungi_compatibility", f"'{strain}' is not in the Strains sheet; row ignored.")
            elif symbol.replace(" ", "") not in syncom.COMPAT_SCORE:
                problem("error", "Fungi_compatibility", f"{strain} vs {fungus}: symbol '{symbol}' must be one of +++, ++, +, -+, -.")
            else:
                fungi.setdefault(strain, {})[fungus] = symbol
    pairs: dict[frozenset, str] = {}
    if "Bacteria_compatibility" in sheets:
        b = sheets["Bacteria_compatibility"].copy()
        b.columns = [_clean(c).lower() for c in b.columns]
        for _, row in b.iterrows():
            a, c, outcome = _clean(row.get("strain_a")), _clean(row.get("strain_b")), _clean(row.get("outcome")).lower()
            if not (a and c and outcome) or a.upper().startswith("EXAMPLE"):
                continue
            if a not in traits or c not in traits:
                problem("warning", "Bacteria_compatibility", f"'{a}' or '{c}' is not in the Strains sheet; row ignored.")
            elif outcome not in {"mixes", "inhibits"}:
                problem("error", "Bacteria_compatibility", f"{a} / {c}: outcome '{outcome}' must be mixes or inhibits.")
            else:
                pairs[frozenset((a, c))] = outcome
    if "Fungus_compatibility" in sheets:
        ff = sheets["Fungus_compatibility"].copy()
        ff.columns = [_clean(c).lower() for c in ff.columns]
        for _, row in ff.iterrows():
            a, c, symbol = _clean(row.get("fungus_a")), _clean(row.get("fungus_b")), _clean(row.get("symbol")).replace(" ", "")
            if not (a and c and symbol) or a.upper().startswith("EXAMPLE"):
                continue
            if a not in traits or c not in traits:
                problem("warning", "Fungus_compatibility", f"'{a}' or '{c}' is not in the Strains sheet; row ignored.")
            elif symbol not in syncom.COMPAT_SCORE:
                problem("error", "Fungus_compatibility", f"{a} / {c}: symbol '{symbol}' must be one of +++, ++, +, -+, -.")
            else:
                pairs[frozenset((a, c))] = syncom.COMPAT_SCORE[symbol]
    if not fungi and not pairs:
        problem("warning", "Compatibility", "No compatibility data supplied: every fit score will be 'not tested'.")

    functions = {c: c.replace("_", " ") for c in trait_cols}
    bank = {"traits": traits, "meta": meta, "fungi": fungi, "pairs": pairs, "antagonism": {},
            "functions": functions, "threshold": syncom.FUNCTION_THRESHOLD, "kind": kinds}
    if any(i["level"] == "error" for i in issues):
        return None, issues
    return bank, issues


def community_from_sheet(sheets: dict[str, pd.DataFrame], bank: dict) -> tuple[list[str], list[str]]:
    """Members/fungi listed in the optional Community sheet (empty lists if absent)."""
    if "Community" not in sheets:
        return [], []
    c = sheets["Community"].copy()
    c.columns = [_clean(x).lower() for x in c.columns]
    members, fungi = [], []
    for _, row in c.iterrows():
        name, kind = _clean(row.get("member")), _clean(row.get("type")).lower()
        if not name or name.upper().startswith("EXAMPLE"):
            continue
        if kind.startswith("fung"):
            fungi.append(name)
        elif name in bank["traits"]:
            members.append(name)
    return members, fungi


def analyze(bank: dict, members: list[str] | None = None, fungi: list[str] | None = None,
            size: int = 4, threshold: int | None = None) -> dict:
    """Assemble (if needed), find gaps, rank chassis, and check robustness.

    Fungi listed in `fungi` that are scored members of the bank (kind = fungus) join the community and their
    functions count; other listed fungi are compatibility partners only.
    """
    fungi = fungi or []
    if threshold is not None:
        bank = {**bank, "threshold": threshold}
    anchors = [f for f in fungi if f in bank["traits"]]
    partners = [f for f in fungi if f not in bank["traits"]]
    auto = not members
    if auto:
        community = syncom.assemble(bank, fungi, size=size, count_anchors=True)
    else:
        community = anchors + [m for m in members if m not in anchors]
    warnings = []
    for a, b in combinations(community, 2):
        pair = syncom.pair_score(a, b, bank)
        if pair is None:
            warnings.append(f"{a} / {b}: interaction not measured; compatibility unresolved.")
        elif pair < 0.5:
            warnings.append(f"{a} / {b}: measured inhibitory or low compatibility; revise the community before editing.")
    for i, a in enumerate(anchors):
        for b in anchors[i + 1:]:
            score = syncom.pair_score(a, b, bank)
            if score is not None and score < 0.5:
                warnings.append(f"{a} and {b} inhibit each other (compatibility {score:.2f}): they should not share a formulation.")
    cov = syncom.coverage(community, bank)
    delivered = syncom.delivered(cov, bank)
    labels = syncom.bank_functions(bank)
    ranking = syncom.rank_chassis(community, partners, bank) if community else []
    for row in ranking:
        row["kind"] = syncom.kind_of(row["strain"], bank)
    return {
        "members": community, "fungi": fungi, "kinds": {m: syncom.kind_of(m, bank) for m in community},
        "auto_assembled": auto, "warnings": warnings,
        "coverage": {labels[f]: s for f, s in cov.items()},
        "delivered": sorted(labels[f] for f in delivered),
        "missing": sorted(labels[f] for f in cov if f not in delivered),
        "ranking": ranking,
        "gaps": syncom.gap_analysis(community, partners, bank) if community else [],
        "sensitivity": syncom.weight_sensitivity(community, partners, bank) if community else {},
        "threshold": bank.get("threshold", syncom.FUNCTION_THRESHOLD),
    }


def results_to_xlsx(result: dict, issues: list[dict]) -> bytes:
    """Package the analysis as a workbook the user can keep or share."""
    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
        pd.DataFrame(result["ranking"]).to_excel(writer, sheet_name="Chassis_ranking", index=False)
        pd.DataFrame(result["gaps"]).to_excel(writer, sheet_name="Functional_gaps", index=False)
        pd.DataFrame({"function": list(result["coverage"]), "best_member_score": list(result["coverage"].values())}
                     ).to_excel(writer, sheet_name="Coverage", index=False)
        pd.DataFrame({"member": result["members"]}).to_excel(writer, sheet_name="Community", index=False)
        notes = [
            ("Community source", "assembled automatically" if result["auto_assembled"] else "supplied by the user"),
            ("Fungal partners", ", ".join(result["fungi"]) or "none"),
            ("Fungi scored as members", ", ".join(m for m, k in result["kinds"].items() if k == "fungus") or "none"),
            ("Delivery threshold (score >=)", result["threshold"]),
            ("Delivered functions", "; ".join(result["delivered"]) or "none"),
            ("Missing functions", "; ".join(result["missing"]) or "none"),
            ("Top choice under 2,000 random weightings (seed 42)",
             "; ".join(f"{k} {v:.0%}" for k, v in result["sensitivity"].items())),
            ("Interpretation", "Decision support only: not measured editing efficiency, safety clearance or field performance."),
        ] + [("warning: community", w) for w in result["warnings"]] + [(f"{i['level']}: {i['where']}", i["message"]) for i in issues]
        pd.DataFrame(notes, columns=["item", "value"]).to_excel(writer, sheet_name="Notes", index=False)
    return buffer.getvalue()
