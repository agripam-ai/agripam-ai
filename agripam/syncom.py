"""Community-aware chassis selection for synthetic microbial communities (SynComs).

Given a bank of tested strains, assemble a compatible community, find its
functional gaps, and rank which member is the best genome-editing chassis.
Scores are transparent decision-support heuristics, not predictions of editing
efficiency or field performance. Blank laboratory cells are missing data.
"""
from __future__ import annotations

import csv
import random
from pathlib import Path

FUNCTIONS = {
    "P_solubilization": "phosphorus solubilization",
    "K_solubilization": "potassium solubilization",
    "N_free_growth": "growth on N-free medium",
    "glucanase": "beta-glucanase",
    "chitinase": "chitinase",
    "catalase": "catalase / stress tolerance",
}
FUNCTION_THRESHOLD = 2  # '++' or stronger counts as a delivered function
COMPAT_SCORE = {"+++": 1.0, "++": 1.0, "+": 0.5, "-+": 0.25, "-": 0.0}
TOOL_PRECEDENCE = {"species": 100.0, "genus": 60.0, "none": 20.0}
WEIGHTS = {"fit": 0.40, "dispensability": 0.35, "editability": 0.25}


def read_tsv(path: Path) -> list[dict]:
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_bank(data_dir: Path) -> dict:
    """Load the strain bank. Decisions use codes and derived fields only, never organism names."""
    traits = {}
    for row in read_tsv(data_dir / "strain_traits.tsv"):
        traits[row["strain"]] = {
            **{f: (int(row[f]) if row[f] != "" else None) for f in FUNCTIONS if f in row},
            "oxidase": int(row["oxidase"]) if row["oxidase"] != "" else None,
        }
    meta = {r["strain"]: r for r in read_tsv(data_dir / "strain_metadata.tsv")}
    fungi = {}
    for row in read_tsv(data_dir / "bacteria_fungi_compatibility.tsv"):
        fungi.setdefault(row["strain"], {})[row["fungus"]] = row["symbol"]
    pairs = {}
    for row in read_tsv(data_dir / "bacteria_pairwise_compatibility.tsv"):
        if row["outcome"] == "mixes":
            pairs.setdefault(frozenset((row["strain"], row["base_strain"])), "mixes")
    antag = {}
    for row in read_tsv(data_dir / "bacteria_pathogen_antagonism.tsv"):
        antag.setdefault(row["strain"], {})[row["pathogen"]] = int(row["score_0_4"])
    return {"traits": traits, "meta": meta, "fungi": fungi, "pairs": pairs, "antagonism": antag,
            "functions": dict(FUNCTIONS), "threshold": FUNCTION_THRESHOLD, "kind": {}}


def safety_status(strain: str, bank: dict) -> tuple[bool, str]:
    """Return (eligible, reason) from the precomputed biosafety_hold field."""
    row = bank["meta"][strain]
    return row["biosafety_hold"] != "yes", row["biosafety_reason"]


def function_score(strain: str, function: str, bank: dict) -> int:
    value = bank["traits"][strain].get(function)
    return 0 if value is None else value


def bank_functions(bank: dict) -> dict[str, str]:
    return bank.get("functions", FUNCTIONS)


def coverage(members: list[str], bank: dict) -> dict[str, int]:
    return {f: max((function_score(m, f, bank) for m in members), default=0) for f in bank_functions(bank)}


def delivered(cov: dict[str, int], bank: dict | None = None) -> set[str]:
    threshold = bank.get("threshold", FUNCTION_THRESHOLD) if bank else FUNCTION_THRESHOLD
    return {f for f, s in cov.items() if s >= threshold}


def dispensability(strain: str, members: list[str], bank: dict) -> tuple[float, list[str]]:
    """Share of delivered community functions that survive if this member's own activity is lost."""
    full = delivered(coverage(members, bank), bank)
    if not full:
        return 100.0, []
    rest = delivered(coverage([m for m in members if m != strain], bank), bank)
    lost = sorted(full - rest)
    return round(100.0 * (1 - len(lost) / len(full)), 1), lost


def kind_of(strain: str, bank: dict) -> str:
    return bank.get("kind", {}).get(strain, "bacterium")


def pair_score(a: str, b: str, bank: dict) -> float | None:
    """Compatibility 0..1 between two strains from any evidence (pair table or bacterium-fungus symbol); None if untested."""
    value = bank["pairs"].get(frozenset((a, b)))
    if value == "mixes":
        return 1.0
    if value == "inhibits":
        return 0.0
    if isinstance(value, (int, float)):
        return float(value)
    for x, y in ((a, b), (b, a)):
        symbol = bank["fungi"].get(x, {}).get(y)
        if symbol is not None:
            return COMPAT_SCORE.get(symbol.replace(" ", ""))
    return None


def fit(strain: str, members: list[str], fungi: list[str], bank: dict) -> tuple[float, int, list[str]]:
    """Mean compatibility with the other members over the tests that exist (no imputation).

    `fungi` are partner fungi that are not scored members; fungi that are members are counted as members.
    """
    values, notes = [], []
    for other in members:
        if other == strain:
            continue
        score = pair_score(strain, other, bank)
        if score is None:
            continue
        values.append(score)
        if score == 0:
            outcome = bank["pairs"].get(frozenset((strain, other)))
            notes.append(f"inhibited by/inhibits {other}" if outcome == "inhibits" else f"inhibits/incompatible with {other}")
    for fungus in fungi:
        if fungus in members:
            continue
        symbol = bank["fungi"].get(strain, {}).get(fungus)
        if symbol is not None:
            score = COMPAT_SCORE.get(symbol.replace(" ", ""))
            if score is not None:
                values.append(score)
                if score == 0:
                    notes.append(f"inhibits/incompatible with {fungus}")
    return (round(100 * sum(values) / len(values), 1) if values else float("nan")), len(values), notes


def editability(strain: str, bank: dict) -> tuple[float, str]:
    row = bank["meta"][strain]
    return TOOL_PRECEDENCE[row["editing_precedent_tier"]], row["editing_precedent_note"]


def rank_chassis(members: list[str], fungi: list[str], bank: dict) -> list[dict]:
    rows = []
    for strain in members:
        ok, reason = safety_status(strain, bank)
        disp, lost = dispensability(strain, members, bank)
        f, n, notes = fit(strain, members, fungi, bank)
        edit, edit_note = editability(strain, bank)
        fit_value = 0.0 if f != f else f  # missing evidence scores 0 and is flagged, never imputed
        composite = (WEIGHTS["fit"] * fit_value + WEIGHTS["dispensability"] * disp
                     + WEIGHTS["editability"] * edit)
        rows.append({"strain": strain, "safety_eligible": ok, "safety_note": reason,
                     "fit_score": f, "fit_evidence_n": n, "fit_notes": "; ".join(notes),
                     "dispensability_score": disp, "functions_lost_if_removed": ";".join(lost),
                     "editability_precedent_score": edit, "editability_note": edit_note,
                     "composite": round(composite if ok else 0.0, 1)})
    eligible = sorted([r for r in rows if r["safety_eligible"]], key=lambda r: -r["composite"])
    blocked = [r for r in rows if not r["safety_eligible"]]
    for i, r in enumerate(eligible, 1):
        r["rank"] = i
    for r in blocked:
        r["rank"] = ""
    return eligible + blocked


def weight_sensitivity(members: list[str], fungi: list[str], bank: dict,
                       draws: int = 2000, seed: int = 42) -> dict[str, float]:
    """Share of random weightings (Dirichlet(1,1,1)) for which each eligible strain ranks first."""
    rng = random.Random(seed)
    base = [r for r in rank_chassis(members, fungi, bank) if r["safety_eligible"]]
    wins = {r["strain"]: 0 for r in base}
    for _ in range(draws):
        w = [rng.expovariate(1.0) for _ in range(3)]
        s = sum(w)
        w = [x / s for x in w]
        best = max(base, key=lambda r: w[0] * (0 if r["fit_score"] != r["fit_score"] else r["fit_score"])
                   + w[1] * r["dispensability_score"] + w[2] * r["editability_precedent_score"])
        wins[best["strain"]] += 1
    return {k: round(v / draws, 3) for k, v in wins.items()}


def gap_analysis(members: list[str], fungi: list[str], bank: dict) -> list[dict]:
    """Functions the community lacks, with the bank strains that carry them and why they can't just join."""
    cov = coverage(members, bank)
    rows = []
    for function, label in bank_functions(bank).items():
        if cov[function] >= bank.get("threshold", FUNCTION_THRESHOLD):
            continue
        donors = []
        for strain in bank["traits"]:
            if strain in members or function_score(strain, function, bank) < 3:
                continue
            f, n, notes = fit(strain, members, fungi, bank)
            ok, _ = safety_status(strain, bank)
            barrier = notes or ([] if f == f else ["compatibility not tested"])
            donors.append((function_score(strain, function, bank), strain, ok, barrier))
        donors.sort(key=lambda d: (-d[0], d[1]))
        blocked = [f"{s} ({sc}/5): " + ("; ".join(b) if b else "compatible") + ("" if ok else " [flagged by project biosafety screen]")
                   for sc, s, ok, b in donors if b or not ok]
        joinable = [f"{s} ({sc}/5)" for sc, s, ok, b in donors if not b and ok]
        rows.append({"missing_function": label, "community_best_score": cov[function],
                     "donors_that_can_join": "; ".join(joinable[:6]),
                     "donors_blocked_by_incompatibility_or_safety": "; ".join(blocked[:6]),
                     "engineering_case": "transplant function into a compatible chassis" if blocked and not joinable
                     else "add a compatible donor" if joinable else "no donor in bank"})
    return rows


def assemble(bank: dict, fungi: list[str], size: int = 4, candidates: list[str] | None = None,
             count_anchors: bool = False) -> list[str]:
    """Greedy set cover: safe, compatible strains that add the most delivered functions.

    With `count_anchors`, the listed fungi that are scored members of the bank start the community, so their
    functions count and they are returned first; `size` limits the bacteria added.
    """
    anchors = [f for f in fungi if f in bank["traits"]] if count_anchors else []
    default_pool = [s for s in bank["traits"] if kind_of(s, bank) == "bacterium"]
    pool = []
    for strain in (candidates or default_pool):
        if strain in anchors or not safety_status(strain, bank)[0]:
            continue
        if any((c := pair_score(strain, fn, bank)) is None or c < 0.5 for fn in fungi):
            continue  # untested or inhibitory against a community fungus
        pool.append(strain)
    chosen: list[str] = list(anchors)
    added = 0
    while added < size and pool:
        have = delivered(coverage(chosen, bank), bank)

        def gain(s):
            new = delivered(coverage(chosen + [s], bank), bank) - have
            return (len(new), sum(function_score(s, f, bank) for f in bank_functions(bank)), s)
        best = max(pool, key=gain)
        if gain(best)[0] == 0 and chosen:
            break
        chosen.append(best)
        pool.remove(best)
        added += 1
        pool = [s for s in pool if (pair_score(s, best, bank) is None or pair_score(s, best, bank) >= 0.5)]
    return chosen
