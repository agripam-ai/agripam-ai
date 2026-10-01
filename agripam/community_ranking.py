"""Bounded, deterministic ranking of community alternatives from measured bank data."""
from itertools import combinations
from math import comb

from . import syncom


def rank_communities(bank, objectives, anchors=(), max_size=4, limit=10000):
    """Enumerate small spaces; use a bounded beam for larger spaces. Never impute measurements."""
    objectives = list(dict.fromkeys(objectives))
    if not objectives or any(f not in syncom.bank_functions(bank) for f in objectives):
        raise ValueError("Choose at least one function from the bank.")
    anchors = tuple(sorted(set(anchors)))
    eligible = sorted(s for s in bank['traits'] if syncom.safety_status(s, bank)[0])
    if any(s not in eligible for s in anchors):
        return {'rows': [], 'evaluated': 0, 'exhaustive': True, 'reason': 'An anchor is missing from the bank or has a biosafety hold.'}
    pool = [s for s in eligible if s not in anchors]
    threshold = bank.get('threshold', 2)
    max_size = min(max_size, len(eligible))
    examined = 0

    def evaluate(members):
        nonlocal examined
        examined += 1
        pairs = [syncom.pair_score(a, b, bank) for a, b in combinations(members, 2)]
        if any(v is not None and v < 0.5 for v in pairs):
            return None
        known = [v for v in pairs if v is not None]
        observed = {f: [bank['traits'][s].get(f) for s in members
                        if bank['traits'][s].get(f) is not None] for f in objectives}
        covered = [f for f, values in observed.items() if values and max(values) >= threshold]
        unknown = [f for f, values in observed.items() if not values]
        gap = [f for f in objectives if f not in covered and f not in unknown]
        coverage = 100 * len(covered) / len(objectives)
        pair_coverage = 100 * len(known) / len(pairs) if pairs else 0
        compatibility = 100 * sum(known) / len(pairs) if pairs else 0
        strength = 100 * sum(max(v) if v else 0 for v in observed.values()) / (5 * len(objectives))
        score = 0.60 * coverage + 0.25 * compatibility + 0.15 * strength
        labels = syncom.bank_functions(bank)
        return {'members': list(members), 'score': round(score, 2),
                'function_coverage_percent': round(coverage, 1),
                'interaction_evidence_percent': round(pair_coverage, 1),
                'observed_compatibility_percent': round(100 * sum(known) / len(known), 1) if known else None,
                'covered': '; '.join(labels[f] for f in covered),
                'below_threshold': '; '.join(labels[f] for f in gap),
                'not_measured': '; '.join(labels[f] for f in unknown),
                'untested_pairs': len(pairs) - len(known)}

    possible = sum(comb(len(pool), n - len(anchors)) for n in range(max(2, len(anchors)), max_size + 1))
    exhaustive = possible <= limit
    rows = []
    if exhaustive:
        for size in range(max(2, len(anchors)), max_size + 1):
            for rest in combinations(pool, size - len(anchors)):
                row = evaluate(tuple(sorted(anchors + rest)))
                if row:
                    rows.append(row)
    else:
        frontier = [anchors]
        for size in range(len(anchors) + 1, max_size + 1):
            candidates = sorted({tuple(sorted(members + (s,))) for members in frontier for s in pool if s not in members})
            layer = []
            for members in candidates:
                if examined >= limit:
                    break
                row = evaluate(members)
                if row:
                    layer.append(row)
                    if size >= 2:
                        rows.append(row)
            layer.sort(key=lambda r: (-r['score'], tuple(r['members'])))
            frontier = [tuple(r['members']) for r in layer[:100]]
            if examined >= limit or not frontier:
                break
    rows.sort(key=lambda r: (-r['score'], len(r['members']), tuple(r['members'])))
    for i, row in enumerate(rows[:100], 1):
        row['rank'] = i
    return {'rows': rows[:100], 'evaluated': examined, 'exhaustive': exhaustive,
            'reason': 'No qualifying combination found.' if not rows else ''}
