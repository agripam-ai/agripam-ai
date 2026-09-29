#!/usr/bin/env python3
"""Leave-one-group-out check of the inferred Type I-C PAM.

For each held-out group (spacer group, or target genome assembly) the PAM is
re-inferred from the remaining observations only, by maximum enrichment over the
matched-background triplet frequencies. The held-out observations are then scored
against that prediction. Deterministic; no random component.
"""
import json
from pathlib import Path

import pandas as pd
from scipy.stats import binomtest

PROJECT = Path(__file__).resolve().parents[1]
EVIDENCE = PROJECT / "results/type_ic_exact_pam_evidence.tsv"
BACKGROUND = PROJECT / "results/pam3_enrichment_vs_matched_controls.tsv"
OUT = PROJECT / "results/heldout_pam_validation"
PSEUDOCOUNT = 0.5

obs = pd.read_csv(EVIDENCE, sep="\t").drop_duplicates(
    ["query_id", "target_accession", "sstart", "send"]
)
bg = pd.read_csv(BACKGROUND, sep="\t").set_index("PAM3")
bg_freq = (bg["control_count"] + PSEUDOCOUNT) / (bg["control_total"] + 64 * PSEUDOCOUNT)


def infer(train):
    counts = train["normalized_PAM3"].value_counts().reindex(bg.index, fill_value=0)
    freq = (counts + PSEUDOCOUNT) / (len(train) + 64 * PSEUDOCOUNT)
    return (freq / bg_freq).sort_values(ascending=False).index[0]


def leave_one_out(column):
    rows = []
    for group in sorted(obs[column].astype(str).unique()):
        held = obs[obs[column].astype(str) == group]
        train = obs[obs[column].astype(str) != group]
        pred = infer(train)
        hits = int((held["normalized_PAM3"] == pred).sum())
        rows.append({"held_out": group, "n_train": len(train), "n_heldout": len(held),
                     "predicted_PAM": pred, "heldout_matching_prediction": hits})
    df = pd.DataFrame(rows)
    hits, total = int(df.heldout_matching_prediction.sum()), int(df.n_heldout.sum())
    p_bg = float(bg_freq[df.predicted_PAM.mode()[0]])
    test = binomtest(hits, total, p_bg, alternative="greater")
    ci = binomtest(hits, total).proportion_ci(method="wilson")
    return df, {"scheme": f"leave-one-{column}-out", "groups": len(df),
                "pooled_heldout_hits": hits, "pooled_heldout_total": total,
                "predicted_PAM_modal": df.predicted_PAM.mode()[0],
                "PAM_stable_across_folds": bool(df.predicted_PAM.nunique() == 1),
                "background_frequency": round(p_bg, 5),
                "binomial_p_greater": test.pvalue,
                "wilson95_low": round(ci.low, 4), "wilson95_high": round(ci.high, 4)}


OUT.mkdir(exist_ok=True)
summary = []
for col, name in (("query_id", "spacer_group"), ("target_assembly", "genome")):
    df, s = leave_one_out(col)
    df.to_csv(OUT / f"leave_one_{name}_out.tsv", sep="\t", index=False)
    summary.append(s)
(OUT / "summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
