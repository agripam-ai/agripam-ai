#!/usr/bin/env python3
"""Build the competition presentation and poster (isolates named by species and strain code).

Usage: python scripts/build_deliverables.py [preview_dir]
Requires python-pptx and Pillow. Numbers are copied from files in this repository.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from deliverable_kit import Slide, build_pptx, render_preview  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "deliverables/assets"
OUT = ROOT / "deliverables/final"

BLUE = "2B7BB9"
DEEP, PANEL, LEAF, MINT, AMBER, INK, MUTED, WHITE = (
    "0E3B36", "155A52", "2E9E6B", "E6F3EC", "E8A33D", "17322E", "5E7A74", "FFFFFF")
PALE_AMBER, PALE_LEAF, BG = "FCF0DA", "D5EDE0", "F7FAF8"
W, H = 13.333, 7.5


def base(kicker, title, notes=""):
    s = Slide(W, H, BG)
    s.text(0.6, 0.32, 9, 0.3, kicker.upper(), size=12, color=LEAF, bold=True)
    s.text(0.6, 0.62, 12.1, 0.95, title, size=30, color=INK, bold=True, valign="t", min_size=24)
    s.notes = notes
    return s


def footer(s, text):
    s.text(0.6, 6.95, 12.1, 0.35, text, size=11, color=MUTED, min_size=9)


def number_circle(s, x, y, n, d=0.5, fill=LEAF):
    s.oval(x, y, d, d, fill=fill)
    s.text(x, y, d, d, str(n), size=16, color=WHITE, bold=True, align="c", valign="m")


def card(s, x, y, w, h, fill=WHITE, line="D5E4DD"):
    s.rect(x, y, w, h, fill=fill, line=line, radius=0.12)


def slides_presentation():
    out = []

    # 1 title
    s = Slide(W, H, DEEP)
    s.text(0.8, 0.9, 7.6, 0.4, "COMPETITION PROJECT · SYNTHETIC BIOLOGY FOR AGRICULTURE", size=13, color="9ED9BC", bold=True)
    s.text(0.8, 1.55, 7.6, 1.1, "AgriPAM-AI", size=60, color=WHITE, bold=True)
    s.text(0.8, 2.85, 7.6, 1.5, "Which microbe in a community should we edit?", size=32, color=WHITE, bold=True)
    s.text(0.8, 4.55, 7.4, 1.2, "Community-aware, chassis-aware CRISPR design for agricultural microbial communities",
           size=18, color="BFE3D0")
    s.text(0.8, 6.75, 8, 0.35, "Isolates are shown by species and strain code · September 2026", size=12, color="9ED9BC")
    tiles = [("23/23", "observations carry the inferred 5′-TTC PAM in the P. polymyxa reference panel"),
             ("60", "automated tests pass; every random step uses seed 42"),
             ("0", "wet-lab claims: every result is a computational hypothesis")]
    for i, (big, label) in enumerate(tiles):
        y = 0.9 + i * 2.05
        s.rect(9.0, y, 3.8, 1.8, fill=PANEL, radius=0.15)
        s.text(9.25, y + 0.12, 3.3, 0.8, big, size=40, color=WHITE, bold=True)
        s.text(9.25, y + 0.95, 3.3, 0.8, label, size=13, color="BFE3D0")
    s.notes = ("Open with the question, not the tool. A tested microbial bank already exists; the open problem is which "
               "member of a community to edit so the community gains a function without losing one.")
    out.append(s)

    # 2 problem
    s = base("The problem", "A tested microbial bank is not yet a designed, editable community",
             "Three separate difficulties: function gaps, chassis choice, and the cost of editing the wrong member.")
    cards = [("Function gaps", "A compatible community rarely covers every function a crop needs. In our worked example the planned community lacks β-glucanase."),
             ("Chassis choice is hard", "PAM availability, restriction–modification and other defenses, and fit with the other members all differ from strain to strain."),
             ("Wrong edits cost twice", "Editing the only provider of a function can remove it, and some strains raise biosafety questions that must be settled first.")]
    for i, (head, body) in enumerate(cards):
        x = 0.6 + i * 4.1
        card(s, x, 1.85, 3.9, 3.7)
        number_circle(s, x + 0.3, 2.15, i + 1)
        s.text(x + 0.3, 2.85, 3.3, 0.7, head, size=20, color=INK, bold=True)
        s.text(x + 0.3, 3.55, 3.3, 1.9, body, size=15, color=INK)
    s.rect(0.6, 5.85, 12.1, 0.85, fill=MINT, radius=0.12)
    s.text(0.9, 5.85, 11.5, 0.85, "Our question: which member should carry a new function, and where in its genome?",
           size=20, color=DEEP, bold=True, valign="m")
    out.append(s)

    # 3 answer
    s = base("Our answer", "One loop from a tested bank to a ranked, auditable edit plan",
             "Emphasise that the first four steps work on the bank tables alone; the genome layer adds site-level design.")
    steps = [("Tested bank", "30 isolates, 7 biochemical traits, antagonism and compatibility tables."),
             ("Assemble", "Compatible, biosafety-eligible strains that cover the most functions."),
             ("Find gaps", "Functions the community lacks, and why a donor cannot simply join."),
             ("Rank chassis", "Fit, dispensability and editing precedent, behind a biosafety gate."),
             ("Design sites", "Native PAM inferred from spacers, plus introduced editors, scored against genes."),
             ("Build and track", "BioBrick parts, checked constructs and a tracker; frozen benchmarks now.")]
    for i, (head, body) in enumerate(steps):
        row, col = divmod(i, 3)
        x, y = 0.6 + col * 4.15, 1.95 + row * 2.35
        card(s, x, y, 3.65, 1.95)
        number_circle(s, x + 0.25, y + 0.22, i + 1, d=0.46)
        s.text(x + 0.85, y + 0.2, 2.6, 0.5, head, size=18, color=INK, bold=True, valign="m")
        s.text(x + 0.25, y + 0.85, 3.2, 1.0, body, size=14, color=INK)
        if col < 2:
            s.arrow(x + 3.68, y + 0.75, 0.42, 0.4, LEAF)
    footer(s, "Steps 1–4 use the bank tables alone; steps 5–6 add the genome.")
    out.append(s)

    # 4 software
    s = base("The software", "How it works, step by step",
             "Walk through the numbered steps; each one writes a table that appears in the final ZIP.")
    items = [("Load a genome", "FASTA, NCBI accession, reads or built-in demo"),
             ("Check quality", "Length, N50, GC, ambiguity, checksum"),
             ("Annotate and screen", "Genes, plasmids, CRISPR-Cas, defenses, mobile elements"),
             ("Find targets", "Native PAM inferred per isolate (64 triplets) plus SpCas9, Cas12a, dCas9; off-target screen"),
             ("Link to genes", "Overlapping gene, neutral regions, function flags"),
             ("Build a portfolio", "Four objectives, Pareto front, blank results ledger"),
             ("Package", "ZIP of tables, logs, versions and checksums")]
    for i, (head, body) in enumerate(items):
        y = 1.7 + i * 0.73
        number_circle(s, 0.6, y + 0.05, i + 1, d=0.45)
        s.text(1.25, y, 6.1, 0.65, [[(head + ": ", {"bold": True}), body]], size=14, color=INK, valign="m", min_size=11)
    card(s, 7.7, 1.75, 5.0, 3.6)
    s.image(7.85, 1.9, 4.7, 2.6, ASSETS / "workflow.png")
    s.text(7.9, 4.55, 4.6, 0.75, "Reference study workflow, P. polymyxa panel: 24 genomes, 249 unique spacers, 23 target observations.",
           size=12, color=MUTED)
    s.rect(7.7, 5.55, 5.0, 1.15, fill=MINT, radius=0.12)
    s.text(7.95, 5.55, 4.5, 1.15, "Missing tools are reported as “not analysed”, never as “not present”.",
           size=14, color=DEEP, bold=True, valign="m")
    out.append(s)

    # 5 coverage grid
    s = base("Community layer", "Functions delivered by the planned community",
             "B1 (B. metallica) is shown because it is in the laboratory design, but it is excluded as an editing chassis by the biosafety gate.")
    cols = ["Phosphorus", "Potassium", "N-free growth", "β-glucanase", "Chitinase", "Catalase"]
    rows = [("B1 · B. metallica", [4, 4, 3, 0, 1, 3]), ("B39 · E. mendocina", [4, 3, 3, 1, 3, 0]), ("B56 · P. taiwanensis", [5, 3, None, 0, 1, 3]),
            ("Community best", [5, 4, 3, 1, 3, 3])]
    heat = {0: "F0F4F2", 1: "DCEFE4", 2: "BFE3D0", 3: "8FD0AE", 4: "5DBA8A", 5: "2E9E6B"}
    x0, y0, cw, ch, lw = 0.5, 1.95, 1.05, 0.72, 2.6
    for j, name in enumerate(cols):
        s.text(x0 + lw + j * cw, y0, cw, 0.6, name, size=12, color=MUTED, bold=True, align="c", valign="m")
    for i, (label, vals) in enumerate(rows):
        y = y0 + 0.7 + i * (ch + 0.08)
        last = i == len(rows) - 1
        s.text(x0, y, lw, ch, label, size=13.5, color=INK, bold=True, valign="m", min_size=10)
        for j, v in enumerate(vals):
            x = x0 + lw + j * cw
            fill = "ECEFEE" if v is None else heat[v]
            gap = last and v is not None and v < 2
            s.rect(x + 0.04, y, cw - 0.08, ch, fill=fill, line=AMBER if gap else None, lw=3)
            s.text(x, y, cw, ch, "n/t" if v is None else str(v), size=16, color=WHITE if (v or 0) >= 4 else INK,
                   bold=True, align="c", valign="m")
    s.text(0.6, 6.0, 8.9, 0.8, "Scores 0–5 from 7-day plate readings; a function counts as delivered at 2 or more. "
           "n/t = not tested. The amber outline marks the community’s one gap.", size=12, color=MUTED)
    call = [("1 gap", "β-glucanase is below the delivery threshold.", AMBER),
            ("B1", "B. metallica: excluded as a chassis by the biosafety hold in the project’s own screen.", "C0533F"),
            ("B39", "E. mendocina: the only strong chitinase source; editing it puts that function at risk.", LEAF)]
    for i, (big, txt, colr) in enumerate(call):
        y = 1.95 + i * 1.6
        card(s, 9.8, y, 2.95, 1.45)
        s.text(10.0, y + 0.1, 2.6, 0.55, big, size=24, color=colr, bold=True)
        s.text(10.0, y + 0.65, 2.6, 0.75, txt, size=12, color=INK)
    out.append(s)

    # 6 ranking
    s = base("Chassis ranking", "Edit the member the community can spare",
             "B56 (P. taiwanensis): its functions are all covered by another member; it scores at or above B39 under every random weighting, and the fungal anchor FI20 ties with it at 90.")
    card(s, 0.6, 1.8, 6.2, 4.7)
    s.bar_chart(0.75, 1.9, 5.9, 3.7, ["FI20 (fungus)\nT. yunnanense", "B56\nP. taiwanensis", "B39\nE. mendocina"], [90.0, 90.0, 75.0], LEAF, title="Composite score, planned community",
                number_format="0.0", vmin=0, vmax=110, point_colors=["9FB3AD", LEAF, "7CC7A0"])
    s.text(0.85, 5.65, 5.7, 0.75, "Composite = 0.40 × fit + 0.35 × dispensability + 0.25 × editing precedent. B1 (B. metallica) is not ranked.",
           size=12, color=MUTED)
    s.rect(7.1, 1.8, 5.6, 1.5, fill=DEEP, radius=0.15)
    s.text(7.35, 1.85, 2.0, 1.4, "100%", size=40, color=WHITE, bold=True, valign="m")
    s.text(9.4, 1.85, 3.2, 1.4, "of 2,000 random weightings keep B56 (P. taiwanensis) at or above B39 (seed 42)", size=14, color="BFE3D0", valign="m")
    s.text(7.1, 3.5, 5.6, 1.6, [
        [("B56, P. taiwanensis: ", {"bold": True}), "every function it provides is covered by another member; its genus has editing precedent. The fungal anchor FI20 ties at 90."],
        [("B39, E. mendocina: ", {"bold": True}), "the only strong chitinase source (score 75), so a poor editing target."]],
        size=14, color=INK, para_after=8)
    card(s, 7.1, 5.15, 5.6, 1.35, fill=PALE_AMBER, line="F2D9A8")
    s.text(7.3, 5.15, 5.2, 1.35, "A safety-screened alternative (E. oleovorans B55 + B. subtilis B3 with T. yunnanense FI20) covers all six functions. B3 ranks first, "
           "but it is the only β-glucanase source, so any edit must preserve that activity.", size=12, color=INK, valign="m")
    out.append(s)

    # 7 design
    s = base("Edit-site design", "From a chosen chassis to a scored target",
             "The worked example uses the built-in demonstration sequence; every number was produced by the code.")
    card(s, 0.6, 1.8, 7.0, 4.9)
    s.text(0.9, 1.95, 6.4, 0.4, "Design score for a 5′-TTC target", size=16, color=LEAF, bold=True)
    s.text(0.9, 2.45, 6.4, 1.5, [
        "0.60 × PAM model probability",
        "+ 0.20 × GC quality",
        "+ 0.10 ÷ exact genome copies",
        "+ 0.10 × sequence quality (homopolymers)"], size=18, color=INK, bold=True, para_after=4)
    s.text(0.9, 4.25, 6.4, 0.4, "Introduced editors scanned on both strands", size=16, color=LEAF, bold=True)
    s.text(0.9, 4.75, 6.4, 1.9, [
        {"runs": [("SpCas9  ", {"bold": True}), "NGG, 20 nt"], "bullet": True},
        {"runs": [("Cas12a  ", {"bold": True}), "TTTV, 23 nt"], "bullet": True},
        {"runs": [("dCas9/CRISPRi  ", {"bold": True}), "NGG, 20 nt"], "bullet": True},
        {"runs": ["Top 2,000 guides per editor screened for 1–2 mismatches; bulges stay unresolved."], "bullet": True}],
        size=14, color=INK, para_after=4)
    s.rect(7.9, 1.8, 4.8, 2.25, fill=DEEP, radius=0.15)
    s.text(8.15, 1.9, 4.3, 1.0, "82.47", size=54, color=WHITE, bold=True)
    s.text(8.15, 3.0, 4.3, 0.95, "design score, demonstration target (PAM probability 81.27%, GC 34.3%, 1 copy)", size=13, color="BFE3D0")
    s.rect(7.9, 4.25, 4.8, 2.45, fill=MINT, radius=0.15)
    s.text(8.15, 4.35, 4.3, 1.0, "92.03", size=54, color=DEEP, bold=True)
    s.text(8.15, 5.45, 4.3, 1.15, "same target if independent evidence places it in a mobile-element region (PAM probability 97.20%)", size=13, color=INK)
    out.append(s)

    # 8 parts to bench
    s = base("Parts to bench", "From a chosen edit to a checked, tracked construct",
             "Parts come from a starter set, an Excel sheet or the iGEM Registry. The checks are screens: a compatible part is "
             "not a working part, so activity in the chassis still has to be measured.")
    module = [("Choose parts", "Starter set, your Excel sheet, or fetch by BBa_ ID from the iGEM Registry."),
              ("Assemble", "Reporter insertion between homology arms, or a CRISPRi guide cassette. The guide RNA stays seamless."),
              ("Screen", "Forbidden sites (RFC10, RFC25, Golden Gate), chassis restriction motifs, host range, repeats."),
              ("Track", "Stage tracker with parental, negative and low-ranked controls; predicted vs measured ρ.")]
    for i, (head, body) in enumerate(module):
        row, col = divmod(i, 2)
        x, y = 0.6 + col * 3.85, 1.85 + row * 2.45
        card(s, x, y, 3.65, 2.25)
        number_circle(s, x + 0.25, y + 0.22, i + 1, d=0.46)
        s.text(x + 0.85, y + 0.2, 2.6, 0.5, head, size=18, color=INK, bold=True, valign="m")
        s.text(x + 0.25, y + 0.85, 3.2, 1.35, body, size=13, color=INK, min_size=10.5)
    card(s, 8.5, 1.85, 4.2, 4.7)
    s.text(8.75, 1.95, 3.7, 0.7, "Example cassette from four Registry parts", size=15, color=LEAF, bold=True)
    s.text(8.75, 2.65, 3.7, 0.9, "920 bp", size=40, color=DEEP, bold=True)
    parts_draw = [("J23100 promoter", 35, LEAF), ("B0034 RBS", 12, AMBER), ("E0040 GFP", 720, "5DBA8A"), ("B0015 terminator", 129, "2B7BB9")]
    bx, bw_total, cur = 8.75, 3.7, 8.75
    for label, bp, colr in parts_draw:
        width = max(0.07, bw_total * bp / 920)
        s.rect(cur, 3.6, width - 0.01, 0.4, fill=colr)
        cur += width
    for i, (label, bp, colr) in enumerate(parts_draw):
        yy = 4.2 + i * 0.4
        s.rect(8.75, yy + 0.07, 0.22, 0.22, fill=colr)
        s.text(9.1, yy, 3.4, 0.36, f"{label}, {bp} bp", size=12.5, color=INK, valign="m")
    s.text(8.75, 5.85, 3.7, 0.6, "plus three 8-bp scars. Sequences fetched from the iGEM Registry.", size=11, color=MUTED, min_size=9)
    out.append(s)

    # 8 reference discovery
    s = base("Reference discovery", "5′-TTC is the inferred Type I-C PAM in P. polymyxa",
             "These numbers belong to the P. polymyxa reference panel only; they say nothing about any isolate. The software tests all 64 three-base motifs against each isolate's own spacer evidence, corrects for multiple testing, reports only motifs with enough observations and never falls back to a default PAM. TTC is the result for this panel. An isolate with no CRISPR-Cas system gets only the introduced-editor scans (SpCas9 NGG, Cas12a TTTV, dCas9 NGG).")
    card(s, 0.6, 1.75, 12.1, 4.15)
    s.image(0.75, 1.85, 11.8, 3.95, ASSETS / "type_ic_pam_summary_v2.png")
    chips = [("q = 6.8e-34", "BH-adjusted TTC enrichment vs 2,100 matched loci"),
             ("AUC 0.9946", "grouped cross-validation; average precision 0.5200 at 0.99% positives"),
             ("Pending", "functional PAM validation in the lab")]
    for i, (big, txt) in enumerate(chips):
        x = 0.6 + i * 4.1
        s.rect(x, 6.05, 3.9, 0.85, fill=PALE_AMBER if i == 2 else MINT, radius=0.12)
        s.text(x + 0.2, 6.05, 1.75, 0.85, big, size=15, color=INK, bold=True, valign="m", min_size=11)
        s.text(x + 1.95, 6.05, 1.85, 0.85, txt, size=10.5, color=INK, valign="m", min_size=8)
    footer(s, "Method: all 64 three-base motifs are tested per isolate and no default PAM is assumed; TTC is this panel's result.")
    out.append(s)

    # 9 stability
    s = base("Stability check", "The PAM survives leaving out any spacer group or genome",
             "Same pipeline, same data source: a stability check, not independent validation.")
    s.rect(0.6, 1.8, 4.6, 3.35, fill=DEEP, radius=0.15)
    s.text(0.85, 1.9, 4.1, 1.2, "23/23", size=60, color=WHITE, bold=True)
    s.text(0.85, 3.2, 4.1, 1.8, "held-out observations matched in both schemes (Wilson 95% interval 0.857–1.000; TTC is 1.7% of background)",
           size=14, color="BFE3D0")
    s.text(0.6, 5.35, 4.6, 1.4, [
        {"runs": ["Leave-one-spacer-group-out: 11 folds"], "bullet": True},
        {"runs": ["Leave-one-genome-out: 15 folds"], "bullet": True},
        {"runs": ["TTC ranked first in every fold"], "bullet": True}], size=14, color=INK, para_after=3)
    card(s, 5.5, 1.8, 7.2, 3.35)
    s.text(5.8, 1.95, 6.6, 0.4, "Consistent with the literature", size=18, color=LEAF, bold=True)
    s.text(5.8, 2.5, 6.6, 2.6, [
        {"runs": ["Functional PAM screening of a B. halodurans Type I-C Cascade gave an NTTC consensus (Leenay et al. 2016)."], "bullet": True},
        {"runs": ["Structural work on Type I-C Cascade describes the Cas8c subunit reading a 5′-TTC PAM."], "bullet": True},
        {"runs": ["Not yet measured in P. polymyxa."], "bullet": True}], size=18, color=INK, para_after=12)
    s.rect(5.5, 5.35, 7.2, 1.35, fill=PALE_AMBER, radius=0.12)
    s.text(5.75, 5.35, 6.7, 1.35, "Caution: the 23 observations are not independent (the 11 spacer groups are the conservative sample), "
           "and every fold draws on the same pipeline.", size=13, color=INK, valign="m")
    out.append(s)

    # 10 external validation
    s = base("External validation", "Frozen predictions tested against deposited experiments",
             "Predictions were frozen before the outcomes were joined. These test guide ranking, not the P. polymyxa PAM.")
    card(s, 0.6, 1.8, 7.7, 4.85)
    s.bar_chart(0.75, 1.9, 7.4, 4.65, ["Yu 2024\n742 guides", "Cross-chassis\nGFP, 33", "crisprHAL\nholdout", "GC-only\nbaseline", "Permutation\ncontrol"],
                [0.486, 0.4525, 0.6325, 0.3534, -0.0086], LEAF, title="Spearman ρ",
                number_format="0.000", vmin=-0.1, vmax=0.8,
                point_colors=[LEAF, LEAF, LEAF, "9FB3AD", "9FB3AD"])
    s.text(8.6, 1.85, 4.1, 3.4, [
        [("Yu 2024, GSE196911: ", {"bold": True}), "pooled ρ 0.486; median within-gene ρ 0.546."],
        [("Cross-chassis GFP: ", {"bold": True}), "ρ 0.4525, bootstrap 95% interval 0.088–0.722; partial transfer only."],
        [("crisprHAL holdout: ", {"bold": True}), "ρ 0.6325 vs 0.3534 for a GC-only baseline; zero train/test overlap."]],
        size=14, color=INK, para_after=9)
    s.rect(8.6, 5.3, 4.1, 1.35, fill=PALE_AMBER, radius=0.12)
    s.text(8.8, 5.3, 3.7, 1.35, "Guide-ranking benchmarks in other genomes. Not tests of any isolate or of the TTC PAM.", size=13, color=INK, valign="m")
    out.append(s)

    # 11 claims and limits
    s = base("Claims and limits", "What we claim, and what we do not",
             "State the limits before the judges do.")
    card(s, 0.6, 1.8, 5.95, 4.85)
    s.text(0.9, 1.95, 5.4, 0.45, "We show", size=20, color=LEAF, bold=True)
    s.text(0.9, 2.55, 5.4, 4.0, [
        {"runs": ["A reproducible, ranked way to choose an editing chassis inside a community"], "bullet": True},
        {"runs": ["5′-TTC as a stable, literature-consistent PAM hypothesis in the P. polymyxa reference panel"], "bullet": True},
        {"runs": ["Guide-ranking transfer on three deposited benchmarks"], "bullet": True},
        {"runs": ["Checked, trackable constructs built from BioBrick parts"], "bullet": True},
        {"runs": ["Every module states when it could not run"], "bullet": True}], size=19, color=INK, para_after=14)
    s.rect(6.75, 1.8, 5.95, 4.85, fill=PALE_AMBER, line="F2D9A8", radius=0.12)
    s.text(7.05, 1.95, 5.4, 0.45, "We do not claim", size=20, color="B5741A", bold=True)
    s.text(7.05, 2.55, 5.4, 4.0, [
        {"runs": ["Measured editing efficiency or field performance"], "bullet": True},
        {"runs": ["That TTC is the PAM of any isolate"], "bullet": True},
        {"runs": ["That a missing analysis means absent biology"], "bullet": True},
        {"runs": ["Biosafety clearance: B1 (B. metallica) needs review before any use"], "bullet": True},
        {"runs": ["Isolate genomes: bank accessions are same-species references"], "bullet": True}], size=19, color=INK, para_after=14)
    out.append(s)

    # 12 next steps
    s = Slide(W, H, DEEP)
    s.text(0.8, 0.5, 9, 0.35, "NEXT", size=13, color="9ED9BC", bold=True)
    s.text(0.8, 0.85, 11.5, 1.0, "From design to the next build-test-learn turn", size=36, color=WHITE, bold=True)
    phases = [("Design", "Done in the software: community, chassis, editing targets and a checked construct, all ranked."),
              ("Build", "Sequence the isolates, replace the reference accessions, re-run on real genomes, then assemble the construct from BioBrick parts."),
              ("Test", "In P. taiwanensis B56, compare the two top-ranked guides with a low-ranked guide and a no-guide control; read out by PCR and amplicon sequencing."),
              ("Learn", "Feed measured outcomes into the results ledger and the tracker, then retest the community for lost or gained function.")]
    for i, (head, body) in enumerate(phases):
        x = 0.8 + i * 3.05
        s.rect(x, 2.3, 2.75, 3.6, fill=PANEL, radius=0.15)
        number_circle(s, x + 0.25, 2.5, i + 1, fill=LEAF if i else "F2B84B")
        s.text(x + 0.25, 3.15, 2.3, 0.5, head, size=20, color=WHITE, bold=True)
        s.text(x + 0.25, 3.7, 2.3, 2.1, body, size=12.5, color="BFE3D0", min_size=10.5)
        if i < 3:
            s.arrow(x + 2.78, 3.95, 0.24, 0.3, "9ED9BC")
    s.text(0.8, 1.85, 11.8, 0.4, "The synthetic-biology cycle: design, build, test, learn. The software covers the first turn of design; the laboratory closes the loop.", size=13, color="9ED9BC", min_size=11)
    s.text(0.8, 6.15, 11.8, 0.7, "Ranked hypotheses today, measured evidence next.", size=22, color=WHITE, bold=True)
    s.text(0.8, 6.85, 11.8, 0.35, "Isolates are shown by species and strain code. Open the software: streamlit run app/streamlit_app.py", size=12, color="9ED9BC")
    s.notes = ("Close on the synthetic-biology cycle. Design is what the software has done. Build: sequence the isolates and assemble the construct. Test: two top-ranked guides against a low-ranked and a no-guide "
               "control. Learn: record everything in the tracker and feed the measured outcomes back into the ranking.")
    out.append(s)
    return out


PW, PH = 17.54, 12.4


def poster():
    s = Slide(PW, PH, BG)

    def T(*a, **k):
        k.setdefault("min_size", 11)
        s.text(*a, **k)
    s.rect(0, 0, PW, 1.85, fill=DEEP)
    T(0.6, 0.25, 16.3, 0.35, "SYNTHETIC BIOLOGY FOR AGRICULTURE · COMPETITION POSTER", size=13, color="9ED9BC", bold=True)
    T(0.6, 0.62, 16.3, 0.75, "AgriPAM-AI: which microbe in a community should we edit?", size=38, color=WHITE, bold=True)
    T(0.6, 1.4, 16.3, 0.35, "Community-aware, chassis-aware CRISPR design · isolates named by species and strain code", size=16, color="BFE3D0")
    cw, gap, x0, y0 = 5.4, 0.25, 0.5, 2.1

    def col(i):
        return x0 + i * (cw + gap)

    # column 1
    x = col(0)
    card(s, x, y0, cw, 2.35)
    T(x + 0.25, y0 + 0.12, cw - 0.5, 0.4, "The problem", size=21, color=LEAF, bold=True)
    T(x + 0.25, y0 + 0.6, cw - 0.5, 1.7, "A tested microbial bank is not yet a designed community. Compatible communities leave function gaps, "
           "and editing the wrong member can remove a function or raise a biosafety question.", size=16, color=INK)
    card(s, x, y0 + 2.55, cw, 3.05)
    T(x + 0.25, y0 + 2.65, cw - 0.5, 0.4, "How it works", size=21, color=LEAF, bold=True)
    T(x + 0.25, y0 + 3.1, cw - 0.5, 2.45, [
        {"runs": [("1  Assemble ", {"bold": True}), "a compatible, safe community"], "after": 3},
        {"runs": [("2  Rank chassis ", {"bold": True}), "by fit, dispensability, precedent"], "after": 3},
        {"runs": [("3  Design sites ", {"bold": True}), "5′-TTC, Cas9, Cas12a, dCas9 targets"], "after": 3},
        {"runs": [("4  Build ", {"bold": True}), "BioBrick parts into a checked construct"], "after": 3},
        {"runs": [("5  Track ", {"bold": True}), "stages, controls, predicted vs measured"], "after": 0}], size=14.5, color=INK)
    card(s, x, y0 + 5.8, cw, 3.7)
    s.image(x + 0.1, y0 + 5.9, cw - 0.2, 2.9, ASSETS / "workflow.png")
    T(x + 0.25, y0 + 8.85, cw - 0.5, 0.6, "Reference study workflow (P. polymyxa, 24 genomes).", size=13, color=MUTED)

    # column 2
    x = col(1)
    card(s, x, y0, cw, 4.55)
    T(x + 0.25, y0 + 0.12, cw - 0.5, 0.4, "Community layer: chassis ranking", size=21, color=LEAF, bold=True)
    s.bar_chart(x + 0.15, y0 + 0.55, cw - 0.3, 2.6, ["FI20\nT. yunnanense", "B56\nP. taiwanensis", "B39\nE. mendocina"], [90.0, 90.0, 75.0], LEAF, title="Composite score",
                number_format="0.0", vmin=0, vmax=110, point_colors=[LEAF, "7CC7A0"])
    T(x + 0.25, y0 + 3.2, cw - 0.5, 1.3, [
        {"runs": [("B56 (P. taiwanensis) at 90, above B39 at 75 ", {"bold": True}), "in 100% of 2,000 weightings; ties with fungus FI20."], "after": 3},
        {"runs": [("B1 (B. metallica) excluded ", {"bold": True}), "by the biosafety gate; ", ("B39 (E. mendocina) ", {"bold": True}), "is the sole chitinase source."], "after": 0}],
        size=14.5, color=INK)
    card(s, x, y0 + 4.75, cw, 2.15)
    T(x + 0.25, y0 + 4.85, cw - 0.5, 0.4, "Genome layer: worked example", size=21, color=LEAF, bold=True)
    T(x + 0.25, y0 + 5.3, cw - 0.5, 1.55, [
        {"runs": [("82.47 ", {"bold": True, "color": DEEP}), "design score for the demonstration 5′-TTC target"], "after": 4},
        {"runs": [("92.03 ", {"bold": True, "color": DEEP}), "with independent mobile-element evidence"], "after": 4},
        {"runs": ["0.60 PAM + 0.20 GC + 0.10 uniqueness + 0.10 sequence quality"], "after": 0}], size=14.5, color=INK)
    card(s, x, y0 + 7.1, cw, 2.4)
    s.bar_chart(x + 0.1, y0 + 7.15, cw - 0.2, 1.95, ["Observ.\n23/23", "Contexts\n16/16", "Spacer groups\n11/11", "Random\n35/2100"],
                [100, 100, 100, 1.7], BLUE, title="Sequences with 5′-TTC (%)", number_format="0.0", vmin=0, vmax=125,
                point_colors=[BLUE, BLUE, BLUE, "9FB3AD"])
    T(x + 0.25, y0 + 9.1, cw - 0.5, 0.4, "Reference panel: TTC in 23/23 observations, q = 6.8e-34", size=13, color=MUTED)

    # column 3
    x = col(2)
    card(s, x, y0, cw, 3.9)
    s.bar_chart(x + 0.15, y0 + 0.1, cw - 0.3, 3.7, ["Yu\n2024", "GFP\ntransfer", "crisprHAL", "GC\nbaseline", "Permuted"],
                [0.486, 0.4525, 0.6325, 0.3534, -0.0086], LEAF, title="Frozen-prediction validation (Spearman ρ)",
                number_format="0.000", vmin=-0.1, vmax=0.8, point_colors=[LEAF, LEAF, LEAF, "9FB3AD", "9FB3AD"])
    s.rect(x, y0 + 4.1, cw, 2.2, fill=DEEP, radius=0.12)
    T(x + 0.25, y0 + 4.15, 2.4, 1.0, "23/23", size=40, color=WHITE, bold=True, valign="m")
    T(x + 2.55, y0 + 4.15, cw - 2.75, 1.0, "held-out observations recover TTC in every leave-one-group-out fold", size=14, color="BFE3D0", valign="m")
    T(x + 0.25, y0 + 5.2, cw - 0.5, 1.05, "Consistent with published Type I-C work (NTTC consensus in B. halodurans; Cas8c reads 5′-TTC). "
           "Not yet measured in P. polymyxa.", size=13.5, color="BFE3D0")
    card(s, x, y0 + 6.5, cw, 1.6, fill=PALE_AMBER, line="F2D9A8")
    T(x + 0.25, y0 + 6.58, cw - 0.5, 0.4, "Limits", size=19, color="B5741A", bold=True)
    T(x + 0.25, y0 + 6.98, cw - 0.5, 1.1, "Computational hypotheses, not measured efficiency. TTC is this panel's result; every isolate is analysed on its own spacers (64 motifs tested). "
           "B1 (B. metallica) needs biosafety review. Bank genomes are same-species references.", size=13.5, color=INK)
    card(s, x, y0 + 8.3, cw, 1.2)
    T(x + 0.25, y0 + 8.35, cw - 0.5, 1.1, [[("Next: ", {"bold": True, "color": LEAF}),
           "sequence the isolates; test two top guides against a low-ranked and a no-guide control in P. taiwanensis B56; feed results back."]],
           size=14.5, color=INK, valign="m")
    s.rect(0, PH - 0.55, PW, 0.55, fill=DEEP)
    T(0.6, PH - 0.55, 16.3, 0.55, "Software: streamlit run app/streamlit_app.py · 60 automated tests · fixed seed 42 · every module reports when it could not run",
           size=14, color="BFE3D0", valign="m")
    return s


def main(preview_dir=None):
    OUT.mkdir(parents=True, exist_ok=True)
    slides = slides_presentation()
    from build_graphical_abstract import build as build_merged_poster  # lazy: it imports this module
    p = build_merged_poster(False)
    build_pptx(slides, OUT / "AgriPAM-AI_competition_presentation.pptx", W, H)
    build_pptx([p], OUT / "AgriPAM-AI_competition_poster.pptx", PW, p.h)
    problems = []
    if preview_dir:
        prev = Path(preview_dir)
        prev.mkdir(parents=True, exist_ok=True)
        for i, s in enumerate(slides, 1):
            problems += [(i, e.get("text")) for e in render_preview(s, prev / f"slide-{i:02d}.png")]
        problems += [("poster", e.get("text")) for e in render_preview(p, prev / "poster.png", px_per_in=70)]
    for where, text in problems:
        print("OVERFLOW", where, str(text)[:80])
    print("built", len(slides), "slides + poster; overflow items:", len(problems))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)
