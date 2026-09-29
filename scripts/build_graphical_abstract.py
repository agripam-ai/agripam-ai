#!/usr/bin/env python3
"""Competition poster built as a graphical abstract: the chain bank -> community -> chassis -> target -> construct on top, evidence and limits below.

Usage: python scripts/build_graphical_abstract.py [preview_dir]
Writes two files:
  deliverables/final/AgriPAM-AI_competition_poster.pptx            (anonymous, for the blind review)
  deliverables/final_nonblind/AgriPAM-AI_competition_poster_with_banner.pptx  (institution banner; NOT for round 1)
"""
import struct
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_deliverables import (ASSETS, BG, BLUE, DEEP, INK, LEAF, AMBER, MINT, MUTED, PALE_AMBER, PANEL, PW, PH,  # noqa: E402
                                WHITE, card, number_circle)
from deliverable_kit import Slide, build_pptx, render_preview  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "manuscript/figures"
BANNER = ROOT / "deliverables/assets_nonblind/institution_banner.png"
PATH_RED = "C0504D"


def png_size(path):
    return struct.unpack(">II", Path(path).read_bytes()[16:24])


def bank_counts():
    inv = pd.read_csv(ROOT / "manuscript/appendices/tsv/S1_inventory.tsv", sep="\t")
    return int((inv.kind == "bacterium").sum()), int((inv.role == "beneficial fungus").sum()), int((inv.role == "pathogen fungus").sum())


def build(banner=False):
    nb, nf, npth = bank_counts()
    rk = pd.read_excel(ROOT / "manuscript/analysis/analysis_results.xlsx", "A_ranking")
    rk = rk[rk.composite.notna()].sort_values("composite", ascending=False)
    bh_ = png_size(BANNER)[1] * PW / png_size(BANNER)[0] if banner else 0.0
    PHH = 2.0 + 6.15 + 0.3 + 4.55 + 0.3 + 0.5 + bh_
    s = Slide(PW, PHH, BG)

    def T(*a, **k):
        k.setdefault("min_size", 10.5)
        s.text(*a, **k)

    dy = 0.0
    if banner:
        bw, bh = png_size(BANNER)
        dy = PW * bh / bw
        s.image(0, 0, PW, dy, BANNER)
    # header
    s.rect(0, dy, PW, 1.75, fill=DEEP)
    T(0.6, dy + 0.2, 16.3, 0.35, "SYNTHETIC BIOLOGY FOR AGRICULTURE · COMPETITION POSTER", size=13, color="9ED9BC", bold=True)
    T(0.6, dy + 0.55, 16.3, 0.7, "AgriPAM-AI: from a tested microbial bank to a designed community and its first edit", size=34, color=WHITE, bold=True)
    T(0.6, dy + 1.25, 16.3, 0.4, "One open workflow: combine, choose what to edit, find the target, build the construct, follow the experiment", size=16, color="BFE3D0")

    n, gap = 5, 0.32
    pw_ = (PW - 0.8 - gap * (n - 1)) / n
    x0, y0, ph = 0.4, dy + 2.0, 6.15
    heads = ["Tested bank", "Design the community", "Choose the chassis", "Find the target", "Build and follow"]
    for i in range(n):
        x = x0 + i * (pw_ + gap)
        card(s, x, y0, pw_, ph)
        number_circle(s, x + 0.2, y0 + 0.2, i + 1, d=0.55, fill=LEAF)
        T(x + 0.9, y0 + 0.2, pw_ - 1.0, 0.6, heads[i], size=19, color=DEEP, bold=True, valign="m")
        if i < n - 1:
            s.arrow(x + pw_ + 0.03, y0 + ph / 2 - 0.2, gap - 0.06, 0.4, LEAF)
    vy, vh = y0 + 1.0, 2.6  # visual zone

    # 1 bank: dots
    x = x0
    dots = [LEAF] * nb + [AMBER] * nf + [PATH_RED] * npth
    cols, d, g = 10, 0.22, 0.05
    gx = x + (pw_ - (cols * d + (cols - 1) * g)) / 2
    for k, c in enumerate(dots):
        s.oval(gx + (k % cols) * (d + g), vy + 0.1 + (k // cols) * (d + g), d, d, fill=c)
    ly = vy + 0.1 + ((len(dots) - 1) // cols + 1) * (d + g) + 0.1
    T(x + 0.25, ly, pw_ - 0.5, 0.9, [
        {"runs": [(f"{nb} bacteria", {"bold": True, "color": LEAF})], "after": 0},
        {"runs": [(f"{nf} beneficial fungi", {"bold": True, "color": "B5741A"})], "after": 0},
        {"runs": [(f"{npth} fungal pathogens", {"bold": True, "color": PATH_RED})], "after": 0}], size=13, color=INK)
    T(x + 0.25, y0 + 3.6, pw_ - 0.5, 2.7, [
        {"runs": ["Seven biochemical traits, antagonism and compatibility, entered as one Excel sheet."], "after": 6},
        {"runs": [("Blank = not tested; never guessed.", {"bold": True})]}], size=15, color=INK)

    # 2 community: functions + members
    x = x0 + (pw_ + gap)
    funcs = ["P", "K", "N-free", "catalase", "glucanase", "chitinase", "biocontrol"]
    for k, f in enumerate(funcs):
        cx = x + 0.3 + (k % 3) * ((pw_ - 0.6) / 3)
        cy = vy + 0.05 + (k // 3) * 0.5
        s.rect(cx, cy, (pw_ - 0.6) / 3 - 0.08, 0.4, fill=MINT, line="BFE3D0", radius=0.2)
        T(cx, cy, (pw_ - 0.6) / 3 - 0.08, 0.4, f, size=11, color=DEEP, bold=True, align="c", valign="m", min_size=8)
    my = vy + 1.75
    for k, (c, lab) in enumerate([(LEAF, "bact."), (LEAF, "bact."), (AMBER, "fungus")]):
        mx = x + 0.55 + k * ((pw_ - 1.5) / 2)
        s.oval(mx, my, 0.6, 0.6, fill=c)
        T(mx - 0.2, my + 0.62, 1.0, 0.3, lab, size=10.5, color=MUTED, align="c", min_size=8)
    T(x + 0.25, y0 + 3.6, pw_ - 0.5, 2.7, [
        {"runs": ["Two or more members, bacteria and fungi, cover every function and get along."], "after": 6},
        {"runs": [("Biosafety gate", {"bold": True}), " holds flagged strains out of any recommendation."]}], size=15, color=INK)

    # 3 chassis: bars
    x = x0 + 2 * (pw_ + gap)
    s.bar_chart(x + 0.15, vy - 0.15, pw_ - 0.3, vh, [str(v) for v in rk.strain], [float(v) for v in rk.composite], LEAF,
                title="Chassis score", number_format="0", vmin=0, vmax=110)
    T(x + 0.25, y0 + 3.6, pw_ - 0.5, 2.7, [
        {"runs": ["Edit the member the others can spare: fit, dispensability, editing precedent."], "after": 6},
        {"runs": [("Robust", {"bold": True}), " under 2,000 random weightings of the criteria."]}], size=15, color=INK)

    # 4 target: 64 triplets
    x = x0 + 3 * (pw_ + gap)
    cols, d, g = 8, 0.23, 0.06
    gx = x + (pw_ - (cols * d + (cols - 1) * g)) / 2
    for k in range(64):
        hit = k == 45
        s.rect(gx + (k % cols) * (d + g), vy - 0.05 + (k // cols) * (d + g), d, d, fill=AMBER if hit else "D5E4DD", radius=0.03)
    T(x + 0.25, vy + 2.3, pw_ - 0.5, 0.35, "64 three-base motifs tested per isolate", size=11.5, color=MUTED, align="c", min_size=9)
    T(x + 0.25, y0 + 3.6, pw_ - 0.5, 2.7, [
        {"runs": [("No default PAM. ", {"bold": True}), "The native PAM is inferred from the isolate's own spacers."], "after": 6},
        {"runs": ["Introduced editors scanned too: SpCas9 NGG, Cas12a TTTV, dCas9."], "after": 6},
        {"runs": [("TTC", {"bold": True, "color": "B5741A"}), " is the result for the P. polymyxa reference panel."]}], size=14.5, color=INK)

    # 5 build: cassette + DBTL
    x = x0 + 4 * (pw_ + gap)
    segs = [(35, "D0A52F" if False else LEAF), (8, "A3B5AE"), (12, AMBER), (8, "A3B5AE"), (720, "5DBA8A"), (8, "A3B5AE"), (129, BLUE)]
    tot, iw = sum(v for v, _ in segs), pw_ - 0.5
    cx = x + 0.25
    for v, c in segs:
        w = max(iw * v / tot, 0.04)
        s.rect(cx, vy + 0.25, w, 0.5, fill=c)
        cx += w
    T(x + 0.15, vy + 0.85, pw_ - 0.3, 0.3, "promoter · RBS · GFP reporter · terminator", size=10.5, color=MUTED, align="c", min_size=8)
    for k, lab in enumerate(["Design", "Build", "Test", "Learn"]):
        bx = x + 0.2 + (k % 2) * ((pw_ - 0.4) / 2)
        by = vy + 1.25 + (k // 2) * 0.6
        s.rect(bx, by, (pw_ - 0.4) / 2 - 0.1, 0.5, fill=PALE_AMBER if k == 0 else MINT, radius=0.2, line="D5E4DD")
        T(bx, by, (pw_ - 0.4) / 2 - 0.1, 0.5, f"{k + 1}  {lab}", size=12.5, color=DEEP, bold=True, align="c", valign="m", min_size=9)
    T(x + 0.25, y0 + 3.6, pw_ - 0.5, 2.7, [
        {"runs": ["BioBrick parts become a checked reporter or CRISPRi construct."], "after": 6},
        {"runs": [("Tracker", {"bold": True}), " keeps controls and compares predicted with measured."]}], size=15, color=INK)

    # evidence row: what supports it, what it does not show, what is next
    ey = y0 + ph + 0.3
    eh = 4.55
    ew = (PW - 0.8 - 2 * 0.3) / 3
    ex = [0.4 + k * (ew + 0.3) for k in range(3)]
    card(s, ex[0], ey, ew, eh)
    T(ex[0] + 0.25, ey + 0.12, ew - 0.5, 0.4, "Native PAM evidence (reference panel)", size=19, color=LEAF, bold=True)
    s.bar_chart(ex[0] + 0.15, ey + 0.55, ew - 0.3, 2.45, ["Observ.\n23/23", "Contexts\n16/16", "Spacer groups\n11/11", "Random\n35/2100"],
                [100, 100, 100, 1.7], BLUE, title="Sequences with 5′-TTC (%)", number_format="0.0", vmin=0, vmax=125,
                point_colors=[BLUE, BLUE, BLUE, "9FB3AD"])
    T(ex[0] + 0.25, ey + 3.05, ew - 0.5, 1.45, [
        {"runs": [("q = 6.8e-34; ", {"bold": True}), "recovered in every leave-one-out fold (23/23)."], "after": 3},
        {"runs": ["Consistent with published Type I-C work. Not yet measured in P. polymyxa."]}], size=14, color=INK)
    card(s, ex[1], ey, ew, eh)
    s.bar_chart(ex[1] + 0.15, ey + 0.1, ew - 0.3, 3.0, ["Yu\n2024", "GFP\ntransfer", "crisprHAL", "GC\nbaseline", "Permuted"],
                [0.486, 0.4525, 0.6325, 0.3534, -0.0086], LEAF, title="Frozen-prediction validation (Spearman ρ)",
                number_format="0.000", vmin=-0.1, vmax=0.8, point_colors=[LEAF, LEAF, LEAF, "9FB3AD", "9FB3AD"])
    T(ex[1] + 0.25, ey + 3.2, ew - 0.5, 1.3, "Guide ranking tested on three other laboratories' data, predictions frozen before outcomes were joined. "
      "These test guide ranking in other genomes, not any isolate or the TTC PAM.", size=14, color=INK)
    card(s, ex[2], ey, ew, 2.55, fill=PALE_AMBER, line="F2D9A8")
    T(ex[2] + 0.25, ey + 0.1, ew - 0.5, 0.4, "Limits", size=19, color="B5741A", bold=True)
    T(ex[2] + 0.25, ey + 0.55, ew - 0.5, 1.95, "Computational hypotheses, not measured editing efficiency. Every isolate is analysed on its own spacers. "
      "B1 (B. metallica) needs biosafety review. Bank genomes are same-species references. No wet-lab result yet.", size=14, color=INK)
    card(s, ex[2], ey + 2.75, ew, 1.8)
    T(ex[2] + 0.25, ey + 2.8, ew - 0.5, 1.7, [[("Next: ", {"bold": True, "color": LEAF}),
      "sequence the isolates; in P. taiwanensis B56 test two top guides against a low-ranked and a no-guide control; feed results back."]],
      size=14.5, color=INK, valign="m")
    s.rect(0, PHH - 0.5, PW, 0.5, fill=DEEP)
    T(0.6, PHH - 0.5, PW - 1.2, 0.5, "Open source · streamlit run app/streamlit_app.py · 60 automated tests · fixed seed 42 · a 5 Mb genome scanned in 4.8 s · isolates named by species and strain code",
      size=13, color="BFE3D0", valign="m")
    return s


def main(preview_dir=None):
    a = ROOT / "deliverables/final"
    b = ROOT / "deliverables/final_nonblind"
    a.mkdir(parents=True, exist_ok=True); b.mkdir(parents=True, exist_ok=True)
    blind, ban = build(False), build(True)
    build_pptx([blind], a / "AgriPAM-AI_competition_poster.pptx", PW, PH)
    build_pptx([ban], b / "AgriPAM-AI_competition_poster_with_banner.pptx", PW, ban.h)
    problems = []
    if preview_dir:
        p = Path(preview_dir); p.mkdir(parents=True, exist_ok=True)
        problems += render_preview(blind, p / "ga_blind.png", px_per_in=70)
        problems += render_preview(ban, p / "ga_banner.png", px_per_in=70)
    for e in problems:
        print("OVERFLOW", str(e.get("text"))[:80])
    print("built merged poster (blind + banner); overflow items:", len(problems))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)
