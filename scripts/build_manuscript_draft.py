#!/usr/bin/env python3
"""Draft the manuscript from the project's own result files (numbers are read, not typed).

Usage: python scripts/build_manuscript_draft.py
Writes manuscript/Manuscript_draft.md (organism names abbreviated as G. species).
Bracketed items [TO CONFIRM ...] and [EXPERIMENT ...] mark what only the laboratory can supply.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from final_check import source_facts  # noqa: E402

TSV = ROOT / "manuscript/appendices/tsv"
XLS = ROOT / "manuscript/analysis/analysis_results.xlsx"
CUR = ROOT / "manuscript/curation/fungal_curation.tsv"

REFS = {
    "review_engineering": "Engineering the plant microbiome: synthetic community approaches to enhance crop protection. https://pubmed.ncbi.nlm.nih.gov/41704789/",
    "review_strategies": "Strategies for tailoring functional microbial synthetic communities. ISME J 2024;18:wrae049. https://academic.oup.com/ismej/article/18/1/wrae049/7636249",
    "metabolic_design": "Designing a synthetic microbial community through genome metabolic modeling to enhance plant–microbe interaction. Environ Microbiome 2023. https://link.springer.com/article/10.1186/s40793-023-00536-3",
    "smetana": "Zelezniak A, et al. Metabolic dependencies drive species co-occurrence in diverse microbial communities. PNAS 2015. https://www.pnas.org/doi/10.1073/pnas.1421834112",
    "micom": "Diener C, Gibbons SM, Resendis-Antonio O. MICOM: metagenome-scale modeling to infer metabolic interactions in the gut microbiota. mSystems 2020. https://journals.asm.org/doi/10.1128/msystems.00606-19",
    "crispor": "Concordet JP, Haeussler M. CRISPOR: intuitive guide selection for CRISPR/Cas9 genome editing experiments and screens. Nucleic Acids Res 2018. https://pubmed.ncbi.nlm.nih.gov/29762716/",
    "chopchop": "Labun K, et al. CHOPCHOP v3: expanding the CRISPR web toolbox beyond genome editing. Nucleic Acids Res 2019;47:W171–W174. https://academic.oup.com/nar/article/47/W1/W171/5491735",
    "spacer2pam": "Rybnicky GA, et al. Spacer2PAM: a computational framework to guide experimental determination of functional CRISPR-Cas system PAM sequences. Nucleic Acids Res 2022;50:3523. https://academic.oup.com/nar/article/50/6/3523/6544632",
    "leenay": "Leenay RT, et al. Identifying and visualizing functional PAM diversity across CRISPR-Cas systems. Mol Cell 2016. https://www.sciencedirect.com/science/article/pii/S1097276516001751",
    "rubin": "Rubin BE, et al. Species- and site-specific genome editing in complex bacterial communities. Nat Microbiol 2022;7:34–47. https://www.nature.com/articles/s41564-021-01014-7",
    "seva": "SEVA 3.0: an update of the Standard European Vector Architecture for enabling portability of genetic constructs among diverse bacterial hosts. Nucleic Acids Res 2020;48:D1164. https://academic.oup.com/nar/article/48/D1/D1164/5628919",
    "igem": "iGEM Registry of Standard Biological Parts. https://parts.igem.org",
    "typeic_structure": "Structural snapshots of R-loop formation by a type I-C CRISPR Cascade. https://pmc.ncbi.nlm.nih.gov/articles/PMC10026943/",
    "hawkins": "Hawkins JS, et al. Mismatch-CRISPRi reveals the co-varying expression-fitness relationships of essential genes in Escherichia coli and Bacillus subtilis. Cell Syst 2020. https://doi.org/10.1016/j.cels.2020.09.009",
    "yu": "Yu et al. CRISPRi guide-efficiency saturating screen in Escherichia coli. Genome Biol 2024. https://doi.org/10.1186/s13059-023-03153-y (data: GEO GSE196911; Zenodo 10.5281/zenodo.10262866)",
    "crisprhal": "crisprHAL: bacterial Cas9 activity data and model. Nat Commun 2023. https://doi.org/10.1038/s41467-023-41143-7 (data: GitHub tbrowne5/crisprHAL; PRJNA939699)",
    "cctyper": "Russel J, et al. CRISPRCasTyper: automated identification, annotation, and classification of CRISPR-Cas loci. CRISPR J 2020;3:462–469. https://doi.org/10.1089/crispr.2020.0059",
    "defensefinder": "Tesson F, et al. Systematic and quantitative view of the antiviral arsenal of prokaryotes. Nat Commun 2022;13:2561. https://doi.org/10.1038/s41467-022-30269-9",
    "fastani": "Jain C, et al. High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. Nat Commun 2018;9:5114. https://doi.org/10.1038/s41467-018-07641-9",
    "pandas": "McKinney W. Data structures for statistical computing in Python. Proceedings of the 9th Python in Science Conference 2010:56–61.",
    "sklearn": "Pedregosa F, et al. Scikit-learn: machine learning in Python. J Mach Learn Res 2011;12:2825–2830.",
    "iqtree": "Nguyen LT, et al. IQ-TREE: a fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. Mol Biol Evol 2015;32:268–274. https://doi.org/10.1093/molbev/msu300",
}
_order: list[str] = []


def c(*keys: str) -> str:
    nums = []
    for k in keys:
        if k not in _order:
            _order.append(k)
        nums.append(str(_order.index(k) + 1))
    return "[" + ", ".join(nums) + "]"


def abbr(name: str) -> str:
    words = re.sub(r"[^A-Za-z ]", " ", name or "").split()
    return f"{words[0][0]}. {words[1]}" if len(words) >= 2 else name


def facts() -> dict:
    f = dict(source_facts())
    inv = pd.read_csv(TSV / "S1_inventory.tsv", sep="\t")
    f["names"] = {r.acronym: abbr(str(r.identification)) for r in inv.itertuples()}
    f["n_bact"] = int((inv.kind == "bacterium").sum())
    f["n_bfun"] = int((inv.role == "beneficial fungus").sum())
    f["n_pfun"] = int((inv.role == "pathogen fungus").sum())
    s2 = pd.read_csv(TSV / "S2b_biochem_summary.tsv", sep="\t")
    f["s2_top"] = s2.iloc[0]
    f["s2"] = s2
    s4 = pd.read_csv(TSV / "S4_bact_antagonism_summary.tsv", sep="\t")
    f["s4"] = s4
    f["n_tested_path"] = len(s4)
    f["never"] = sorted(set(inv[inv.kind == "bacterium"].acronym) - set(s4.strain))
    s5 = pd.read_csv(TSV / "S5b_fungi_summary.tsv", sep="\t")
    f["s5"] = s5
    s6 = pd.read_csv(TSV / "S6a_bact_vs_fungi.tsv", sep="\t")
    w2 = s6[(s6.reading == "week 2") & (s6.fungus == "FI20")]
    f["fi20_good"] = int(w2.symbol.isin(["++", "+++"]).sum())
    f["fi20_bad"] = sorted(w2[w2.symbol == "-"].strain)
    qc = pd.read_csv(TSV / "S10_completeness_QC.tsv", sep="\t")
    f["qc"] = {r.iloc[0]: r.iloc[1] for _, r in qc.iterrows()}
    f["cur"] = pd.read_csv(CUR, sep="\t").set_index("acronym")
    x = pd.ExcelFile(XLS)
    f["summary"] = x.parse("summary").set_index("scenario")
    for tag in "ABCD":
        rk = x.parse(f"{tag}_ranking")
        f[f"{tag}_members"] = rk.strain.tolist() if tag != "B" else rk.strain.tolist()
        f[f"{tag}_top"] = rk[rk["rank"].notna() & (rk["rank"] != "")].sort_values("rank").iloc[0]
        cov = x.parse(f"{tag}_coverage")
        f[f"{tag}_cov"] = cov
        f[f"{tag}_ranking"] = rk
    return f


def build() -> str:
    f = facts()
    n = f["names"]
    b1, top4 = f["s2_top"], f["s4"].head(4)
    s5 = f["s5"].set_index("fungus")
    A, B, C, D = (f[f"{t}_ranking"] for t in "ABCD")
    flagged = ["B1", "BI1", "CC25", "B34", "B54", "B57"] + [a for a in f["cur"].index if f["cur"].loc[a, "biosafety_hold"] == "yes"]
    ttc, ctx, grp = f["ttc_observations"], f["contexts"], f["spacer_groups"]
    bm = json.loads((ROOT / "manuscript/analysis/software_benchmark.json").read_text())
    ops = {o["operation"]: o for o in bm["operations"]}
    sb = {
        "lib": bm["code_lines"]["library (agripam/)"], "app": bm["code_lines"]["application (app/)"],
        "scripts": bm["code_lines"]["scripts"],
        "scan1": ops["5'-TTC target scan and scoring, 1 Mb genome"]["seconds"], "scan5": ops["5'-TTC target scan and scoring, 5 Mb genome"]["seconds"],
        "cand5": ops["5'-TTC target scan and scoring, 5 Mb genome"]["note"], "spcas": ops["SpCas9 (NGG) target scan, 1 Mb genome"],
        "machine": bm["machine"], "suite": ops["Full test suite"]["note"],
    }
    md = []
    add = md.append

    add("# AgriPAM-AI: open software that takes a microbial bank to a designed community, an editing target and a screened construct, with an application to an agricultural bank\n")
    add("Working draft (software article with an application). Organism names follow the manuscript convention (G. species). Items in brackets [TO CONFIRM …] and [EXPERIMENT …] can only be supplied by the laboratory and must be resolved before submission.\n")

    add("# Abstract\n")
    add(f"**Background.** Laboratories that hold banks of tested isolates lack a tool that links the design of a community to the genome edit that would improve it: a tested bank is not yet a designed community, and adding a function by editing raises the question of which member to edit. "
        f"**Results.** We present AgriPAM-AI, open software (a Python library and a seven-window web application driven by an Excel workbook; {sb['lib']:,} lines of library code and {sb['app']:,} of application code, with an automated test suite) that (i) assembles a community of two or more components, bacteria and fungi, from scored agronomic traits and pairwise compatibility, behind a biosafety gate; "
        f"(ii) ranks each member as an editing chassis by compatibility, dispensability and editing precedent; (iii) identifies editing targets, by inferring a native protospacer-adjacent motif (PAM) from spacer evidence across all 64 three-base motifs and by scanning introduced editors; "
        f"and (iv) assembles a screened standard-parts construct and tracks the experiment, with every output versioned and checksummed. On a laptop the community analysis of a 35-strain bank takes milliseconds and a 5 Mb genome is scanned for candidate targets in {sb['scan5']:.1f} s. "
        f"We applied the software to a bank of {f['n_bact']} bacterial isolates, {f['n_bfun']} beneficial fungi and {f['n_pfun']} fungal pathogens. It showed that the strains leading the functional rankings carry biosafety flags, that the strongest bacterial antagonists inhibit the strongest fungal antagonist, and that a fungal anchor supplies the biocontrol that no eligible bacterium provides. "
        f"For the Type I-C system of P. polymyxa the inference returned 5′-TTC ({ttc} observations; q = {f['q_value']:.1e}), recovered in every leave-one-out fold. "
        f"The guide-ranking components reproduced published outcomes on three deposited data sets (Spearman ρ = {f['yu_rho']}, {f['hk_rho']} and {f['ch_rho']}). "
        f"**Conclusions.** The software turns a laboratory bank into a ranked, auditable plan for what to combine, what to edit and how to test it. All outputs are computational hypotheses. [EXPERIMENT: one sentence on the wet-lab confirmation, if available at submission.]\n")

    add("# Introduction\n")
    add(f"Synthetic microbial communities (SynComs) are increasingly proposed to protect crops and improve their nutrition, and reviews describe disease suppression in several crops together with design principles based on assays, genomes and ecology {c('review_engineering', 'review_strategies')}. "
        f"Computational approaches complement these principles: genome-scale metabolic models predict cross-feeding and select minimal or key-species communities {c('metabolic_design', 'smetana', 'micom')}. "
        f"In practice, however, a laboratory usually holds a bank of isolates characterised by plate assays, not by curated metabolic models, and the questions it faces are practical: which isolates can be combined, which functions are still missing, and which member could carry a new function without disturbing the others.\n")
    add(f"Genome editing offers a way to add a function, but tools for it are built around a single organism that has already been chosen. Guide-design servers rank guides on one genome {c('crispor', 'chopchop')}; "
        f"PAM tools infer or measure the motif of a native system {c('spacer2pam', 'leenay')}; part catalogues and portable vectors supply standard parts for many bacterial hosts {c('igem', 'seva')}; "
        f"and experimental strategies can edit chosen members inside a community {c('rubin')}. One recent review anticipates that integrating machine learning and gene editing will improve SynCom formulation {c('review_engineering')}, which presents that integration as future work. "
        f"In the sources we examined, we did not find a workflow that connects the choice of community, the choice of the member to edit, the design of the edit, and the construct and its follow-up in one auditable chain.\n")
    add(f"Here we describe such a workflow as software, AgriPAM-AI, and apply it to a laboratory bank. Our contributions are: (i) a single open application and library, run from one Excel workbook, that carries the chain from bank to tracker; a community design step that scores bacteria and fungi and respects a biosafety gate; a dispensability-based rule for choosing which member to edit; "
        f"native PAM inference with stability checks and frozen-prediction validation; and a standard-parts construct designer with an experiment tracker; (ii) a reproducibility layer (fixed seeds, checksums, automated tests, a documented workflow for every window); and (iii) an application to a real bank. We state throughout what the results do and do not show.\n")

    add("# Results\n")
    add("## The software\n")
    add(f"AgriPAM-AI is a Python library (agripam/, {sb['lib']:,} lines) with a web application (app/, {sb['app']:,} lines) built on Streamlit, and command-line scripts ({sb['scripts']:,} lines) that regenerate every table and figure. Tabular data use pandas {c('pandas')} and Excel input and output use openpyxl; the learned PAM model uses scikit-learn {c('sklearn')}. "
        f"Figure 1 shows the architecture. The user supplies up to four inputs (a bank workbook, a genome, a parts list and, for validation, a held-out table) and obtains tables, a results workbook, construct files and a provenance record. The application has seven windows, ordered by the sequence of work (Table 1); the interface offers a selectable text size and a collapsible in-app summary of how the software works.\n")
    add("![**Figure 1.** Architecture of AgriPAM-AI: four inputs, seven application windows, four outputs, and the library modules beneath.](figures/Fig1_architecture.png)\n")
    add("**Table 1.** The seven windows: purpose, input and output.\n")
    add("| Window | Purpose | Input | Output |")
    add("|---|---|---|---|")
    for row in [
        ("1 Start here", "Orient the user; show why a PAM alone is not enough; a ten-step project plan", "Choices on screen", "On-screen plan"),
        ("2 Genome evaluation", "One genome from assembly to ranked edit sites and a multi-objective design portfolio", "FASTA, reads or NCBI accession; optional mobile-element reference", "Tables, logs and version record in one ZIP"),
        ("3 SynCom candidate bank", "Community design, gap analysis, chassis ranking; native-bank evidence and PAM reference panels", "Bank workbook (strains, traits, compatibility)", "Community, gaps, chassis ranking, robustness, results workbook"),
        ("4 Agricultural editing knowledgebase", "What is known about editing agricultural microbes; protocol planner", "Search terms; optional reviewed records", "Filtered evidence, protocol"),
        ("5 Parts & constructs", "Assemble and screen a reporter or CRISPRi construct; follow the experiment", "Parts (starter set, Excel, iGEM Registry); neutral region or guide", "FASTA, GenBank, PCR sizes, tracker workbook"),
        ("6 External validation", "Check predictions against other laboratories' data; evaluate a frozen held-out table", "Predictions with outcomes joined afterwards", "Correlations, controls, report"),
        ("7 Software & reproducibility", "Editing-system atlas; reproducibility contract", "None", "Versions, seeds, checksums, commands"),
    ]:
        add("| " + " | ".join(row) + " |")
    add("")
    add("### Inputs and outputs\n")
    add("The bank workbook has one row per strain with its kind (bacterium or fungus), biosafety fields, editing precedent and function scores from 0 to 5, where a blank means not tested and is never imputed; two further sheets hold fungus and bacterium compatibility. A downloadable template and a filled synthetic example are provided. "
        "The results workbook contains the strains, the community, function coverage, gaps, the chassis ranking and the compatibility sheets. Genome runs return one ZIP with every table, log and version record; construct runs return FASTA and GenBank files and a tracker.\n")
    add("### Which PAMs the software searches\n")
    add("The software does not assume a PAM. For a native system, the discovery step tests all 64 three-base motifs against the protospacers matched by the organism's own spacers, corrects for multiple testing, and reports only motifs with enough observations; it never falls back to a default motif. For introduced editors it scans fixed, published motifs: SpCas9 (NGG), Cas12a (TTTV) and dCas9 for CRISPRi (NGG). "
        "The 5′-TTC motif reported below is the output of the native inference for the P. polymyxa reference panel, not a property assumed for every isolate; the dedicated 5′-TTC designer and scan are conveniences for that result, and a different isolate must be analysed on its own spacer evidence. Isolates with no detected CRISPR-Cas system receive only the introduced-editor scans.\n")
    add("### Reproducibility, testing and performance\n")
    add(f"Random procedures use a fixed seed (42); outputs carry SHA-256 checksums and version records; a reproducibility contract in window 7 lists the commands that rebuild every result. The automated suite ({sb['suite'].lower()}) covers the scoring rules, file parsers, construct screens and target scanning; the interface layer is tested less thoroughly, and one optional interface test is known to fail under the current Streamlit test harness. "
        f"Timings on {sb['machine']} were: reading the 35-strain bank and assembling a community with a 2,000-weighting robustness check, about 0.01 s; assembling and screening a reporter construct with 750-bp arms, about 0.002 s; the leave-one-group-out PAM stability check, about 1 s; scanning a 1 Mb genome for 5′-TTC targets, {sb['scan1']:.1f} s, and a 5 Mb genome, {sb['scan5']:.1f} s ({sb['cand5']}); and a SpCas9 (NGG) scan of 1 Mb, {sb['spcas']['seconds']:.1f} s ({sb['spcas']['note']}). "
        f"Time grows roughly linearly with genome length, so the interactive scans are practical for bacterial genomes.\n")
    add("### Design choices\n")
    add("Three choices shape the software. It never imputes missing assays, so a gap is visible as a gap. It puts a biosafety gate before any ranking, so a flagged strain can be a community member but never a recommended chassis. "
        "And every reference panel is labelled with what it does and does not show, so that a fixed reference study is not mistaken for a result on an uploaded strain.\n")
    add("## Application to an agricultural bank\n")
    add("The remaining sections apply the software to a laboratory bank, following the order of the windows.\n")
    add("## The bank and its assay layers\n")
    add(f"The bank comprises {f['n_bact']} bacterial isolates, {f['n_bfun']} beneficial fungi (FI1, FI2, FI20, FI62) and {f['n_pfun']} fungal pathogens used as plate targets (Appendix S1). "
        f"Four assay layers were recorded: seven biochemical tests for the bacteria (Appendix S2), antagonism of bacteria and of beneficial fungi against the pathogens (S3 to S5), "
        f"compatibility of bacteria with the beneficial fungi and with each other (S6, S7), and fungus–fungus interactions (S8). Symbols were converted to 0–5 scores; blank cells were treated as not evaluated. "
        f"Coverage is uneven (Appendix S10): {f['qc']['Bacteria vs pathogenic fungi']}; the bacterium–bacterium matrix uses only five base strains; and fungus–fungus compatibility was tested only against FI20. "
        f"{len(f['never'])} of {f['n_bact']} bacteria were never tested against pathogens ({', '.join(f['never'])}).\n")
    add("## Functional and antagonism profiles\n")
    add(f"By total biochemical score, {b1.strain} ({n.get(b1.strain, '')}) ranked first (total {int(b1.total_score)}, mean intensity {b1.mean_intensity}, positive in {int(b1['positive_function_groups (of 7)'])} of 7 function groups), followed by "
        + ", ".join(f"{r.strain} ({int(r.total_score)})" for r in f['s2'].iloc[1:4].itertuples()) + " (Table 2). "
        f"Against fungal pathogens, {top4.iloc[0].strain} was the broadest bacterial antagonist (mean {top4.iloc[0].mean_intensity_0_4} on a 0–4 scale over {int(top4.iloc[0].tests_evaluated)} tests), followed by "
        + ", ".join(f"{r.strain} ({r.mean_intensity_0_4})" for r in top4.iloc[1:].itertuples()) + ". "
        f"Among fungi, FI20 gave the maximum score against every pathogen in every reading (mean {s5.loc['FI20', 'mean_intensity_0_4']}; {int(s5.loc['FI20', 'tests_evaluated'])} tests), ahead of FI62 ({s5.loc['FI62', 'mean_intensity_0_4']}) and FI2 ({s5.loc['FI2', 'mean_intensity_0_4']}); FI1 has one reading only and is not comparable.\n")
    add("## Compatibility constrains what can be combined\n")
    add(f"With FI20, {f['fi20_good']} bacteria grew together or potentiated the fungus, whereas {', '.join(f['fi20_bad'])} inhibited it (Appendix S6). FI62 inhibited FI20 with a large, well-defined halo (S8). "
        f"Plate remarks recorded further bacterium–bacterium inhibitions, for example of B3 by B2, B34, B54 and CC26 (S7). "
        f"Thus the best individuals do not make a single community: the two strongest enzyme-linked bacterial antagonists inhibit the strongest fungal antagonist, and the two best fungal antagonists inhibit each other.\n")
    add("## Community design by the software\n")
    add(f"We entered the bank in the software's input format and ran four scenarios (Table 3). Strains carrying a biosafety flag from the project's screen ({', '.join(flagged)}) are never ranked as chassis. "
        f"With FI20 as a scored member (scenario A), the software assembled {', '.join(f['A_members'])}, which delivered every function of the scale; without a fungus (scenario B) the best eligible bacteria, {', '.join(f['B_members'])}, left biocontrol breadth undelivered. "
        f"The lab-planned community (B1, B39, B56 with FI20; scenario C) lacked β-glucanase, and B1 was excluded as a chassis by the gate. Offered two fungal anchors that inhibit each other (scenario D, FI20 and FI62), the software raised a warning.\n")
    add("## Choosing the member to edit\n")
    a_top, c_top = f["A_top"], f["C_top"]

    def why(row):
        return (f"composite {row.composite}, dispensability {row.dispensability_score}, editing precedent {int(row.editability_precedent_score)}"
                + ("; the rest of the community covers all of its functions" if float(row.dispensability_score) == 100.0 else ""))
    cur = f["cur"]
    add(f"For each community the software ranks members by compatibility with the others, dispensability (the share of the community's functions that survive if the member's own activity is lost) and editing precedent. "
        f"In scenario A the top-ranked chassis was {a_top.strain} ({why(a_top)}); in scenario C it was {c_top.strain} ({why(c_top)}). "
        f"Fungi are ranked with the same rules. Their editing precedent and biosafety were curated from the literature: genus-level precedent exists for FI1, FI2 and FI20 (CRISPR/Cas9 systems published for other species of [taxon withheld], [taxon withheld] and [taxon withheld]), none for the unidentified FI62; "
        f"FI2 and FI62 are held (clinical reports for the species; unidentified isolate), and FI1 and FI20 carry a caution because their genera contain mycotoxigenic or opportunistic species. This curation is decision support and needs biosafety-committee confirmation. "
        f"Under 2,000 random weightings of the three criteria, the outcome was: scenario A, {f['summary'].loc['A', 'robustness']}; scenario C, {f['summary'].loc['C', 'robustness']}. "
        f"Scenario C contains B1, which carries a biosafety flag: the coverage it provides describes the planned community, not a recommendation to use it.\n")
    add("![**Figure 2.** Community design and chassis ranking in scenarios A and C. Heat map: function scores (0 to 5); bars: composite chassis score.](figures/Fig3_community_design.png)\n")
    add("## Editing targets: an inferred native PAM and introduced editors\n")
    add(f"In the P. polymyxa reference panel (24 genomes, 10 with CRISPR-Cas), the orientation-normalised motif 5′-TTC was found next to protospacers in {ttc} accession-coordinate observations, {ctx}/{ctx} unique contexts and {grp}/{grp} spacer groups, against {f['background']} matched random loci "
        f"(Benjamini–Hochberg q = {f['q_value']:.2e}). A PAM-plus-mobile-context model reached ROC AUC {f['auc']} and average precision {f['ap']:.2f} under grouped cross-validation at 0.99% positive prevalence; because only 15 positive contexts exist, this is a ranking model, not an efficiency estimate. "
        f"Re-inferring the PAM with each of {f['heldout_folds'][0]} spacer groups, and separately each of {f['heldout_folds'][1]} genomes, left out recovered TTC in every fold ({f['heldout'][0]} held-out observations in both schemes). "
        f"The motif agrees with a Type I-C system that gave an NTTC consensus in a functional screen {c('leenay')} and with structural work describing recognition of a 5′-TTC PAM by Cas8c {c('typeic_structure')}. "
        f"The 23 observations are not independent, and the check uses the same pipeline; it establishes stability, not validation. Introduced editors (SpCas9, Cas12a, dCas9) are scanned on both strands and scored on editability, deliverability, agronomic value and preservation safety, with a bounded one-to-two-mismatch off-target screen.\n")
    add("![**Figure 3.** Evidence for the 5′-TTC PAM.](figures/Fig4_pam_evidence.png)\n")
    add("## Independent validation of the guide-ranking components\n")
    add(f"Predictions were frozen before outcomes were joined. On a published CRISPRi screen in E. coli ({f['yu_n']} guides) {c('yu')} the pooled Spearman correlation was {f['yu_rho']} (within-gene median {f['yu_median']}; permutation control {f['yu_null']}). "
        f"Transfer of E. coli activity to B. subtilis for {f['hk_n']} sequence-matched guides gave ρ = {f['hk_rho']} (95% interval {f['hk_lo']}–{f['hk_hi']}), i.e. partial transfer {c('hawkins')}. "
        f"A sequence model fitted on the released training guides of a bacterial Cas9 data set reached ρ = {f['ch_rho']} on {f['ch_n']:,} held-out guides, against {f['ch_base']} for a GC-only baseline {c('crisprhal')}. "
        f"These benchmarks test guide ranking in other genomes; they do not test the TTC PAM or any member of the bank.\n")
    add("![**Figure 4.** Frozen-prediction benchmarks with their controls.](figures/Fig5_validation.png)\n")
    add("## From parts to bench\n")
    add(f"The construct designer assembles either a reporter insertion between genome-derived homology arms or a CRISPRi guide cassette from standard parts (starter set, user sheet or the iGEM Registry {c('igem')}). "
        f"Junctions use the BioBrick scar, except inside the guide RNA. Each part and the whole construct are screened for forbidden sites (RFC10, RFC25, Golden Gate), for the chassis's own restriction motifs, homopolymers, GC extremes, repeats and declared host range. "
        f"As an example, a reporter cassette built from four Registry parts (J23100, B0034, E0040, B0015) is 920 bp including three 8-bp scars; registry parts are characterised mainly in E. coli, so activity in the chosen chassis must be measured. "
        f"A tracker records six stages per construct with parental, negative and low-ranked controls and compares predicted with measured effects.\n")
    add("![**Figure 5.** A reporter cassette assembled from four Registry parts.](figures/Fig6_cassette.png)\n")
    add("## Experimental confirmation\n")
    add("[EXPERIMENT: replace this section with the results of the verification plan: (A) PAM dependence in the native Type I-C strain, (B) reporter insertion in the top-ranked member, and (C) CRISPRi guide ranking, with replicates, controls and effect sizes. If not available, delete this heading and state in the Discussion that the results are computational.]\n")

    add("# Discussion\n")
    add("The workflow changes the order in which decisions are made. Instead of choosing an organism and then designing guides, it starts from what the community needs and lets compatibility, biosafety and dispensability decide what to edit. "
        "In our bank this mattered: the strains with the highest functional and antagonism scores were flagged, and the best antagonists could not be combined, so a naive ranking by strength would have proposed a community that fails on compatibility or safety.\n")
    add("Compared with metabolic-model approaches, ours does not predict cross-feeding and cannot anticipate interactions the plates did not test; compared with guide-design servers, it adds the choice of organism and the construct. "
        "The two are complementary: metabolic models could supply interaction scores where plate data are missing, and guide servers could refine the per-target scores.\n")
    add(f"Software limitations are that it runs locally for one user, that the full genome workflow depends on external tools the user must install, that the interface is less thoroughly tested than the library, and that the community scores depend on the quality of the input workbook. Scientific limitations are substantial. Compatibility rests on semi-quantitative plate readings with sparse bacterium–bacterium coverage; the inferred PAM rests on {ttc.split('/')[0]} non-independent observations and has not been measured; "
        "genome-based editability of the isolates is not yet computed and the reference accessions are same-species genomes; fungal precedent and biosafety are literature-based and no species-level editing report exists for the fungal isolates; and registry parts are characterised mainly in E. coli. "
        "Several top-ranked isolates belong to a genus that includes opportunistic pathogens, and any use requires a biosafety review.\n")

    add("# Methods\n")
    add("## Bank and assays\n")
    add("Bacterial biochemistry (catalase, oxidase, phosphate and potassium solubilisation, growth on N-free media, β-glucanase and chitinase) was read at 24 h, 72 h and 7 d; antagonism and compatibility were read on dual-culture plates at one to three weeks. "
        "[TO CONFIRM: exact media, incubation temperature and time, inoculum, and number of replicates for each assay.] Symbols were converted to numbers ('-' = 0 to '+++++' = 5; antagonism 0–4), and blank cells were not scored. "
        "The scoring reproduces the earlier layered analysis exactly on the strains compared (Appendix S2, S4, S5).\n")
    add("## Community scoring\n")
    add("A function is delivered at a score of 2 or more. Fit is the mean compatibility with the other members over the tests that exist (no imputation). Dispensability is the share of the community's delivered functions that survive if the member's own activity is lost. "
        "Editing precedent is 100, 60 or 20 for species-level, genus-level or none. The composite is 0.40 fit + 0.35 dispensability + 0.25 precedent; strains with a biosafety flag are not ranked. "
        "Robustness was checked under 2,000 random weightings (seed 42). Communities are assembled by greedy set cover among eligible strains, excluding strains inhibited by a chosen member and requiring evidence of compatibility with the fungal anchor.\n")
    add("## Native PAM inference and model\n")
    add(f"CRISPR arrays and Cas operons were annotated with CRISPRCasTyper {c('cctyper')} and MinCED; defense systems with DefenseFinder {c('defensefinder')}; species boundaries were checked by FastANI {c('fastani')} and the phylogeny by Parsnp and IQ-TREE {c('iqtree')}. "
        "Spacers were matched exactly to phage, plasmid and genome sequences, orientation-normalised and tested against 2,100 matched random loci with Benjamini–Hochberg correction over the 64 triplets. "
        "The model is a class-weighted logistic regression on the three PAM bases and a mobile-element flag (liblinear, seed 42) evaluated by grouped five-fold cross-validation. Held-out stability re-infers the PAM with each spacer group or genome withheld.\n")
    add("## Validation protocol\n")
    add("Held-out predictions were generated without the outcome column, frozen, and only then joined to outcomes by guide identifier; permutation controls and 2,000-resample bootstraps accompany each benchmark (docs/EXTERNAL_VALIDATION.md).\n")
    add("## Constructs and screens\n")
    add("Parts are assembled with the 8-bp BioBrick scar (or none), homology arms default to 750 bp, and diagnostic-PCR sizes are computed for primers 150 bp outside the arms. Forbidden-site sets are RFC10 (EcoRI, XbaI, SpeI, PstI, NotI), RFC25 (adds AgeI, NgoMIV) and Golden Gate (BsaI, BsmBI).\n")
    add("## Implementation\n")
    add("The library is written in Python (3.9 or later) and organised by role: syncom and syncom_io (community scoring and workbook input/output), core and workflow (genome analysis), editor_targeting and pam_discovery (targets), multiobjective (portfolio), validation (held-out evaluation), and biobricks and parts_io (parts and constructs). The web application (Streamlit) imports the library and adds no scoring logic of its own, so every number on screen can be reproduced from the command line. "
        "External specialist tools (CRISPRCasTyper, DefenseFinder, FastANI, Parsnp, IQ-TREE) are optional: when one is not installed, its module reports 'not analysed' rather than returning a partial result. A container recipe is provided for a fixed environment. The optional iGEM Registry fetch needs internet access.\n")
    add("## Software, data and code availability\n")
    add("The software (Python, Streamlit) and all tables are released with a versioned archive and SHA-256 checksums (the application runs locally, single-user) [TO CONFIRM: licence and repository]. Appendix S1 to S10 and the software input workbook accompany the paper. Random procedures use seed 42.\n")

    add("# Tables and figure legends\n")
    add("**Table 2.** Top ten bacterial isolates by total biochemical score (Appendix S2).\n")
    add("| Rank | Strain | Species | Total score | Mean intensity | Function groups (of 7) |")
    add("|---|---|---|---|---|---|")
    for _, r in f["s2"].head(10).iterrows():
        add(f"| {int(r['rank'])} | {r['strain']} | {n.get(r['strain'], '')} | {int(r['total_score'])} | {r['mean_intensity']} | {int(r['positive_function_groups (of 7)'])} |")
    add("")
    add("**Table 3.** Communities assembled by the software (Appendix workbook 'analysis results').\n")
    add("| Scenario | Members | Functions missing | Top chassis | Not ranked (biosafety flag) |")
    add("|---|---|---|---|---|")
    for tag, row in f["summary"].iterrows():
        add(f"| {tag}. {row.description} | {row.members} | {row.functions_missing} | {row.top_chassis} | {row.not_ranked_biosafety} |")
    add("")
    add("Screenshots of each window at full resolution should be added by the team as Supplementary Figures S1 to S7 [TO CONFIRM]. [EXPERIMENT: Figure 6, wet-lab results.]\n")

    add("# Declarations\n")
    add("**Ethics and biosafety.** [TO CONFIRM: institutional biosafety review for the [taxon withheld] cepacia complex strains and any regulatory requirements.] **Author contributions (CRediT).** [TO CONFIRM.] **Competing interests.** [TO CONFIRM.] **Funding.** [TO CONFIRM.]\n")

    add("# References\n")
    for i, key in enumerate(_order, 1):
        add(f"{i}. {REFS[key]}")
    return "\n".join(md) + "\n"


if __name__ == "__main__":
    out = ROOT / "manuscript/Manuscript_draft.md"
    text = build()
    out.write_text(text)
    open_items = len(re.findall(r"\[(TO CONFIRM|EXPERIMENT)", text))
    print(f"wrote {out} ({len(text.split())} words, {open_items} open items, {len(_order)} references)")
