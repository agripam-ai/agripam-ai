# How AgriPAM-AI / RhizoForge-Select works, step by step

This guide follows the software from installation to a final ranked recommendation. Every number in the worked examples was produced by the code in this repository. Isolates are identified by strain code and species (for example B56, P. taiwanensis; FI20, T. yunnanense).

## 0. What the software does, and what it does not

RhizoForge-Select answers three practical questions for an agricultural microbial community:

1. **Which member of a tested community should be genome-edited?** (SynCom chassis selector)
2. **Which editing sites in that member's genome are compatible with its defenses and its likely PAM?** (Genome workspace)
3. **How strong is the evidence behind each answer?** (Evidence, validation and reproducibility tabs)

It produces **ranked, auditable computational hypotheses**. It does not measure editing efficiency, and it never turns a missing analysis into a negative biological finding. Every score is a decision-support heuristic. Physical validation of the inferred PAM and of any edit is still required.

## 1. Setup (2 minutes)

```bash
./scripts/install_app.zsh      # lightweight app environment, project-local
./scripts/run_app.zsh          # opens http://localhost:8501
```

The lightweight environment runs every module that does not need specialist software and marks the rest as "unavailable" in the module-status table. For the complete workflow (gene annotation, CRISPRCasTyper, DefenseFinder, MOB-suite, geNomad), create the full environment with `./scripts/install_full_pipeline.zsh` and start it with `./scripts/run_full_app.zsh --server.port 8502`.

Check the installation at any time:

```bash
python -m unittest discover -s tests     # 60 tests, 1 optional app test skipped
```

## 2. The app at a glance

| Tab | Purpose | Uses |
|---|---|---|
| **Start here** | Interactive learning model and a project planner that lists what a defensible editing project needs | No data required |
| **Genome evaluation** | Analyse one genome end to end (steps 3 to 9 below) | Uploaded FASTA, NCBI accession, raw reads or the built-in demo |
| **SynCom candidate bank** | Native-bank evidence, the SynCom chassis selector (step 10), and the fixed *P. polymyxa* reference results | Laboratory bank tables in `data/syncom/`, precomputed results |
| **Agricultural editing knowledgebase** | Searchable catalogue of strains with published editing routes, PAM requirements and validation status | `data/agricultural_editing_knowledgebase.tsv` |
| **Parts & constructs** | Choose BioBrick parts, build a checked construct for the chosen edit, and follow the experiment | Starter parts, an Excel sheet or the iGEM Registry; the genome from the Genome evaluation tab or an uploaded FASTA |
| **External validation** | Frozen-prediction benchmarks against deposited experiments, and an upload box for your own held-out table | `data/external_validation/` |
| **Software & reproducibility** | Editing-system atlas, provenance and the reproducibility contract | Project configuration |

The sidebar has a **Text size** selector (Normal, Large, Extra large; Large by default). Tables follow the browser zoom (Ctrl or Cmd and +).

Reference results (the *P. polymyxa* PAM model, TTC designer and evidence tables) are always labelled as precomputed reference results. They are never presented as output for an uploaded genome.

## 3. Genome workspace, step 1: provide a genome

Choose one input on the **Genome evaluation** tab:

- **Upload FASTA:** a multi-contig assembly (`.fa`, `.fasta`, `.fna`).
- **NCBI assembly accession:** a versioned `GCF_...` or `GCA_...` accession; the genome and its organism name are retrieved with the NCBI Datasets command (internet needed).
- **Built-in demonstration:** a tiny offline sequence for presentations.
- **Raw reads:** paired Illumina, long-read or hybrid FASTQ, assembled with SPAdes, Flye or Unicycler when installed.

Accessions are validated before any command runs, so only well-formed assembly identifiers reach the shell.

## 4. Step 2: genome quality control

The sequence is normalised (upper case, IUPAC symbols validated) and summarised: contig count, total length, GC content, ambiguous-base burden, N50 and L50, longest and shortest contig, and a SHA-256 checksum that ties every later result to this exact sequence. A genome that fails a quality threshold is marked "review" rather than silently analysed.

Worked example, built-in demonstration sequence: 1 contig, 53 bp, GC 39.6%, no ambiguous bases, SHA-256 `274e29c4…`.

## 5. Step 3: annotation and screens

The workflow then runs, in order, and records each as a stage with a status of `complete`, `limited`, `unavailable`, `failed` or `not_run`:

| Stage | Engine | If the engine is missing |
|---|---|---|
| Coding sequences, proteins, GFF | Bakta, then Prokka | Prodigal genes without product names (status `limited`) |
| Plasmid prediction | MOB-suite | labelled contig-level heuristic |
| CRISPR-Cas systems | CRISPRCasTyper | `unavailable`; absence is never reported |
| Defense systems | DefenseFinder | `unavailable`; absence is never reported |
| Restriction-modification, mobile elements, agricultural functions, alternative editing systems | Screens of the annotated products (geNomad for mobile elements) | `limited` without product names |
| Candidate neutral regions | Exclusion windows around genes and functional, defense and mobile features | `unavailable` without coordinates |

Agricultural-function hits (antifungal compounds, lipopeptide synthesis, hydrolytic enzymes, phosphate solubilisation, nitrogen metabolism, siderophores, root colonisation) are candidate associations. Neutral regions are not proven safe harbors.

## 6. Step 4: find target sites

Two independent scans run on both strands:

- **Inferred native Type I-C targets (5′-TTC).** Each protospacer is scored (worked example below).
- **Introduced editors:** SpCas9 (`NGG`, 3′ PAM, 20 nt), Cas12a (`TTTV`, 5′ PAM, 23 nt) and dCas9/CRISPRi (`NGG`, 20 nt). Sites use IUPAC matching and are reported in 1-based genome coordinates.

For each introduced editor the software counts exact copies genome-wide and runs a **bounded off-target screen**: for the top 2,000 unique guides, every site with one or two substitutions is counted with a three-seed guarantee. Bulges, non-canonical PAMs and nuclease-specific mismatch tolerance are outside this screen and are labelled unresolved.

### Worked example: the 5′-TTC design score

For the demonstration sequence the scanner finds one candidate, `AGP000001`, at positions 10 to 44 on the plus strand (protospacer `AACTAAATAAACAACAAAGGACTCCATACTGGTTG`, PAM `TTC`).

```
design_score = 0.60 × PAM model probability
             + 0.20 × GC quality        (1 − |GC% − 50| / 50)
             + 0.10 / exact genome copies
             + 0.10 × sequence quality  (penalises homopolymers longer than 4)
```

| Quantity | Value |
|---|---|
| PAM model probability for TTC, no mobile context | 81.27% |
| GC content of the protospacer | 34.29% |
| Exact genome copies | 1 |
| Longest homopolymer | 3 |
| **Design score** | **82.47** |
| Design score if independent evidence puts the target in a mobile-element region (PAM probability 97.20%) | 92.03 |

The PAM probability comes from the 64-triplet table `machine_learning/type_ic_all_64_pam_scores.tsv`. It is a ranking score, not a cleavage probability.

## 7. Step 5: link each target to its gene

`pam_gene_context.tsv` maps every target to its overlapping or nearest coding sequence and classifies it as within a CDS, promoter-proximal or intergenic. Conservative flags mark annotations related to growth or to beneficial agronomic functions. These flags are screening aids, not essentiality calls.

## 8. Step 6: strain-specific PAM discovery (optional)

If the user supplies a versioned phage, plasmid or mobile-element FASTA together with the host genome, the software:

1. extracts the host's CRISPR spacers from the CRISPRCasTyper output;
2. searches both orientations of the reference for exact full-length protospacer matches;
3. reads the normalised 3-base sequence immediately 5′ of each match;
4. ranks all 64 triplets by count, background enrichment, Wilson interval, nominal P and Benjamini-Hochberg Q.

Only a motif supported by at least three observations and Q < 0.10 is passed on to the gene-context step. If evidence is absent, the software says so; it never substitutes TTC silently.

## 9. Step 7: multi-objective portfolio

Candidates are scored on four objectives, each between 0 and 100:

| Objective | Built from |
|---|---|
| Editability | PAM support (50%), guide uniqueness (35%), repair support (15%) |
| Deliverability | Penalties for restriction-modification hits and for detected defenses |
| Agronomic value | Whether the target relates to a beneficial function or promoter |
| Preservation safety | Penalties for growth-related or beneficial-function annotations, CDS disruption and repeated protospacers |

A candidate is rejected if safety is below 35 or editability is below 25. The rest are compared on the **Pareto front** (no other candidate is at least as good on all four objectives and better on one). The output names a best predicted design, a biologically safer alternative and an easier experimental alternative, and writes a blank experimental-results ledger so measured outcomes can be joined later.

## 10. Step 8: which member of the community should be edited?

This is the community-level layer (`agripam/syncom.py`, `scripts/run_syncom_selector.py`), shown on the **SynCom candidate bank** tab.

**Bring your own data.** On the **SynCom candidate bank** tab, open "Analyze your own bank". Download the Excel template (or the filled synthetic example), enter one row per strain, upload it, choose the fungal partner(s), and press **Run SynCom analysis**. The workbook has four sheets: `Strains` (code, optional identification, `biosafety_hold` yes/no, `pathogen_screen_score`, `editing_precedent` species/genus/none, then any number of function columns scored 0 to 5, blank = not tested), `Fungi_compatibility` (`+++`, `++`, `+`, `-+`, `-`), `Bacteria_compatibility` (`mixes` or `inhibits`) and an optional `Community` sheet. Rows starting with EXAMPLE are ignored. The software validates every cell, reports errors and warnings in plain language, and returns the community, function coverage, gaps, the chassis ranking and a robustness check, with a downloadable results workbook. Blank biosafety cells are flagged and treated as eligible, so assess biosafety before any use.

**Input for the laboratory example.** `data/syncom/` holds the laboratory tables built from the tested bank: seven biochemical traits for 30 bacterial isolates (0 to 5), antagonism against fungal pathogens, compatibility with beneficial fungi, and pairwise bacterial compatibility. Blank cells are missing data, never negatives.

**Fungi as members.** A strain with `kind = fungus` and function scores (for example biocontrol breadth) is a scored community member: its functions count toward the community, it is ranked as a possible chassis like any bacterium, and its compatibility is checked against the bacteria (`Fungi_compatibility`) and against other fungi (`Fungus_compatibility`). Two fungi that inhibit each other raise a warning; a bacterium inhibited by a fungal anchor is left out of the assembled community. Fungi listed as partners but absent from the Strains sheet are compatibility partners only, as before.

**Procedure.**

1. **Assemble** (`assemble`): greedy set cover. Add the biosafety-eligible strain, compatible with the community fungus, that delivers the most new functions. A function counts as delivered at a score of 2 or more.
2. **Find gaps** (`gap_analysis`): list functions the community lacks, the bank strains that carry them, and why a strain cannot simply join (inhibits the community fungus, or fails the safety gate).
3. **Rank chassis** (`rank_chassis`): for each member:
   - *fit*: mean compatibility with the other members over the tests that exist;
   - *dispensability*: the share of the community's delivered functions that survive if that member's own activity is lost;
   - *editing precedent*: 100 for species-level, 60 for genus-level, 20 for none;
   - composite = 0.40 × fit + 0.35 × dispensability + 0.25 × precedent.
4. **Safety gate**: a strain marked as a biosafety hold in the project cohort table, or with a pathogen-screen score of 30 or more, is never ranked.
5. **Check robustness** (`weight_sensitivity`): re-rank under 2,000 random weightings (seed 42) and report how often each strain stays first.

**Worked result (planned community FI20 with B1, B39, B56).** B1 is excluded by the safety gate. B56 is the top bacterium (composite 90, tied with the fungal anchor FI20) because the other members already cover its functions and its genus has editing precedent. B39 scores 75 (dispensability 85.7): it is the only strong chitinase producer, so editing it puts that function at risk. B56 stays at or above B39 under all 2,000 random weightings. The community lacks β-glucanase; the strongest producers inhibit FI20 or are flagged, and compatible producers (B21, B28, B3, B20) exist, so the tool recommends adding one of those.

## 11. Step 9: from parts to bench (Parts & constructs tab)

`agripam/biobricks.py`, `agripam/parts_io.py` and `agripam/constructs_ui.py` connect the chosen edit to a DNA construct and to the experiment that follows.

1. **Parts library.** Load the starter set (four iGEM Registry parts with their provenance, plus an sgRNA scaffold to verify), upload the parts workbook (`data/biobricks/Parts_input_template.xlsx`: `part_id`, name, role, sequence, hosts, source), or fetch parts by ID (for example `BBa_J23100`) from the iGEM Registry. IDs are validated before any request is made, and a failed request is reported without stopping the app. Edit each part's role and declared host range in the table.
2. **Choose the edit.** *Reporter insertion at a neutral site:* pick a candidate neutral region from the genome workflow (or type a contig and position), the arm length and the orientation, then a promoter, optional RBS, reporter and optional terminator. The construct is left arm, cassette, right arm. *CRISPRi guide cassette:* pick a dCas9 site from the workflow or paste a 20-nt spacer, then a promoter, an sgRNA scaffold and an optional terminator. Promoter, spacer and scaffold are always joined without a scar, because a scar inside the guide RNA would break it.
3. **Screen.** Each part and the whole construct are checked for forbidden internal sites (RFC10: EcoRI, XbaI, SpeI, PstI, NotI; RFC25; Golden Gate BsaI/BsmBI), for the chassis restriction motifs you enter (IUPAC, both strands), for long homopolymers, extreme GC and repeated 24-mers, and for whether the part's declared host range includes your chassis.
4. **Outputs.** An annotated map, FASTA and GenBank files (they open in common sequence editors), and, for insertions, the expected diagnostic PCR sizes for wild type and edited alleles.
5. **Follow the output.** The tracker holds one row per construct plus the controls the experiment needs (parental, negative, and a low-ranked guide for CRISPRi). Fill the stage columns (designed, assembled, sequence confirmed, delivered, edit confirmed, phenotype measured), the assay, replicates, measured value, units and normalised effect. Upload the file again to see the progress chart, warnings about missing controls or fewer than three replicates, and the Spearman agreement between the predicted score and the measured effect (five or more constructs needed).

A part that passes every screen is compatible, not proven: activity in the chosen chassis has to be measured, and the module says so in its output. Registry data come mostly from *E. coli*.

## 12. Step 10: read the evidence and validation tabs

- **Reference study.** In 24 *P. polymyxa* genomes, 5′-TTC appears in 23/23 accession-coordinate observations, 16/16 unique contexts and 11/11 spacer groups, against 35/2,100 matched background loci (Q = 6.755e-34). The grouped PAM-plus-mobile model reaches ROC AUC 0.9946 and average precision 0.5200 on 0.99% positive prevalence.
- **Held-out consistency check** (`scripts/heldout_pam_validation.py`): re-inferring the PAM with each spacer group (11) or each target genome (15) withheld recovers TTC in every fold and matches 23/23 held-out observations. This is a stability check in the same pipeline, not independent validation.
- **External validation.** Predictions are frozen before outcomes are joined. GSE196911: 742 guides, pooled Spearman ρ 0.4860, within-gene median 0.5459, permutation control −0.0086. Cross-chassis GFP benchmark: 33 matched guides, ρ 0.4525, guide-bootstrap 95% interval 0.0876 to 0.7223. crisprHAL released holdout: ρ 0.6325 against 0.3534 for a GC-only baseline. An unseen-genome challenge deliberately abstains outside the reference model's scope.
- **You can test your own model:** upload a held-out table on the External validation tab; the software computes Spearman ρ, ranking utility, calibration, bootstrap intervals and a permutation control.

## 12. What you get out

The genome workflow exports one ZIP: normalised FASTA, GFF and protein FASTA, module folders, logs, target and gene-context tables, the off-target screen, the portfolio and results ledger, `pipeline_status.tsv`, an uncertainty report and a JSON provenance manifest with software versions, parameters and checksums. The SynCom selector writes `chassis_ranking.tsv`, `gap_analysis.tsv` and `summary.json` to `results/syncom_selector/`.

## 13. Reproduce the headline results

```bash
python -m unittest discover -s tests                 # deterministic checks
python scripts/build_competition_demo_report.py      # results/competition_demo/
python scripts/heldout_pam_validation.py             # results/heldout_pam_validation/
python scripts/annotate_syncom_metadata.py           # derived safety/precedent fields
python scripts/run_syncom_selector.py                # results/syncom_selector/
```

`scripts/build_syncom_bank.py` rebuilds the bank tables from the laboratory workbook (needs `openpyxl`). All random procedures use the fixed seed 42.

## 14. Interpretation boundaries

- The inferred TTC PAM, every target, and every chassis ranking are computational hypotheses until measured.
- Scores are not editing efficiencies, off-target absence, safety certificates or field-performance predictions.
- SynCom compatibility rests on few plate tests; the number of tests behind each fit score is reported (`fit_evidence_n`).
- Genome accessions used for the bank are same-species references, not the isolates. Sequencing the isolates is the next step.
- A missing specialist tool means "not analysed", never "not present".
