# AgriPAM-AI: engine details

AI-guided, chassis-aware genome-editing design for agricultural microbial synthetic communities.

The current release contains the **AgriPAM-AI evidence engine**. It ranks host-compatible CRISPR designs from comparative defense genomics and provides an interactive genome-analysis workflow. The biocontrol claim for a named strain such as B34 remains a prospective application hypothesis until activity is measured.

## Competition package

The complete competition materials are in `deliverables/`:

- `final/AgriPAM-AI_competition_presentation.pptx` (12 slides with speaker notes)
- `final/AgriPAM-AI_competition_poster.pptx`
- `NARRATION_SCRIPT.md` (eight-minute script)
- `COMPETITION_PROJECT_REPORT.md`
- `COMPETITION_COMPLIANCE_AUDIT.md`

The step-by-step guide to how the software works is `docs/HOW_IT_WORKS.md`.

The presentation and poster distinguish computational predictions from pending physical validation.

Computational discovery and ranking of native CRISPR-Cas systems from agricultural bacterial genomes.

## SynCom chassis selector

`scripts/run_syncom_selector.py` reads the tested strain bank in `data/syncom/` and ranks which member of a compatible community is the best editing chassis (biosafety-gated, dispensability-aware). Your own bank goes in the Excel template (`data/syncom/SynCom_input_template.xlsx`, rebuilt by `scripts/make_syncom_template.py`); upload it in the app's SynCom tab to get the ranking. Rebuild the bank tables from the laboratory workbook with `scripts/build_syncom_bank.py` (needs `openpyxl`). Results are in `results/syncom_selector/` and shown in the app's SynCom tab.

## Parts and constructs

The **Parts & constructs** tab builds a checked DNA construct for the chosen edit from BioBrick parts (starter set, Excel sheet, or the iGEM Registry by `BBa_` ID), supports a reporter insertion at a neutral site or a CRISPRi guide cassette, screens it against RFC10/RFC25/Golden Gate rules, chassis restriction motifs and host range, exports FASTA and GenBank, and follows the experiment in a tracker with controls and predicted-vs-measured agreement. See `docs/HOW_IT_WORKS.md`, step 9.

## Reproducibility rules

- Keep every NCBI assembly accession attached to its files and results.
- Record every excluded genome in `metadata/exclusion_log.tsv`.
- Record software and database versions in `metadata/software_versions.tsv`.
- Use the fixed seed in `config/project.yaml` for modelling, sampling, and bootstrap analysis.
- Treat inferred PAMs and CRISPR activity as predictions until experimentally validated.
- Strain-specific PAM discovery extracts CRISPRCasTyper spacers, searches an explicitly
  supplied versioned mobile-sequence FASTA, normalizes exact protospacer matches on both
  strands, and ranks all 64 triplets with counts, background enrichment, Wilson intervals
  and multiple-testing correction. A reference checksum and every observation are exported.
- Do not change thresholds after viewing results without documenting the change and rerunning the full analysis.

## Current evidence

- A 5'-TTC Type I-C PAM is supported by 23/23 accession-coordinate observations,
  16/16 unique target contexts, and 11/11 spacer groups in *P. polymyxa*.
- TTC occurs in 35/2100 matched-background observations (BH-adjusted
  enrichment q = 6.755e-34).
- The grouped PAM-plus-mobile model achieved ROC AUC 0.9946 and average
  precision 0.5200 on a dataset with 0.99% positive prevalence.
- The PAM remains computationally inferred; functional validation is pending.

## Competition demonstration

AgriPAM-AI includes a judge-facing Streamlit application that explicitly
separates the **Genome workspace** for an uploaded strain from the fixed
**Reference discovery: comparative P. polymyxa study**. The reference
PAM model, TTC designer and evidence tabs are labelled as precomputed reference
results and are never presented as outputs from the uploaded genome.

The first page is now the **Genome workspace**. It accepts:

- an uploaded multi-contig FASTA file;
- a versioned NCBI assembly accession such as `GCF_000597985.1`; or
- a built-in demonstration sequence for offline presentations.
- raw paired-end, long-read, or hybrid FASTQ data when the corresponding
  assembler is installed.

The full workflow exports normalized FASTA, assembly QC (including N50/L50 and ambiguity burden), coding-sequence annotation,
protein FASTA, GFF, plasmid predictions, CRISPR-Cas output, DefenseFinder
status/results, restriction-modification and mobile-element screens,
agricultural-function candidates, candidate neutral genomic regions, logs,
pipeline status, an uncertainty report, effector–chassis compatibility, and a provenance manifest in one ZIP package. See
`docs/FULL_PIPELINE.md` for evidence levels and limitations.

Every inferred 5′-TTC target is also linked to its overlapping or nearest gene
in `pam_gene_context.tsv`. The report distinguishes direct CDS overlap,
promoter-proximal targets, and intergenic targets, and adds conservative growth
and agronomic-function screening flags. These flags support candidate review;
they do not by themselves prove essentiality or phenotypic effect.

For each genome, the app validates IUPAC DNA symbols, reports genome size,
contig count, GC content, ambiguous bases and a sequence checksum. It scans
both strands for unambiguous 5'-TTC Type I-C candidates, calculates GC and
homopolymer checks, counts exact within-genome copies, ranks candidates and
exports TSV and JSON results. The NCBI accession route requires an internet
connection and the bundled or system NCBI Datasets command. The FASTA route
works offline.

The lightweight environment runs every available module and marks unavailable
specialist modules explicitly. Install `environment-full.yml` to enable the
complete tool-aware workflow. No unavailable module is interpreted as a
negative biological finding.

Introduced SpCas9, Cas12a and dCas9/CRISPRi candidates receive exact whole-genome
copy counts and a bounded substitution off-target screen. For the top 2,000
unique guides per editor, every compatible site with one or two mismatches is
counted using a three-seed guarantee. Bulges, noncanonical PAMs and empirical
nuclease-specific mismatch tolerance remain outside this screen and are labelled
as unresolved. Scores are ranking heuristics, not activity probabilities.

Install the lightweight application dependencies in an isolated project-local
directory (the scientific environments are not modified):

```bash
./scripts/install_app.zsh
```

Run the application:

```bash
./scripts/run_app.zsh
```

Open the local address shown in the terminal, normally `http://localhost:8501`.

For the full environment:

```bash
./scripts/install_full_pipeline.zsh
./scripts/run_full_app.zsh --server.port 8502
```

## Deterministic checks

Run the test suite and reconstruct the competition summary:

```bash
make test
make demo-report
```

The report is written to `results/competition_demo/` as TSV and JSON. Candidate
exports from the application are also available as TSV and JSON.

The current unit-test suite covers FASTA parsing, genome summary statistics,
candidate ranking, strand handling, PAM lookup and safe NCBI accession
validation.

## Independent external validation

`data/external_validation_datasets.tsv` is the auditable dataset registry and
`docs/EXTERNAL_VALIDATION.md` defines the leakage-resistant protocol. For GSE196911,
742 predictions were generated without the outcome file and frozen before joining
the deposited OD1 results. The reproduced pooled Spearman rho is 0.4860, the median
within-gene rho is 0.5459, and the deterministic permutation control is -0.0086.
This is external experimental validation in one *E. coli* genome, not evidence of
cross-genome generalization. A second, distinct Hawkins et al. benchmark freezes
the experimental activity of 33 matched GFP guides in *E. coli* before comparing
them with held-out *B. subtilis* measurements (Spearman rho 0.4525; guide-bootstrap
95% interval 0.0876–0.7223; permutation p 0.0090). This supports partial
cross-chassis rank transfer for one reporter; it is not an AgriPAM model test or
genome-wide validation. A third benchmark fits an interpretable sequence model
to 20,195 released crisprHAL TevSpCas9 training guides and freezes predictions
for 5,049 released test guides before outcomes are joined. It reaches rho 0.6325
versus 0.3534 for a GC-only baseline, with zero exact train/test sequence overlap.
This is a predefined same-study, same-chassis guide holdout—not independent
external-study or held-out-genome validation. Reproduce it with:

```bash
python scripts/fetch_crisprhal_released_data.py
python scripts/prepare_crisprhal_released_holdout.py
```

The frozen unseen-genome challenge analyzes GCF_withheld and deliberately
abstains from a strain-specific native-system claim because the genome is outside
the *P. polymyxa* reference-model scope. GSE74926 is retained as biological context,
not misused as a guide-efficiency benchmark. B26 is
a prospective WGS/validation case and has no completed wet-lab result.

## Interpretation boundary

The target scanner identifies sequences compatible with the inferred 5'-TTC
PAM. Its model output is a prioritization score, not a measured cleavage or
editing efficiency. The bounded mismatch screen does not replace a complete,
nuclease-specific off-target analysis or experimental validation before a
candidate can be described as a functional editing tool.
