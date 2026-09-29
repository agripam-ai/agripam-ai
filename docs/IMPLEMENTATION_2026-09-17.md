# AgriPAM-AI implementation record — 2026-09-17

## Changed

- Added `agripam/validation.py`, a model-agnostic evaluator that requires deposited outcomes to be joined only after predictions are frozen.
- Added Spearman correlation, threshold-gated ROC AUC and average precision, top-k enrichment, probability-only calibration/Brier score, deterministic outcome permutation, genome-grouped split manifests, ablation comparison and checksummed result bundles.
- Added `data/external_validation_datasets.tsv` with evidence roles and eligibility decisions for GSE196911, GSE74926 and prospective B26.
- Added an **External validation** app page showing provenance, eligibility, split rules, real-table upload, metrics, controls, limitations and JSON export.
- Added `docs/EXTERNAL_VALIDATION.md` with the prespecified train/test and outcome-blinding protocol.
- Added `tests/test_validation.py` covering metric calculation, metric gating, deterministic negative controls and disjoint genome holdout.
- Updated the README to expose the validation boundary and B26 status.
- Downloaded the processed GSE196911 purine-screen deposit, verified its SHA-256 checksum and normalized 742 unique guides across nine genes without adding a prediction column.
- Added `scripts/prepare_gse196911_validation.py`; the primary endpoint is declared as negative OD1 log2 fold-change so that higher values mean stronger depletion.
- Retrieved and checksummed the article supplement. Table S13 mapped all 742 GEO records to 742 unique guide sequences and genomic contexts.
- Exported a blinded prediction-input table containing no logFC or p-value columns, plus a separate mapping manifest and checksum.
- Added the authors' Supplementary Table S14 values as a clearly labelled published reference, not an AgriPAM-AI result.
- Added within-gene Spearman reporting and its median across genes so gene-level depletion differences cannot masquerade as guide-ranking performance.
- Retrieved the authors' released 129-feature random-forest model and exact purine-screen predictor features from the versioned Zenodo archive.
- Added a narrow, documented scikit-learn 0.24-to-1.3 tree compatibility loader; it preserves all tree nodes and adds only the newer missing-value routing byte as zero.
- Added a range-based Zenodo retrieval script with pinned uncompressed sizes and SHA-256 checks for the model, feature headers and purine-screen feature table.
- Expanded the FASTA/accession workflow with taxon-agnostic assembly QC: N50, L50, contig span, ambiguity burden and explicit review flags. The report is included in every result ZIP and manifest.
- Added exact whole-assembly guide-copy counting for introduced Cas9/Cas12a/CRISPRi targets, plus GC, homopolymer and transparent priority scores. Approximate mismatch/bulge analysis remains explicitly unresolved.
- Added a downloadable module-by-module uncertainty report separating specialist calls, transparent screens, unavailable analyses and required next actions.
- Added an exact three-seed off-target search for the top 2,000 unique guides per introduced editor. It guarantees retrieval of PAM-compatible sites at one or two substitutions and explicitly excludes bulges/noncanonical PAM tolerance.
- Added off-target-adjusted guide priorities while retaining the unadjusted sequence/uniqueness score and every component needed to audit the adjustment.
- Added an effector–chassis compatibility matrix for SpCas9, Cas12a and dCas9/CRISPRi that reports targetability, exact uniqueness, bounded off-target status, defense/delivery evidence, repair evidence and uncertainty without claiming efficiency.
- Added 2,000-resample gene-cluster bootstrap intervals and a within-gene permutation null, keeping the biological gene—not each guide—as the uncertainty unit. The reproduced pooled-rho 95% interval is 0.4022–0.5522; the within-gene median interval is 0.3085–0.6155.
- Added outcome-blind, fixed-seed feature-group permutation sensitivity for sequence context, thermodynamics, position, composition and promoter status. Sequence context shows the largest reliance signal (pooled rho change -0.3732); this is explicitly not described as a retrained or causal ablation.
- Generated and checksummed 742 outcome-blind predictions, then joined the separately stored OD1 outcomes and exported a complete validation bundle.
- Reproduced the Supplementary Table S14 per-gene OD1 correlations. The primary all-nine-gene OD1 median rho is 0.5459, pooled rho is 0.4860, and the outcome-permutation control is -0.0086.
- Added continuous top-k ranking utility. The top 10% has mean depletion 1.677 log2 units above the full-set mean. ROC/PR and calibration remain gated because no binary threshold was pre-specified.
- Added leave-one-gene-out influence analysis: pooled rho remains 0.4669–0.5063 after removing any one of the nine genes.
- Added an exploratory random-effects synthesis across genes: rho 0.4842 (approximate 95% CI 0.3564–0.5942), with I-squared 75.2%. The dashboard exposes this heterogeneity and the negative *purE* result rather than hiding it.
- Added within-gene permutation significance for top-k utility. The top 1%, 5% and 10% lifts each attain the minimum one-sided empirical p of 1/2,001 at 2,000 resamples.

## Verified

- The released model generated 742/742 predictions without access to the outcome file, and the prediction/output hashes are recorded.
- The reproduced per-gene OD1 correlations agree with the authors' Supplementary Table S14 to displayed precision.
- 31 project tests passed and one optional app test was skipped in the general run; the optional Streamlit application test passed separately.

## Remains

- Define the AgriPAM-AI guide-efficiency predictor and its training sources. The existing TTC PAM model is a different prediction task and must not be presented as a CRISPRi guide-efficiency model.
- Pre-register any active-guide threshold before adding ROC/PR. Do not derive a threshold from this test set.
- Add genuinely independent genomes/experiments before claiming held-out-genome generalization. GSE196911 contains one genome.
- Implement retrained feature-group ablations on development data; input zeroing is not a substitute for retraining and is therefore not reported as ablation evidence.
- Add uncertainty intervals by genome-aware bootstrap once the real held-out dataset is available.
- Sequence B26, run genome QC and lock predictions before any physical validation. B26 remains prospective until then.

## Cross-chassis extension

- Added the authors' deposited Hawkins et al. 2020 GFP CRISPRi replicate files and variant-to-parent map with SHA-256 provenance.
- Added an outcome-separated *E. coli* to *B. subtilis* transfer pipeline for 33 fully matched guides. The *E. coli* scores are frozen before the held-out *B. subtilis* table is joined.
- Reproduced a cross-chassis Spearman rho of 0.4525, guide-bootstrap 95% interval 0.0876–0.7223 and 2,000-permutation empirical p = 0.0090.
- Added replicate-reliability, continuous top-k results, downloadable tables, dashboard provenance and explicit limitations.
- This addition closes part of the cross-organism evidence gap, but does not provide a learned-model held-out-genome validation because it uses one shared GFP reporter and an experimental transfer score.

## Predefined crisprHAL bacterial Cas9 holdout

- Added a third outcome-separated benchmark using the released crisprHAL TevSpCas9 split: 20,195 training guides and 5,049 testing guides, with zero exact 37-nt sequence overlap.
- Added an interpretable positional-sequence ridge model with regularization selected by five-fold cross-validation on the training table only. Test outcomes are written separately before prediction and joined only after the prediction file is frozen and checksummed.
- The sequence model reaches test Spearman rho 0.6325 versus 0.3534 for a GC-only baseline (delta rho 0.2791). The 2,000-resample guide-bootstrap interval is 0.6134–0.6502 and the permutation p-value is 1/2,001.
- The top 10% activity lift is 2.5479 units with one-sided permutation p = 1/2,001. The dashboard labels this as a same-study, same-chassis guide-level holdout—not an independent external study or held-out-genome validation.
- Added an exhaustive sequence-similarity leakage audit against all 20,195 training contexts. Only 18 test guides have at least 18/20 nearest-guide identity; excluding them leaves 5,031 guides with rho 0.6339, showing that the result is not driven by the small near-identical subset.
