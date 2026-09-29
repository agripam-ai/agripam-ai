# External validation protocol

AgriPAM-AI separates model development, retrospective external validation and prospective validation.

1. Register the deposited dataset, publication, version and checksum.
2. Normalize guide sequence, organism, genome/strain, experimental condition and measured outcome without generating predictions.
3. Partition by complete genome/study, never random guide alone. The held-out genome identifiers are written to a manifest.
4. Train or select the predictor using the training partition only. Feature transforms and thresholds are fixed there.
5. Generate and freeze held-out predictions without the outcome column present.
6. Join outcomes by `guide_id` only after freezing predictions, then calculate Spearman correlation. Calculate ROC AUC, average precision and top-k enrichment only for a pre-specified, biologically justified activity threshold. Calculate calibration only for probability-valued predictions.
7. Run an outcome-permutation negative control and compare feature-ablation models using the identical split.
8. Export the normalized table, split manifest, configuration, checksums, software version and metrics.

For GSE196911, Table S13 of the peer-reviewed supplement is joined to GEO by gene and all three deposited time-point measurements. The resulting blinded prediction input contains sequence/context fields but no logFC or p-value fields. The released Yu et al. forest is then run on those inputs and its predictions are written and checksummed before the outcome file is joined. Primary reporting includes within-gene Spearman correlations and their median because a pooled correlation can be inflated by between-gene depletion differences.

The reproduced primary OD1 result is 742 guides across nine genes: pooled Spearman rho = 0.4860 and median within-gene rho = 0.5459. A deterministic outcome-permutation control gives rho = -0.0086. A 2,000-resample gene-cluster bootstrap gives a pooled-rho 95% interval of 0.4022–0.5522 and a median-within-gene interval of 0.3085–0.6155. A 2,000-resample within-gene permutation test gives empirical two-sided p = 1/2,001 for both statistics; this is the minimum attainable value at that resampling depth. The top 10% by frozen prediction has mean depletion 1.677 log2 units above the all-guide mean. ROC/PR and calibration are not reported because no binary threshold was pre-registered and the model output is not a calibrated probability. This is a held-out experiment in one *E. coli* genome, not a held-out-genome validation.

Robustness analyses make the result more informative without changing the frozen predictions. Leave-one-gene-out pooled rho ranges from 0.4669 to 0.5063, so no single gene produces the overall association. An exploratory DerSimonian-Laird synthesis of the nine within-gene Spearman correlations gives rho = 0.4842 (approximate 95% CI 0.3564–0.5942) and I-squared = 75.2%, demonstrating substantial biological heterogeneity rather than uniform performance. One gene, *purE*, has a negative within-gene correlation; it is retained and displayed as a failure case. Within-gene outcome permutations give one-sided empirical p = 1/2,001 for the observed top 1%, 5% and 10% mean lifts. The random-effects calculation uses approximate Fisher-z variances for Spearman correlations and is therefore labelled exploratory.

The forest was released under scikit-learn 0.24.2. Reproduction under 1.3 uses a narrow compatibility loader that adds the later `missing_go_to_left` node byte as zero while preserving all released nodes, thresholds and leaf values. No retraining is performed. The per-gene OD1 correlations reproduce Supplementary Table S14 values to displayed precision. The paper's headline 0.588838 is a different aggregate: the median excluding purE/purK across all three time points.

Feature-group reliance is assessed without outcomes by jointly permuting each feature group across the 742 blinded inputs with seed 42, freezing those sensitivity predictions, and evaluating them only after the outcome join. Sequence-context permutation reduces pooled rho from 0.4860 to 0.1128 and median within-gene rho from 0.5459 to 0.1981. Position permutation reduces them to 0.3914 and 0.4124. Thermodynamic, composition and promoter perturbations cause smaller changes in this dataset. These are model-reliance sensitivity results under distributional perturbation—not retrained ablations, causal effects or proof that a feature group is biologically unimportant.

Reproduction starts with `scripts/fetch_yu2024_released_model.py`, which retrieves only the three pinned ZIP members required from Zenodo and rejects any size or SHA-256 mismatch. `scripts/generate_yu2024_frozen_predictions.py` then accepts the blinded table, released feature table, model and headers and emits a checksummed prediction file without reading the held-out outcome file.

GSE196911 / Zenodo 10.5281/zenodo.10262866 is registered as the primary public retrospective benchmark. It is not displayed as completed until the deposited guide-level table is imported, normalized and checksummed. GSE74926 supports the feasibility of *B. subtilis* CRISPRi but is excluded from guide-efficiency metrics because its GEO deposit is not a clean guide-level outcome benchmark.

## Cross-chassis transfer benchmark

Hawkins et al. (2020; DOI 10.1016/j.cels.2020.09.009; SRA PRJNA574461) measured the same GFP-targeting CRISPRi library in *E. coli* and *B. subtilis*. The authors' public `mismatch_crispri/gfpdata` release contains 33 fully matched parent guides and their single- and double-mismatch derivatives, with two biological replicate files per chassis. AgriPAM-AI reconstructs each fully matched parent measurement from the released variant-to-parent map, averages the two *E. coli* replicates, writes and checksums those values as frozen transfer scores, and only then joins the separately written *B. subtilis* replicate means by guide ID.

Across the 33 sequence-matched guides, the *E. coli* transfer score has Spearman rho = 0.4525 against the held-out *B. subtilis* activity. A 2,000-resample guide bootstrap gives a 95% interval of 0.0876–0.7223, and a 2,000-permutation two-sided test gives empirical p = 0.0090. Replicate agreement is rho = 0.8222 in *E. coli* and rho = 0.9873 in *B. subtilis*. The top 20% ranked in *E. coli* has mean *B. subtilis* activity 0.0394 units above the all-guide mean.

This benchmark is deliberately labelled a sequence-matched experimental transfer baseline. The source value is an *E. coli* experiment, not an AgriPAM-trained prediction. All guides target one shared GFP reporter, overlapping targets are not independent, and the result does not constitute genome-wide, native-locus or B26 validation. ROC/PR and calibration are omitted because no binary activity threshold or probability output was pre-specified.

Reproduction uses `scripts/prepare_hawkins2020_cross_species_validation.py`. Raw releases, the original-guide map, frozen scores, held-out outcomes, joined table, metrics and SHA-256 values are retained under `data/external_validation/`.

## crisprHAL TevSpCas9 predefined released holdout

The crisprHAL repository provides separate TevSpCas9 training and testing tables for a bacterial activity model. AgriPAM-AI uses this split as a third benchmark with deliberately narrower language: it is a predefined guide-level holdout from one study and one *C. rodentium* chassis, not an external-study or held-out-genome validation.

`scripts/prepare_crisprhal_released_holdout.py` validates each 37-nt context, writes the 5,049 test sequences and outcomes to separate files, fits an interpretable ridge model using only the 20,195 released training guides, selects regularization by five-fold training-only cross-validation, freezes and checksums test predictions, and only then joins the held-out activities. There are zero exact 37-nt sequence overlaps between the released training and test tables.

The frozen sequence model gives Spearman rho = 0.6325 on the released test set, compared with rho = 0.3534 for a GC-only baseline (delta rho = 0.2791). A 2,000-resample guide bootstrap gives a 95% interval of 0.6134–0.6502, and the outcome-permutation test gives empirical two-sided p = 1/2,001. The top 10% has a mean activity lift of 2.5479 relative to the entire test set, with permutation p = 1/2,001. ROC, precision-recall and calibration are omitted because no binary activity threshold or calibrated probability was pre-specified. Guide-level resampling cannot account for unreported gene or target-cluster dependence.

A complete nearest-sequence audit compares each test context with all 20,195 training contexts. Eighteen of 5,049 test guides have a nearest training guide with at least 18 of 20 identical positions; seven contexts have at least 34 of 37 identical positions. After excluding all test guides above 17 of 20 nearest-guide identity, 5,031 test guides remain and Spearman rho is 0.6339. Thus the headline association is not driven by the small near-identical subset, although the released split remains guide-level rather than genome-level.

B26 remains a prospective case. Its WGS, native-system confirmation, PAM inference, guide predictions and physical outcomes are pending.
