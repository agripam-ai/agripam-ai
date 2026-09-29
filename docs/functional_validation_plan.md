# Native-strain functional validation plan

## Claim to test

AgriPAM-AI predicts that a native Type I-C system in *P. polymyxa*
recognizes a 5'-TTC PAM. Computational evidence is not treated as functional
validation.

## Decision gates

1. **Genome identity and quality** — confirm the isolate's taxonomic placement,
   assembly quality, and contamination status.
2. **System presence** — confirm that the assembly contains an intact Type I-C
   Cas operon and an associated CRISPR array. Confirm key loci and array sequence
   by PCR/Sanger sequencing where feasible.
3. **Candidate eligibility** — select a native spacer with a unique target and
   design otherwise identical TTC and non-TTC target constructs. Complete an
   off-target review before testing.
4. **Functional comparison** — compare the TTC target, PAM-mutant target, and
   non-targeting control using a predefined interference or reporter readout.
5. **Interpretation** — claim PAM-dependent activity only if the TTC condition is
   reproducibly different from both controls across independent biological
   replicates. Report effect size, uncertainty, exclusions, and all raw data.

## Minimum controls

- Matched protospacer with the predicted 5'-TTC PAM.
- Same protospacer with a prespecified non-TTC PAM substitution.
- Non-targeting sequence control.
- Transformation/assay viability control.
- No-template and positive controls for diagnostic PCR.

## Predefined outputs

- Sequencing and assembly QC report.
- Type I-C locus map and PCR/Sanger confirmation table.
- Candidate and off-target report exported by AgriPAM-AI.
- Raw replicate-level experimental measurements.
- Effect sizes with confidence intervals and an explicitly defined exclusion log.
- Final status: supported, not supported, or inconclusive.

## Competition fallback

If functional testing cannot be completed before submission, report PCR/Sanger
confirmation as physical verification of system presence and label PAM function
as pending. Never represent locus confirmation as cleavage or editing validation.
