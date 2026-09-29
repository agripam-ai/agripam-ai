**Aim.** From a tested microbial bank to an edited, tracked synthetic community (SynCom): which strains to combine, which member to edit, where and how to edit it, and how to follow the result.

**What it does**
1. **Designs a community of two or more components** (bacteria and fungi) from agronomic traits (growth promotion, biocontrol, enzymes) and pairwise compatibility, and lists the functions it still lacks.
2. **Chooses the member to edit**: the one the rest of the community can spare, behind a biosafety gate.
3. **Finds editing targets**: inferred native systems (for example a Type I-C PAM) and introduced editors (SpCas9, Cas12a, dCas9), scored against genes, defenses and off-target risk.
4. **Connects to synthetic biology**: BioBrick parts (starter set, Excel sheet or iGEM Registry) are assembled into a checked construct, either a reporter insertion or a CRISPRi guide cassette.
5. **Follows the output**: a tracker with controls compares predicted and measured results.

**How it works, in order**
Bank → community → chassis → genome targets → construct → tracker. The first three steps need only bank tables; the genome and construct steps add sequence.

**Evidence status**
- The inferred 5′-TTC PAM in the *P. polymyxa* reference panel is supported by 23/23 observations and recovered in every leave-one-out fold; it is not yet measured.
- Guide ranking transfers on three deposited benchmarks (Spearman ρ 0.49, 0.45, 0.63).
- No project wet-lab result exists yet; a verification plan and tracker are provided.

**Reading the results.** Scores are decision support, not editing efficiency or safety clearance. A missing analysis means "not analysed", never "not present". Blank laboratory cells mean "not tested".
