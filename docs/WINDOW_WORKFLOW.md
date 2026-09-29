# The software window by window

This guide walks through every window (tab) in the order they appear, with what to enter, what to expect and how to read it. Labels are quoted from the interface.

## What the software can do

1. **Analyse a microbial bank.** Score tested strains (bacteria and fungi) for growth-promotion traits, enzymes and biocontrol, and record how they behave together.
2. **Design a community of two or more components** from those traits and from pairwise compatibility, with a biosafety gate, and find the functions the community lacks.
3. **Choose what to edit.** Rank the members as editing chassis: prefer the member the others can spare.
4. **Find editing targets** beyond natural systems: infer a native PAM from spacer evidence, and scan introduced editors on the strain's genome, scored against genes, defenses and off-targets.
5. **Integrate with synthetic biology.** Assemble BioBrick parts into a checked reporter or CRISPRi construct.
6. **Follow the result.** Track each construct and its controls, and compare prediction with measurement.

## The order of work

**Bank → community → chassis → genome targets → construct → tracker → evidence.**

The tabs are: Start here · Genome evaluation · SynCom candidate bank · Agricultural editing knowledgebase · Parts & constructs · External validation · Software & reproducibility. The sidebar holds the **Text size** selector, **How this software works**, and the **Analysis contract**.

---

## 1. Start here

**Purpose:** orient a new user and show why a PAM alone is not enough.

- **Interactive CRISPR learning model.** Choose the *Editing system*, *Educational objective*, *Target-site evidence* and *Guide specificity*, and switch *Compatible delivery demonstrated* and *Compatible repair strategy*. The panel shows what is and is not established. Two closed panels explain the simulation and the IUPAC notation. It is educational: it designs no guide and predicts no efficiency.
- **Plan your target-organism project.** Ten numbered steps (define the objective, verify the strain and genome, choose the editing system, identify a target, define repair, establish delivery, evaluate defense barriers, define controls, confirm the genotype, validate function and safety). Each carries a status telling you whether the software addresses it (for example "Run Genome evaluation") or whether it needs laboratory work. Set *Intended modification* and *Primary measurable readout* to tailor the plan.

**Output:** a plan on screen; nothing is saved.

## 2. Genome evaluation

**Purpose:** analyse one genome from assembly to ranked edit sites.

1. Fill *Target organism name* and *Strain or isolate name*, and pick the *PAM-inference reference collection*.
2. Supply the genome in one of three ways: *Genome assembly (FASTA)*, *Raw reads (FASTQ)* with *Sequencing technology*, or a *Versioned NCBI assembly accession*. Optionally add a *Separate mobile-sequence reference* for strain-specific PAM discovery, and tick *Run installed specialist analyses*.
3. Press **Run genome workflow**.

**What appears:** automatic mobile-reference evidence; genome quality control; module status; an editing-readiness assessment; introduced-effector compatibility with the chassis; an explorer of introduced-editor targets (SpCas9 NGG, Cas12a TTTV, dCas9/CRISPRi NGG) with a bounded off-target screen; and the **Multi-objective design portfolio** (editability, deliverability, agronomic value and preservation safety, with the Pareto set and a blank results ledger). One ZIP holds every table, log and version record.

**Reading it:** a module marked unavailable or limited means "not analysed". Candidate neutral regions are candidates, not safe harbors.

## 3. SynCom candidate bank

**Purpose:** design the community and choose the member to edit.

- **Analyze your own bank** (first section). Download the Excel template (or the filled synthetic example), enter one row per strain with *kind* (bacterium or fungus), biosafety fields, editing precedent and function scores 0 to 5 (blank = not tested), plus fungus and bacterium compatibility. Upload it, choose the *Fungal anchor(s) of the community*, *Assemble one for me* or *I choose the members*, and the delivery threshold, then press **Run SynCom analysis**. You get the community, its function coverage, the gaps and possible donors, the chassis ranking (fit, dispensability, editing precedent; biosafety-flagged strains are not ranked), a robustness check, and a results workbook.
- **Community-aware chassis selection.** A precomputed worked example, selectable by scenario.
- **Native-bank evidence and candidate selection.** Phenotypes and proposed roles, reference discovery across the candidate collection, the public-reference screen, the fungal editing-readiness reference, live all-accession CRISPR-Cas results, a completed pilot, and a table of which isolates to sequence first (downloadable).
- **Reference panels (closed by default).** The 64-triplet PAM model; the 5′-TTC candidate designer (type a candidate 3-nt PAM and say whether the target lies near a mobile element); the auditable model and control evidence; the native-bank defense and editing-system atlas; and a quick reference-TTC scan on an uploaded genome, an accession, or a built-in demonstration.

**Reading it:** reference panels describe the fixed *P. polymyxa* reference study; they are not results for an uploaded strain.

## 4. Agricultural editing knowledgebase

**Purpose:** what is already known about editing agricultural microbes.

- **Thirty-isolate native-system discovery cohort:** choose priorities and download the table.
- **Search and filter the evidence:** search by organism, strain or role; filter by evidence level, source type, agricultural status, editing route and PAM evidence; download the filtered records. You can add reviewed records by upload.
- **Evidence interpretation and update policy:** what "experimentally demonstrated" and "computationally inferred" mean here.
- **Target-strain editing protocol planner:** pick a target organism and editing route to get a step-by-step protocol (downloadable).

## 5. Parts & constructs

**Purpose:** turn the chosen edit into DNA and follow it.

1. **Parts library:** load the starter library, upload a parts workbook, or fetch parts by BBa_ ID from the iGEM Registry; edit each part's role and declared host range.
2. **Choose the edit and build the construct:** set the *Chassis (species)*, the *Assembly standards to screen*, optional *Chassis restriction motifs*, the *Junction between parts* and the *Edit type*. For a reporter insertion, pick a neutral region (from the genome workflow or typed), the arm length, the orientation and the parts. For a CRISPRi cassette, pick a guide (from the workflow or a pasted 20-nt spacer), a promoter, an sgRNA scaffold and an optional terminator.
3. **Constructs, checks and downloads:** length, GC, errors and warnings, an annotated map, FASTA and GenBank files, and expected diagnostic-PCR sizes.
4. **Follow the output:** create tracker rows with the controls (parental, negative, low-ranked guide), fill stages and measurements, download, upload again to see progress, warnings and the predicted-vs-measured correlation.

**Reading it:** a construct that passes the screens is compatible with the standard, not proven to work in your chassis.

## 6. External validation

**Purpose:** show how far the predictions are supported by other laboratories' experiments.

Sections: dataset provenance and eligibility; the deposited held-out outcome package (742 guides) with gene-cluster uncertainty, robustness and ranking utility; the published reference results; the predefined bacterial Cas9 released holdout; the cross-chassis transfer benchmark; the frozen unseen-genome challenge; the leakage-control contract; a tool to **evaluate a frozen held-out table** you upload (predictions plus subsequently joined outcomes, with an activity threshold if your outcome is binary); and the B26 prospective case. Every benchmark states what it does not show.

## 7. Software & reproducibility

**Purpose:** what modules exist and how to reproduce results.

- **Editing-system modules and capabilities:** an atlas of the editing systems the software screens and what evidence each needs.
- **Reproducibility contract:** versions, seeds (42), checksums, and the commands that rebuild the results.

---

## Four typical journeys

| I have | Go to | Then |
|---|---|---|
| A tested bank | SynCom candidate bank → Analyze your own bank | choose the chassis → Genome evaluation for its genome |
| A genome to edit | Genome evaluation | Parts & constructs for the construct |
| A construct to test | Parts & constructs (sections 3 and 4) | fill the tracker as the experiment proceeds |
| A claim to check | External validation | Software & reproducibility for the commands |
