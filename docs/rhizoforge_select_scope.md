# RhizoForge-Select scope

## Product architecture

**RhizoForge-Select** is the competition-facing product. **AgriPAM-AI** is its current evidence engine. The engine reads a FASTA or NCBI accession, summarizes the chassis, identifies candidate Type I-C 5′-TTC targets, scores sequence quality and exports auditable designs.

The next application layer can use the same ranking engine to design a regulatory edit in a validated beneficial strain such as B34, with a defined *[taxon withheld]* biocontrol assay as the endpoint.

## Core hypothesis

An AI model that accounts for chassis defense systems, PAM compatibility and sequence context can prioritize regulatory edits that improve a beneficial strain's activity against *[taxon withheld]* while preserving nutrient-transformation functions and avoiding unnecessary effects on beneficial fungi.

This is a testable hypothesis, not a result established by the current genome survey.

## What can be claimed now

- The three representative bacterial panels have distinct defense landscapes.
- *P. polymyxa* contains validated Type I-C and Type III-B systems in the analyzed panel.
- A 5′-TTC motif is a high-priority computational candidate for *P. polymyxa* Type I-C targeting.
- The software converts those observations into ranked, reviewable candidate sequences.

## What requires new evidence

- B34 genome sequence and annotation.
- Identification of a regulatory target linked to antifungal activity.
- Demonstration that the edit changes antifungal activity against *[taxon withheld]*.
- Measurements showing nutrient-transformation functions remain intact.
- Non-target fungal compatibility testing.

## Recommended competition framing

Present RhizoForge-Select as a **closed-loop design and validation system**:

1. Learn chassis constraints from genomes.
2. Rank an editor and guide for the selected chassis.
3. Predict the functional trade-off before editing.
4. Test the highest-ranked design with a measurable reporter or biocontrol assay.
5. Feed the measured result back into the model.

Do not present the current *P. polymyxa* PAM discovery as proof that an isolate has the same PAM or that a proposed edit improves biocontrol.
