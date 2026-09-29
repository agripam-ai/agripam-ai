# Interim discovery screen: *Paenibacillus polymyxa* CRISPR–mobilome interactions

## Scope

This report ranks hypotheses for follow-up. It does **not** claim novelty or
function. Results are derived from the curated 24-genome collection, confident
CRISPRCasTyper operons, exact spacer matches, NCBI replicon metadata, and the
existing PAM-orientation workflow.

## Audited findings

1. The collection contains 12 confident Cas operons: seven Type I-C and five
   Type III-B.
2. All seven Type I-C operons occur on chromosomal replicons. Two of the five
   Type III-B operons occur on plasmids: pT1-8_1 and pl29.
3. Type I-C gene content is highly conserved. Six operons share the same
   forward gene profile and the seventh has the corresponding reverse
   orientation.
4. Repeat-only subtype predictions disagree with the adjacent Cas-operon
   subtype for 26 of 28 array–operon links. In genomes with only Type I-C, the
   immediately adjacent arrays are commonly labelled III-B by the repeat
   classifier. In the genome carrying both systems, proximity links the
   apparently reciprocal labels to the opposite operons. Therefore, array
   ownership must be assigned from genomic context and distance, not from the
   repeat classifier alone.
5. Two distinct exact-match spacers, ASP000006 and ASP000086, target plasmid
   pT1-8_1 with an orientation-normalized TTC flank. Each spacer occurs in the
   Type I-C-linked arrays of the same six other genomes.
6. pT1-8_1 encodes a confident Type III-B operon. The two protospacers do not
   overlap its Cas genes; they occur in a distant region containing a
   discoidin-domain protein, a helix–turn–helix protein, and nearby phage-holin
   annotations.
7. Strain T1-8 does not contain either of the two matching spacer sequences in
   the audited membership table. This is inter-strain targeting of a resident
   plasmid, not self-targeting within T1-8.

## Leading hypothesis

The strongest current hypothesis is that a Type I-C-bearing *P. polymyxa*
lineage acquired multiple spacers against a circulating mobile element related
to pT1-8_1, a plasmid that itself carries Type III-B defence machinery. This
could represent a lineage-specific barrier to the horizontal spread of a
defence-bearing plasmid.

The general phenomenon of plasmid-encoded CRISPR systems and inter-mobile-
element competition is already known. What remains unresolved is whether the
specific Type I-C/pT1-8_1 interaction, its repeated spacer acquisition, and its
effect on defence-system circulation in *P. polymyxa* are new and functional.

## Alternative explanations that must be excluded

- The plasmid may contain conserved sequence shared with unrelated mobile
  elements, so the historical invader may not have been pT1-8_1 itself.
- The exact matches may reflect vertical inheritance of spacers among closely
  related genomes rather than six independent acquisition events.
- TTC may be an acquisition motif but not a functional interference PAM.
- Assembly or annotation errors may affect plasmid assignment.
- The Type I-C systems may be transcriptionally silent or otherwise inactive.
- Repeat-subtype misclassification could have propagated into earlier labels;
  all PAM evidence must be reissued with context-based array ownership.

## Decisive analyses

1. Search a much broader *Paenibacillus* genome and plasmid collection for the
   two spacers, their protospacers, pT1-8_1 homologues, and Type III-B plasmids.
2. Build host and plasmid phylogenies to distinguish vertical inheritance from
   independent acquisition and horizontal transfer.
3. Compare the pT1-8_1 target region across related plasmids and test whether
   the two protospacers are conserved backbone sites or phage-derived cargo.
4. Reassign every array to a Cas operon using distance, orientation, subtype,
   and single-system genomes; repeat the PAM analysis after that curation.
5. Confirm the target flanks directly from the deposited sequence and test
   alternative PAM orientation conventions.
6. Experimentally compare transformation or maintenance of target plasmids
   carrying TTC versus PAM-mutant controls in an active Type I-C host.

## Current claim boundary

Permitted: “Comparative analysis identified two recurrent Type I-C-linked
spacers that match a phage-associated region of a Type III-B-bearing plasmid,
suggesting a testable barrier to defence-plasmid circulation.”

Not yet permitted: “Type I-C prevents transfer of Type III-B plasmids,” “TTC is
the validated interference PAM,” or “a new CRISPR mechanism was discovered.”

