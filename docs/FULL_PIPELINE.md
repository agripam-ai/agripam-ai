# RhizoForge-Select full genome workflow

## Inputs

- assembled bacterial genome in FASTA format;
- versioned NCBI assembly accession (`GCF_...` or `GCA_...`); or
- raw FASTQ reads: paired Illumina, long reads, or hybrid data.

Raw reads are assembled with SPAdes, Flye, or Unicycler. An assembled FASTA is
the fastest and most reproducible competition demonstration input.

## Outputs and evidence levels

| Output | Preferred engine | Fallback |
|---|---|---|
| CDS, proteins and GFF | Bakta, then Prokka | Prodigal CDS/proteins without functional products |
| Plasmids | MOB-suite | labelled contig-level heuristic |
| CRISPR-Cas | CRISPRCasTyper | unavailable, never reported as absence |
| Defense systems | DefenseFinder | unavailable, never reported as absence |
| Restriction-modification | annotated-product screen | limited if product annotation is absent |
| Mobile elements | geNomad plus annotated products | annotated-product screen |
| Agricultural functions | annotated-product screening | limited if product annotation is absent |
| Candidate neutral regions | exclusion-window screen | unavailable without coordinates |
| PAM-to-gene consequences | target/CDS coordinate overlap and nearest-gene mapping | limited without CDS annotation |

The agricultural-function table searches annotations for antifungal compounds,
lipopeptide synthesis, secretion, hydrolytic enzymes, phosphate solubilization,
nitrogen metabolism, siderophore production, and root-colonization traits. These
are candidate associations. Publication-grade claims require domain/cluster
confirmation (for example antiSMASH for biosynthetic clusters), synteny review,
and experimental validation.

Candidate neutral regions are regions separated from annotated genes and all
identified functional, defense, and mobile features. They are **not proven safe
harbors**. Essentiality, conservation/synteny, transcription, growth effects,
and experimental stability must be evaluated before integration.

`pam_gene_context.tsv` reports every 5′-TTC-compatible target, its coordinates,
overlapping or nearest CDS, gene/product annotation, coding or regulatory
context, and screening flags for growth-related and agronomic functions.
Keyword-based growth flags are not essential-gene calls. A selected edit must
still be checked against dedicated essentiality evidence, conservation,
transcript context, and experimental phenotype.

## Installation

```bash
chmod +x scripts/install_full_pipeline.zsh scripts/run_full_app.zsh
./scripts/install_full_pipeline.zsh
```

Initialize the databases required by the versions of Bakta, DefenseFinder,
CRISPRCasTyper, MOB-suite, and geNomad installed in the environment. Database
downloads are intentionally not hidden inside the application: their versions
must be recorded for reproducibility and their licenses/redistribution terms
must be respected.

Run:

```bash
./scripts/run_full_app.zsh --server.port 8502
```

By default the environment is installed at `.envs/rhizoforge-select` inside the project. Set `RHIZOFORGE_ENV_PREFIX` to use another location.
This no-space prefix is required because several Bioconda command launchers do
not support environment paths containing spaces. The project itself remains in
its original directory.

The downloadable ZIP contains the normalized genome, GFF, protein FASTA,
module-specific folders, logs, tabular screens, pipeline status, and a JSON run
manifest. Missing modules remain visible in `pipeline_status.tsv`.

## Strain-specific PAM discovery

Upload a versioned phage, plasmid or mobile-element FASTA collection together
with the host genome and enable specialist analyses. RhizoForge-Select extracts
spacers from CRISPRCasTyper, searches both reference orientations for exact
full-length protospacer matches, retrieves the normalized 5′ adjacent triplet,
and exports:

- `crispr_arrays.tsv` and `crispr_spacers.tsv`;
- `pam_protospacer_observations.tsv`;
- `pam_64_triplet_ranking.tsv`, with counts, background frequency, enrichment,
  Wilson intervals, nominal P values and BH-adjusted Q values;
- `pam_logo_frequencies.tsv` and `pam_discovery_summary.json`.

Only motifs supported by at least three observations and BH Q < 0.10 are passed
to the gene-context mapper. The program never silently substitutes TTC when
strain-specific evidence is absent. The current exact-match policy is recorded
in the manifest. PAM calls remain computational hypotheses until experimental
PAM-library validation.
