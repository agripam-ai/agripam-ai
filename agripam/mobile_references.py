"""Build taxonomically tiered mobile-element references from local collections."""
from __future__ import annotations

import csv
import re
from pathlib import Path


SPECIES_COMPLEX = {
    "bacillus velezensis": ("bacillus velezensis", "bacillus amyloliquefaciens", "bacillus siamensis", "bacillus subtilis"),
    "bacillus subtilis": ("bacillus subtilis", "bacillus velezensis", "bacillus amyloliquefaciens", "bacillus siamensis"),
}

# Genome databases may retain former genus names after taxonomic revisions.
# These aliases prevent a nomenclatural change from hiding relevant references.
TAXONOMIC_ALIASES = {}
GENUS_ALIASES = {}


def _tier(label: str, organism: str) -> tuple[int, str] | None:
    label, organism = label.lower(), organism.lower().strip()
    genus = organism.split()[0] if organism else ""
    species_labels = TAXONOMIC_ALIASES.get(organism, (organism,))
    if organism and any(species and species in label for species in species_labels):
        return 2, "same species"
    relatives = SPECIES_COMPLEX.get(organism, ())
    if relatives and any(relative in label for relative in relatives):
        return 3, "related species complex"
    genus_labels = GENUS_ALIASES.get(genus, (genus,))
    if any(re.search(rf"\b{re.escape(candidate)}\b", label) for candidate in genus_labels):
        return 4, "same genus"
    return None


def _stream_fasta(path: Path):
    header = None
    chunks: list[str] = []
    with path.open(errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:].strip(), []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield header, "".join(chunks)


def build_tiered_mobile_reference(root: Path, organism: str) -> tuple[str | None, dict]:
    """Extract species-complex/genus phages and plasmids without loading broad FASTAs at once."""
    root = Path(root)
    manifest = root / "metadata" / "caudoviricetes_refseq_sequence_manifest.tsv"
    phage_fasta = root / "data" / "phages" / "processed" / "caudoviricetes_refseq_exact_deduplicated.fna"
    accession_tiers: dict[str, dict] = {}
    if manifest.exists():
        with manifest.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                classification = _tier(f"{row.get('host_name', '')} {row.get('virus_name', '')}", organism)
                if classification:
                    tier, scope = classification
                    accession_tiers[str(row.get("accession", ""))] = {
                        "tier": tier, "scope": scope, "source": "phage",
                        "record_name": row.get("virus_name", ""),
                        "reported_host": row.get("host_name", ""),
                    }

    records: list[str] = []
    inventory: list[dict] = []
    counts = {2: 0, 3: 0, 4: 0}
    if phage_fasta.exists() and accession_tiers:
        for header, sequence in _stream_fasta(phage_fasta):
            accession = header.split()[0]
            if accession in accession_tiers:
                details = accession_tiers[accession]
                tier, scope = details["tier"], details["scope"]
                records.append(f">{accession} source=phage evidence_tier={tier} scope={scope.replace(' ', '_')}\n{sequence}")
                inventory.append({
                    "accession": accession, "source_type": "phage", "evidence_tier": tier,
                    "taxonomic_scope": scope, "record_name": details["record_name"],
                    "reported_host": details["reported_host"], "sequence_length_bp": len(sequence),
                    "ncbi_record": f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}",
                })
                counts[tier] += 1

    for path in sorted((root / "data" / "plasmids" / "ncbi_raw_chunks").glob("*.fasta")):
        for header, sequence in _stream_fasta(path):
            classification = _tier(header, organism)
            if classification:
                tier, scope = classification
                accession = header.split()[0]
                records.append(f">{accession} source=plasmid evidence_tier={tier} scope={scope.replace(' ', '_')}\n{sequence}")
                inventory.append({
                    "accession": accession, "source_type": "plasmid", "evidence_tier": tier,
                    "taxonomic_scope": scope, "record_name": header,
                    "reported_host": "parsed from NCBI FASTA definition", "sequence_length_bp": len(sequence),
                    "ncbi_record": f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}",
                })
                counts[tier] += 1

    metadata = {
        "organism": organism, "records": len(records),
        "tier_2_same_species": counts[2], "tier_3_species_complex": counts[3],
        "tier_4_same_genus": counts[4],
        "reference_inventory": inventory,
        "confidence_policy": "Tier 2 species > Tier 3 species complex > Tier 4 genus; all require experimental validation",
    }
    return ("\n".join(records) if records else None), metadata
