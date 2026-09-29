"""NCBI Datasets adapter for the interactive application."""
import json
import re
import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path

ACCESSION_PATTERN = re.compile(r"^(GC[AF]_\d{9}\.\d+)$", re.I)

def validate_assembly_accession(accession: str) -> str:
    accession = accession.strip().upper()
    if not ACCESSION_PATTERN.fullmatch(accession):
        raise ValueError("Enter a versioned assembly accession such as GCF_000597985.1")
    return accession

def find_datasets_executable(project_root: Path) -> str:
    bundled = project_root / ".tools" / "bin" / "datasets"
    if bundled.is_file():
        return str(bundled)
    executable = shutil.which("datasets")
    if executable:
        return executable
    raise RuntimeError("NCBI Datasets is unavailable. Upload a FASTA file instead.")

def download_genome_fasta(accession: str, project_root: Path, timeout: int = 240) -> str:
    accession = validate_assembly_accession(accession)
    with tempfile.TemporaryDirectory(prefix="agripam_ncbi_") as tmp:
        archive = Path(tmp) / "genome.zip"
        command = [find_datasets_executable(project_root), "download", "genome", "accession",
                   accession, "--include", "genome", "--filename", str(archive), "--no-progressbar"]
        result = subprocess.run(command, capture_output=True, text=True, timeout=timeout, check=False)
        if result.returncode != 0:
            detail = result.stderr.strip().splitlines()[-1] if result.stderr.strip() else "download failed"
            raise RuntimeError(f"NCBI download failed: {detail}")
        if not zipfile.is_zipfile(archive):
            raise RuntimeError("NCBI returned an incomplete archive; retry or upload FASTA")
        with zipfile.ZipFile(archive) as package:
            names = sorted(n for n in package.namelist() if n.endswith((".fna", ".fa", ".fasta")))
            if not names:
                raise RuntimeError("The NCBI package contained no genome FASTA")
            return package.read(names[0]).decode()


def get_assembly_identity(accession: str, project_root: Path, timeout: int = 120) -> dict[str, str]:
    """Return organism and strain labels from an NCBI assembly summary."""
    accession = validate_assembly_accession(accession)
    command = [find_datasets_executable(project_root), "summary", "genome", "accession",
               accession, "--as-json-lines"]
    result = subprocess.run(command, capture_output=True, text=True, timeout=timeout, check=False)
    if result.returncode != 0 or not result.stdout.strip():
        return {"accession": accession, "organism": "NCBI organism (metadata unavailable)", "strain": accession}
    try:
        record = json.loads(result.stdout.splitlines()[0])
    except (json.JSONDecodeError, IndexError):
        return {"accession": accession, "organism": "NCBI organism (metadata unavailable)", "strain": accession}
    organism = (record.get("organism", {}).get("organism_name")
                or record.get("organism", {}).get("organismName")
                or "NCBI organism")
    infraspecific = record.get("organism", {}).get("infraspecific_names", {}) or record.get("organism", {}).get("infraspecificNames", {})
    strain = infraspecific.get("strain") or record.get("assembly_info", {}).get("assembly_name") or accession
    return {"accession": accession, "organism": str(organism), "strain": str(strain)}
