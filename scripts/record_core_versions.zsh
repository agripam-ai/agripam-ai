#!/bin/zsh
set -euo pipefail

PROJECT_DIR="${0:A:h:h}"
MAMBA_BIN="$PROJECT_DIR/.tools/bin/micromamba"
MAMBA_ROOT="$PROJECT_DIR/.tools/micromamba-root"
ENV_PREFIX="$PROJECT_DIR/.envs/agripam-core"
VERSION_FILE="$PROJECT_DIR/metadata/software_versions.tsv"
ACCESS_DATE="$(date +%F)"

export MAMBA_ROOT_PREFIX="$MAMBA_ROOT"

if [[ ! -x "$MAMBA_BIN" || ! -d "$ENV_PREFIX" ]]; then
  print -u2 "Core environment is missing. Run scripts/bootstrap_tools.zsh first."
  exit 1
fi

print "software\tversion\tdatabase_name\tdatabase_version_or_date\tcommand_or_source\taccess_date" > "$VERSION_FILE"

datasets_version="$("$MAMBA_BIN" run --prefix "$ENV_PREFIX" datasets version 2>&1 | tail -n 1 | tr '\t' ' ')"
dataformat_version="$("$MAMBA_BIN" run --prefix "$ENV_PREFIX" dataformat version 2>&1 | tail -n 1 | tr '\t' ' ')"
seqkit_version="$("$MAMBA_BIN" run --prefix "$ENV_PREFIX" seqkit version 2>&1 | tail -n 1 | tr '\t' ' ')"
blast_version="$("$MAMBA_BIN" run --prefix "$ENV_PREFIX" blastn -version 2>&1 | head -n 1 | tr '\t' ' ')"
python_version="$("$MAMBA_BIN" run --prefix "$ENV_PREFIX" python --version 2>&1 | tail -n 1 | tr '\t' ' ')"

print "NCBI Datasets\t${datasets_version}\tNA\tNA\tconda-forge/bioconda\t${ACCESS_DATE}" >> "$VERSION_FILE"
print "NCBI dataformat\t${dataformat_version}\tNA\tNA\tconda-forge/bioconda\t${ACCESS_DATE}" >> "$VERSION_FILE"
print "SeqKit\t${seqkit_version}\tNA\tNA\tbioconda\t${ACCESS_DATE}" >> "$VERSION_FILE"
print "BLAST+\t${blast_version}\tNA\tNA\tbioconda\t${ACCESS_DATE}" >> "$VERSION_FILE"
print "Python\t${python_version}\tNA\tNA\tconda-forge\t${ACCESS_DATE}" >> "$VERSION_FILE"

print "Recorded core program versions in $VERSION_FILE"
