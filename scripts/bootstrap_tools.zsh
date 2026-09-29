#!/bin/zsh
set -euo pipefail

PROJECT_DIR="${0:A:h:h}"
TOOLS_DIR="$PROJECT_DIR/.tools"
BIN_DIR="$TOOLS_DIR/bin"
MAMBA_ROOT="$TOOLS_DIR/micromamba-root"
MAMBA_BIN="$BIN_DIR/micromamba"
ENV_PREFIX="$PROJECT_DIR/.envs/agripam-core"
DOWNLOAD_URL="https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-osx-arm64"

mkdir -p "$BIN_DIR" "$MAMBA_ROOT" "$PROJECT_DIR/.envs"

if [[ "$(uname -m)" != "arm64" ]]; then
  print -u2 "This bootstrap file is configured for an Apple Silicon Mac (arm64)."
  exit 1
fi

if [[ ! -x "$MAMBA_BIN" ]]; then
  print "Downloading the official Micromamba Apple Silicon binary..."
  curl -L --fail --retry 3 "$DOWNLOAD_URL" -o "$MAMBA_BIN"
  chmod 755 "$MAMBA_BIN"
fi

export MAMBA_ROOT_PREFIX="$MAMBA_ROOT"

print "Creating the isolated core bioinformatics environment..."
"$MAMBA_BIN" create --yes --prefix "$ENV_PREFIX" --file "$PROJECT_DIR/config/core-environment.yaml"

print "Verifying installed programs..."
"$MAMBA_BIN" run --prefix "$ENV_PREFIX" datasets version
"$MAMBA_BIN" run --prefix "$ENV_PREFIX" dataformat version
"$MAMBA_BIN" run --prefix "$ENV_PREFIX" seqkit version
"$MAMBA_BIN" run --prefix "$ENV_PREFIX" blastn -version
"$MAMBA_BIN" run --prefix "$ENV_PREFIX" python --version

print ""
print "Core environment created at:"
print "$ENV_PREFIX"
print ""
print "Run a command inside it with:"
print "$MAMBA_BIN run --prefix $ENV_PREFIX COMMAND"

