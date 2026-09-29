#!/bin/zsh
set -euo pipefail

PROJECT="${0:A:h:h}"
# Bioconda generates launchers that are not safe when their prefix contains
# spaces. Keep the scientific environment in a stable, no-space location.
ENV_PREFIX="${RHIZOFORGE_ENV_PREFIX:-$PROJECT/.envs/rhizoforge-select}"
LOCAL_MICROMAMBA="$PROJECT/.tools/bin/micromamba"
export MAMBA_ROOT_PREFIX="$PROJECT/.tools/rhizoforge-mamba-root"
export XDG_CACHE_HOME="$PROJECT/.tools/rhizoforge-cache"
mkdir -p "$MAMBA_ROOT_PREFIX" "$XDG_CACHE_HOME"

if [[ -x "$LOCAL_MICROMAMBA" ]]; then
  SOLVER="$LOCAL_MICROMAMBA"
elif command -v micromamba >/dev/null 2>&1; then
  SOLVER="$(command -v micromamba)"
elif command -v mamba >/dev/null 2>&1; then
  SOLVER="$(command -v mamba)"
elif command -v conda >/dev/null 2>&1; then
  SOLVER="$(command -v conda)"
else
  echo "Install micromamba, mamba, or conda first." >&2
  exit 1
fi

"$SOLVER" create --yes --no-rc --override-channels \
  --channel conda-forge --channel bioconda \
  --prefix "$ENV_PREFIX" --file "$PROJECT/environment-full.yml"
echo "Full RhizoForge-Select environment: $ENV_PREFIX"
echo "Some tools require their licensed/current reference databases; see docs/FULL_PIPELINE.md."
