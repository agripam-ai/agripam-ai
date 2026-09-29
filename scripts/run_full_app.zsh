#!/bin/zsh
set -euo pipefail

PROJECT="${0:A:h:h}"
ENV_PREFIX="${RHIZOFORGE_ENV_PREFIX:-$PROJECT/.envs/rhizoforge-select}"
PYTHON="$ENV_PREFIX/bin/python"
if [[ ! -x "$PYTHON" ]]; then
  echo "Full environment missing. Run: $PROJECT/scripts/install_full_pipeline.zsh" >&2
  exit 1
fi
export PATH="$ENV_PREFIX/bin:$PATH"
export CCTYPER_DB="$PROJECT/data/databases/cctyper_1.8.0"
exec "$PYTHON" -m streamlit run "$PROJECT/app/streamlit_app.py" --global.developmentMode false "$@"
