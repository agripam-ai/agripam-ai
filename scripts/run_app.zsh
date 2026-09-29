#!/bin/zsh
set -euo pipefail

PROJECT="${0:A:h:h}"
PYTHON="$PROJECT/.envs/cctyper-local/bin/python"
VENDOR="$PROJECT/.app_vendor"

if [[ ! -d "$VENDOR/streamlit" ]]; then
  echo "Streamlit is not installed. Run: $PROJECT/scripts/install_app.zsh" >&2
  exit 1
fi

export PYTHONPATH="$VENDOR${PYTHONPATH:+:$PYTHONPATH}"
exec "$PYTHON" -m streamlit run "$PROJECT/app/streamlit_app.py" \
  --global.developmentMode false \
  "$@"
