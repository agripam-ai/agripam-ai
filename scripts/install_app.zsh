#!/bin/zsh
set -euo pipefail

PROJECT="${0:A:h:h}"
PYTHON="$PROJECT/.envs/cctyper-local/bin/python"
VENDOR="$PROJECT/.app_vendor"

mkdir -p "$VENDOR"
"$PYTHON" -m pip install \
  --disable-pip-version-check \
  --upgrade \
  --target "$VENDOR" \
  -r "$PROJECT/requirements-app.txt"

echo "Application dependencies installed in: $VENDOR"
