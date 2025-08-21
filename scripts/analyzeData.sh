#!/bin/bash

export SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Backward-compatible wrapper to the new evtDisplay.sh script
echo "[analyzeData.sh] Deprecated: please use scripts/evtDisplay.sh instead." >&2
exec "$SCRIPTS_DIR/evtDisplay.sh" "$@"
