#!/bin/bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Install the analysis program where it is expected, then run it on the example.
cp "$SCRIPT_DIR/solve.py" /app/btt_analyze.py
mkdir -p /app/output
python3 /app/btt_analyze.py /app/data/example /app/output/example_results.json
