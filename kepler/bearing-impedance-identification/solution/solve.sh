#!/bin/bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p /app/output
python3 "$SCRIPT_DIR/solve.py" --data /app/data --out /app/output/bearing_coefficients.json
