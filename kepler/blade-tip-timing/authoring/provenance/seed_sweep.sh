#!/bin/bash
# usage (from the bundle root): SCRATCH=/some/dir authoring/provenance/seed_sweep.sh SEED
# Regenerates the dataset on a new seed into SCRATCH and scores the reference,
# the steady-state fit and the no-Southwell variant against the stated limits.
set -e
W=${SCRATCH:-/tmp/btt_seeds}/s$1
mkdir -p "$W"
BTT_SEED=$1 BTT_ROOT="$W" python3 authoring/provenance/generate.py > "$W/gen.log"
python3 authoring/provenance/check.py "$W" | sed "s/^/seed $1  /"
