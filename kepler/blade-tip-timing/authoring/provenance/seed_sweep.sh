#!/bin/bash
# usage (from the bundle root): SCRATCH=/some/dir authoring/provenance/seed_sweep.sh OFFSET
# Regenerates every run configuration (example + held-out rotors) with the
# seeds shifted by OFFSET into SCRATCH, then scores the reference, the
# steady-state fit and the no-Southwell variant against the stated limits.
set -e
W=${SCRATCH:-/tmp/btt_seeds}/o$1
mkdir -p "$W"
for n in example rotor_a rotor_b rotor_c; do
  BTT_ROOT="$W" BTT_SEED_OFFSET=$1 BTT_ONLY=$n python3 authoring/provenance/generate.py > "$W/gen_$n.log" &
done
wait
python3 authoring/provenance/check.py "$W" | sed "s/^/offset $1  /"
