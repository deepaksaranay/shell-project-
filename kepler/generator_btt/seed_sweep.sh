#!/bin/bash
# usage: SCRATCH=/some/dir generator_btt/seed_sweep.sh SEED   (run from kepler/)
# Regenerates the dataset on a new seed and scores the reference, the
# steady-state fit and the no-Southwell variant against the stated limits.
set -e
SP=${SCRATCH:-/tmp/btt_seeds}; S=$1; W=$SP/s$S
mkdir -p $W/blade-tip-timing/environment/data $W/blade-tip-timing/tests $W/gen
sed -e "s/default_rng(7710)/default_rng($S)/" generator_btt/generate.py > $W/gen/generate.py
python3 $W/gen/generate.py > $W/gen.log
python3 generator_btt/check.py $W/blade-tip-timing | sed "s/^/seed $S  /"
