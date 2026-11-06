#!/usr/bin/env bash
# Surge-vessel design check: simulate the submitted design and record 0/1.
set -u
LOGDIR=/logs/verifier
mkdir -p "$LOGDIR"
reward=0
if python3 -m pytest -q -rA -p no:cacheprovider --ctrf "$LOGDIR/ctrf.json" /tests/test_outputs.py; then
  reward=1
fi
printf '%s\n' "$reward" > "$LOGDIR/reward.txt"
