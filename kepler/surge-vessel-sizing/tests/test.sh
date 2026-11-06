#!/usr/bin/env bash
# Surge-vessel design check: simulate the submitted design and record 0/1.
set -u
mkdir -p /logs/verifier
echo 0 > /logs/verifier/reward.txt
if python3 -m pytest -q -rA -p no:cacheprovider --ctrf /logs/verifier/ctrf.json /tests/test_outputs.py; then
  echo 1 > /logs/verifier/reward.txt
else
  echo 0 > /logs/verifier/reward.txt
fi
