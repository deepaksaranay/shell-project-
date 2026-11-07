#!/bin/sh
# Compare the reported two-point contact forces with the reference values.
mkdir -p /logs/verifier
cd /tests || { echo 0 > /logs/verifier/reward.txt; exit 0; }
if python3 -m pytest /tests/test_outputs.py -rA --ctrf /logs/verifier/ctrf.json; then
    echo 1 > /logs/verifier/reward.txt
else
    echo 0 > /logs/verifier/reward.txt
fi
