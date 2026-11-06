#!/bin/bash
# Writes exactly 0 or 1 to /logs/verifier/reward.txt; 1 only when every test passes.
mkdir -p /logs/verifier
echo 0 > /logs/verifier/reward.txt
export PYTHONDONTWRITEBYTECODE=1
if ! python3 -c "import pytest, ctrf" >/dev/null 2>&1; then
  pip install --no-cache-dir -q pytest==8.4.1 pytest-json-ctrf==0.3.5 >/dev/null 2>&1
fi
python3 -m pytest /tests/test_outputs.py -rA -p no:cacheprovider --ctrf /logs/verifier/ctrf.json
rc=$?
if [ "$rc" -eq 0 ]; then
  echo 1 > /logs/verifier/reward.txt
else
  echo 0 > /logs/verifier/reward.txt
fi
exit 0
