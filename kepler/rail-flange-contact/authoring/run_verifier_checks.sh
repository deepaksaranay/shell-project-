#!/bin/sh
# Runs the real tests/test_outputs.py (unchanged except for the result path)
# on every file in authoring/evidence/verifier_inputs/, plus a missing file,
# a symlink to the answer key and a copy of the answer key, and records the
# pytest output in authoring/evidence/verifier_runs.log.
set -u
HERE=$(cd "$(dirname "$0")/.." && pwd)
W=$(mktemp -d)
cp -r "$HERE/tests" "$W/tests"
mkdir -p "$W/out"
sed -i "s#/app/output/flange_contact.json#$W/out/flange_contact.json#" "$W/tests/test_outputs.py"
LOG="$HERE/authoring/evidence/verifier_runs.log"
: > "$LOG"
run() {
    echo "=== $1" >> "$LOG"
    (cd "$W/tests" && python3 -m pytest -q -rA test_outputs.py 2>&1 | grep -E "^(PASSED|FAILED|ERROR)|passed|failed|error" ) >> "$LOG"
}
rm -f "$W/out/flange_contact.json"; run "no output file"
ln -s "$W/tests/answer_key.json" "$W/out/flange_contact.json"; run "symlink to the answer key"; rm "$W/out/flange_contact.json"
cp "$HERE/tests/answer_key.json" "$W/out/flange_contact.json"; run "copy of the answer key"
for f in "$HERE"/authoring/evidence/verifier_inputs/*.json; do
    cp "$f" "$W/out/flange_contact.json"; run "$(basename "$f")"
done
rm -rf "$W"
cat "$LOG"
