#!/usr/bin/env bash
# Oracle: size the vessel with the independent transient model in solve.py.
set -eu
here=$(dirname "$(readlink -f "$0")")
exec python3 "$here/solve.py"
