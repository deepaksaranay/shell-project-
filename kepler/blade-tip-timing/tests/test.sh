#!/bin/bash
mkdir -p /logs/verifier

run_tests() {
  # Preferred: isolated test environment via uv. Fallback: the image's own
  # pytest (installed in the Dockerfile) if uv cannot be fetched.
  if command -v curl >/dev/null 2>&1 || (apt-get update >/dev/null 2>&1 && apt-get install -y curl >/dev/null 2>&1); then
    curl -LsSf https://astral.sh/uv/0.9.5/install.sh 2>/dev/null | sh >/dev/null 2>&1
    [ -f "$HOME/.local/bin/env" ] && source "$HOME/.local/bin/env"
    if command -v uvx >/dev/null 2>&1; then
      uvx -p 3.13 -w pytest==8.4.1 -w numpy==2.3.3 pytest /tests/test_outputs.py -rA
      return $?
    fi
  fi
  python3 -m pytest /tests/test_outputs.py -rA
}

run_tests
if [ $? -eq 0 ]; then
  echo 1 > /logs/verifier/reward.txt
else
  echo 0 > /logs/verifier/reward.txt
fi
