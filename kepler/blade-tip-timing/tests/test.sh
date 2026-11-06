#!/bin/bash
mkdir -p /logs/verifier

apt-get update >/dev/null 2>&1
apt-get install -y curl >/dev/null 2>&1
curl -LsSf https://astral.sh/uv/0.9.5/install.sh | sh >/dev/null 2>&1
source "$HOME/.local/bin/env"

uvx -p 3.13 -w pytest==8.4.1 -w numpy==2.3.3 pytest /tests/test_outputs.py -rA

if [ $? -eq 0 ]; then
  echo 1 > /logs/verifier/reward.txt
else
  echo 0 > /logs/verifier/reward.txt
fi
