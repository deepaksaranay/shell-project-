#!/bin/bash
# usage: run_seed.sh SEED
set -e
SP=${SCRATCH:-/tmp/btt_seeds}; mkdir -p $SP; K=/home/user/shell-project-/kepler; S=$1; W=$SP/s$S
mkdir -p $W/blade-tip-timing/environment/data $W/blade-tip-timing/tests $W/gen
sed -e "s/default_rng(7710)/default_rng($S)/" $K/generator_btt/generate.py > $W/gen/generate.py
python3 $W/gen/generate.py > $W/gen.log
python3 - "$W" <<'PY'
import sys, json, numpy as np
W = sys.argv[1]
sys.path.insert(0, "/home/user/shell-project-/kepler/blade-tip-timing/solution"); import solve as S
T = json.load(open(f"{W}/blade-tip-timing/tests/expected.json"))
for tag, sw in [("reference", True), ("no_southwell", False)]:
    o = S.analyse(f"{W}/blade-tip-timing/environment/data", use_southwell=sw)
    ef = max(abs(o["blades"][k]["resonance_frequency_Hz"]/v["resonance_frequency_Hz"]-1) for k,v in T["blades"].items())
    ez = [o["blades"][k]["damping_ratio"]/v["damping_ratio"]-1 for k,v in T["blades"].items()]
    ok = o["engine_order"]==T["engine_order"] and ef<=0.0025 and max(map(abs,ez))<=0.25 and abs(np.mean(ez))<=0.07
    print(f"seed {W[-4:]} {tag:12s} EO={o['engine_order']} df={ef:.5f} dzmax={max(map(abs,ez)):.3f} dzmean={np.mean(ez):+.3f} PASS={ok}")
PY
