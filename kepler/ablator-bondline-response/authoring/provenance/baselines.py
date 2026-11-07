#!/usr/bin/env python3
"""Shortcut / flawed models, to show the tolerances reject them."""
import json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ablation
name = sys.argv[1]
if name == "coarse":
    out = ablation.run(Na=12, beta=0.01, verbose=False)
else:
    out = ablation.run(Na=120, verbose=False, opts={name: True})
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), f"baseline_{name}.json"), "w"), indent=1)
print(name, "done")
