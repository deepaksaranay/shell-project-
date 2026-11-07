#!/usr/bin/env python3
"""Shortcut models, to show the tolerances reject them."""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import twopoint
name = sys.argv[1]
el = 1.0 if name == "coarse" else 0.4
m = twopoint.Model(el_mm=el, opts={} if name == "coarse" else {name: True})
res = []
for cs in m.case["load_cases"]:
    ys = m.case["initial_guess"]["lateral_shift_mm"]
    g = m.vgap(ys, 0.0); zc = -np.min(g[np.isfinite(g)]) - 0.06
    res.append(twopoint.solve_case(m, cs, ys, zc, verbose=False))
json.dump({"element_mm": el, "cases": res}, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), f"baseline_{name}.json"), "w"), indent=1)
print(name, "done")
