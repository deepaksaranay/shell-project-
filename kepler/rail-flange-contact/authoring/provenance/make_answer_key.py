#!/usr/bin/env python3
"""tests/answer_key.json from the 0.3 mm generator run (layout differs from the output file)."""
import json, os
here = os.path.dirname(os.path.abspath(__file__))
ref = json.load(open(os.path.join(here, "truth_el0.3.json")))
case = json.load(open(os.path.join(here, "..", "..", "environment", "data", "case.json")))
loads = {c["id"]: c["vertical_load_N"] for c in case["load_cases"]}
key = {"cases": {}, "limits": {"normal_force_vs_Q": 0.01, "creep_force_vs_Q": 0.01, "contact_angle_deg": 0.3, "lateral_shift_mm": 0.01},
       "provenance": "authoring/provenance/twopoint.py, 0.3 mm elements; the 0.4 mm run differs by at most 0.16% of Q"}
for c in ref["cases"]:
    e = {"Q": loads[c["id"]], "shift": round(c["lateral_shift_mm"], 6), "patches": {}}
    for k in ("tread", "flange"):
        p = c[k]
        e["patches"][k] = None if p is None else [round(p["contact_angle_deg"], 4), round(p["normal_force_N"], 2), round(p["Fx_N"], 2), round(p["Fs_N"], 2)]
    key["cases"][c["id"]] = e
json.dump(key, open(os.path.join(here, "..", "..", "tests", "answer_key.json"), "w"), indent=1)
print("ok")
