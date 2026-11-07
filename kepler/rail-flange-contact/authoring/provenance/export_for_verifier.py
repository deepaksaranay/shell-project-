#!/usr/bin/env python3
"""Write generator / baseline results in exactly the task's output schema
(/app/output/flange_contact.json layout), so the real verifier can be run on them.
Usage: python3 export_for_verifier.py <generator_or_baseline_json> <out_json>"""
import json, sys
d = json.load(open(sys.argv[1]))
cases = []
for c in d["cases"]:
    e = {"id": c["id"], "lateral_shift_mm": float(c["lateral_shift_mm"])}
    for k in ("tread", "flange"):
        p = c.get(k)
        e[k] = None if p is None else {f: float(p[f]) for f in ("contact_angle_deg", "normal_force_N", "Fx_N", "Fs_N")}
    cases.append(e)
json.dump({"cases": cases}, open(sys.argv[2], "w"), indent=2)
