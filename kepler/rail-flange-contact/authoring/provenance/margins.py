#!/usr/bin/env python3
"""Per-case errors of a result file against tests/answer_key.json, as share of each limit."""
import json, math, os, sys
key = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tests", "answer_key.json")))
L = key["limits"]
d = json.load(open(sys.argv[1]))
worst = {}
for c in d["cases"]:
    r = key["cases"][c["id"]]; Q = r["Q"]
    row = [c["id"], "shift %+.4f" % (c["lateral_shift_mm"] - r["shift"])]
    worst["shift"] = max(worst.get("shift", 0), abs(c["lateral_shift_mm"] - r["shift"]) / L["lateral_shift_mm"])
    for k in ("tread", "flange"):
        p, rp = c[k], r["patches"][k]
        if p is None or rp is None:
            row.append(f"{k}: {p is None}/{rp is None}"); continue
        da = p["contact_angle_deg"] - rp[0]; dN = p["normal_force_N"] - rp[1]
        dF = math.hypot(p["Fx_N"] - rp[2], p["Fs_N"] - rp[3])
        worst["angle"] = max(worst.get("angle", 0), abs(da) / L["contact_angle_deg"])
        worst["N"] = max(worst.get("N", 0), abs(dN) / (L["normal_force_vs_Q"] * Q))
        worst["F"] = max(worst.get("F", 0), dF / (L["creep_force_vs_Q"] * Q))
        row.append(f"{k}: ang {da:+.3f} N {dN:+.0f} F {dF:.0f}")
    print("  ".join(row))
print("worst share of limit:", {k: f"{v:.0%}" for k, v in worst.items()})
