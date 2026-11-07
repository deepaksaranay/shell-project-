#!/usr/bin/env python3
"""Worst error per quantity for a result file, as value and as share of the limit."""
import json, sys, os
key = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tests", "answer_key.json")))
d = json.load(open(sys.argv[1]))
snaps = d["snapshots"]; peak = d["peak"]["peak_bondline_temperature_K"] if "peak" in d else d["peak_bondline_temperature_K"]
for q, lim in key["limits"].items():
    if q == "peak_bondline_K":
        e = peak - key["peak_bondline_K"]; tol = lim["abs"]
        print(f"{q:28s} err {e:+.3f}  ({abs(e)/tol:.0%} of limit)"); continue
    worst = max(((s[q] - key["series"][q][str(int(round(s["t_s"])))], lim["abs"] + lim["rel"] * abs(key["series"][q][str(int(round(s["t_s"])))]), s["t_s"]) for s in snaps), key=lambda w: abs(w[0]) / w[1])
    print(f"{q:28s} err {worst[0]:+.4f} at {worst[2]} s ({abs(worst[0])/worst[1]:.0%} of limit)")
