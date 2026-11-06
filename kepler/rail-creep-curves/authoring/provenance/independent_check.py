"""Independent correct implementation, scored through the verifier.

Differs from both the reference-value generator (transient scheme) and the
shipped solution: the steady-rolling equations s = w + K q are solved
directly with an active-set method (stick elements by a linear solve, slip
elements with iterated slip directions), on a 0.45 mm grid whose origin is
shifted by a quarter element relative to the reference grids.
Writes creep_forces.json in the task's output format.
"""
import json, math, sys
import numpy as np
import problem as P
from contact import Grid, solve_normal, solve_tangential

el = 0.45
yr, g = P.gap_function()
x0, x1, y0, y1 = -8.0 + el / 4, 8.0 + el / 4, -18.0 + el / 4, 18.0 + el / 4
nx, ny = int(round((x1 - x0) / el)), int(round((y1 - y0) / el))
gr = Grid(x0 * 1e-3, x1 * 1e-3, nx, y0 * 1e-3, y1 * 1e-3, ny)
h = (np.interp(gr.Y * 1e3, yr, g) + (gr.X * 1e3) ** 2 / (2 * P.RW_MM)) * 1e-3
p, delta, _ = solve_normal(gr, h, P.N, P.E, P.NU)
out = {"contact_area_mm2": float((p > 0).sum() * gr.dA * 1e6), "max_pressure_MPa": float(p.max() / 1e6), "cases": []}
K = None
for cid, xi, eta, phi in P.CASES:
    r = solve_tangential(gr, p, P.E, P.NU, P.MU, xi, eta, phi, K=K)
    K = r["K"]
    out["cases"].append({"id": cid, "Fx_N": float(r["Fx"]), "Fy_N": float(r["Fy"])})
    print(cid, round(r["Fx"], 1), round(r["Fy"], 1), "iterations", r["iterations"], "converged", r["converged"], flush=True)
json.dump(out, open(sys.argv[1] if len(sys.argv) > 1 else "independent_output.json", "w"), indent=2)
print(json.dumps({k: v for k, v in out.items() if k != "cases"}))
