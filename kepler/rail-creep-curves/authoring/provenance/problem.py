"""Task problem definition + truth computation (run from this directory)."""
import json, math, sys, time
import numpy as np
from geometry import rail_profile, wheel_profile
from contact import Grid, solve_normal, solve_tangential_rolling

E, NU, N, MU, RW_MM = 210e9, 0.28, 110e3, 0.35, 460.0
CASES = [  # (id, xi, eta, phi [1/m])
    ("A", 0.0005, 0.0, 0.0),
    ("B", 0.0020, 0.0, 0.0),
    ("C", 0.0, 0.0010, 0.0),
    ("D", 0.0, 0.0030, 0.0),
    ("E", 0.0010, 0.0010, 0.0),
    ("F", 0.0, 0.0, 0.50),
    ("G", 0.0010, 0.0005, 0.30),
    ("H", -0.0015, 0.0, -0.60),
]


def gap_function():
    yr, zr = rail_profile(); yw, zw = wheel_profile()
    g = zw - zr
    return yr, g - g.min()


def solve(el_mm, cases=CASES, verbose=True):
    yr, g = gap_function()
    x0, x1, y0, y1 = -8.0, 8.0, -18.0, 18.0
    nx, ny = int(round((x1 - x0) / el_mm)), int(round((y1 - y0) / el_mm))
    gr = Grid(x0 * 1e-3, x1 * 1e-3, nx, y0 * 1e-3, y1 * 1e-3, ny)
    h = (np.interp(gr.Y * 1e3, yr, g) + (gr.X * 1e3) ** 2 / (2 * RW_MM)) * 1e-3
    p, delta, _ = solve_normal(gr, h, N, E, NU)
    out = {"element_mm": el_mm, "area_mm2": float((p > 0).sum() * gr.dA * 1e6),
           "p_max_MPa": float(p.max() / 1e6), "approach_um": float(delta * 1e6), "cases": {}}
    for cid, xi, eta, phi in cases:
        t = time.time()
        r = solve_tangential_rolling(gr, p, E, NU, MU, xi, eta, phi)
        out["cases"][cid] = {"Fx_N": float(r["Fx"]), "Fy_N": float(r["Fy"]),
                             "adhesion_fraction": r["adhesion_fraction"], "steps": r["steps"],
                             "converged": bool(r["converged"])}
        if verbose:
            print(f"  el={el_mm} {cid}: Fx={r['Fx']:9.1f} Fy={r['Fy']:9.1f} N  adh={r['adhesion_fraction']:.2f} "
                  f"steps={r['steps']} conv={r['converged']} ({time.time()-t:.0f}s)", flush=True)
    return out


if __name__ == "__main__":
    el = float(sys.argv[1])
    res = solve(el)
    json.dump(res, open(f"truth_el{el}.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != "cases"}))
