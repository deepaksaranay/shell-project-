"""Shortcut methods, to check that the tolerances reject them.
1. 'hertz': replace the measured gap by its local quadratic (Hertz ellipse
   at the minimum-gap point), then solve exactly.
2. 'coarse': exact method on a 1.0 mm grid.
3. 'linear': Kalker linear theory (full adhesion) on the exact patch,
   capped at mu*N (no slip region).
"""
import json, sys, numpy as np
import problem as P
from contact import Grid, solve_normal, solve_tangential_rolling, _kmatrix

which = sys.argv[1]
if which == "coarse":
    res = P.solve(1.0, verbose=False)
elif which == "hertz":
    yr, g = P.gap_function()
    i = np.argmin(g); y0 = yr[i]
    sel = np.abs(yr - y0) <= 2.0
    c2 = np.polyfit(yr[sel] - y0, g[sel], 2)[0]       # g ~ c2 (y-y0)^2
    gq = c2 * (yr - y0) ** 2
    P.gap_function = lambda: (yr, gq)
    res = P.solve(0.5, verbose=False)
elif which == "linear":
    yr, g = P.gap_function()
    el = 0.5; x0, x1, y0_, y1 = -8.0, 8.0, -18.0, 18.0
    nx, ny = int(round((x1 - x0) / el)), int(round((y1 - y0_) / el))
    gr = Grid(x0 * 1e-3, x1 * 1e-3, nx, y0_ * 1e-3, y1 * 1e-3, ny)
    h = (np.interp(gr.Y * 1e3, yr, g) + (gr.X * 1e3) ** 2 / (2 * P.RW_MM)) * 1e-3
    p, d, _ = solve_normal(gr, h, P.N, P.E, P.NU)
    C = np.nonzero(p > 0)[0]; m = len(C)
    K = _kmatrix(gr, P.E, P.NU, C)
    res = {"cases": {}}
    for cid, xi, eta, phi in P.CASES:
        w = np.concatenate([xi - phi * gr.Y[C], eta + phi * gr.X[C]])
        q = np.linalg.solve(K, -w)
        F = np.array([q[:m].sum(), q[m:].sum()]) * gr.dA
        nF = np.hypot(*F)
        if nF > P.MU * P.N:
            F *= P.MU * P.N / nF
        res["cases"][cid] = {"Fx_N": float(F[0]), "Fy_N": float(F[1])}
json.dump(res, open(f"baseline_{which}.json", "w"), indent=1)
print(which, json.dumps(res)[:600])
