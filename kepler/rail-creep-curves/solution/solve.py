#!/usr/bin/env python3
"""Creep forces for a measured wheel/rail profile pair with Kalker's exact
(half-space, boundary-element) rolling-contact theory.

1. Undeformed separation h(x, y) = (wheel z - rail z)(y) + x^2 / (2 R_wheel).
2. Normal problem on a uniform grid of rectangular pressure elements
   (closed-form Love/Boussinesq influence numbers), active-set solve with
   the total load as a side condition.
3. Tangential problem for each creepage case: steady rolling obtained as the
   fixed point of the discrete transient rolling scheme with rolling step =
   element length. Each step is the convex problem
      min 1/2 q.Uq + q.(w dx - u_prev(upstream))  subject to |q_I| <= mu p_I,
   solved by accelerated projected gradient. Quasi-identical materials, so
   the normal and tangential problems are uncoupled.
"""
import json
import math
import os

import numpy as np

DATA = "/app/data"
OUT = "/app/output/creep_forces.json"
ELEMENT_MM = 0.4


def corner_sum(X, Y, hx, hy, kind):
    """Rectangle integrals of 1/r ('n'), x^2/r^3 ('xx'), y^2/r^3 ('yy'), xy/r^3 ('xy')."""
    def xa(a, b):
        r = np.zeros_like(a)
        nz = a != 0
        r[nz] = a[nz] * np.arcsinh(b[nz] / np.abs(a[nz]))
        return r
    tot = np.zeros_like(X)
    for sx in (1, -1):
        for sy in (1, -1):
            a = X + sx * hx / 2
            b = Y + sy * hy / 2
            sg = sx * sy
            if kind == "n":
                tot += sg * (xa(a, b) + xa(b, a))
            elif kind == "xx":
                tot += sg * xa(b, a)
            elif kind == "yy":
                tot += sg * xa(a, b)
            else:
                tot -= sg * np.hypot(a, b)
    return tot


def main():
    with open(os.path.join(DATA, "contact_case.json")) as f:
        case = json.load(f)
    rail = np.loadtxt(os.path.join(DATA, "rail_profile.csv"), delimiter=",", skiprows=1)
    wheel = np.loadtxt(os.path.join(DATA, "wheel_profile.csv"), delimiter=",", skiprows=1)
    E = case["materials"]["wheel"]["youngs_modulus_Pa"]
    nu = case["materials"]["wheel"]["poisson_ratio"]
    G = E / (2 * (1 + nu))
    N = case["normal_load_N"]
    mu = case["coefficient_of_friction"]
    Rw = case["wheel_rolling_radius_mm"] * 1e-3

    yy = rail[:, 0]
    sep = np.interp(yy, wheel[:, 0], wheel[:, 1]) - rail[:, 1]
    sep -= sep.min()
    y_c = yy[np.argmin(sep)]

    el = ELEMENT_MM * 1e-3
    # candidate region: separation below 0.3 mm laterally, +-9 mm longitudinally
    band = yy[sep < 0.3]
    ylo, yhi = (band.min() - 1) * 1e-3, (band.max() + 1) * 1e-3
    nxe = int(math.ceil(0.018 / el)); nye = int(math.ceil((yhi - ylo) / el))
    xs = -nxe * el / 2 + el * (np.arange(nxe) + 0.5)
    ys = ylo + el * (np.arange(nye) + 0.5)
    XX, YY = np.meshgrid(xs, ys, indexing="ij")
    Xf, Yf = XX.ravel(), YY.ravel()
    h = np.interp(Yf * 1e3, yy, sep) * 1e-3 + Xf ** 2 / (2 * Rw)

    def infl(rows, cols, kind):
        return corner_sum(Xf[rows][:, None] - Xf[cols][None, :], Yf[rows][:, None] - Yf[cols][None, :], el, el, kind)

    # ---- normal problem
    cand = np.nonzero(h < 2.5e-4)[0]
    Azz = (1 - nu) / (math.pi * G) * infl(cand, cand, "n")
    hc = h[cand]
    act = np.argsort(hc)[: max(8, len(hc) // 10)]
    for _ in range(300):
        k = len(act)
        M = np.zeros((k + 1, k + 1))
        M[:k, :k] = Azz[np.ix_(act, act)]
        M[:k, k] = -1.0
        M[k, :k] = el * el
        sol = np.linalg.solve(M, np.concatenate([-hc[act], [N]]))
        pa, delta = sol[:k], sol[k]
        if (pa < 0).any():
            act = act[pa >= 0]
            continue
        pc = np.zeros(len(cand)); pc[act] = pa
        gapc = hc + Azz @ pc - delta
        add = np.setdiff1d(np.nonzero(gapc < -1e-14)[0], act)
        if len(add) == 0:
            break
        act = np.sort(np.concatenate([act, add]))
    p = np.zeros(len(h)); p[cand] = pc

    # ---- tangential problem
    C = np.nonzero(p > 0)[0]
    m = len(C)
    ix, iy = C // nye, C % nye
    ok = ix + 1 < nxe
    upst = np.where(ok, (ix + 1) * nye + iy, C)
    c0 = 1 / (math.pi * G)
    def tang(rows):
        I1, Ixx, Iyy, Ixy = (infl(rows, C, k) for k in ("n", "xx", "yy", "xy"))
        return np.block([[c0 * ((1 - nu) * I1 + nu * Ixx), c0 * nu * Ixy],
                         [c0 * nu * Ixy, c0 * ((1 - nu) * I1 + nu * Iyy)]])
    U = tang(C)
    V = tang(upst)
    V[np.concatenate([~ok, ~ok])] = 0.0
    L = np.linalg.eigvalsh(U)[-1]
    g = mu * p[C]

    def project(v):
        vx, vy = v[:m], v[m:]
        t = np.hypot(vx, vy)
        f = np.where(t > g, g / np.maximum(t, 1e-300), 1.0)
        return np.concatenate([vx * f, vy * f])

    results = []
    for cs in case["creepage_cases"]:
        xi, eta, phi = cs["xi"], cs["eta"], cs["phi_per_m"]
        w = np.concatenate([xi - phi * Yf[C], eta + phi * Xf[C]]) * el
        q = np.zeros(2 * m)
        for step in range(3000):
            b = w - V @ q
            z = q.copy(); y = q.copy(); t = 1.0
            for _ in range(4000):
                znew = project(y - (U @ y + b) / L)
                tn = 0.5 * (1 + math.sqrt(1 + 4 * t * t))
                y = znew + (t - 1) / tn * (znew - z)
                done = np.max(np.abs(znew - z)) <= 1e-10 * g.max()
                z, t = znew, tn
                if done:
                    break
            change = np.max(np.abs(z - q))
            q = z
            if change <= 1e-9 * g.max() and step > 2:
                break
        results.append({"id": cs["id"], "Fx_N": float(q[:m].sum() * el * el), "Fy_N": float(q[m:].sum() * el * el)})
        print(results[-1], flush=True)

    out = {"contact_area_mm2": float(m * el * el * 1e6),
           "max_pressure_MPa": float(p.max() / 1e6),
           "cases": results}
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps({k: v for k, v in out.items() if k != "cases"}))


if __name__ == "__main__":
    main()
