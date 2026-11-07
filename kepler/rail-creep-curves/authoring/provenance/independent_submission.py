#!/usr/bin/env python3
"""Structurally independent full submission (written separately from the
answer-key generator and from solution/solve.py), scored through the real
verifier. Every stage uses a different algorithm or discretisation:

  * elements: non-square rectangles, 0.30 mm (rolling direction) x 0.40 mm,
    on a grid not aligned with any other implementation;
  * influence products: zero-padded FFT convolution of the closed-form
    rectangle kernels (DC-FFT) instead of dense matrices;
  * normal problem: Polonsky-Keer constrained conjugate gradient with the
    load as a side condition (instead of an active-set direct solve);
  * tangential problem: each rolling step of Kalker's transient scheme is
    solved by ADMM (operator splitting with one Cholesky factorisation,
    disk projection) instead of projected gradient or a steady active set;
    upstream displacements come from the FFT convolution of the previous
    traction field.

Usage: python3 independent_submission.py <data_dir> <output_json> [case_ids]
"""
import glob
import json
import math
import os
import sys

import numpy as np
from scipy.linalg import cho_factor, cho_solve

DX, DY = 0.30e-3, 0.40e-3


def rect_kernel(nx, ny, kind):
    """Kernel of offsets (i*DX, j*DY), i in [-(nx-1), nx-1], j in [-(ny-1), ny-1]."""
    i = np.arange(-(nx - 1), nx) * DX
    j = np.arange(-(ny - 1), ny) * DY
    X, Y = np.meshgrid(i, j, indexing="ij")
    def xa(a, b):
        out = np.zeros_like(a); nz = a != 0
        out[nz] = a[nz] * np.arcsinh(b[nz] / np.abs(a[nz]))
        return out
    k = np.zeros_like(X)
    for sx in (1, -1):
        for sy in (1, -1):
            a, b = X + sx * DX / 2, Y + sy * DY / 2
            s = sx * sy
            if kind == "1/r":
                k += s * (xa(a, b) + xa(b, a))
            elif kind == "xx":
                k += s * xa(b, a)
            elif kind == "yy":
                k += s * xa(a, b)
            else:
                k -= s * np.hypot(a, b)
    return k


class Conv:
    """Linear convolution of an (nx, ny) field with a kernel via padded FFT."""
    def __init__(self, kernel, nx, ny):
        self.nx, self.ny = nx, ny
        self.shape = (2 * nx - 1 + nx, 2 * ny - 1 + ny)
        self.kf = np.fft.rfft2(kernel, s=self.shape)

    def __call__(self, f):
        full = np.fft.irfft2(self.kf * np.fft.rfft2(f, s=self.shape), s=self.shape)
        return full[self.nx - 1: 2 * self.nx - 1, self.ny - 1: 2 * self.ny - 1]


def main(data_dir, out_path):
    case = json.load(open(os.path.join(data_dir, "contact_case.json")))
    rail = np.loadtxt(os.path.join(data_dir, "rail_profile.csv"), delimiter=",", skiprows=1)
    wheel = np.loadtxt(os.path.join(data_dir, "wheel_profile.csv"), delimiter=",", skiprows=1)
    E = case["materials"]["rail"]["youngs_modulus_Pa"]; nu = case["materials"]["rail"]["poisson_ratio"]
    G = E / (2 * (1 + nu)); N = case["normal_load_N"]; mu = case["coefficient_of_friction"]
    R = case["wheel_rolling_radius_mm"] * 1e-3

    y_mm = wheel[:, 0]
    sep = wheel[:, 1] - np.interp(y_mm, rail[:, 0], rail[:, 1])
    sep -= sep.min()
    yb = y_mm[sep < 0.25]
    y0, y1 = (yb.min() - 1.3) * 1e-3, (yb.max() + 1.3) * 1e-3
    nx = int(math.ceil(0.0172 / DX)); ny = int(math.ceil((y1 - y0) / DY))
    xc = (np.arange(nx) - (nx - 1) / 2) * DX + 0.07e-3
    yc = y0 + (np.arange(ny) + 0.5) * DY
    Xg, Yg = np.meshgrid(xc, yc, indexing="ij")
    h = np.interp(Yg * 1e3, y_mm, sep) * 1e-3 + Xg ** 2 / (2 * R)
    dA = DX * DY

    # ---- normal problem: Polonsky-Keer CG
    Kn = Conv((1 - nu) / (math.pi * G) * rect_kernel(nx, ny, "1/r"), nx, ny)
    p = np.where(h < 5e-5, N / dA / max(1, (h < 5e-5).sum()), 0.0)
    t = np.zeros_like(p); Gold = 1.0; delta = 0.0
    for it in range(5000):
        g = Kn(p) + h
        I = p > 0
        g = g - g[I].mean()
        Gc = np.sum(g[I] ** 2)
        t = np.where(I, g + delta * (Gc / Gold) * t, 0.0)
        Gold = Gc
        r = Kn(t)
        r = r - r[I].mean()
        tau = np.sum(g[I] * t[I]) / np.sum(r[I] * t[I])
        pold = p.copy()
        p = np.where(I, p - tau * t, p)
        p[p < 0] = 0.0
        Iol = (p == 0) & (g < 0)
        if Iol.any():
            delta = 0.0
            p[Iol] -= tau * g[Iol]
        else:
            delta = 1.0
        p *= N / (p.sum() * dA)
        if np.abs(p - pold).sum() * dA / N < 1e-11:
            break
    area = (p > 0).sum() * dA

    # ---- tangential problem: transient scheme, ADMM per step
    c = 1 / (math.pi * G)
    k1, kxx, kyy, kxy = (rect_kernel(nx, ny, k) for k in ("1/r", "xx", "yy", "xy"))
    kernels = {"xx": c * ((1 - nu) * k1 + nu * kxx), "yy": c * ((1 - nu) * k1 + nu * kyy), "xy": c * nu * kxy}
    convs = {k: Conv(v, nx, ny) for k, v in kernels.items()}
    I = np.nonzero(p > 0)
    ci, cj = I
    m = len(ci)
    di = ci[:, None] - ci[None, :] + nx - 1
    dj = cj[:, None] - cj[None, :] + ny - 1
    U = np.block([[kernels["xx"][di, dj], kernels["xy"][di, dj]],
                  [kernels["xy"][di, dj], kernels["yy"][di, dj]]])
    rho = float(np.mean(np.diag(U)))
    fac = cho_factor(U + rho * np.eye(2 * m))
    gbound = mu * p[I]
    up_ok = ci + 1 < nx

    def upstream_disp(qx, qy):
        Qx = np.zeros((nx, ny)); Qy = np.zeros((nx, ny))
        Qx[I] = qx; Qy[I] = qy
        ux = convs["xx"](Qx) + convs["xy"](Qy)
        uy = convs["xy"](Qx) + convs["yy"](Qy)
        iu = np.minimum(ci + 1, nx - 1)
        return np.concatenate([np.where(up_ok, ux[iu, cj], 0.0), np.where(up_ok, uy[iu, cj], 0.0)])

    def project(v):
        vx, vy = v[:m], v[m:]
        n = np.hypot(vx, vy)
        f = np.where(n > gbound, gbound / np.maximum(n, 1e-300), 1.0)
        return np.concatenate([vx * f, vy * f])

    # checkpoints: one file per worker; a run with a case subset (argv[3],
    # e.g. "DE") only solves those cases, so cases can run in parallel
    subset = sys.argv[3] if len(sys.argv) > 3 else None
    ckpt = out_path + (f".partial.{subset}.json" if subset else ".partial.json")
    done = {}
    for f in sorted(glob.glob(out_path + ".partial*.json")):
        done.update(json.load(open(f)))
    results = []
    for cs in case["creepage_cases"]:
        if cs["id"] in done:  # resume after an interruption
            results.append(done[cs["id"]])
            continue
        if subset and cs["id"] not in subset:
            continue
        w = np.concatenate([cs["xi"] - cs["phi_per_m"] * Yg[I], cs["eta"] + cs["phi_per_m"] * Xg[I]]) * DX
        q = np.zeros(2 * m); z = np.zeros(2 * m); lam = np.zeros(2 * m)
        for step in range(4000):
            b = w - upstream_disp(q[:m], q[m:])
            for _ in range(3000):
                qn = cho_solve(fac, rho * (z - lam) - b)
                zn = project(qn + lam)
                lam = lam + qn - zn
                rp = np.max(np.abs(qn - zn)); rd = np.max(np.abs(zn - z))
                z = zn
                if max(rp, rd) <= 1e-10 * gbound.max():
                    break
            change = np.max(np.abs(z - q))
            q = z.copy()
            if change <= 1e-9 * gbound.max() and step > 2:
                break
        results.append({"id": cs["id"], "Fx_N": float(q[:m].sum() * dA), "Fy_N": float(q[m:].sum() * dA)})
        done[cs["id"]] = results[-1]
        mine = json.load(open(ckpt)) if os.path.exists(ckpt) else {}
        mine[cs["id"]] = results[-1]
        json.dump(mine, open(ckpt, "w"))
        print(results[-1], "steps", step + 1, flush=True)
    if len(results) < len(case["creepage_cases"]):
        return
    out = {"contact_area_mm2": float(area * 1e6), "max_pressure_MPa": float(p.max() / 1e6), "cases": results}
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps({k: v for k, v in out.items() if k != "cases"}))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
