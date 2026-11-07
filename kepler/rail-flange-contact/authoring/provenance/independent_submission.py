#!/usr/bin/env python3
"""Structurally independent submission (written separately from the
answer-key generator and from solution/solve.py), scored through the real
verifier. Different at every numerical stage:

  * contact points: minima of the vertical separation found by bounded
    scalar minimisation inside each bracketed dip (no fixed resampling);
  * elements: non-square 0.25 mm (rolling) x 0.35 mm (lateral) rectangles,
    patch grid offset from the contact point by a fraction of an element;
  * influence products: zero-padded FFT convolution of closed-form kernels;
  * normal problem: Polonsky-Keer constrained conjugate gradient
    (displacement control, no load constraint);
  * tangential problem: Kalker's transient scheme, each rolling step solved
    by ADMM (one Cholesky factorisation per patch, disk projection);
  * equilibrium: alternate (i) Newton on lateral shift / vertical position
    with the patch creep forces held fixed and (ii) re-solving the creep
    forces at the new position, until both stop changing.

Usage: python3 independent_submission.py <data_dir> <output_json>
"""
import json
import math
import os
import sys

import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize_scalar

DX, DS = 0.25e-3, 0.35e-3
XH = 9.0e-3


def kernel(nx, ny, kind):
    i = np.arange(-(nx - 1), nx) * DX
    j = np.arange(-(ny - 1), ny) * DS
    X, Y = np.meshgrid(i, j, indexing="ij")
    def xa(a, b):
        o = np.zeros_like(a); k = a != 0
        o[k] = a[k] * np.arcsinh(b[k] / np.abs(a[k]))
        return o
    out = np.zeros_like(X)
    for sx in (1, -1):
        for sy in (1, -1):
            a, b = X + sx * DX / 2, Y + sy * DS / 2
            sg = sx * sy
            if kind == "1/r":
                out += sg * (xa(a, b) + xa(b, a))
            elif kind == "xx":
                out += sg * xa(b, a)
            elif kind == "yy":
                out += sg * xa(a, b)
            else:
                out -= sg * np.hypot(a, b)
    return out


class Conv:
    def __init__(self, k, nx, ny):
        self.nx, self.ny = nx, ny
        self.shape = (3 * nx, 3 * ny)
        self.kf = np.fft.rfft2(k, s=self.shape)

    def __call__(self, f):
        r = np.fft.irfft2(self.kf * np.fft.rfft2(f, s=self.shape), s=self.shape)
        return r[self.nx - 1:2 * self.nx - 1, self.ny - 1:2 * self.ny - 1]


def main(d, out_path):
    case = json.load(open(os.path.join(d, "case.json")))
    R = np.loadtxt(os.path.join(d, "rail_profile.csv"), delimiter=",", skiprows=1)
    W = np.loadtxt(os.path.join(d, "wheel_profile.csv"), delimiter=",", skiprows=1)
    yr, zr, yw, zw = R[:, 0], R[:, 1], W[:, 0], W[:, 1]
    sr_slope = np.gradient(zr, yr)
    E = case["materials"]["youngs_modulus_Pa"]; nu = case["materials"]["poisson_ratio"]; G = E / (2 * (1 + nu))
    mu = case["coefficient_of_friction"]; r0 = case["wheel"]["tape_circle_radius_mm"]
    zt = float(np.interp(case["wheel"]["tape_circle_y_mm"], yw, zw))

    def gfun(y, ys, zc):
        return np.interp(y - ys, yw, zw) + zc - np.interp(y, yr, zr)

    def contact_points(ys, zc):
        split = case["frames"]["region_split_y_w_mm"]
        pts = []
        for lo, hi in ((yw[0], split), (split, yw[-1])):
            a, b = max(yr[0], lo + ys), min(yr[-1], hi + ys)
            if b - a < 0.05:
                continue
            y = np.union1d(yr[(yr >= a) & (yr <= b)], yw[(yw + ys >= a) & (yw + ys <= b)] + ys)
            g = gfun(y, ys, zc)
            k = int(np.argmin(g))
            r = minimize_scalar(lambda t: gfun(t, ys, zc), bounds=(y[max(k - 1, 0)], y[min(k + 1, len(y) - 1)]),
                                method="bounded", options={"xatol": 1e-7})
            pts.append((r.x, math.atan(np.interp(r.x, yr, sr_slope)), a, b))
        return pts

    def build_patch(pt, ys, zc):
        y0, gm, ylo, yhi = pt
        c, s = math.cos(gm), math.sin(gm); z0 = np.interp(y0, yr, zr)
        def rot(y, z):
            return (y - y0) * c + (z - z0) * s, -(y - y0) * s + (z - z0) * c
        mr = (yr >= ylo) & (yr <= yhi); ywt = yw + ys; mw = (ywt >= ylo) & (ywt <= yhi)
        a, b = rot(yr[mr], zr[mr]); e, f = rot(ywt[mw], zw[mw] + zc)
        oa, oe = np.argsort(a), np.argsort(e)
        lo, hi = max(a.min(), e.min()), min(a.max(), e.max())
        sgrid = np.arange(lo, hi, 0.05)
        hs = (np.interp(sgrid, e[oe], f[oe]) - np.interp(sgrid, a[oa], b[oa])) * 1e-3
        k = np.nonzero(hs < 2.5e-4)[0]
        if len(k) == 0:
            return None
        s_lo, s_hi = sgrid[k[0]] * 1e-3, sgrid[k[-1]] * 1e-3
        ri = r0 - (np.interp(y0 - ys, yw, zw) - zt)
        Rx = ri * 1e-3 / c
        nx = int(round(2 * XH / DX)); ny = int(math.ceil((s_hi - s_lo) / DS)) + 2
        xc = -XH + DX * (np.arange(nx) + 0.5) + 0.31 * DX
        sc = s_lo - DS + DS * (np.arange(ny) + 0.5)
        Xg, Sg = np.meshgrid(xc, sc, indexing="ij")
        H = np.interp(Sg, sgrid * 1e-3, hs, left=1.0, right=1.0) + Xg ** 2 / (2 * Rx)
        return dict(X=Xg, S=Sg, H=H, nx=nx, ny=ny, r=ri, gamma=gm)

    kcache = {}
    def convs(nx, ny):
        if (nx, ny) not in kcache:
            k1, kxx, kyy, kxy = (kernel(nx, ny, t) for t in ("1/r", "xx", "yy", "xy"))
            cc = 1 / (math.pi * G)
            kcache[(nx, ny)] = dict(n=Conv((1 - nu) * cc * k1, nx, ny),
                                    raw={"xx": cc * ((1 - nu) * k1 + nu * kxx), "yy": cc * ((1 - nu) * k1 + nu * kyy), "xy": cc * nu * kxy})
            kcache[(nx, ny)]["c"] = {k: Conv(v, nx, ny) for k, v in kcache[(nx, ny)]["raw"].items()}
        return kcache[(nx, ny)]

    def normal_cg(P):
        H = P["H"]; K = convs(P["nx"], P["ny"])["n"]
        if H.min() >= 0:
            return np.zeros_like(H)
        Estar = E / (2 * (1 - nu ** 2))
        p = np.where(H < 0, -H * Estar / math.sqrt(max((H < 0).sum(), 1) * DX * DS), 0.0)
        t = np.zeros_like(p); Gold = 1.0; delta = 0.0
        hscale = np.abs(H[H < 0]).max()
        for it in range(50000):
            g = K(p) + H
            I = p > 0
            Gc = np.sum(g[I] ** 2)
            t = np.where(I, g + delta * (Gc / Gold) * t, 0.0); Gold = Gc
            r = K(t)
            tau = np.sum(g[I] * t[I]) / max(np.sum(r[I] * t[I]), 1e-300)
            pold = p.copy()
            p = np.where(I, p - tau * t, p); p[p < 0] = 0.0
            ol = (p == 0) & (g < 0)
            if ol.any():
                delta = 0.0; p[ol] -= tau * g[ol]
            else:
                delta = 1.0
            if it % 10 == 0:
                gg = K(p) + H; In = p > 0
                if In.any() and np.abs(gg[In]).max() < 1e-9 * hscale and (In.all() or gg[~In].min() > -1e-9 * hscale):
                    break
        return p

    def tangential(P, p, creep, warm):
        nx, ny = P["nx"], P["ny"]
        cv = convs(nx, ny); raw = cv["raw"]; cc = cv["c"]
        I = np.nonzero(p > 0); ci, cj = I; m = len(ci)
        di = ci[:, None] - ci[None, :] + nx - 1; dj = cj[:, None] - cj[None, :] + ny - 1
        U = np.block([[raw["xx"][di, dj], raw["xy"][di, dj]], [raw["xy"][di, dj], raw["yy"][di, dj]]])
        rho = float(np.mean(np.diag(U)))
        fac = cho_factor(U + rho * np.eye(2 * m))
        gb = mu * p[I]
        upok = ci + 1 < nx
        xi, eta, phi = creep
        w = np.concatenate([xi - phi * P["S"][I], eta + phi * P["X"][I]]) * DX
        def upstream(q):
            Qx = np.zeros((nx, ny)); Qy = np.zeros((nx, ny)); Qx[I] = q[:m]; Qy[I] = q[m:]
            ux = cc["xx"](Qx) + cc["xy"](Qy); uy = cc["xy"](Qx) + cc["yy"](Qy)
            iu = np.minimum(ci + 1, nx - 1)
            return np.concatenate([np.where(upok, ux[iu, cj], 0.0), np.where(upok, uy[iu, cj], 0.0)])
        def proj(v):
            n = np.hypot(v[:m], v[m:]); f = np.where(n > gb, gb / np.maximum(n, 1e-300), 1.0)
            return np.concatenate([v[:m] * f, v[m:] * f])
        if warm is not None and warm[0].shape == p.shape:
            Qx0, Qy0 = warm
            q = proj(np.concatenate([Qx0[I], Qy0[I]]))
        else:
            q = np.zeros(2 * m)
        z = q.copy(); lam = np.zeros(2 * m)
        for step in range(4000):
            b = w - upstream(q)
            for _ in range(4000):
                qn = cho_solve(fac, rho * (z - lam) - b)
                zn = proj(qn + lam)
                lam += qn - zn
                res = max(np.max(np.abs(qn - zn)), np.max(np.abs(zn - z)))
                z = zn
                if res <= 1e-10 * gb.max():
                    break
            ch = np.max(np.abs(z - q)); q = z.copy()
            if ch <= 1e-9 * gb.max() and step > 2:
                break
        Qx = np.zeros((nx, ny)); Qy = np.zeros((nx, ny)); Qx[I] = q[:m]; Qy[I] = q[m:]
        dA = DX * DS
        return q[:m].sum() * dA, q[m:].sum() * dA, (Qx, Qy)

    def normal_state(ys, zc):
        pats = []
        for pt in contact_points(ys, zc):
            P = build_patch(pt, ys, zc)
            if P is None or P["H"].min() >= 0:
                continue
            p = normal_cg(P)
            N = p.sum() * DX * DS
            if N > 0:
                pats.append((P, p, N))
        return pats

    def forces(pats, tang):
        Fy = Fz = 0.0
        for P, p, N in pats:
            kind = "flange" if abs(P["gamma"]) > math.radians(20) else "tread"
            Fx_, Fs_ = tang.get(kind, (0.0, 0.0))
            Fy += Fs_ * math.cos(P["gamma"]) - N * math.sin(P["gamma"])
            Fz += Fs_ * math.sin(P["gamma"]) + N * math.cos(P["gamma"])
        return Fy, Fz

    results = []
    ckpt = out_path + ".partial.json"
    done = json.load(open(ckpt)) if os.path.exists(ckpt) else {}
    for cs in case["load_cases"]:
        if cs["id"] in done:
            results.append(done[cs["id"]]); continue
        Q, Y = cs["vertical_load_N"], cs["lateral_force_N"]
        V = cs["speed_m_per_s"]; Om = cs["wheel_angular_velocity_rad_per_s"]; psi = cs["yaw_angle_rad"]
        ys = case["initial_guess"]["lateral_shift_mm"]
        y = np.linspace(yr[0], yr[-1], 20001)
        zc = -np.nanmin(np.where((y - ys > yw[0]) & (y - ys < yw[-1]), gfun(y, ys, 0.0), np.nan)) - 0.06
        tang = {}; warm = {}
        for outer in range(40):
            def res_(a, b):
                Fy, Fz = forces(normal_state(a, b), tang)
                return np.array([Fz - Q, Fy - Y])
            for _ in range(60):
                r = res_(ys, zc)
                if np.max(np.abs(r)) < 1e-6 * Q:
                    break
                J = np.column_stack([(res_(ys + 2e-4, zc) - r) / 2e-4, (res_(ys, zc + 2e-5) - r) / 2e-5])
                dlt = np.clip(np.linalg.solve(J, -r), [-0.05, -0.01], [0.05, 0.01])
                ys, zc = ys + dlt[0], zc + dlt[1]
            pats = normal_state(ys, zc)
            new = {}; info = {}
            for P, p, N in pats:
                kind = "flange" if abs(P["gamma"]) > math.radians(20) else "tread"
                k = Om * P["r"] * 1e-3 / V
                creep = (1 - k * math.cos(psi), -k * math.sin(psi) * math.cos(P["gamma"]), -(Om / V) * math.cos(psi) * math.sin(P["gamma"]))
                Fx_, Fs_, wq = tangential(P, p, creep, warm.get(kind))
                warm[kind] = wq
                new[kind] = (Fx_, Fs_)
                info[kind] = {"contact_angle_deg": math.degrees(P["gamma"]), "normal_force_N": N,
                              "area_mm2": float((p > 0).sum() * DX * DS * 1e6), "max_pressure_MPa": float(p.max() / 1e6),
                              "Fx_N": Fx_, "Fs_N": Fs_}
            change = max([abs(new[k][j] - tang.get(k, (0, 0))[j]) for k in new for j in (0, 1)] + [0.0])
            tang = new
            print(cs["id"], outer, round(ys, 5), round(zc, 6), {k: round(v["normal_force_N"]) for k, v in info.items()}, round(change, 1), flush=True)
            if change < 2e-4 * Q:
                break
        rec = {"id": cs["id"], "lateral_shift_mm": float(ys), "tread": info.get("tread"), "flange": info.get("flange"),
               "total_Fx_N": float(sum(v[0] for v in tang.values()))}
        results.append(rec); done[cs["id"]] = rec
        json.dump(done, open(ckpt, "w"))
    json.dump({"cases": results}, open(out_path, "w"), indent=2)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
