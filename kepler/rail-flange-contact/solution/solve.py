#!/usr/bin/env python3
"""Two-point (tread + flange) wheel-rail contact with Kalker's exact theory
and the wheel's force equilibrium.

1. Contact search on the vertical separation g(y) (resampled at 0.01 mm):
   local minima are contact points, contact angle = rail slope there, the
   potential patch is the basin between neighbouring maxima of g.
2. Per patch: rotate both profiles into the tangent-plane axes (s, n),
   normal separation n_w(s) - n_r(s) + x^2 cos(gamma) / (2 r_i); half-space
   boundary elements (closed-form Boussinesq/Cerruti rectangle integrals);
   displacement-controlled active-set normal solve.
3. Local creepages from the wheelset kinematics; steady rolling as the fixed
   point of Kalker's transient scheme, each step a convex problem solved by
   accelerated projected gradient (FISTA).
4. Equilibrium sum F_z = Q, sum F_y = Y in the unknown lateral shift and
   vertical position, solved by Broyden's method on the full model (creep
   forces included), starting from a Jacobian of the normal-force model.
"""
import json
import math
import os

import numpy as np

DATA = os.environ.get("RFC_DATA", "/app/data")
OUT = os.environ.get("RFC_OUT", "/app/output/flange_contact.json")
EL = 0.4e-3          # element size (m)
XHALF = 9.0e-3       # half length of the patch grid in the rolling direction
TOL_IN, TOL_STEP = float(os.environ.get("RFC_TIN", 1e-8)), float(os.environ.get("RFC_TSTEP", 1e-7))


def corner_sum(X, Y, hx, hy, kind):
    def xa(a, b):
        r = np.zeros_like(a); nz = a != 0
        r[nz] = a[nz] * np.arcsinh(b[nz] / np.abs(a[nz]))
        return r
    tot = np.zeros_like(X)
    for sx in (1, -1):
        for sy in (1, -1):
            a = X + sx * hx / 2; b = Y + sy * hy / 2; sg = sx * sy
            if kind == "n":
                tot += sg * (xa(a, b) + xa(b, a))
            elif kind == "xx":
                tot += sg * xa(b, a)
            elif kind == "yy":
                tot += sg * xa(a, b)
            else:
                tot -= sg * np.hypot(a, b)
    return tot


class Setup:
    def __init__(self):
        c = json.load(open(os.path.join(DATA, "case.json")))
        self.c = c
        r = np.loadtxt(os.path.join(DATA, "rail_profile.csv"), delimiter=",", skiprows=1)
        w = np.loadtxt(os.path.join(DATA, "wheel_profile.csv"), delimiter=",", skiprows=1)
        self.yr, self.zr, self.yw, self.zw = r[:, 0], r[:, 1], w[:, 0], w[:, 1]
        self.dzr = np.gradient(self.zr, self.yr)
        self.E = c["materials"]["youngs_modulus_Pa"]; self.nu = c["materials"]["poisson_ratio"]
        self.G = self.E / (2 * (1 + self.nu)); self.mu = c["coefficient_of_friction"]
        self.r0 = c["wheel"]["tape_circle_radius_mm"]
        self.zw_tape = float(np.interp(c["wheel"]["tape_circle_y_mm"], self.yw, self.zw))
        self.cache = {}
        self.kcache = {}

    def points(self, ys, zc):
        """Contact point of each wheel region = minimum of the vertical separation
        over the region (evaluated on a 0.01 mm resampling)."""
        split = self.c["frames"]["region_split_y_w_mm"]
        out = []
        for lo, hi in ((self.yw[0], split), (split, self.yw[-1])):
            a, b = max(self.yr[0], lo + ys), min(self.yr[-1], hi + ys)
            if b - a < 0.05:
                continue
            yf = np.arange(a, b, 0.01)
            g = np.interp(yf - ys, self.yw, self.zw) + zc - np.interp(yf, self.yr, self.zr)
            k = int(np.argmin(g))
            out.append(dict(y=yf[k], g=g[k], gamma=math.atan(np.interp(yf[k], self.yr, self.dzr)), ylo=a, yhi=b))
        return out

    def patch(self, pt, ys, zc):
        gm = pt["gamma"]; c, s = math.cos(gm), math.sin(gm)
        y0 = pt["y"]; z0 = np.interp(y0, self.yr, self.zr)
        def rot(y, z):
            return (y - y0) * c + (z - z0) * s, -(y - y0) * s + (z - z0) * c
        mr = (self.yr >= pt["ylo"]) & (self.yr <= pt["yhi"])
        ywt = self.yw + ys
        mw = (ywt >= pt["ylo"]) & (ywt <= pt["yhi"])
        if mr.sum() < 3 or mw.sum() < 3:
            return None
        sr, nr = rot(self.yr[mr], self.zr[mr]); sw, nw = rot(ywt[mw], self.zw[mw] + zc)
        lo, hi = max(sr.min(), sw.min()), min(sr.max(), sw.max())
        el = EL * 1e3
        sc = np.arange(math.ceil(lo / el) * el, hi, el) + el / 2
        sc = sc[sc < hi]
        h = (np.interp(sc, np.sort(sw), nw[np.argsort(sw)]) - np.interp(sc, np.sort(sr), nr[np.argsort(sr)])) * 1e-3
        k = h < 2.5e-4
        sc, h = sc[k] * 1e-3, h[k]
        if len(sc) == 0:
            return None
        ri = self.r0 - (np.interp(y0 - ys, self.yw, self.zw) - self.zw_tape)
        R = ri * 1e-3 / c
        nx = int(round(2 * XHALF / EL))
        xs = -XHALF + EL * (np.arange(nx) + 0.5)
        ny = int(round((sc.max() - sc.min()) / EL)) + 1
        ss = sc.min() + EL * np.arange(ny)
        XX, SS = np.meshgrid(xs, ss, indexing="ij")
        H = np.interp(SS, sc, h, left=1.0, right=1.0) + XX ** 2 / (2 * R)
        return dict(X=XX.ravel(), S=SS.ravel(), H=H.ravel(), nx=nx, ny=ny, r=ri, gamma=gm, pt=pt)

    def kernels(self, nx, ny):
        """Influence numbers for every grid offset (the grid is regular, so
        each matrix is a lookup into these tables)."""
        key = (nx, ny)
        if key not in self.kcache:
            i = np.arange(-(nx - 1), nx) * EL; j = np.arange(-(ny - 1), ny) * EL
            DXo, DSo = np.meshgrid(i, j, indexing="ij")
            c0 = 1 / (math.pi * self.G); nu = self.nu
            I1, Ixx, Iyy, Ixy = (corner_sum(DXo, DSo, EL, EL, k) for k in ("n", "xx", "yy", "xy"))
            self.kcache[key] = {"n": (1 - nu) * c0 * I1, "xx": c0 * ((1 - nu) * I1 + nu * Ixx),
                                "yy": c0 * ((1 - nu) * I1 + nu * Iyy), "xy": c0 * nu * Ixy}
        return self.kcache[key]

    def normal(self, P):
        H = P["H"]
        cand = np.nonzero(H < 2e-5)[0]
        if not (H[cand] < 0).any():
            return np.zeros_like(H)
        K = self.kernels(P["nx"], P["ny"])["n"]
        ix, iy = cand // P["ny"], cand % P["ny"]
        A = K[ix[:, None] - ix[None, :] + P["nx"] - 1, iy[:, None] - iy[None, :] + P["ny"] - 1]
        h = H[cand]
        act = np.nonzero(h < 0)[0]
        for _ in range(400):
            pa = np.linalg.solve(A[np.ix_(act, act)], -h[act])
            if (pa < 0).any():
                act = act[pa >= 0]; continue
            pc = np.zeros(len(cand)); pc[act] = pa
            add = np.setdiff1d(np.nonzero(h + A @ pc < -1e-13)[0], act)
            if len(add) == 0:
                break
            act = np.sort(np.concatenate([act, add]))
        p = np.zeros_like(H); p[cand] = pc
        return p

    def tangential(self, P, p, creep, q0=None):
        C = np.nonzero(p > 0)[0]
        m = len(C)
        X, S = P["X"], P["S"]
        ix = C // P["ny"]; iy = C % P["ny"]
        ok = ix + 1 < P["nx"]
        up = np.where(ok, (ix + 1) * P["ny"] + iy, C)
        K = self.kernels(P["nx"], P["ny"])
        def infl(rows):
            a = rows // P["ny"]; b = rows % P["ny"]
            di = a[:, None] - ix[None, :] + P["nx"] - 1; dj = b[:, None] - iy[None, :] + P["ny"] - 1
            return np.block([[K["xx"][di, dj], K["xy"][di, dj]], [K["xy"][di, dj], K["yy"][di, dj]]])
        U = infl(C); Vm = infl(up); Vm[np.concatenate([~ok, ~ok])] = 0.0
        L = np.linalg.eigvalsh(U)[-1]
        g = self.mu * p[C]
        xi, eta, phi = creep
        w = np.concatenate([xi - phi * S[C], eta + phi * X[C]]) * EL
        def proj(v):
            t = np.hypot(v[:m], v[m:]); f = np.where(t > g, g / np.maximum(t, 1e-300), 1.0)
            return np.concatenate([v[:m] * f, v[m:] * f])
        q = np.zeros(2 * m) if q0 is None or len(q0) != 2 * m else proj(q0)
        for step in range(4000):
            b = w - Vm @ q
            z = q.copy(); yv = q.copy(); t = 1.0
            for _ in range(5000):
                zn = proj(yv - (U @ yv + b) / L)
                if np.dot(yv - zn, zn - z) > 0:      # adaptive restart of the momentum
                    t = 1.0
                tn = 0.5 * (1 + math.sqrt(1 + 4 * t * t))
                yv = zn + (t - 1) / tn * (zn - z)
                done = np.max(np.abs(zn - z)) <= TOL_IN * g.max()
                z, t = zn, tn
                if done:
                    break
            ch = np.max(np.abs(z - q)); q = z
            if ch <= TOL_STEP * g.max() and step > 2:
                break
        dA = EL * EL
        return q[:m].sum() * dA, q[m:].sum() * dA, q

    def creep(self, cs, ri, gm):
        V = cs["speed_m_per_s"]; Om = cs["wheel_angular_velocity_rad_per_s"]; psi = cs["yaw_angle_rad"]
        k = Om * ri * 1e-3 / V
        return (1 - k * math.cos(psi), -k * math.sin(psi) * math.cos(gm), -(Om / V) * math.cos(psi) * math.sin(gm))

    def evaluate(self, ys, zc, cs, with_tang=True):
        res = []
        for pt in self.points(ys, zc):
            P = self.patch(pt, ys, zc)
            if P is None:
                continue
            p = self.normal(P)
            N = p.sum() * EL * EL
            if N <= 0:
                continue
            item = dict(P=P, p=p, N=N, Fx=0.0, Fs=0.0)
            if with_tang:
                cr = self.creep(cs, P["r"], P["gamma"])
                key = round(math.degrees(P["gamma"]))
                Fx, Fs, q = self.tangential(P, p, cr, self.cache.get(key))
                self.cache[key] = q
                item.update(Fx=Fx, Fs=Fs, creep=cr)
            res.append(item)
        Fy = sum(i["Fs"] * math.cos(i["P"]["gamma"]) - i["N"] * math.sin(i["P"]["gamma"]) for i in res)
        Fz = sum(i["Fs"] * math.sin(i["P"]["gamma"]) + i["N"] * math.cos(i["P"]["gamma"]) for i in res)
        return res, Fy, Fz


def solve_case(st, cs, ys, zc):
    Q, Y = cs["vertical_load_N"], cs["lateral_force_N"]
    target = np.array([Q, Y])
    def F(v, tang):
        res, Fy, Fz = st.evaluate(v[0], v[1], cs, tang)
        return np.array([Fz, Fy]) - target, res
    v = np.array([ys, zc])
    # stage 1: Newton on the normal-force-only model gets close and gives a Jacobian
    for _ in range(40):
        r, _ = F(v, False)
        J = np.column_stack([(F(v + [1e-4, 0], False)[0] - r) / 1e-4, (F(v + [0, 1e-5], False)[0] - r) / 1e-5])
        d = np.clip(np.linalg.solve(J, -r), [-0.05, -0.01], [0.05, 0.01])
        v = v + d
        if np.max(np.abs(r)) < 1e-6 * Q:
            break
    # stage 2: Broyden on the full model (creep forces included)
    r, res = F(v, True)
    for it in range(60):
        if np.max(np.abs(r)) < 1e-4 * Q:
            break
        d = np.clip(np.linalg.solve(J, -r), [-0.05, -0.01], [0.05, 0.01])
        v2 = v + d
        r2, res2 = F(v2, True)
        dv, dr = v2 - v, r2 - r
        J = J + np.outer(dr - J @ dv, dv) / (dv @ dv)
        v, r, res = v2, r2, res2
    print(cs["id"], "residual", r, flush=True)
    out = {"id": cs["id"], "lateral_shift_mm": float(v[0])}
    tread = [i for i in res if abs(i["P"]["gamma"]) < math.radians(20)]
    flange = [i for i in res if abs(i["P"]["gamma"]) >= math.radians(20)]
    for name, lst in (("tread", tread), ("flange", flange)):
        if not lst:
            out[name] = None; continue
        i = lst[0]
        out[name] = {"contact_angle_deg": math.degrees(i["P"]["gamma"]), "normal_force_N": float(i["N"]),
                     "Fx_N": float(i["Fx"]), "Fs_N": float(i["Fs"])}
    return out


def main():
    st = Setup()
    c = st.c
    results = []
    for cs in c["load_cases"]:
        ys = c["initial_guess"]["lateral_shift_mm"]
        pts = st.points(ys, 0.0)
        zc = -min(p["g"] for p in pts) - 0.06
        st.cache = {}
        results.append(solve_case(st, cs, ys, zc))
        print(json.dumps(results[-1]), flush=True)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump({"cases": results}, f, indent=2)


if __name__ == "__main__":
    main()
