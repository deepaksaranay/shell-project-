#!/usr/bin/env python3
"""Answer-key generator: two-point (tread + flange) wheel-rail contact with
Kalker's exact theory on each patch and the wheel's vertical/lateral force
equilibrium (including the creep forces) solved for the wheel position.

Model (as stated in environment/data/case.json):
  * vertical separation g(y) = z_w(y - y_shift) + z_c - z_r(y); each local
    minimum of g is a contact point, its basin (between the neighbouring local
    maxima) the potential patch region; contact angle = rail slope there;
  * each patch is planar in the plane tangent to the rail at the contact
    point; normal separation from the profiles rotated by the contact angle,
    plus x^2/(2 r_i / cos(gamma_i)) for the wheel's longitudinal curvature;
  * normal problem displacement-controlled; tangential problem steady
    rolling (transient-scheme fixed point, convex step problems);
  * local creepages from the wheelset kinematics; forces transformed to the
    track frame; sum F_z = Q, sum F_y = Y.
"""
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from contact import Grid, _rect_integrals, solve_tangential_rolling  # noqa: E402

DATA = os.environ.get("RFC_DATA") or os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "environment", "data")


def load():
    case = json.load(open(os.path.join(DATA, "case.json")))
    r = np.loadtxt(os.path.join(DATA, "rail_profile.csv"), delimiter=",", skiprows=1)
    w = np.loadtxt(os.path.join(DATA, "wheel_profile.csv"), delimiter=",", skiprows=1)
    return case, r, w


class Model:
    def __init__(self, el_mm=0.4, xwin_mm=(9.0, 9.0), opts=None):
        self.opts = opts or {}
        self.case, r, w = load()
        self.yr, self.zr = r[:, 0], r[:, 1]
        self.yw, self.zw = w[:, 0], w[:, 1]
        c = self.case
        self.E = c["materials"]["youngs_modulus_Pa"]; self.nu = c["materials"]["poisson_ratio"]
        self.mu = c["coefficient_of_friction"]
        self.r0 = c["wheel"]["tape_circle_radius_mm"]
        self.zw0 = float(np.interp(0.0, self.yw, self.zw))
        self.el = el_mm
        self.xwin = xwin_mm
        self.slope_r = np.gradient(self.zr, self.yr)
        self.slope_w = np.gradient(self.zw, self.yw)

    # ---------------- geometry
    def vgap(self, ys, zc):
        y = self.yr
        ok = (y - ys >= self.yw[0]) & (y - ys <= self.yw[-1])
        g = np.full_like(y, np.inf)
        g[ok] = np.interp(y[ok] - ys, self.yw, self.zw) + zc - self.zr[ok]
        return g

    def contact_points(self, ys, zc):
        g = self.vgap(ys, zc)
        fin = np.isfinite(g)
        idx = np.nonzero(fin[1:-1] & fin[:-2] & fin[2:] & (g[1:-1] < g[:-2]) & (g[1:-1] <= g[2:]))[0] + 1
        maxima = np.nonzero(fin[1:-1] & fin[:-2] & fin[2:] & (g[1:-1] > g[:-2]) & (g[1:-1] >= g[2:]))[0] + 1
        lo_edge = np.nonzero(fin)[0][0]; hi_edge = np.nonzero(fin)[0][-1]
        pts = []
        for k in idx:
            # parabolic refinement of the minimum
            a, b, c3 = g[k - 1], g[k], g[k + 1]
            den = a - 2 * b + c3
            d = 0.5 * (a - c3) / den if den > 0 else 0.0
            yc = self.yr[k] + d * (self.yr[k + 1] - self.yr[k])
            left = maxima[maxima < k]; right = maxima[maxima > k]
            ylo = self.yr[left[-1]] if len(left) else self.yr[lo_edge]
            yhi = self.yr[right[0]] if len(right) else self.yr[hi_edge]
            gam = math.atan(float(np.interp(yc, self.yr, self.slope_r)))
            pts.append({"y": yc, "gmin": b - 0.25 * (a - c3) * d, "gamma": gam, "ylo": ylo, "yhi": yhi})
        return pts

    def patch_gap(self, pt, ys, zc):
        """Normal separation on the patch grid (m), plus the grid and rolling radius."""
        gam = pt["gamma"]; cg, sg = math.cos(gam), math.sin(gam)
        if self.opts.get("horizontal_patches"):
            cg, sg = 1.0, 0.0
        y0 = pt["y"]; z0 = float(np.interp(y0, self.yr, self.zr))
        # rotate rail polyline and wheel polyline (track frame) into (s, n)
        def rot(y, z):
            return (y - y0) * cg + (z - z0) * sg, -(y - y0) * sg + (z - z0) * cg
        sr, nr = rot(self.yr, self.zr)
        ywt = self.yw + ys; zwt = self.zw + zc
        sw, nw = rot(ywt, zwt)
        # basin limits in s
        s_lo = ((pt["ylo"] - y0) * cg + (np.interp(pt["ylo"], self.yr, self.zr) - z0) * sg)
        s_hi = ((pt["yhi"] - y0) * cg + (np.interp(pt["yhi"], self.yr, self.zr) - z0) * sg)
        el = self.el
        # s grid covering the basin part where separation is small
        s_all = np.arange(math.floor(s_lo / el) * el, s_hi, el)
        sc = s_all + el / 2
        sc = sc[(sc > s_lo) & (sc < s_hi)]
        def interp_mono(s, n, q):
            o = np.argsort(s)
            return np.interp(q, s[o], n[o], left=np.nan, right=np.nan)
        # restrict to the local branch: rail/wheel points within the basin y-range
        mr = (self.yr >= pt["ylo"] - 1) & (self.yr <= pt["yhi"] + 1)
        mw = (ywt >= pt["ylo"] - 1) & (ywt <= pt["yhi"] + 1)
        hs = (interp_mono(sw[mw], nw[mw], sc) - interp_mono(sr[mr], nr[mr], sc)) * 1e-3
        keep = np.isfinite(hs) & (hs < 0.25e-3)
        sc, hs = sc[keep], hs[keep]
        r_i = self.r0 - (float(np.interp(y0 - ys, self.yw, self.zw)) - self.zw0)
        if self.opts.get("same_radius"):
            r_i = self.r0
        Rx = r_i * 1e-3 / cg
        if len(sc) == 0:
            return None
        x0, x1 = -self.xwin[0], self.xwin[1]
        nx = int(round((x1 - x0) / el))
        smin, smax = sc.min() - el / 2, sc.max() + el / 2
        ny = int(round((smax - smin) / el))
        gr = Grid(x0 * 1e-3, x1 * 1e-3, nx, smin * 1e-3, smax * 1e-3, ny)
        hgrid = np.interp(gr.Y, sc * 1e-3, hs, left=1.0, right=1.0) + gr.X ** 2 / (2 * Rx)
        return gr, hgrid, r_i

    # ---------------- normal problem (displacement controlled)
    def normal(self, gr, h):
        G = self.E / (2 * (1 + self.nu))
        C = np.nonzero(h < 0)[0]
        if len(C) == 0:
            return np.zeros_like(h)
        I1, _, _, _ = _rect_integrals(gr.X[C], gr.Y[C], gr.X[C], gr.Y[C], gr.dx, gr.dy)
        cand = np.nonzero(h < 2e-5)[0]
        I1c, _, _, _ = _rect_integrals(gr.X[cand], gr.Y[cand], gr.X[cand], gr.Y[cand], gr.dx, gr.dy)
        A = (1 - self.nu) / (math.pi * G) * I1c
        hc = h[cand]
        act = np.nonzero(hc < 0)[0]
        for _ in range(500):
            pa = np.linalg.solve(A[np.ix_(act, act)], -hc[act])
            if (pa < 0).any():
                act = act[pa >= 0]
                continue
            pc = np.zeros(len(cand)); pc[act] = pa
            gap = hc + A @ pc
            add = np.setdiff1d(np.nonzero(gap < -1e-13)[0], act)
            if len(add) == 0:
                break
            act = np.sort(np.concatenate([act, add]))
        p = np.zeros_like(h); p[cand] = pc
        return p

    # ---------------- creepages
    def creepages(self, cs, r_i, gam):
        V = cs["speed_m_per_s"]; Om = cs["wheel_angular_velocity_rad_per_s"]; psi = cs["yaw_angle_rad"]
        rr = r_i * 1e-3
        xi = 1 - Om * rr * math.cos(psi) / V
        eta = -(Om * rr / V) * math.sin(psi) * math.cos(gam)
        phi = -(Om / V) * math.cos(psi) * math.sin(gam)
        return xi, eta, phi

    # ---------------- state at a wheel position
    def state(self, ys, zc, cs, tang=None, need_tang=True):
        pts = self.contact_points(ys, zc)
        patches = []
        for pt in pts:
            pg = self.patch_gap(pt, ys, zc)
            if pg is None:
                continue
            gr, h, r_i = pg
            if h.min() >= 0:
                continue
            p = self.normal(gr, h)
            N = p.sum() * gr.dA
            if N <= 0:
                continue
            patches.append({"pt": pt, "grid": gr, "p": p, "N": N, "r": r_i})
        for pa in patches:
            gam = pa["pt"]["gamma"]
            pa["creep"] = self.creepages(cs, pa["r"], gam)
            if need_tang:
                res = solve_tangential_rolling(pa["grid"], pa["p"], self.E, self.nu, self.mu, *pa["creep"])
                pa["Fx"], pa["Fs"] = res["Fx"], res["Fy"]
                pa["steps"] = res["steps"]; pa["conv"] = res["converged"]
        return patches

    @staticmethod
    def totals(patches, tang_scale=None):
        Fy = Fz = Fx = 0.0
        for i, pa in enumerate(patches):
            g = pa["pt"]["gamma"]
            Fs, Fxx = pa.get("Fs", 0.0), pa.get("Fx", 0.0)
            if tang_scale is not None:
                Fxx, Fs = tang_scale[i]
            Fy += Fs * math.cos(g) - pa["N"] * math.sin(g)
            Fz += Fs * math.sin(g) + pa["N"] * math.cos(g)
            Fx += Fxx
        return Fx, Fy, Fz


def classify(patches):
    """tread = contact point with the smallest |gamma|, flange = the other (if any)."""
    out = {"tread": None, "flange": None}
    if not patches:
        return out
    srt = sorted(patches, key=lambda p: abs(p["pt"]["gamma"]))
    out["tread"] = srt[0] if abs(srt[0]["pt"]["gamma"]) < math.radians(20) else None
    fl = [p for p in patches if abs(p["pt"]["gamma"]) >= math.radians(20)]
    out["flange"] = max(fl, key=lambda p: abs(p["pt"]["gamma"])) if fl else None
    if out["tread"] is None and len(srt) == 1:
        out["flange"] = srt[0]
    return out


def solve_case(m, cs, ys0, zc0, verbose=True):
    Q, Y = cs["vertical_load_N"], cs["lateral_force_N"]
    ys, zc = ys0, zc0
    ratios = None      # tangential force per unit normal force from the last full solve, per patch class
    for outer in range(30):
        # inner: Newton on (ys, zc) with tangential forces modelled as ratio * N (frozen ratios)
        def resid(ys_, zc_):
            pats = m.state(ys_, zc_, cs, need_tang=False)
            cl = classify(pats)
            sc = []
            for pa in pats:
                key = "tread" if pa is cl["tread"] else "flange"
                rx, rs = ratios[key] if (ratios and key in ratios) else (0.0, 0.0)
                if m.opts.get("no_tangential_equilibrium"):
                    rx, rs = 0.0, 0.0
                sc.append((rx * pa["N"], rs * pa["N"]))
            Fx, Fy, Fz = Model.totals(pats, sc)
            return np.array([Fz - Q, Fy - Y]), pats
        for it in range(40):
            R, pats = resid(ys, zc)
            if np.max(np.abs(R)) < 1e-6 * Q:
                break
            J = np.zeros((2, 2)); hy, hz = 1e-4, 1e-5
            J[:, 0] = (resid(ys + hy, zc)[0] - R) / hy
            J[:, 1] = (resid(ys, zc + hz)[0] - R) / hz
            d = np.linalg.solve(J, -R)
            lim = np.array([0.05, 0.01])                 # mm per Newton step
            d = np.clip(d, -lim, lim)
            ys, zc = ys + d[0], zc + d[1]
        # full solve with tangential problems at this position
        pats = m.state(ys, zc, cs, need_tang=True)
        cl = classify(pats)
        new = {}
        for key in ("tread", "flange"):
            pa = cl[key]
            if pa is not None:
                new[key] = (pa["Fx"] / pa["N"], pa["Fs"] / pa["N"])
        Fx, Fy, Fz = Model.totals(pats)
        err = max(abs(Fz - Q), abs(Fy - Y))
        if m.opts.get("no_tangential_equilibrium"):
            err = 0.0
        if verbose:
            print(f"  outer {outer}: ys={ys:.5f} zc={zc:.6f} Fz-Q={Fz-Q:9.2f} Fy-Y={Fy-Y:9.2f} "
                  + " ".join(f"{k}:N={cl[k]['N']:.0f}" for k in cl if cl[k] is not None), flush=True)
        if err < 2e-4 * Q:
            break
        ratios = new
    out = {"id": cs["id"], "lateral_shift_mm": ys, "vertical_position_mm": zc, "Fx_N": Fx, "Fy_N": Fy, "Fz_N": Fz}
    for key in ("tread", "flange"):
        pa = cl[key]
        if pa is None:
            out[key] = None
            continue
        g = pa["pt"]["gamma"]
        out[key] = {"contact_angle_deg": math.degrees(g), "contact_y_mm": pa["pt"]["y"], "rolling_radius_mm": pa["r"],
                    "normal_force_N": pa["N"], "area_mm2": float((pa["p"] > 0).sum() * pa["grid"].dA * 1e6),
                    "max_pressure_MPa": float(pa["p"].max() / 1e6), "Fx_N": pa["Fx"], "Fs_N": pa["Fs"],
                    "creepages": pa["creep"], "steps": pa["steps"], "converged": bool(pa["conv"])}
    return out


if __name__ == "__main__":
    el = float(sys.argv[1]) if len(sys.argv) > 1 else 0.4
    m = Model(el_mm=el)
    res = []
    ys, zc = m.case["initial_guess"]["lateral_shift_mm"], m.case["initial_guess"]["vertical_position_mm"]
    only = sys.argv[2] if len(sys.argv) > 2 else None
    for cs in m.case["load_cases"]:
        if only and cs["id"] not in only:
            continue
        print("case", cs["id"], flush=True)
        g = m.vgap(ys, 0.0)
        zc = -np.min(g[np.isfinite(g)]) - 0.06            # start slightly interpenetrated
        r = solve_case(m, cs, ys, zc)
        res.append(r)
        print(json.dumps({k: v for k, v in r.items() if k not in ("tread", "flange")}), flush=True)
    json.dump({"element_mm": el, "cases": res}, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), f"truth_el{el}.json"), "w"), indent=1)
