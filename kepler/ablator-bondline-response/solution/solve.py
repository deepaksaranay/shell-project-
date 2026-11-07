#!/usr/bin/env python3
"""Charring-ablator response of a stagnation-point heat shield.

Formulation (non-conservative CMA form of the energy equation):
    rho cp dT/dt|x = d/dx(k dT/dx) + mdot dh_g/dx + (h_g - hbar) drho/dt|x
    hbar = (rho_v h_v - rho_c h_c) / (rho_v - rho_c)
    d mdot/dx = drho/dt|x   (gas flows to the surface without storage)
written in the Landau coordinate eta = (x - s)/(L - s) that follows the
receding surface, discretised with finite differences on a graded grid.
The surface temperature has no thermal mass: it is the root of the surface
energy balance with a second-order one-sided conduction flux, and the
blowing-corrected film coefficient is solved at the same time. The
ablator/aluminium interface uses a two-material half-cell balance. The
semi-discrete system is integrated with an implicit Runge-Kutta (Radau)
method.
"""
import json
import math
import os

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

DATA = os.environ.get("ABL_DATA", "/app/data")
OUT = os.environ.get("ABL_OUT", "/app/output/heatshield_response.json")
SIG = 5.670374419e-8
NA, NS = 150, 8


def csv(name):
    return np.loadtxt(os.path.join(DATA, name), delimiter=",", skiprows=1)


class Solid:
    def __init__(self, tab, hf):
        self.T, self.cp = tab[:, 0], tab[:, 1]
        self.k = tab[:, 2]
        self.eps = tab[:, 3] if tab.shape[1] > 3 else None
        H = np.concatenate([[0.0], np.cumsum(0.5 * (self.cp[1:] + self.cp[:-1]) * np.diff(self.T))])
        self.Hn = H
        self.hf = hf
        self.off = self._int(np.array([298.0]))[0]

    def _int(self, t):
        T, cp, H = self.T, self.cp, self.Hn
        k = np.clip(np.searchsorted(T, t) - 1, 0, len(T) - 2)
        tc = np.clip(t, T[0], T[-1])
        d = tc - T[k]
        sl = (cp[k + 1] - cp[k]) / (T[k + 1] - T[k])
        base = H[k] + cp[k] * d + 0.5 * sl * d * d
        return base + np.where(t < T[0], cp[0] * (t - T[0]), 0.0) + np.where(t > T[-1], cp[-1] * (t - T[-1]), 0.0)

    def h(self, t):
        return self.hf + self._int(np.asarray(t, float)) - self.off

    def c(self, t):
        return np.interp(t, self.T, self.cp)

    def kk(self, t):
        return np.interp(t, self.T, self.k)


def main():
    case = json.load(open(os.path.join(DATA, "case.json")))
    ab = case["ablator"]; dec = ab["decomposition"]; G = dec["Gamma"]
    cA, cB = dec["components"]["A"], dec["components"]["B"]; rC = dec["components"]["C"]["rho_kg_per_m3"]
    rv = G * (cA["rho_v_kg_per_m3"] + cB["rho_v_kg_per_m3"]) + (1 - G) * rC
    rc = G * (cA["rho_r_kg_per_m3"] + cB["rho_r_kg_per_m3"]) + (1 - G) * rC
    V = Solid(csv(ab["virgin"]["properties_file"]), ab["virgin"]["formation_enthalpy_J_per_kg_at_298K"])
    Cc = Solid(csv(ab["char"]["properties_file"]), ab["char"]["formation_enthalpy_J_per_kg_at_298K"])
    gas = csv(ab["pyrolysis_gas"]["enthalpy_file"])
    hg = lambda t: np.interp(t, gas[:, 0], gas[:, 1])
    sb = case["substrate"]; Al = Solid(csv(sb["properties_file"]), 0.0); rAl = sb["density_kg_per_m3"]
    th = csv(case["surface"]["thermochemistry_file"])
    bT = np.unique(th[:, 0]); bG = np.unique(th[:, 1])
    BC = th[:, 2].reshape(len(bT), len(bG)); HW = th[:, 3].reshape(len(bT), len(bG))
    tr = csv("trajectory.csv")
    L = ab["thickness_mm"] * 1e-3; Ls = sb["thickness_mm"] * 1e-3
    Tsink = case["surface"]["radiation_sink_temperature_K"]; T0 = case["initial_temperature_K"]

    def interp2(tab, T, g):
        T = min(max(T, bT[0]), bT[-1]); g = min(max(g, bG[0]), bG[-1])
        i = min(np.searchsorted(bT, T, side="right") - 1, len(bT) - 2)
        j = min(np.searchsorted(bG, g, side="right") - 1, len(bG) - 2)
        a = (T - bT[i]) / (bT[i + 1] - bT[i]); b = (g - bG[j]) / (bG[j + 1] - bG[j])
        return (1 - a) * (1 - b) * tab[i, j] + a * (1 - b) * tab[i + 1, j] + (1 - a) * b * tab[i, j + 1] + a * b * tab[i + 1, j + 1]

    def yv_of(rho):
        return np.clip(rv * (rho - rc) / (rho * (rv - rc)), 0.0, 1.0)

    def react(T, ra, rb):
        out = []
        for P, r in ((cA, ra), (cB, rb)):
            f = np.clip((r - P["rho_r_kg_per_m3"]) / P["rho_v_kg_per_m3"], 0, None)
            out.append(-P["k_per_s"] * P["rho_v_kg_per_m3"] * f ** P["psi"] * np.exp(-P["E_over_R_K"] / T))
        return out

    # graded grid in eta (fine near the surface), nodes 0..NA
    eta = (np.exp(3.5 * np.arange(NA + 1) / NA) - 1) / (np.exp(3.5) - 1)
    h1, h2 = eta[1] - eta[0], eta[2] - eta[0]
    # one-sided second-order first-derivative weights at eta=0 on nonuniform spacing
    w0 = -(h1 + h2) / (h1 * h2); w1 = h2 / (h1 * (h2 - h1)); w2 = -h1 / (h2 * (h2 - h1))
    dzs = Ls / NS

    def surface(t, ra, rb, T1, T2, mg, s):
        """Solve SEB for T_w; return T_w, mdot_c, sdot."""
        C0 = np.interp(t, tr[:, 0], tr[:, 1]); hr = np.interp(t, tr[:, 0], tr[:, 2]); qr = np.interp(t, tr[:, 0], tr[:, 3])
        rho_w = G * (ra + rb) + (1 - G) * rC
        y = yv_of(np.array([rho_w]))[0]
        dx = L - s

        def film(Tw):
            def f(C):
                g = mg / C
                B = g + interp2(BC, Tw, g)
                return C - C0 * (math.log1p(B) / B if B > 1e-12 else 1.0)
            C = C0 if f(C0) <= 0 else brentq(f, 1e-14, C0, xtol=1e-16, rtol=1e-13)
            g = mg / C
            return C, interp2(BC, Tw, g) * C, interp2(HW, Tw, g)

        def resid(Tw):
            C, mc, hw = film(Tw)
            eps = y * np.interp(Tw, V.T, V.eps) + (1 - y) * np.interp(Tw, Cc.T, Cc.eps)
            hs = y * V.h(Tw) + (1 - y) * Cc.h(Tw)
            qin = (C * (hr - hw) + mg * hg(Tw) + mc * hs - (mg + mc) * hw
                   + eps * qr - eps * SIG * (Tw ** 4 - Tsink ** 4))
            k = y * V.kk(Tw) + (1 - y) * Cc.kk(Tw)
            dTdx = (w0 * Tw + w1 * T1 + w2 * T2) / dx
            return qin + k * dTdx, mc
        Tw = brentq(lambda T: resid(T)[0], 200.0, 5000.0, xtol=1e-9)
        mc = resid(Tw)[1]
        return Tw, mc, mc / rho_w

    n = NA + 1
    def unpack(y):
        return y[:NA], y[NA:NA + NS], y[NA + NS:NA + NS + n], y[NA + NS + n:NA + NS + 2 * n], y[-2], y[-1]

    def rhs(t, y, want=False):
        Ti, Ts, ra, rb, s, _ = unpack(y)        # Ti: ablator nodes 1..NA (NA = interface)
        dx = L - s
        rho = G * (ra + rb) + (1 - G) * rC
        # gas flux: mdot(eta) = integral from eta to 1 of -drho/dt dx (trapezoid)
        def gasflux(Tall):
            a, b = react(Tall, ra, rb)
            dr = G * (a + b)
            seg = -0.5 * (dr[1:] + dr[:-1]) * np.diff(eta) * dx
            mf = np.concatenate([np.cumsum(seg[::-1])[::-1], [0.0]])
            return mf, dr, a, b
        # T_w depends on mdot_g which depends weakly on T_w via node 0: two fixed-point passes
        Tw = Ti[0]
        for _ in range(3):
            Tall = np.concatenate([[Tw], Ti])
            mf, dr, a, b = gasflux(Tall)
            Tw, mc, sdot = surface(t, ra[0], rb[0], Ti[0], Ti[1], mf[0], s)
        Tall = np.concatenate([[Tw], Ti])
        mf, dr, a, b = gasflux(Tall)
        yv = yv_of(rho)
        k = yv * V.kk(Tall) + (1 - yv) * Cc.kk(Tall)
        cp = yv * V.c(Tall) + (1 - yv) * Cc.c(Tall)
        hbar = (rv * V.h(Tall) - rc * Cc.h(Tall)) / (rv - rc)
        hgas = hg(Tall)
        x = eta * dx
        vel = sdot * (1 - eta)                   # dx/dt of a grid node
        # interior nodes 1..NA-1: central differences on nonuniform grid
        i = np.arange(1, NA)
        hm = x[i] - x[i - 1]; hp = x[i + 1] - x[i]
        kf_p = 0.5 * (k[i] + k[i + 1]); kf_m = 0.5 * (k[i] + k[i - 1])
        cond = 2 * (kf_p * (Tall[i + 1] - Tall[i]) / hp - kf_m * (Tall[i] - Tall[i - 1]) / hm) / (hm + hp)
        dhg = (hgas[i + 1] - hgas[i]) / hp       # upwind: gas arrives from deeper nodes
        dTdx = (Tall[i + 1] - Tall[i]) / hp      # upwind for grid motion (material enters from deeper side)
        src = cond + mf[i] * dhg + (hgas[i] - hbar[i]) * dr[i]
        dT = src / (rho[i] * cp[i]) + vel[i] * dTdx
        # interface node NA: half ablator cell + half aluminium cell
        ha = x[NA] - x[NA - 1]
        Tsub = np.concatenate([[Tall[NA]], Ts])
        kal = Al.kk(Tsub)
        q_in = 0.5 * (k[NA] + k[NA - 1]) * (Tall[NA - 1] - Tall[NA]) / ha
        q_out = 0.5 * (kal[0] + kal[1]) * (Tsub[0] - Tsub[1]) / dzs
        capI = rho[NA] * cp[NA] * ha / 2 + rAl * Al.c(Tall[NA]) * dzs / 2
        srcI = (hgas[NA] - hbar[NA]) * dr[NA] * ha / 2   # mdot = 0 at the bondline
        dTI = (q_in - q_out + srcI) / capI
        # substrate nodes
        fl = 0.5 * (kal[1:] + kal[:-1]) * np.diff(Tsub) / dzs   # +x conduction flux reversed sign
        net = np.concatenate([fl[1:], [0.0]]) - fl
        wid = np.full(NS, dzs); wid[-1] = dzs / 2
        dTs = net / (rAl * Al.c(Ts) * wid)
        # densities in the moving frame
        gra = np.concatenate([(ra[1:] - ra[:-1]) / np.diff(x), [0.0]])
        grb = np.concatenate([(rb[1:] - rb[:-1]) / np.diff(x), [0.0]])
        dra = a + vel * gra; drb = b + vel * grb
        out = np.concatenate([dT, [dTI], dTs, dra, drb, [sdot, mf[0] + mc]])
        if want:
            return out, Tw
        return out

    y0 = np.concatenate([np.full(NA, T0), np.full(NS, T0), np.full(n, cA["rho_v_kg_per_m3"]), np.full(n, cB["rho_v_kg_per_m3"]), [0.0, 0.0]])
    times = case["output_times_s"]
    atol = np.concatenate([np.full(NA + NS, 1e-3), np.full(2 * n, 1e-5), [1e-10, 1e-7]])
    sol = solve_ivp(rhs, (0, times[-1]), y0, method="Radau", rtol=1e-6, atol=atol, dense_output=True, max_step=2.0)
    print(sol.message, sol.nfev, flush=True)

    def snap(t):
        y = sol.sol(t)
        _, mc_Tw = rhs(t, y, want=True)
        Ti, Ts, ra, rb, s, ml = unpack(y)
        rho = G * (ra + rb) + (1 - G) * rC
        tau = (rho - rc) / (rv - rc)
        depth = eta * (L - s)
        def front(lv):
            k = np.nonzero(tau >= lv)[0]
            if len(k) == 0:
                return depth[-1]
            k = k[0]
            if k == 0:
                return 0.0
            return depth[k - 1] + (lv - tau[k - 1]) / (tau[k] - tau[k - 1]) * (depth[k] - depth[k - 1])
        return {"t_s": float(t), "recession_mm": float(s * 1e3), "surface_temperature_K": float(mc_Tw),
                "bondline_temperature_K": float(Ti[-1]), "pyrolysis_front_depth_mm": float(front(0.98) * 1e3),
                "char_front_depth_mm": float(front(0.02) * 1e3), "mass_loss_kg_per_m2": float(ml)}

    tt = np.linspace(0, times[-1], 3001)
    Tb = sol.sol(tt)[NA - 1]
    j = int(np.argmax(Tb))
    fine = np.linspace(tt[max(j - 1, 0)], tt[min(j + 1, len(tt) - 1)], 401)
    Tf = sol.sol(fine)[NA - 1]
    out = {"snapshots": [snap(t) for t in times],
           "peak_bondline_temperature_K": float(Tf.max())}
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
