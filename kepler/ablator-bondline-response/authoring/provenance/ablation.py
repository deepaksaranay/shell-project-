#!/usr/bin/env python3
"""Reference generator: 1-D charring-ablator response (answer-key source).

Method of lines on a contracting, surface-attached vertex-centred
finite-volume grid (nodes keep their fractional position xi in the remaining
ablator thickness), conservative energy/mass balances with Reynolds-transport
terms for the moving control volumes, quasi-steady pyrolysis-gas flow,
surface energy balance with the blowing-corrected film coefficient solved
by a scalar root find, integrated with a stiff BDF solver.
"""
import json
import math
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

DATA = os.environ.get("ABL_DATA") or os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "environment", "data")
SIG = 5.670374419e-8


def load_csv(name):
    with open(os.path.join(DATA, name)) as f:
        head = f.readline().strip().split(",")
    return head, np.loadtxt(os.path.join(DATA, name), delimiter=",", skiprows=1)


class CpTable:
    """cp piecewise linear in T (held constant outside), h exact integral."""
    def __init__(self, T, cp, hf):
        self.T, self.cp = T, cp
        seg = 0.5 * (cp[1:] + cp[:-1]) * np.diff(T)
        self.H = np.concatenate([[0.0], np.cumsum(seg)])
        self.hf = hf
        self.h298 = self._H(np.array([298.0]))[0]

    def _H(self, t):
        T, cp, H = self.T, self.cp, self.H
        t = np.asarray(t, float)
        out = np.empty_like(t)
        lo = t <= T[0]; hi = t >= T[-1]; mid = ~(lo | hi)
        out[lo] = cp[0] * (t[lo] - T[0])
        out[hi] = H[-1] + cp[-1] * (t[hi] - T[-1])
        k = np.clip(np.searchsorted(T, t[mid]) - 1, 0, len(T) - 2)
        d = t[mid] - T[k]
        slope = (cp[k + 1] - cp[k]) / (T[k + 1] - T[k])
        out[mid] = H[k] + cp[k] * d + 0.5 * slope * d * d
        return out

    def h(self, t):
        return self.hf + self._H(t) - self.h298

    def c(self, t):
        return np.interp(t, self.T, self.cp)


class Model:
    def __init__(self, opts=None):
        self.opts = opts or {}
        case = json.load(open(os.path.join(DATA, "case.json")))
        self.case = case
        ab = case["ablator"]; dec = ab["decomposition"]; comp = dec["components"]
        self.G = dec["Gamma"]
        self.A, self.B, self.rhoC = comp["A"], comp["B"], comp["C"]["rho_kg_per_m3"]
        self.rho_v = self.G * (self.A["rho_v_kg_per_m3"] + self.B["rho_v_kg_per_m3"]) + (1 - self.G) * self.rhoC
        self.rho_c = self.G * (self.A["rho_r_kg_per_m3"] + self.B["rho_r_kg_per_m3"]) + (1 - self.G) * self.rhoC
        _, v = load_csv(ab["virgin"]["properties_file"])
        _, c = load_csv(ab["char"]["properties_file"])
        self.vT, self.vk, self.veps = v[:, 0], v[:, 2], v[:, 3]
        self.cT, self.ck, self.ceps = c[:, 0], c[:, 2], c[:, 3]
        self.hv = CpTable(v[:, 0], v[:, 1], ab["virgin"]["formation_enthalpy_J_per_kg_at_298K"])
        self.hc = CpTable(c[:, 0], c[:, 1], ab["char"]["formation_enthalpy_J_per_kg_at_298K"])
        _, g = load_csv(ab["pyrolysis_gas"]["enthalpy_file"])
        self.gT, self.gh = g[:, 0], g[:, 1]
        sb = case["substrate"]
        _, s = load_csv(sb["properties_file"])
        self.al = CpTable(s[:, 0], s[:, 1], 0.0)
        self.alT, self.alk = s[:, 0], s[:, 2]
        self.rho_al = sb["density_kg_per_m3"]
        _, b = load_csv(case["surface"]["thermochemistry_file"])
        self.bT = np.unique(b[:, 0]); self.bG = np.unique(b[:, 1])
        self.Bc = b[:, 2].reshape(len(self.bT), len(self.bG))
        self.Hw = b[:, 3].reshape(len(self.bT), len(self.bG))
        _, tr = load_csv("trajectory.csv")
        self.tr = tr
        self.L = ab["thickness_mm"] * 1e-3
        self.Ls = sb["thickness_mm"] * 1e-3
        self.T0 = case["initial_temperature_K"]
        self.Tsink = case["surface"]["radiation_sink_temperature_K"]

    def hg(self, T):
        return np.interp(T, self.gT, self.gh)

    def bilinear(self, tab, T, bg):
        T = min(max(T, self.bT[0]), self.bT[-1]); bg = min(max(bg, self.bG[0]), self.bG[-1])
        i = min(np.searchsorted(self.bT, T, side="right") - 1, len(self.bT) - 2)
        j = min(np.searchsorted(self.bG, bg, side="right") - 1, len(self.bG) - 2)
        a = (T - self.bT[i]) / (self.bT[i + 1] - self.bT[i]); c = (bg - self.bG[j]) / (self.bG[j + 1] - self.bG[j])
        return ((1 - a) * (1 - c) * tab[i, j] + a * (1 - c) * tab[i + 1, j]
                + (1 - a) * c * tab[i, j + 1] + a * c * tab[i + 1, j + 1])

    def env(self, t):
        tr = self.tr
        return (np.interp(t, tr[:, 0], tr[:, 1]), np.interp(t, tr[:, 0], tr[:, 2]), np.interp(t, tr[:, 0], tr[:, 3]))

    def yv(self, rho):
        return np.clip(self.rho_v * (rho - self.rho_c) / (rho * (self.rho_v - self.rho_c)), 0.0, 1.0)

    def rates(self, T, rA, rB):
        out = []
        for P, r in ((self.A, rA), (self.B, rB)):
            frac = np.clip((r - P["rho_r_kg_per_m3"]) / P["rho_v_kg_per_m3"], 0.0, None)
            out.append(-P["k_per_s"] * P["rho_v_kg_per_m3"] * frac ** P["psi"] * np.exp(-P["E_over_R_K"] / T))
        return out

    def surface(self, t, Tw, mg, yv_w, rho_w):
        """Returns (C_H blowing corrected, mdot_c, h_w, q_net_without_hs_term)."""
        C0, hr, qr = self.env(t)
        def bc_of(C):
            bg = mg / C
            return bg, self.bilinear(self.Bc, Tw, 0.0 if self.opts.get("bc_at_bg0") else bg)
        def f(C):
            bg, bc = bc_of(C)
            B = bg + bc
            om = 1.0 if (B < 1e-12 or self.opts.get("no_blowing")) else math.log1p(B) / B   # 2*lambda = 1
            return C - C0 * om
        if f(C0) <= 0:
            C = C0
        else:
            C = brentq(f, 1e-14, C0, xtol=1e-16, rtol=1e-13)
        bg, bc = bc_of(C)
        hw = self.bilinear(self.Hw, Tw, bg)
        mc = 0.0 if self.opts.get("no_recession") else bc * C
        eps = yv_w * np.interp(Tw, self.vT, self.veps) + (1 - yv_w) * np.interp(Tw, self.cT, self.ceps)
        q = C * (hr - hw) - (mg + mc) * hw + eps * qr - eps * SIG * (Tw ** 4 - self.Tsink ** 4)
        return C, mc, hw, q


def run(Na=160, Ns=8, beta=4.0, rtol=1e-7, atol_T=1e-4, verbose=True, opts=None, t_end=None):
    m = Model(opts)
    gas_conv = 0.0 if (opts or {}).get("no_gas_convection") else 1.0
    xi = (np.exp(beta * np.arange(Na + 1) / Na) - 1) / (np.exp(beta) - 1)
    fxi = 0.5 * (xi[1:] + xi[:-1])                      # interior faces in xi
    wxi = np.diff(np.concatenate([[0.0], fxi, [1.0]]))   # CV widths in xi (ablator part)
    dzs = m.Ls / Ns
    zs = m.L + dzs * np.arange(1, Ns + 1)                # substrate nodes (beyond interface)
    ws_sub = np.full(Ns, dzs); ws_sub[-1] = dzs / 2       # substrate CV widths
    half_sub = dzs / 2                                    # substrate half-cell on interface node
    nT = Na + 1 + Ns
    iT = slice(0, nT); iA = slice(nT, nT + Na + 1); iB = slice(nT + Na + 1, nT + 2 * (Na + 1))
    iS = nT + 2 * (Na + 1); iM = iS + 1
    A, B = m.A, m.B

    def rhs(t, y):
        T = y[iT]; rA = y[iA]; rB = y[iB]; s = y[iS]
        Ta = T[: Na + 1]
        rho = m.G * (rA + rB) + (1 - m.G) * m.rhoC
        yv = m.yv(rho)
        dA, dB = m.rates(Ta, rA, rB)
        dr = m.G * (dA + dB)                              # d rho / dt (fixed frame), <= 0
        Lr = m.L - s
        D = Lr * wxi                                      # ablator CV widths
        prod = -dr * D                                    # gas production per CV
        mf = np.concatenate([np.cumsum(prod[::-1])[::-1], [0.0]])  # mf[j] = flux at face j-1/2 (j=0 surface)
        mg = mf[0]
        C, mc, hw, qs = m.surface(t, Ta[0], mg, yv[0], rho[0])
        sdot = mc / rho[0]
        # face velocities (ablator faces), surface face velocity = sdot
        wf = sdot * (1 - fxi)                             # interior faces
        # conductivities
        ka = yv * np.interp(Ta, m.vT, m.vk) + (1 - yv) * np.interp(Ta, m.cT, m.ck)
        xa = s + xi * Lr
        Gf = np.zeros(Na)                                 # +x energy flux at interior ablator faces
        kf = 0.5 * (ka[1:] + ka[:-1])
        Gf = -kf * (Ta[1:] - Ta[:-1]) / np.diff(xa) - gas_conv * mf[1:Na + 1] * m.hg(Ta[1:])
        # substrate faces: interface node <-> sub node 0, sub i <-> i+1
        Tsub = np.concatenate([[Ta[-1]], T[Na + 1:]])
        ksub = np.interp(Tsub, m.alT, m.alk)
        Gs = -0.5 * (ksub[1:] + ksub[:-1]) * np.diff(Tsub) / dzs
        # energies
        hv = m.hv.h(Ta); hc = m.hc.h(Ta)
        h = yv * hv + (1 - yv) * hc
        e = rho * h
        cp = yv * m.hv.c(Ta) + (1 - yv) * m.hc.c(Ta)
        # ---- ablator mass balances (per CV)
        def mass_rate(ri, di):
            inn = np.concatenate([ri[1:] * wf, [0.0]])     # entering through right face
            out = np.concatenate([[ri[0] * sdot], ri[1:] * wf])   # leaving through left face (upwind: own CV)
            return di * D + inn - out
        dMA = mass_rate(rA, dA); dMB = mass_rate(rB, dB)
        dD = -sdot * wxi
        drA = (dMA - rA * dD) / D; drB = (dMB - rB * dD) / D
        drho = m.G * (drA + drB)
        # ---- energy balance per CV
        # without gas convection the gas still leaves at the surface with h_g(T_w)
        Gl = np.concatenate([[qs + mc * h[0] + (1 - gas_conv) * mg * m.hg(Ta[0])], Gf])        # left-face flux (surface: SEB incl. mc*h_s)
        Gr = np.concatenate([Gf, [Gs[0]]])
        conv_in = np.concatenate([e[1:] * wf, [0.0]])
        conv_out = np.concatenate([[e[0] * sdot], e[1:] * wf])
        dE = Gl - Gr + conv_in - conv_out
        # shortcut variant: gas leaves where it is produced (no heat pick-up on the way out)
        dE = dE + (1 - gas_conv) * dr * D * m.hg(Ta)
        Mab = rho * D
        dMab = m.G * (dMA + dMB) + (1 - m.G) * mass_rate(np.full(Na + 1, m.rhoC), 0.0)  # fibre is carried by the grid too
        dyv_drho = m.rho_v * m.rho_c / (rho ** 2 * (m.rho_v - m.rho_c)) * ((yv > 0) & (yv < 1))
        num = dE - h * dMab - Mab * (hv - hc) * dyv_drho * drho
        cap = Mab * cp
        cap[-1] += m.rho_al * half_sub * m.al.c(np.array([Ta[-1]]))[0]
        dTa = num / cap
        # substrate nodes
        Gsl = Gs
        Gsr = np.concatenate([Gs[1:], [0.0]])
        dTs = (Gsl - Gsr) / (m.rho_al * ws_sub * m.al.c(T[Na + 1:]))
        return np.concatenate([dTa, dTs, drA, drB, [sdot, mg + mc]])

    y0 = np.concatenate([np.full(nT, m.T0), np.full(Na + 1, A["rho_v_kg_per_m3"]), np.full(Na + 1, B["rho_v_kg_per_m3"]), [0.0, 0.0]])
    times = m.case["output_times_s"]
    atol = np.concatenate([np.full(nT, atol_T), np.full(2 * (Na + 1), 1e-6), [1e-11, 1e-8]])
    nb = Na + 1
    sol = solve_ivp(rhs, (0, times[-1]), y0, method="BDF", rtol=rtol, atol=atol, dense_output=True, max_step=1.0)
    if verbose:
        print("nfev", sol.nfev, "njev", sol.njev, sol.message, file=sys.stderr)

    def snapshot(t):
        y = sol.sol(t)
        T = y[iT]; rho = m.G * (y[iA] + y[iB]) + (1 - m.G) * m.rhoC; s = y[iS]
        tau = (rho - m.rho_c) / (m.rho_v - m.rho_c)
        depth = xi * (m.L - s)
        def front(level):
            k = np.nonzero(tau >= level)[0]
            if len(k) == 0:
                return float(depth[-1])
            k = k[0]
            if k == 0:
                return 0.0
            return float(depth[k - 1] + (level - tau[k - 1]) / (tau[k] - tau[k - 1]) * (depth[k] - depth[k - 1]))
        return {"t_s": t, "recession_mm": s * 1e3, "surface_temperature_K": float(T[0]),
                "bondline_temperature_K": float(T[Na]), "pyrolysis_front_depth_mm": front(0.98) * 1e3,
                "char_front_depth_mm": front(0.02) * 1e3, "mass_loss_kg_per_m2": float(y[iM])}
    snaps = [snapshot(t) for t in times]
    tt = np.linspace(0, times[-1], 4001)
    Tb = sol.sol(tt)[Na]
    k = int(np.argmax(Tb))
    from scipy.optimize import minimize_scalar
    lo, hi = tt[max(k - 1, 0)], tt[min(k + 1, len(tt) - 1)]
    r = minimize_scalar(lambda t: -sol.sol(t)[Na], bounds=(lo, hi), method="bounded", options={"xatol": 1e-4})
    peak = {"peak_bondline_temperature_K": float(-r.fun), "time_of_peak_s": float(r.x)}
    Tw = sol.sol(tt)[0]
    peak["peak_surface_temperature_K"] = float(Tw.max())
    return {"snapshots": snaps, "peak": peak, "grid": {"Na": Na, "Ns": Ns, "beta": beta, "rtol": rtol}}


if __name__ == "__main__":
    Na = int(sys.argv[1]) if len(sys.argv) > 1 else 160
    out = run(Na=Na)
    print(json.dumps(out, indent=1))
