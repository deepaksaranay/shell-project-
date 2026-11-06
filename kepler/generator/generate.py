#!/usr/bin/env python3
"""Generate the shipped rig records and the sealed truth for
bearing-impedance-identification.

NOT part of the task bundle. Run from the repository root:
    python3 kepler/generator/generate.py

Physical model (frequency domain, per excitation tone w):
  housing: (-w^2 Ms + Hs + Hb) xh - Hb xs = F
  shaft:   (Zr + Hb) xs - Hb xh = 0
  Hb = Kb - w^2 Mb + i w Cb   (bearing film, relative motion d = xh - xs)
  Hs = Ks - w^2 Ms_res + i w Cs (support flexures, absolute housing motion)
  Zr = Kr - w^2 mr + i w Cr   (rotor on its rolling-element support bearings)
The idle shaker's stinger loads the housing with F_idle = -Zst * xh_idle, and
both load cells record the force actually transmitted to the housing.
"""
import json
import os

import numpy as np

ROOT = os.path.join(os.path.dirname(__file__), "..", "bearing-impedance-identification")
DATA = os.path.join(ROOT, "environment", "data")
TESTS = os.path.join(ROOT, "tests")

rng = np.random.default_rng(20261106)

FS = 4096.0
T = 6.0
N = int(FS * T)
t = np.arange(N) / FS
TONES = np.arange(20.0, 321.0, 10.0)  # Hz, all bin-centred (df = 1/6 Hz)

# sensor sensitivities (stated in the rig sheet)
LC_SENS = 1.124e-3         # V/N   (5 mV/lbf)
ACC_SENS = 0.100 / 9.80665  # V/(m/s^2)  (100 mV/g)
PROX_SENS = 7.874e3        # V/m   (200 mV/mil), output = -sens * gap
PROX_ANG = np.deg2rad([45.0, 135.0])
GAP0 = np.array([1.27e-3, 1.24e-3])  # static gaps, m

# housing + support
MS = 18.5  # kg, floating housing incl. instrumentation (stated)
KS = np.array([[1.20e7, 0.04e7], [0.04e7, 1.45e7]])
CS = np.array([[850.0, 0.0], [0.0, 960.0]])
MS_RES = np.array([[0.35, 0.0], [0.0, 0.40]])  # flexure participating mass

# rotor on its own support bearings (not stated -- solver must not need it)
KR = np.diag([3.1e8, 3.4e8])
MR = np.diag([42.0, 42.0])
CR = np.diag([4.0e3, 4.5e3])

# idle shaker stinger
KST, CST = 3.0e6, 60.0

SPEEDS = [4100, 6250, 8300, 10450]
LOADS = [1000, 3000]


def bearing_truth(rpm, load):
    s = rpm / 10000.0
    w = load / 1000.0
    K = np.array([
        [4.0e7 + 1.5e7 * w + 2.0e7 * s, 2.6e7 * s + 0.3e7 * w],
        [-(1.9e7 * s + 0.6e7 * w), 6.0e7 + 2.4e7 * w + 1.5e7 * s],
    ])
    cxx = 1.1e5 / (0.55 + s) * (1 + 0.10 * w)
    cyy = 1.5e5 / (0.55 + s) * (1 + 0.12 * w)
    C = np.array([[cxx, 0.09 * cxx], [0.07 * cyy, cyy]])
    M = np.array([[2.2 + 0.6 * s, 0.0], [0.0, 3.0 + 0.5 * s]])
    return K, C, M


def schroeder_phases(n):
    k = np.arange(1, n + 1)
    return -np.pi * k * (k - 1) / n


def simulate(Kb, Cb, Mb, axis, rpm_actual, runout, rotating, amp=140.0):
    """Return dict of channel voltages for one shaker test."""
    nt = len(TONES)
    ph = schroeder_phases(nt) + rng.uniform(-0.2, 0.2, nt)
    Fdrive = amp * (1 + 0.15 * rng.standard_normal(nt)) * np.exp(1j * ph)
    F = np.zeros((nt, 2), complex)
    XH = np.zeros((nt, 2), complex)
    D = np.zeros((nt, 2), complex)
    idle = 1 - axis
    for k, f in enumerate(TONES):
        w = 2 * np.pi * f
        Hb = Kb - w**2 * Mb + 1j * w * Cb
        Hs = KS - w**2 * MS_RES + 1j * w * CS
        Zr = KR - w**2 * MR + 1j * w * CR
        A_h = -w**2 * MS * np.eye(2) + Hs + Hb
        # eliminate shaft: xs = (Zr+Hb)^-1 Hb xh
        T_s = np.linalg.solve(Zr + Hb, Hb)
        A = A_h - Hb @ T_s
        A[idle, idle] += KST + 1j * w * CST
        rhs = np.zeros(2, complex)
        rhs[axis] = Fdrive[k]
        xh = np.linalg.solve(A, rhs)
        xs = T_s @ xh
        Fk = rhs.copy()
        Fk[idle] = -(KST + 1j * w * CST) * xh[idle]
        F[k], XH[k], D[k] = Fk, xh, xh - xs

    def synth(C):
        out = np.zeros(N)
        for k, f in enumerate(TONES):
            out += np.real(C[k] * np.exp(2j * np.pi * f * t))
        return out

    fx, fy = synth(F[:, 0]), synth(F[:, 1])
    ax = synth(-(2 * np.pi * TONES) ** 2 * XH[:, 0])
    ay = synth(-(2 * np.pi * TONES) ** 2 * XH[:, 1])
    dx, dy = synth(D[:, 0]), synth(D[:, 1])

    Om = 2 * np.pi * rpm_actual / 60.0
    phi0 = rng.uniform(0, 2 * np.pi)
    phi = Om * t + phi0
    prox = []
    for j, th in enumerate(PROX_ANG):
        dp = dx * np.cos(th) + dy * np.sin(th)
        gap = GAP0[j] + dp  # housing moving toward +theta opens the gap
        if rotating:
            a1, b1, a2, b2 = runout
            gap = gap + a1 * np.cos(phi - th + b1) + a2 * np.cos(2 * (phi - th) + b2)
        gap = gap + 0.25e-6 * rng.standard_normal(N)
        prox.append(-PROX_SENS * gap)

    if rotating:
        frac = np.mod(phi, 2 * np.pi)
        kp = np.where(frac < 0.8, -9.0, -1.2) + 0.02 * rng.standard_normal(N)
    else:
        kp = -1.2 + 0.02 * rng.standard_normal(N)

    frms = np.std(fx if axis == 0 else fy)
    arms = np.std(ax if axis == 0 else ay)
    ch = {
        "Fx_V": LC_SENS * (fx + 0.006 * frms * rng.standard_normal(N)),
        "Fy_V": LC_SENS * (fy + 0.006 * frms * rng.standard_normal(N)),
        "ax_V": ACC_SENS * (ax + 0.006 * arms * rng.standard_normal(N)),
        "ay_V": ACC_SENS * (ay + 0.006 * arms * rng.standard_normal(N)),
        "prox45_V": prox[0],
        "prox135_V": prox[1],
        "keyphasor_V": kp,
    }
    return ch


def write_csv(path, ch):
    cols = ["Fx_V", "Fy_V", "ax_V", "ay_V", "prox45_V", "prox135_V", "keyphasor_V"]
    arr = np.column_stack([t] + [ch[c] for c in cols])
    header = "time_s," + ",".join(cols)
    np.savetxt(path, arr, delimiter=",", header=header, comments="",
               fmt=["%.6f"] + ["%.6e"] * len(cols))


def main():
    os.makedirs(DATA, exist_ok=True)
    zeros = np.zeros((2, 2))
    rows = []
    for axis, name in [(0, "x"), (1, "y")]:
        ch = simulate(zeros, zeros, zeros, axis, 0.0, None, rotating=False)
        write_csv(os.path.join(DATA, f"baseline_dry_{name}.csv"), ch)

    truth = {}
    cid = 0
    for load in LOADS:
        for rpm in SPEEDS:
            cid += 1
            tag = f"C{cid:02d}"
            Kb, Cb, Mb = bearing_truth(rpm, load)
            rpm_act = rpm * (1 + rng.uniform(-0.006, 0.006))
            runout = (rng.uniform(14e-6, 24e-6), rng.uniform(0, 2 * np.pi),
                      rng.uniform(4e-6, 8e-6), rng.uniform(0, 2 * np.pi))
            files = {}
            for axis, name in [(0, "x"), (1, "y")]:
                ch = simulate(Kb, Cb, Mb, axis, rpm_act, runout, rotating=True)
                fn = f"{tag}_{name}.csv"
                write_csv(os.path.join(DATA, fn), ch)
                files[name] = fn
            rows.append((tag, rpm, load, files["x"], files["y"]))
            truth[tag] = {"K": Kb.tolist(), "C": Cb.tolist(), "M": Mb.tolist(),
                          "rpm_actual": rpm_act}

    with open(os.path.join(DATA, "test_log.csv"), "w") as f:
        f.write("condition,nominal_speed_rpm,static_load_N,x_shaker_file,y_shaker_file\n")
        for r in rows:
            f.write(",".join(str(v) for v in r) + "\n")

    with open(os.path.join(TESTS, "expected.json"), "w") as f:
        json.dump(truth, f, indent=1)
    print("wrote", len(rows), "conditions")


if __name__ == "__main__":
    main()
