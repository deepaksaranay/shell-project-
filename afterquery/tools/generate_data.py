"""Deterministic generator for the bearing-impedance-identification task data.

Not part of the agent-visible bundle. Writes CSV data into environment/data and
the ground truth into tests/expected.json.
"""
import json
import pathlib
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1] / "bearing-impedance-identification"
DATA = ROOT / "environment" / "data"
G0 = 9.80665

fs = 10240.0
T = 2.0
N = int(fs * T)
t = np.arange(N) / fs
rpm = 6000.0
f_sync = rpm / 60.0
m_stator = 24.6  # kg, bearing test (stator + bearing insert)
m_tare = 17.9  # kg, tare test (bearing insert removed)

freqs = np.array([f for f in range(20, 341, 20) if f % 100 != 0], dtype=float)

# Bearing truth (rotordynamic coefficients, reaction force on stator = -(K x + C xdot + M xddot))
K = np.array([[3.20e7, 1.10e7], [-1.40e7, 4.50e7]])
C = np.array([[5.50e4, -8.0e3], [-5.0e3, 7.20e4]])
M = np.array([[6.0, 0.4], [-0.3, 7.5]])
# Support structure (pitch stabilizers) in parallel with bearing
Ks = np.array([[1.8e6, 0.0], [0.0, 2.3e6]])
Cs = np.array([[350.0, 0.0], [0.0, 420.0]])


def H(Kx, Cx, Mx, w):
    return Kx + 1j * w * Cx - w**2 * Mx


rng = np.random.default_rng(20261106)


def make_test(direction, H_fun, runout, label, mass):
    """Simulate one shaker test. direction 0 -> x shaker, 1 -> y shaker."""
    nf = len(freqs)
    amp = 180.0 + 60.0 * rng.random(nf)
    phase = 2 * np.pi * rng.random(nf)
    # small leakage into the orthogonal direction from shaker misalignment
    leak = 0.06 * amp * np.exp(1j * 2 * np.pi * rng.random(nf))
    F_t = np.zeros((2, N))
    X_t = np.zeros((2, N))
    A_t = np.zeros((2, N))
    for k, f in enumerate(freqs):
        w = 2 * np.pi * f
        Fv = np.zeros(2, complex)
        Fv[direction] = amp[k] * np.exp(1j * phase[k])
        Fv[1 - direction] = leak[k]
        # m a = F - H x,  a = -w^2 x  ->  (H - m w^2) x = F
        Z = H_fun(w) - mass * w**2 * np.eye(2)
        Xv = np.linalg.solve(Z, Fv)
        Av = -w**2 * Xv
        e = np.exp(1j * w * t)
        F_t += np.real(np.outer(Fv, e))
        X_t += np.real(np.outer(Xv, e))
        A_t += np.real(np.outer(Av, e))
    rel_disp = X_t.copy()
    if runout:
        # synchronous rotor runout / unbalance seen only by the proximity probes
        for h, a in [(1, 3.5e-6), (2, 1.2e-6), (3, 0.6e-6)]:
            w = 2 * np.pi * f_sync * h
            ph = 2 * np.pi * rng.random(2)
            rel_disp[0] += a * np.cos(w * t + ph[0])
            rel_disp[1] += a * np.sin(w * t + ph[1])
        rel_disp += np.array([[12.0e-6], [-35.0e-6]])  # static eccentricity
    # sensor noise
    F_t += rng.normal(0, 2.0, F_t.shape)
    A_t += rng.normal(0, 0.05, A_t.shape)
    rel_disp += rng.normal(0, 0.015e-6, rel_disp.shape)
    hdr = "time_s,force_x_N,force_y_N,accel_x_g,accel_y_g,disp_x_um,disp_y_um"
    arr = np.column_stack([t, F_t[0], F_t[1], A_t[0] / G0, A_t[1] / G0,
                           rel_disp[0] * 1e6, rel_disp[1] * 1e6])
    np.savetxt(DATA / f"{label}.csv", arr, delimiter=",", header=hdr, comments="",
               fmt=["%.7f", "%.4f", "%.4f", "%.6f", "%.6f", "%.5f", "%.5f"])


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    total = lambda w: H(K, C, M, w) + H(Ks, Cs, np.zeros((2, 2)), w)
    tare = lambda w: H(Ks, Cs, np.zeros((2, 2)), w)
    make_test(0, total, True, "bearing_shake_x", m_stator)
    make_test(1, total, True, "bearing_shake_y", m_stator)
    make_test(0, tare, False, "tare_shake_x", m_tare)
    make_test(1, tare, False, "tare_shake_y", m_tare)
    rig = {
        "sample_rate_hz": fs,
        "record_length_s": T,
        "rotor_speed_rpm": rpm,
        "stator_mass_kg": m_stator,
        "tare_stator_mass_kg": m_tare,
        "notes": "stator_mass_kg applies to the bearing_* tests (stator with bearing insert installed); tare_stator_mass_kg applies to the tare_* tests (bearing insert removed).",
        "standard_gravity_m_per_s2": G0,
        "channels": {
            "force_x_N": "shaker force applied to the stator in +x, newtons",
            "force_y_N": "shaker force applied to the stator in +y, newtons",
            "accel_x_g": "absolute stator acceleration in +x, units of g",
            "accel_y_g": "absolute stator acceleration in +y, units of g",
            "disp_x_um": "stator-minus-rotor relative displacement in +x from proximity probes, micrometres",
            "disp_y_um": "stator-minus-rotor relative displacement in +y from proximity probes, micrometres",
        },
    }
    (DATA / "rig.json").write_text(json.dumps(rig, indent=2) + "\n")
    truth = {"frequencies_hz": freqs.tolist()}
    for name, mat in (("K", K), ("C", C), ("M", M)):
        for i, a in enumerate("xy"):
            for j, b in enumerate("xy"):
                truth[f"{name}{a}{b}"] = float(mat[i, j])
    (ROOT / "tests" / "expected.json").write_text(json.dumps(truth, indent=2) + "\n")


if __name__ == "__main__":
    main()
