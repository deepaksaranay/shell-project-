#!/bin/bash
set -euo pipefail
mkdir -p /app/output

cat > /tmp/identify.py <<'PY'
import json
import numpy as np

DATA = "/app/data"
rig = json.load(open(f"{DATA}/rig.json"))
fs = rig["sample_rate_hz"]
MASS = {"bearing": rig["stator_mass_kg"], "tare": rig["tare_stator_mass_kg"]}
g0 = rig.get("standard_gravity_m_per_s2", 9.80665)


def load(name):
    d = np.genfromtxt(f"{DATA}/{name}.csv", delimiter=",", names=True)
    F = np.vstack([d["force_x_N"], d["force_y_N"]])
    A = np.vstack([d["accel_x_g"], d["accel_y_g"]]) * g0
    X = np.vstack([d["disp_x_um"], d["disp_y_um"]]) * 1e-6
    n = F.shape[1]
    f = np.fft.rfftfreq(n, 1.0 / fs)
    return f, np.fft.rfft(F, axis=1), np.fft.rfft(A, axis=1), np.fft.rfft(X, axis=1)


def tone_bins(Fspec):
    """Excitation bins = clear peaks in the shaker force spectrum."""
    mag = np.abs(Fspec).max(axis=0)
    mag[0] = 0.0
    thr = 20.0 * np.median(mag)
    return np.flatnonzero(mag > thr)


def impedance(prefix):
    m_s = MASS[prefix]
    f, Fx, Ax, Xx = load(f"{prefix}_shake_x")
    _, Fy, Ay, Xy = load(f"{prefix}_shake_y")
    bins = np.intersect1d(tone_bins(Fx), tone_bins(Fy))
    Hs = []
    for b in bins:
        # columns = the two tests
        Fm = np.column_stack([Fx[:, b] - m_s * Ax[:, b], Fy[:, b] - m_s * Ay[:, b]])
        Dm = np.column_stack([Xx[:, b], Xy[:, b]])
        Hs.append(Fm @ np.linalg.inv(Dm))
    return f[bins], np.array(Hs)


fb, Hb = impedance("bearing")
ft, Ht = impedance("tare")
assert np.allclose(fb, ft), (fb, ft)
H = Hb - Ht
w = 2 * np.pi * fb

out = {"frequencies_hz": [float(v) for v in fb]}
# Re(H) = K - w^2 M ; Im(H) = w C  (independent least squares per element)
Ar = np.column_stack([np.ones_like(w), -w**2])
for i, a in enumerate("xy"):
    for j, b in enumerate("xy"):
        (k, m), *_ = np.linalg.lstsq(Ar, H[:, i, j].real, rcond=None)
        c = float(np.sum(w * H[:, i, j].imag) / np.sum(w * w))
        out[f"K{a}{b}"] = float(k)
        out[f"C{a}{b}"] = c
        out[f"M{a}{b}"] = float(m)

json.dump(out, open("/app/output/bearing_coefficients.json", "w"), indent=2)
print(json.dumps(out, indent=2))
PY

python3 /tmp/identify.py
