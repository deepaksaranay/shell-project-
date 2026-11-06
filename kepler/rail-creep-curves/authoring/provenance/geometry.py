"""Wheel and rail transverse profiles for the task (track coordinates, mm).

Rail head: crown arc R300 for |y - yc| <= 10 mm, blended tangentially into
R80 shoulder arcs (profile defined by its curvature, integrated twice).
Wheel tread: worn (hollow) tread with concave transverse radius 360 mm over
the worn band, convex 1:40 cone outside it, smoothly blended; the wheel
profile is given as the height of the wheel surface above the rail datum.
"""
import numpy as np

Y = np.round(np.arange(-40.0, 40.0 + 1e-9, 0.1), 1)


def _integrate_curvature(y, kappa, y0=0.0):
    """z(y) with z(y0)=0, z'(y0)=0 and z'' = -kappa (for a crown, kappa>0)."""
    dy = y[1] - y[0]
    i0 = int(np.argmin(np.abs(y - y0)))
    slope = np.zeros_like(y); z = np.zeros_like(y)
    for i in range(i0 + 1, len(y)):
        slope[i] = slope[i - 1] - 0.5 * (kappa[i] + kappa[i - 1]) * dy
        z[i] = z[i - 1] + 0.5 * (slope[i] + slope[i - 1]) * dy
    for i in range(i0 - 1, -1, -1):
        slope[i] = slope[i + 1] + 0.5 * (kappa[i] + kappa[i + 1]) * dy
        z[i] = z[i + 1] - 0.5 * (slope[i] + slope[i + 1]) * dy
    return z


def smoothstep(t):
    t = np.clip(t, 0, 1)
    return t * t * (3 - 2 * t)


def rail_profile():
    y = Y.copy()
    yc = 0.0
    d = np.abs(y - yc)
    k = 1 / 300 + (1 / 80 - 1 / 300) * smoothstep((d - 9.0) / 2.0)
    return y, _integrate_curvature(y, k, yc)


def wheel_profile():
    y = Y.copy()
    # hollow-worn band centred at +1.0 mm: concave (negative curvature) R360,
    # blending to the convex flank/cone (curvature ~ 1/150) outside |d| > 12 mm
    d = y.copy()
    k_w = -1 / 360 + (1 / 150 + 1 / 360) * smoothstep((np.abs(d) - 11.0) / 3.0)
    zw = -_integrate_curvature(y, k_w, 0.0)       # wheel surface bends upward where convex
    zw = zw + 0.0012 * y                         # slight relative tilt (wheelset roll/cone)
    # small asymmetric wear feature: a 4 um deep, 6 mm wide groove at y = +6 mm
    zw = zw + 0.004 * np.exp(-((y - 6.0) / 3.0) ** 2)
    return y, zw
