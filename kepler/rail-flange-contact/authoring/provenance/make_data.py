#!/usr/bin/env python3
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from geometry2 import rail_raw, wheel_raw, sample
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "environment", "data")
os.makedirs(OUT, exist_ok=True)
yr, zr = sample(*rail_raw()); yw, zw = sample(*wheel_raw())
zw = zw - np.interp(0.0, yw, zw)
for name, y, z in (("rail_profile.csv", yr, zr), ("wheel_profile.csv", yw, zw)):
    with open(os.path.join(OUT, name), "w") as f:
        f.write("y_mm,z_mm\n")
        for a, b in zip(y, z):
            f.write(f"{a:.3f},{b:.5f}\n")
V = 25.0; r0 = 460.0; Om0 = V / (r0 * 1e-3)
def lc(i, Q, Y, psi_mrad, brake=0.0):
    return {"id": i, "vertical_load_N": Q, "lateral_force_N": Y, "yaw_angle_rad": psi_mrad * 1e-3,
            "speed_m_per_s": V, "wheel_angular_velocity_rad_per_s": round(Om0 * (1 - brake), 6)}
case = {
    "materials": {"youngs_modulus_Pa": 210e9, "poisson_ratio": 0.28, "note": "wheel and rail are the same steel"},
    "coefficient_of_friction": 0.35,
    "wheel": {"tape_circle_radius_mm": r0, "tape_circle_y_mm": 0.0,
              "profile": "wheel_profile.csv: wheel surface height z_w (mm, upward) against wheel lateral coordinate y_w (mm), section through the wheel axis; z_w = 0 at the tape circle (y_w = 0); the wheel radius at y_w is r(y_w) = tape_circle_radius_mm - z_w(y_w)"},
    "rail": {"profile": "rail_profile.csv: rail head height z_r (mm, upward) against track lateral coordinate y (mm); the rail is straight"},
    "frames": {
        "track": "x along the track in the direction of travel, y lateral pointing toward the track centre (the flange side), z upward",
        "wheel_position": "the wheel profile sits at y = y_w + lateral_shift and its surface at height z_w(y - lateral_shift) + vertical_position; lateral_shift and vertical_position are unknown and follow from equilibrium (no wheel roll)",
        "contact_search": "the vertical separation g(y) = z_w(y - lateral_shift) + vertical_position - z_r(y); each local minimum of g is a contact point; the potential patch of a contact point is the part of the profiles between the neighbouring local maxima of g; a patch is in contact when the bodies interpenetrate there",
        "patch_plane": "each patch is planar, in the plane tangent to the rail profile at its contact point; the contact angle gamma is the rail slope angle there (atan dz_r/dy, negative where the rail falls toward +y); local axes: x (travel), s = (cos gamma, sin gamma) in (y, z), n = (-sin gamma, cos gamma) pointing from the rail into the wheel",
        "patch_gap": "undeformed normal separation in a patch: both profiles (wheel in its current position) rotated into the (s, n) axes, separation n_wheel(s) - n_rail(s), plus x^2 / (2 R) with R = r_i / cos(gamma), r_i the wheel radius at the contact point (the rail is straight); negative values are interpenetration",
        "elastic_independence": "the patches are far apart and are treated as elastically independent half-spaces; materials are quasi-identical so normal and tangential problems decouple"
    },
    "kinematics": {
        "rigid_slip": "rigid slip velocity of wheel material relative to the rail, divided by the speed V, at local patch coordinates (x, s): (xi - phi*s, eta + phi*x)",
        "xi": "xi_i = 1 - (Omega * r_i / V) * cos(psi)",
        "eta": "eta_i = -(Omega * r_i / V) * sin(psi) * cos(gamma_i)",
        "phi": "phi_i = -(Omega / V) * cos(psi) * sin(gamma_i)   (1/m)",
        "rolling": "steady rolling; material passes through each patch from +x (leading edge) to -x"
    },
    "equilibrium": "the wheelset applies to this wheel a downward force Q and a lateral force Y toward -y (pressing the flange against the rail); in equilibrium the contact forces exerted by the rail on the wheel, summed over the patches and resolved in the track frame, give sum F_z = Q and sum F_y = Y. A patch force with normal component N and tangential components (F_x, F_s) contributes F_y = F_s cos(gamma) - N sin(gamma), F_z = F_s sin(gamma) + N cos(gamma)",
    "initial_guess": {"lateral_shift_mm": 7.5, "vertical_position_mm": 0.0},
    "load_cases": [
        lc("A", 80e3, 20e3, 0.0),
        lc("B", 80e3, 30e3, -1.5),
        lc("C", 80e3, 40e3, -3.0),
        lc("D", 80e3, 25e3, 1.5),
        lc("E", 80e3, 30e3, -1.5, brake=0.0012),
        lc("F", 60e3, 36e3, -3.0),
    ],
}
json.dump(case, open(os.path.join(OUT, "case.json"), "w"), indent=2)
print("ok", len(yr), len(yw))
