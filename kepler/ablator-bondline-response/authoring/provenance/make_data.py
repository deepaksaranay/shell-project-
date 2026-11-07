#!/usr/bin/env python3
"""Generate the synthetic (but realistic) material, surface-thermochemistry and
trajectory data for the charring-ablator task. Nothing here is taken from a
published material card; the numbers are chosen to be in the range of a
mid-density carbon/phenolic-type ablator on an aluminium structure."""
import json
import os
import numpy as np

OUT = os.path.join(os.path.dirname(__file__), "..", "..", "environment", "data")
os.makedirs(OUT, exist_ok=True)

T = np.arange(250.0, 4001.0, 50.0)

def wcsv(name, header, cols, fmt):
    with open(os.path.join(OUT, name), "w") as f:
        f.write(",".join(header) + "\n")
        for row in zip(*cols):
            f.write(",".join(fm % v for fm, v in zip(fmt, row)) + "\n")

# ---- solid properties (virgin / char), tabulated vs temperature
x = (T - 250) / 3750
cp_v = 960 + 1150 * (1 - np.exp(-(T - 250) / 520)) + 40 * x
k_v = 0.215 + 0.21 * x + 0.08 * x ** 2
cp_c = 720 + 1380 * (1 - np.exp(-(T - 250) / 900))
k_c = 0.42 + 0.36 * x + 3.1 * x ** 3          # radiative contribution at high T
eps_v = np.full_like(T, 0.86); eps_c = np.full_like(T, 0.90)
wcsv("virgin_properties.csv", ["T_K", "cp_J_per_kgK", "k_W_per_mK", "emissivity"],
     [T, cp_v, k_v, eps_v], ["%.0f", "%.2f", "%.5f", "%.3f"])
wcsv("char_properties.csv", ["T_K", "cp_J_per_kgK", "k_W_per_mK", "emissivity"],
     [T, cp_c, k_c, eps_c], ["%.0f", "%.2f", "%.5f", "%.3f"])

# ---- pyrolysis gas enthalpy (equilibrium mixture, absolute, J/kg)
h_pg = -1.45e6 + 1850 * (T - 298) + 0.62 * (T - 298) ** 2 + 2.2e6 / (1 + np.exp(-(T - 2300) / 260))
wcsv("pyrolysis_gas.csv", ["T_K", "h_J_per_kg"], [T, h_pg], ["%.0f", "%.1f"])

# ---- aluminium substrate
Ts = np.arange(250.0, 901.0, 50.0)
cp_al = 880 + 0.47 * (Ts - 300) - 1.2e-4 * (Ts - 300) ** 2
k_al = 165 + 0.02 * (Ts - 300)
wcsv("substrate_properties.csv", ["T_K", "cp_J_per_kgK", "k_W_per_mK"], [Ts, cp_al, k_al], ["%.0f", "%.2f", "%.3f"])

# ---- surface thermochemistry table B'_c(T, B'_g), h_w(T, B'_g)
Tb = np.arange(300.0, 4001.0, 100.0)
Bg = np.array([0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0])
rows = []
for t in Tb:
    for bg in Bg:
        ox = (0.168 + 0.05 * bg / (1 + bg)) / (1 + np.exp(-(t - 980) / 70))
        sub = 1.6e-3 * np.exp((t - 3000) / 140) * (1 + 0.6 * bg)
        bc = ox + sub
        h_air = 1010 * (t - 300) + 0.16 * (t - 300) ** 2
        hpg = np.interp(t, T, h_pg)
        hw = (h_air + bc * (-2.6e6 + 1900 * (t - 300)) + bg * hpg) / (1 + bc + bg)
        rows.append((t, bg, bc, hw))
r = np.array(rows)
wcsv("surface_thermochemistry.csv", ["T_w_K", "Bprime_g", "Bprime_c", "h_w_J_per_kg"],
     [r[:, 0], r[:, 1], r[:, 2], r[:, 3]], ["%.0f", "%.3f", "%.6f", "%.1f"])

# ---- trajectory (heating environment at the stagnation point)
t = np.arange(0.0, 600.1, 2.0)
pulse = np.where(t <= 110, np.sin(np.pi * np.clip(t, 0, 110) / 110) ** 1.6, 0.0)
CH0 = 0.0035 + 0.17 * pulse * (1 - 0.25 * t / 110)
hr = np.where(t <= 110, 2.7e7 - 2.3e7 * (t / 110) ** 1.3, 3.0e5)
hr = np.where(t > 110, 3.0e5, hr)
qrad = 7.5e5 * np.where(t <= 110, np.sin(np.pi * np.clip(t, 0, 110) / 110) ** 4, 0.0)
wcsv("trajectory.csv", ["t_s", "rhoe_ue_CH0_kg_per_m2s", "h_r_J_per_kg", "q_rad_in_W_per_m2"],
     [t, CH0, hr, qrad], ["%.1f", "%.6f", "%.1f", "%.1f"])

case = {
    "description": "Stagnation-point heat-shield response: charring ablator bonded to an aluminium structure",
    "ablator": {
        "thickness_mm": 22.0,
        "decomposition": {
            "model": "three-component: rho = Gamma*(rho_A + rho_B) + (1 - Gamma)*rho_C; for i in A, B: d(rho_i)/dt = -k_i * rho_v_i * ((rho_i - rho_r_i)/rho_v_i)**psi_i * exp(-E_over_R_i / T); C does not decompose",
            "Gamma": 0.5,
            "components": {
                "A": {"rho_v_kg_per_m3": 230.0, "rho_r_kg_per_m3": 0.0, "k_per_s": 1.6e4, "E_over_R_K": 8700.0, "psi": 3.0},
                "B": {"rho_v_kg_per_m3": 690.0, "rho_r_kg_per_m3": 460.0, "k_per_s": 3.9e9, "E_over_R_K": 20300.0, "psi": 3.0},
                "C": {"rho_kg_per_m3": 300.0}
            }
        },
        "virgin": {"properties_file": "virgin_properties.csv", "formation_enthalpy_J_per_kg_at_298K": -1.05e6},
        "char": {"properties_file": "char_properties.csv", "formation_enthalpy_J_per_kg_at_298K": 0.0},
        "pyrolysis_gas": {"enthalpy_file": "pyrolysis_gas.csv"},
        "property_mixing": "Solid cp, k, emissivity and absorptivity are blended between char and virgin with the virgin mass fraction y_v = rho_v*(rho - rho_c)/(rho*(rho_v - rho_c)): property = y_v*virgin + (1 - y_v)*char. The solid enthalpy per unit mass is h = y_v*h_v(T) + (1 - y_v)*h_c(T), with h_v(T) = h_f,v + integral from 298 K of cp_v dT and likewise for char."
    },
    "substrate": {"material": "aluminium", "thickness_mm": 3.0, "density_kg_per_m3": 2700.0,
                  "properties_file": "substrate_properties.csv",
                  "bond": "perfect thermal contact, impermeable to pyrolysis gas"},
    "back_face": "adiabatic",
    "initial_temperature_K": 300.0,
    "surface": {
        "thermochemistry_file": "surface_thermochemistry.csv",
        "absorptivity_equals_emissivity": True,
        "radiation_sink_temperature_K": 300.0,
        "blowing_correction": "C_H = C_H0 * ln(1 + 2*lambda*B') / (2*lambda*B'), B' = B'_c + B'_g, lambda = 0.5 (C_H = C_H0 at B' = 0)",
        "definitions": "B'_g = mdot_g / (rhoe_ue_CH), B'_c = mdot_c / (rhoe_ue_CH), with rhoe_ue_CH the blowing-corrected value; the char removal rate mdot_c gives the recession rate sdot = mdot_c / rho_w, rho_w the solid density at the surface",
        "energy_balance": "rhoe_ue_CH*(h_r - h_w) + mdot_g*h_g(T_w) + mdot_c*h_s(T_w) - (mdot_g + mdot_c)*h_w + alpha*q_rad_in - eps*sigma*(T_w**4 - T_sink**4) = conducted heat flux into the solid at the surface",
        "pyrolysis_gas": "flows to the heated surface instantaneously without storage, in thermal equilibrium with the solid; mdot_g is the total decomposition rate integrated through the thickness"
    },
    "interpolation": "All tables are interpolated linearly (thermochemistry: bilinear in T_w and B'_g; cp linear in T so that h(T) is its exact integral); values outside a table are held at the end value.",
    "stefan_boltzmann_W_per_m2K4": 5.670374419e-8,
    "output_times_s": [20, 40, 60, 80, 100, 150, 300, 600],
    "definitions_for_output": {
        "recession_mm": "surface recession s(t) from the original surface",
        "surface_temperature_K": "T_w(t)",
        "bondline_temperature_K": "temperature at the ablator/aluminium interface",
        "pyrolysis_front_depth_mm": "depth below the current surface where the density degree tau = (rho - rho_c)/(rho_v - rho_c) first reaches 0.98, linear interpolation between grid values",
        "char_front_depth_mm": "depth below the current surface where tau first reaches 0.02",
        "mass_loss_kg_per_m2": "integral of mdot_g + mdot_c from 0 to t"
    }
}
with open(os.path.join(OUT, "case.json"), "w") as f:
    json.dump(case, f, indent=2)
print("written", sorted(os.listdir(OUT)))
