# Validation and shortcut attempts (shipped data)

Limits: EO exact; f within 0.25%; zeta within 25% per blade and fleet-mean signed error within +/-7%; steady resonant amplitude within 15% per blade.

```
constant f=706, zeta=0.0024, amp=0.4 mm, EO=6: df_max=0.0209 dz_max=0.950 dz_mean=+0.211 da_max=0.887 -> fail
steady-state SDOF fit + measured peak: df_max=0.0032 dz_max=0.353 dz_mean=+0.152 da_max=0.271 -> fail
reference f, zeta + measured (transient) peak amplitude: df_max=0.0003 dz_max=0.075 dz_mean=+0.005 da_max=0.271 -> fail
reference f, amp + fleet-mean zeta: df_max=0.0003 dz_max=0.733 dz_mean=+0.076 da_max=0.034 -> fail
reference but EO=3: df_max=0.0003 dz_max=0.075 dz_mean=+0.005 da_max=0.034 -> fail
reference: df_max=0.0003 dz_max=0.075 dz_mean=+0.005 da_max=0.034 -> PASS
```

- Skipping the OPR repair: blade assignment breaks and the analysis fails (least-squares does not converge).
- Verifier (pytest 8.4.1 + pytest-json-ctrf 0.3.5): nop -> reward 0 (75 failed); oracle -> reward 1 on three repeat runs (75 passed).
- Seed sweep (seeds 11/23/37/41/53 + shipped 7710): reference passes on all; steady-state SDOF + measured peak and no-Southwell fail on all.
