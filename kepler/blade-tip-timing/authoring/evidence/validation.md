# Shortcut / cheat attempts (shipped data)

Scored against the stated requirements (EO exact, f within 0.25%, zeta within 25% per blade, mean zeta error within +/-7%).

```
constant f=706 Hz, zeta=0.0024, EO=6: df_max=0.0209 dz_max=0.950 dz_mean=+0.211 -> fail
steady-state SDOF fit: df_max=0.0032 dz_max=0.351 dz_mean=+0.152 -> fail
steady f + fleet-mean zeta: df_max=0.0032 dz_max=0.953 dz_mean=+0.213 -> fail
reference f + reference fleet-mean zeta (no per-blade damping): df_max=0.0003 dz_max=0.734 dz_mean=+0.077 -> fail
reference but EO=3 (alias): df_max=0.0003 dz_max=0.074 dz_mean=+0.005 -> fail
reference: df_max=0.0003 dz_max=0.074 dz_mean=+0.005 -> PASS
```

Verifier: nop -> reward 0 (51 failed); oracle -> reward 1 on three repeat runs (51 passed), pytest 8.4.1 + pytest-json-ctrf 0.3.5, CTRF report written to /logs/verifier/ctrf.json.
Seed sweep (authoring/provenance/seed_sweep.sh, seeds 11/23/37 + shipped 7710): reference passes on all; steady-state SDOF and no-Southwell fail on all.
