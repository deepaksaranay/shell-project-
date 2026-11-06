# Validation record

Verifier emulation uses the same test.sh with pytest 8.4.1 and pytest-json-ctrf 0.3.5.

| case | reward |
|---|---|
| nop (no design file) | 0 (1 failed, 3 errors) |
| oracle: 17.9 m³, 7.966 m³ gas, 0.12 m orifice | 1 (4 passed). Min abs head 3.90 m, max gauge 102.4 m, gas 0.870 |

## Lazy and mistaken designs

```
budget, half gas, big orifice (18.0 m3, gas 9.0, d 0.35): ========================= 2 failed, 2 passed in 2.52s ========================== reward=0
budget, 30% gas (18.0 m3, gas 5.4, d 0.10): ========================= 1 failed, 3 passed in 2.38s ========================== reward=0
budget, 70% gas (18.0 m3, gas 12.6, d 0.10): ========================= 1 failed, 3 passed in 2.54s ========================== reward=0
pump-node-only sizing (12.0 m3, gas 4.2, d 0.08): ========================= 1 failed, 3 passed in 2.45s ========================== reward=0
isothermal-gas sizing (15.0 m3, gas 6.0, d 0.08): ========================= 1 failed, 3 passed in 2.59s ========================== reward=0
oversized (25.0 m3, gas 11.0, d 0.10): ========================= 1 failed, 3 passed in 2.47s ========================== reward=0
budget, 45% gas, min orifice (18.0 m3, gas 8.1, d 0.05): ============================== 4 passed in 2.35s =============================== reward=1
```

The only lazy design that passes sits inside the narrow feasible gas window (0.42–0.46 at the budget volume). The usual 50%-gas rule of thumb fails.
