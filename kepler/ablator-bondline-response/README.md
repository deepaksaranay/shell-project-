# ablator-bondline-response

The task is to compute the one-dimensional thermal response of a charring
heat-shield ablator, 22 mm thick, bonded to 3 mm of aluminium, over a
stagnation-point entry heating pulse (110 s) and the heat soak after it, out
to 600 s. At eight output times it reports:

- surface recession;
- surface and bondline temperatures;
- pyrolysis-front and char-front depths;
- mass loss.

It also reports the peak bondline temperature.

## Difficulty

**Who does this work.** In practice this analysis is done by
thermal-protection-system (TPS) engineers at spacecraft and entry-vehicle
developers, space agencies and their contractors. They use it to:

- size the heat-shield thickness so that the bondline to the aluminium
  structure stays below the adhesive limit through entry and the heat soak
  after it;
- rebuild flight data with the same models.

The answer needs a full charring-ablation model, with every piece consistent
with the others:

- three-component Arrhenius decomposition (two resin phases and an inert
  fibre);
- pyrolysis gas that picks up heat on its way out through the char;
- absolute, formation-referenced enthalpies, which set the heat of
  pyrolysis;
- virgin/char property blending by virgin mass fraction;
- a receding surface;
- an aluminium substrate;
- a surface energy balance using the B′ table with the blowing correction.

In that surface balance the corrected film coefficient depends on
B′_c(T_w, B′_g), and B′_g depends on the coefficient, so it is a fixed
point at every instant.

The data are synthetic but realistic (`authoring/provenance/make_data.py`):
a mid-density phenolic-type ablator (610 → 380 kg/m³), a B′ table with an
oxidation plateau and the onset of sublimation, and peak cold-wall heating
of a few MW/m². No published result applies.

Each shortcut fails clearly. These were scored through the verifier, with
details in `authoring/evidence/validation.md`:

| model | what goes wrong | verifier |
|---|---|---|
| no blowing correction | recession +23%, surface +225 K | 0/7 |
| pyrolysis-gas heat pick-up ignored | bondline +37 K, char front +1.9 mm | 1/7 |
| no recession | recession and mass loss | 0/7 |
| B′_c read at B′_g = 0 | recession −5.7% | 6/7 (fails) |
| 12-node uniform grid | bondline +15 K, surface −85 K | 1/7 |

The moving surface makes the bookkeeping easy to get subtly wrong. While
building the generator, two bugs each shifted the bondline by several
kelvin while the other outputs still looked plausible:

- an upwind index in the grid-motion transport;
- a missing fibre-mass term.

Both were caught only by cross-checking against independent codes and an
analytic slab solution (below).

## Reference solution

`solution/solve.py` (about 6 minutes on 4 CPUs) writes the energy equation
in the non-conservative form

  ρ c_p ∂T/∂t = ∂/∂x(k ∂T/∂x) + ṁ ∂h_g/∂x + (h_g − h̄) ∂ρ/∂t,
  h̄ = (ρ_v h_v − ρ_c h_c)/(ρ_v − ρ_c).

It then takes these steps:

1. Gets the gas flux by integrating the decomposition rate through the
   thickness.
2. Transforms the equation to the Landau coordinate η = (x − s)/(L − s),
   which follows the receding surface.
3. Discretises it with finite differences on a graded 150-node grid.
4. Takes the surface temperature as the root of the surface energy balance
   against a second-order one-sided conduction flux, with an inner root
   find for the blowing-corrected film coefficient.
5. Treats the bondline as a two-material half-cell and integrates in time
   with Radau.

Run through the verifier, it passes 7/7 and uses at most 21% of any limit:

| check | worst error |
|---|---|
| bondline temperature | +0.62 K |
| peak bondline temperature | +0.53 K |
| pyrolysis front | +0.061 mm |
| surface temperature | +0.72 K |

## Verification

`tests/test_outputs.py` reads `/app/output/heatshield_response.json`. It
must be a regular file (no symlinks) with exactly the documented keys, one
snapshot per output time and finite numbers. At every output time:

| quantity | limit |
|---|---|
| recession | 3% + 0.01 mm |
| surface temperature | 12 K |
| bondline temperature | 3 K |
| pyrolysis-front depth | 3% + 0.1 mm |
| char-front depth | 3% + 0.1 mm |
| mass loss | 2% + 0.01 kg/m² |
| peak bondline temperature | 3 K |

**Where the answer key comes from.** `tests/answer_key.json` is written by
`authoring/provenance/make_answer_key.py` from the generator
`authoring/provenance/ablation.py`:

- conservative finite volumes on a contracting, surface-attached grid with
  Reynolds-transport terms for the moving control volumes;
- BDF time integration;
- 480 ablator nodes. The 240-node run differs by at most 0.30 K and 0.3%.

The shipped solution does not produce the key.

**How the limits were set.** The limits come from how closely structurally
different correct implementations agree with the key. Each one was scored
through the verifier:

| implementation | how it differs | largest share of a limit | result |
|---|---|---|---|
| `authoring/provenance/independent_submission.py` | separately written; material-fixed uniform cell-centred grid with surface-cell consumption and merging (no grid motion), own backward-Euler/Picard stepping, series-resistance bondline | 8% (char front +0.016 mm; bondline 0.01 K) | 7/7 |
| `solution/solve.py` (shipped) | non-conservative Landau-coordinate finite differences, algebraic surface node, Radau | 21% (bondline +0.62 K) | 7/7 |
| generator on modest grids | 60 or 100 uniform nodes, 60 graded nodes | up to 77% (60 uniform: peak bondline +2.3 K) | 7/7 each |
| generator, 240 nodes | — | 10% | — |

Before the comparison, the generator and the shipped solution were both
checked against the analytic solution for a constant-flux slab with an
insulated back. Both matched to within 0.1 K at the surface and 0.06 K at
the bondline. They were also checked against each other on cases with
recession only and with no decomposition.

Correct methods sit well inside the limits, and every shortcut above fails
by a wide margin. A missing file, a symlink, a copy of the answer key, or
extra or renamed keys all score 0.

Outputs are in `authoring/evidence/`.
