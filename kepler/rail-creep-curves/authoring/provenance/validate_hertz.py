# Validation against Kalker's table (circular contact, g = a/b = 1):
#   nu = 0.00: C11 = C22 = 3.40 ;  nu = 0.25: C11 = 4.12, C22 = 3.67
# and full sliding |F| -> mu N at large creepage.
import math, time
import numpy as np
from contact import Grid, solve_normal, solve_tangential
E, R, N, mu = 210e9, 0.3, 50e3, 0.3
for nu in (0.0, 0.25):
    G = E / (2 * (1 + nu))
    for n in (24, 32):
        gr = Grid(-0.0075, 0.0075, n, -0.0075, 0.0075, n)
        h = (gr.X**2 + gr.Y**2) / (2 * R)
        t = time.time()
        p, d, _ = solve_normal(gr, h, N, E, nu)
        c2 = (p > 0).sum() * gr.dA / math.pi
        r = solve_tangential(gr, p, E, nu, mu, 1e-6, 0, 0)
        K = r["K"]
        r2 = solve_tangential(gr, p, E, nu, mu, 0, 1e-6, 0, K=K)
        line = f"nu={nu} n={n}: C11={-r['Fx']/(G*c2*1e-6):.3f} C22={-r2['Fy']/(G*c2*1e-6):.3f}"
        for xi in (5e-4, 1e-3, 2e-3, 1e-2):
            rr = solve_tangential(gr, p, E, nu, mu, xi, 0, 0, K=K)
            line += f" | xi={xi:g}: Fx/muN={rr['Fx']/(mu*N):+.3f} adh={rr['adhesion_fraction']:.2f} it={rr['iterations']}{'' if rr['converged'] else ' NOCONV'}"
        print(line, f"({time.time()-t:.1f}s)", flush=True)
