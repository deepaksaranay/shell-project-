"""Half-space steady rolling contact (Kalker's exact theory, discretised).

Quasi-identical elastic bodies (same E, nu), so the normal and tangential
problems decouple. Uniform-traction rectangular elements on a regular grid;
influence coefficients from closed-form rectangle integrals of the
Boussinesq/Cerruti point solutions (relative displacement of the two bodies).

Normal problem: active-set solution of  h + u_z(p) - delta = 0 on C,
p >= 0, gap >= 0 off C, sum(p) dA = N.

Tangential problem (steady rolling, material flowing in -x, leading edge at
+x): slip s = w + (u_I - u_I')/dx with I' the element one step upstream
(x + dx), rigid slip w = (xi - phi*y, eta + phi*x). In C:
adhesion s = 0 with |q| <= mu*p; slip q = -mu*p*s/|s|. Solved by
element-wise Gauss-Seidel sweeps from the leading edge.
"""
import math

import numpy as np


def _asinh_ratio(a, b):
    """a*asinh(b/|a|), with the a = 0 limit 0."""
    out = np.zeros_like(a)
    m = a != 0
    out[m] = a[m] * np.arcsinh(b[m] / np.abs(a[m]))
    return out


def _rect_integrals(xo, yo, xc, yc, dx, dy):
    """Integrals over source rectangles (centre xc, yc) seen from points (xo, yo)
    of 1/r, x^2/r^3, y^2/r^3 and x*y/r^3 (x, y = observation - source)."""
    X = xo[:, None] - xc[None, :]
    Y = yo[:, None] - yc[None, :]
    I1 = np.zeros_like(X); Ixx = np.zeros_like(X); Iyy = np.zeros_like(X); Ixy = np.zeros_like(X)
    for sx, ax in ((1, X + dx / 2), (-1, X - dx / 2)):
        for sy, ay in ((1, Y + dy / 2), (-1, Y - dy / 2)):
            sg = sx * sy
            fx = _asinh_ratio(ax, ay)   # x asinh(y/|x|)
            fy = _asinh_ratio(ay, ax)   # y asinh(x/|y|)
            I1 += sg * (fx + fy)
            Ixx += sg * fy
            Iyy += sg * fx
            Ixy += sg * (-np.sqrt(ax * ax + ay * ay))
    return I1, Ixx, Iyy, Ixy


class Grid:
    def __init__(self, x0, x1, nx, y0, y1, ny):
        self.dx = (x1 - x0) / nx
        self.dy = (y1 - y0) / ny
        self.xc = x0 + self.dx * (np.arange(nx) + 0.5)
        self.yc = y0 + self.dy * (np.arange(ny) + 0.5)
        self.nx, self.ny = nx, ny
        X, Y = np.meshgrid(self.xc, self.yc, indexing="ij")
        self.X, self.Y = X.ravel(), Y.ravel()      # index = ix*ny + iy
        self.dA = self.dx * self.dy


def solve_normal(grid, h, N, E, nu, max_iter=200):
    G = E / (2 * (1 + nu))
    I1, _, _, _ = _rect_integrals(grid.X, grid.Y, grid.X, grid.Y, grid.dx, grid.dy)
    Azz = (1 - nu) / (math.pi * G) * I1
    n = len(h)
    C = np.argsort(h)[: max(4, n // 20)]
    C = np.sort(C)
    for _ in range(max_iter):
        m = len(C)
        M = np.zeros((m + 1, m + 1))
        M[:m, :m] = Azz[np.ix_(C, C)]
        M[:m, m] = -1.0
        M[m, :m] = grid.dA
        rhs = np.concatenate([-h[C], [N]])
        sol = np.linalg.solve(M, rhs)
        pC, delta = sol[:m], sol[m]
        if np.any(pC < 0):
            C = C[pC >= 0] if np.any(pC >= 0) else C
            continue
        p = np.zeros(n); p[C] = pC
        gap = h + Azz @ p - delta
        outside = np.setdiff1d(np.arange(n), C)
        pen = outside[gap[outside] < -1e-12 * max(1.0, abs(delta))]
        if len(pen) == 0:
            return p, delta, Azz
        C = np.sort(np.concatenate([C, pen]))
    raise RuntimeError("normal problem did not converge")


def tangential_influence(grid, E, nu, rows, cols):
    G = E / (2 * (1 + nu))
    I1, Ixx, Iyy, Ixy = _rect_integrals(grid.X[rows], grid.Y[rows], grid.X[cols], grid.Y[cols], grid.dx, grid.dy)
    c = 1.0 / (math.pi * G)
    return c * ((1 - nu) * I1 + nu * Ixx), c * ((1 - nu) * I1 + nu * Iyy), c * nu * Ixy


def _kmatrix(grid, E, nu, C):
    ix, iy = C // grid.ny, C % grid.ny
    up = (ix + 1) * grid.ny + iy                      # upstream neighbour (x + dx)
    valid_up = ix + 1 < grid.nx
    Uxx, Uyy, Uxy = tangential_influence(grid, E, nu, C, C)
    Ux2, Uy2, Ux2y = tangential_influence(grid, E, nu, np.where(valid_up, up, C), C)
    Ux2[~valid_up] = 0; Uy2[~valid_up] = 0; Ux2y[~valid_up] = 0
    m = len(C)
    K = np.empty((2 * m, 2 * m))
    K[:m, :m] = (Uxx - Ux2) / grid.dx
    K[:m, m:] = (Uxy - Ux2y) / grid.dx
    K[m:, :m] = (Uxy - Ux2y) / grid.dx
    K[m:, m:] = (Uyy - Uy2) / grid.dx
    return K


def solve_tangential(grid, p, E, nu, mu, xi, eta, phi, K=None, max_iter=3000, relax=0.7, q_tol=1e-6):
    """Active-set solution: adhesion elements satisfy s = 0 exactly (direct
    linear solve); slip elements carry q = -mu p e with e the slip direction,
    updated from the computed slip until the stick/slip partition and the
    directions stop changing."""
    C = np.nonzero(p > 0)[0]
    m = len(C)
    if K is None:
        K = _kmatrix(grid, E, nu, C)
    w = np.concatenate([xi - phi * grid.Y[C], eta + phi * grid.X[C]])
    g = mu * p[C]
    slip = np.zeros(m, bool)
    ex = np.zeros(m); ey = np.zeros(m)
    q = np.zeros(2 * m)
    q_prev = None
    for it in range(max_iter):
        A = np.nonzero(~slip)[0]; S = np.nonzero(slip)[0]
        qS = np.concatenate([-g[S] * ex[S], -g[S] * ey[S]])
        idxA = np.concatenate([A, A + m]); idxS = np.concatenate([S, S + m])
        q = np.zeros(2 * m)
        q[idxS] = qS
        if len(A):
            rhs = -(w[idxA] + K[np.ix_(idxA, idxS)] @ qS)
            q[idxA] = np.linalg.solve(K[np.ix_(idxA, idxA)], rhs)
        s = w + K @ q
        qx, qy = q[:m], q[m:]
        sx, sy = s[:m], s[m:]
        changed = False
        tq = np.hypot(qx, qy)
        over = (~slip) & (tq > g * (1 + 1e-9))
        if over.any():
            ratio = np.where(over, tq / g, 0)
            k = np.argsort(-ratio)[: max(1, int(over.sum() * 0.5))]
            k = k[over[k]]
            slip[k] = True
            ex[k] = -qx[k] / tq[k]; ey[k] = -qy[k] / tq[k]
            changed = True
        Sset = np.nonzero(slip)[0]
        if len(Sset):
            # local stick test: traction that would cancel the slip here, all else fixed
            a = K[Sset, Sset]; b = K[Sset, Sset + m]; c = K[Sset + m, Sset]; d = K[Sset + m, Sset + m]
            rx = sx[Sset] - a * qx[Sset] - b * qy[Sset]
            ry = sy[Sset] - c * qx[Sset] - d * qy[Sset]
            det = a * d - b * c
            qsx = -(d * rx - b * ry) / det
            qsy = -(a * ry - c * rx) / det
            stick = np.hypot(qsx, qsy) < g[Sset] * (1 - 1e-6)
            if stick.any() and not over.any():
                slip[Sset[stick]] = False
                changed = True
            keep = Sset[~stick]
            ss = np.hypot(sx[keep], sy[keep])
            ok = ss > 0
            keep = keep[ok]; ss = ss[ok]
            tx = sx[keep] / ss; ty = sy[keep] / ss
            rk = relax if it < 50 else relax * 50.0 / it   # damp late oscillations
            nx_ = ex[keep] + rk * (tx - ex[keep]); ny_ = ey[keep] + rk * (ty - ey[keep])
            nn = np.hypot(nx_, ny_)
            nx_ /= nn; ny_ /= nn
            dmax = np.max(np.hypot(tx - ex[keep], ty - ey[keep])) if len(keep) else 0.0
            ex[keep] = nx_; ey[keep] = ny_
        dq = np.max(np.abs(q - q_prev)) if q_prev is not None else np.inf
        q_prev = q.copy()
        if not changed and dq > q_tol * g.max():
            changed = True
        if not changed:
            break
    dA = grid.dA
    return dict(Fx=qx.sum() * dA, Fy=qy.sum() * dA, qx=qx, qy=qy, C=C, iterations=it + 1,
                adhesion_fraction=float((~slip).sum()) / m, converged=not changed, K=K)


def _proj_disk(qx, qy, g):
    t = np.hypot(qx, qy)
    f = np.where(t > g, g / np.maximum(t, 1e-300), 1.0)
    return qx * f, qy * f


def solve_tangential_rolling(grid, p, E, nu, mu, xi, eta, phi, tol=1e-9, max_steps=2000, inner_tol=1e-10):
    """Steady rolling as the fixed point of Kalker's transient scheme.

    Each rolling step dq = dx (material advances one element) is the convex
    problem  min 1/2 q'Uq + q'b  s.t. |q_I| <= mu p_I,  b = w dx - u_prev(I'),
    whose KKT conditions are exactly stick (S = 0) / slip (q = -mu p S/|S|)
    for the step slip S = w dx + u(q) - u_prev(I'). U is symmetric positive
    definite, so the inner problem is solved reliably (FISTA with disk
    projection). Steps are repeated until the tractions stop changing.
    """
    C = np.nonzero(p > 0)[0]
    m = len(C)
    ix, iy = C // grid.ny, C % grid.ny
    up = (ix + 1) * grid.ny + iy
    valid_up = ix + 1 < grid.nx
    Uxx, Uyy, Uxy = tangential_influence(grid, E, nu, C, C)
    U = np.block([[Uxx, Uxy], [Uxy, Uyy]])
    rows_up = np.where(valid_up, up, C)
    Vxx, Vyy, Vxy = tangential_influence(grid, E, nu, rows_up, C)
    Vxx[~valid_up] = 0; Vyy[~valid_up] = 0; Vxy[~valid_up] = 0
    V = np.block([[Vxx, Vxy], [Vxy.copy(), Vyy]])   # displacement at upstream points
    # Uxy is symmetric for equal rectangles; V's off-diagonal blocks use the same kernel
    Lip = np.linalg.eigvalsh(U)[-1]
    g = mu * p[C]
    w = np.concatenate([xi - phi * grid.Y[C], eta + phi * grid.X[C]]) * grid.dx
    q = np.zeros(2 * m)
    for step in range(max_steps):
        b = w - V @ q
        # FISTA on min 1/2 q'Uq + b'q over disks
        z = q.copy(); y = q.copy(); tk = 1.0
        for _ in range(5000):
            grad = U @ y + b
            zx, zy = _proj_disk(y[:m] - grad[:m] / Lip, y[m:] - grad[m:] / Lip, g)
            znew = np.concatenate([zx, zy])
            tn = 0.5 * (1 + math.sqrt(1 + 4 * tk * tk))
            y = znew + ((tk - 1) / tn) * (znew - z)
            if np.max(np.abs(znew - z)) <= inner_tol * g.max():
                z = znew
                break
            z, tk = znew, tn
        dq = np.max(np.abs(z - q))
        q = z
        if dq <= tol * g.max() and step > 2:
            break
    qx, qy = q[:m], q[m:]
    S = w + U @ q - V @ q
    adh = np.hypot(qx, qy) < g * (1 - 1e-6)
    return dict(Fx=qx.sum() * grid.dA, Fy=qy.sum() * grid.dA, qx=qx, qy=qy, C=C, steps=step + 1,
                adhesion_fraction=float(adh.sum()) / m, converged=dq <= tol * g.max())
