"""Rail head (with gauge corner) and worn wheel (tread, flange root, flange)
transverse profiles, built from curvature-vs-arc-length laws and sampled as
z(y) on a uniform y grid (mm). y points toward the track centre (gauge side),
z upward. Rail profile: rail-fixed. Wheel profile: wheel-fixed coordinates at
the section through the wheel axis; tape circle at y_w = 0."""
import numpy as np

DS = 0.002  # arc step (mm) for integration


def smooth(t):
    t = np.clip(t, 0, 1)
    return t * t * (3 - 2 * t)


def integrate(kappa_fn, length, theta_stop=None, sign=+1):
    """Integrate dtheta/dsigma = kappa(sigma) from (0,0,theta=0) along +y (sign=+1) or -y (sign=-1)."""
    n = int(length / DS)
    y = [0.0]; z = [0.0]; th = 0.0
    for i in range(n):
        s = (i + 0.5) * DS
        k = kappa_fn(s, th)
        th += k * DS
        y.append(y[-1] + sign * np.cos(th) * DS)
        z.append(z[-1] + sign * np.sin(th) * DS)
    return np.array(y), np.array(z)


def rail_raw():
    def k_gauge(s, th):       # toward +y (gauge side); theta is slope angle d z/d y
        if th <= np.radians(-72):
            return 0.0
        return -(1 / 300 + (1 / 80 - 1 / 300) * smooth((s - 9) / 2) + (1 / 13 - 1 / 80) * smooth((s - 21) / 2))
    def k_field(s, th):       # toward -y: going left, a crown curves down => slope angle increases
        if th >= np.radians(60):
            return 0.0
        return +(1 / 300 + (1 / 80 - 1 / 300) * smooth((s - 9) / 2) + (1 / 13 - 1 / 80) * smooth((s - 23) / 2))
    yg, zg = integrate(k_gauge, 46.0, sign=+1)
    yf, zf = integrate(k_field, 44.0, sign=-1)
    y = np.concatenate([yf[::-1], yg[1:]]); z = np.concatenate([zf[::-1], zg[1:]])
    return y, z


def wheel_raw():
    def k_gauge(s, th):       # toward +y: hollow band, then flange root R15 down to the flank, flange tip
        if s < 9.0:
            return 1 / 360
        if th > np.radians(-68) and s < 40:
            if not hasattr(k_gauge, "flank_start"):
                pass
            return 1 / 360 * (1 - smooth((s - 9) / 1.5)) - (1 / 15) * smooth((s - 9) / 1.5)
        k_gauge.reached = getattr(k_gauge, "reached", s)
        if s < k_gauge.reached + 7.0:
            return 0.0                                  # straight flank at about 68 degrees
        return 1 / 10                                    # convex flange tip, turning back up
    def k_field(s, th):       # toward -y: hollow band then convex tread/chamfer (1/150)
        return -(1 / 360) + (1 / 360 + 1 / 150) * smooth((s - 11) / 3)
    yg, zg = integrate(k_gauge, 48.0, sign=+1)
    yf, zf = integrate(k_field, 38.0, sign=-1)
    keep = np.nonzero(np.diff(yg) > 0)[0]
    last = keep[-1] + 1 if len(keep) == len(yg) - 1 else int(np.argmax(np.diff(yg) <= 0))
    yg, zg = yg[:last + 1], zg[:last + 1]
    y = np.concatenate([yf[::-1], yg[1:]]); z = np.concatenate([zf[::-1], zg[1:]])
    z = z + 0.0012 * y                                             # slight tilt
    return y, z


def sample(y, z, step=0.05):
    yy = np.round(np.arange(np.ceil(y.min() / step) * step, np.floor(y.max() / step) * step + 1e-9, step), 3)
    return yy, np.interp(yy, y, z)


if __name__ == "__main__":
    yr, zr = rail_raw(); yw, zw = wheel_raw()
    print("rail y range", yr.min(), yr.max(), "z min", zr.min())
    print("wheel y range", yw.min(), yw.max(), "z min", zw.min())
    th = np.degrees(np.arctan2(np.diff(zw), np.diff(yw)))
    print("wheel slope range", th.min(), th.max())
