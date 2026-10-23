#!/usr/bin/env python3
"""Independent check of the surface-layer solve in src/wind_energy/SurfaceStability.H.

The Monin-Obukhov equations for one cell are

    u_* = kappa U / (ln(z/z0) - psi_m(z/L))
    L   = -u_*^3 theta / (kappa g q)
    q   = kappa u_* (T_s - theta) / (ln(z/z0) - psi_h(z/L))   (given T_s)
    T_s = theta + q (ln(z/z0) - psi_h(z/L)) / (kappa u_*)       (given q)

with the Dyer (1974) stability functions. This script solves them in two
ways:

* ``root``: a bracketed root of the residual in zeta = z/L (scipy brentq),
  written from the equations above and not from the C++ code;
* ``bisection``: a line-by-line port of surface_stability::solve, which
  bisects in zeta = z/L over [max(-5, zeta_h), min(2, zeta_peak)] and holds
  zeta at the bound when there is no root in that range.

It sweeps wind speed (0-20 m/s), reference height (2-30 m), roughness
(0.01-1 m), heat flux and surface temperature, reports the largest
difference from the root where a root exists inside the bounds, checks that
every case is finite with T_s bounded, and prints the pinned values used in
unit_tests/wind_energy/test_surface_stability.cpp.

Usage: python3 tools/surface_stability_sweep.py
"""

import math

import numpy as np
from scipy.optimize import brentq

KAPPA = 0.41
G = 9.81
BETA_M = BETA_H = 16.0
GAMMA_M = GAMMA_H = 5.0
EPS = np.finfo(float).eps
LARGE = np.finfo(float).max * 1.0e-10


def psi_m(zeta):
    if zeta > 0:
        return -GAMMA_M * zeta
    x = (1.0 - BETA_M * zeta) ** 0.25
    return (
        2.0 * math.log(0.5 * (1.0 + x))
        + math.log(0.5 * (1.0 + x * x))
        - 2.0 * math.atan(x)
        + 0.5 * math.pi
    )


def psi_h(zeta):
    if zeta > 0:
        return -GAMMA_H * zeta
    x = math.sqrt(1.0 - BETA_H * zeta)
    return 2.0 * math.log(0.5 * (1.0 + x))


def root(mode, value, wspd, theta, zref, z0):
    """Solve for zeta = zref/L with a bracketed root finder."""
    lr = math.log(zref / z0)

    def state(zeta):
        ustar = KAPPA * wspd / (lr - psi_m(zeta))
        if mode == "heat_flux":
            q = value
        else:
            q = KAPPA * ustar * (value - theta) / (lr - psi_h(zeta))
        return ustar, q

    def residual(zeta):
        ustar, q = state(zeta)
        # zeta L = zref  with  L = -u*^3 theta / (kappa g q)
        return zeta * (-(ustar**3) * theta) - zref * KAPPA * G * q

    if (mode == "heat_flux" and value == 0.0) or (
        mode != "heat_flux" and value == theta
    ):
        zeta = 0.0
    else:
        # Scan the physical branch (both log-law denominators positive) for
        # a sign change, then refine. Stable surfaces (q < 0) have zeta > 0
        # and unstable surfaces zeta < 0.
        stable = (value < 0.0) if mode == "heat_flux" else (value < theta)
        if stable:
            grid = np.logspace(-10, 3, 1300)
        else:
            grid = -np.logspace(-10, 3, 1300)
        grid = np.array(
            [z for z in grid if lr - psi_m(z) > 0.0 and lr - psi_h(z) > 0.0]
        )
        vals = [residual(z) for z in grid]
        zeta = None
        for a, b, fa, fb in zip(grid[:-1], grid[1:], vals[:-1], vals[1:]):
            if fa == 0.0:
                zeta = a
                break
            if fa * fb < 0:
                zeta = brentq(residual, a, b, xtol=1e-14, rtol=1e-14)
                break
        if zeta is None:
            return None
    ustar, q = state(zeta)
    lr_h = lr - psi_h(zeta)
    ts = theta + q * lr_h / (KAPPA * ustar) if mode == "heat_flux" else value
    obl = zref / zeta if zeta != 0.0 else LARGE
    return dict(utau=ustar, L=obl, q=q, Ts=ts, iters=0)


ZETA_MIN = -5.0
ZETA_MAX = 2.0


def bisection(mode, value, wspd, theta, zref, z0):
    """Port of surface_stability::solve for the iterated inputs."""
    lr = math.log(zref / min(z0, zref * math.exp(-1.0)))
    speed = max(wspd, 1.0e-6)
    is_flux = mode == "heat_flux"
    forcing = -value if is_flux else (theta - value)
    zeta, limited = 0.0, False
    if abs(forcing) > EPS:
        c_flux = zref * KAPPA * G * value / ((KAPPA * speed) ** 3 * theta) if is_flux else 0.0
        rib = 0.0 if is_flux else G * zref * (theta - value) / (theta * speed**2)

        def residual(z):
            dm = lr - psi_m(z)
            if is_flux:
                return z + c_flux * dm**3
            return z - rib * dm * dm / (lr - psi_h(z))

        if forcing > 0.0:
            lo, hi = 0.0, ZETA_MAX
            if is_flux:
                hi = min(hi, (1.0 / math.sqrt(-3.0 * c_flux * GAMMA_M) - lr) / GAMMA_M)
            if hi <= 0.0 or residual(hi) < 0.0:
                zeta, limited = ZETA_MAX, True
        else:
            x = 2.0 * math.exp(0.25 * lr) - 1.0
            lo, hi = max(ZETA_MIN, (1.0 - x * x) / BETA_H), 0.0
            if residual(lo) > 0.0:
                zeta, limited = lo, True
        if not limited:
            r_lo = residual(lo)
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                if (residual(mid) > 0.0) == (r_lo > 0.0):
                    lo = mid
                else:
                    hi = mid
            zeta = 0.5 * (lo + hi)
    utau = KAPPA * speed / (lr - psi_m(zeta))
    obl = zref / zeta if zeta != 0.0 else LARGE
    if is_flux:
        q = -zeta * utau**3 * theta / (zref * KAPPA * G) if limited else value
        ts = theta + q * (lr - psi_h(zeta)) / (utau * KAPPA)
    else:
        ts = value
        q = -(theta - ts) * utau * KAPPA / (lr - psi_h(zeta))
    return dict(utau=utau, L=obl, q=q, Ts=ts, zeta=zeta, limited=limited)


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-12)


def sweep():
    """Compare with the bracketed root over a wide grid of conditions."""
    theta = 300.0
    n = n_root = n_limited = 0
    worst = 0.0
    worst_case = None
    max_dt = 0.0
    nonfinite = 0
    for zref in (2.0, 5.0, 18.75, 30.0):
        for z0 in (0.01, 0.1, 0.5, 1.0):
            for wspd in (0.0, 0.5, 1.0, 3.0, 8.0, 20.0):
                cases = [("heat_flux", q) for q in
                         (-0.1, -0.05, -0.02, -0.005, 0.0, 0.01, 0.05, 0.2, 0.5)]
                cases += [("surface_temperature", theta + d) for d in
                          (-10.0, -5.0, -2.0, -0.5, 0.0, 0.5, 2.0, 5.0, 10.0)]
                for mode, value in cases:
                    n += 1
                    b = bisection(mode, value, wspd, theta, zref, z0)
                    vals = (b["utau"], b["q"], b["Ts"], b["L"])
                    if not all(math.isfinite(v) for v in vals) or b["utau"] <= 0.0:
                        nonfinite += 1
                        print("NON-FINITE OR NON-POSITIVE u*:", mode, value, wspd, zref, z0, b)
                        continue
                    max_dt = max(max_dt, abs(b["Ts"] - theta))
                    if b["limited"]:
                        n_limited += 1
                        continue
                    if wspd == 0.0 or z0 > zref * math.exp(-1.0):
                        continue
                    r = root(mode, value, wspd, theta, zref, z0)
                    if r is None:
                        continue
                    n_root += 1
                    du = rel(b["utau"], r["utau"])
                    if du > worst:
                        worst, worst_case = du, (mode, value, wspd, zref, z0)
    print(f"cases: {n}, compared with the bracketed root: {n_root}, "
          f"held at a zeta bound: {n_limited}, non-finite: {nonfinite}")
    print(f"largest relative u* difference from the root: {worst:.3e} at {worst_case}")
    print(f"largest |T_s - theta| over all cases (K): {max_dt:.3f}")


def pinned():
    print("\nPinned values for test_surface_stability.cpp")
    for z in (-1.0, -0.1, 0.0, 0.5):
        print(f"psi_m({z}) = {psi_m(z):.15g}   psi_h({z}) = {psi_h(z):.15g}")
    for mode, value in (("heat_flux", 0.05), ("heat_flux", -0.01), ("surface_temperature", 302.0)):
        r = root(mode, value, 8.0, 300.0, 15.0, 0.1)
        print(
            f"{mode} {value}: utau={r['utau']:.12g} L={r['L']:.12g} "
            f"q={r['q']:.12g} Ts={r['Ts']:.12g}"
        )


if __name__ == "__main__":
    sweep()
    pinned()
