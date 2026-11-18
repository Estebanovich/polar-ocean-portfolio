"""Carbonate chemistry of pkg/dic, ported for the H6b temperature term (reported, not in the sum).

Coefficients: pkg/dic/carbon_chem.F, CARBON_COEFFS (lines 481-777, OCMIP-2 CO2CALC), without
CARBONCHEM_TOTALPHSCALE (undefined in our DIC_OPTIONS.h). Solver: CALC_PCO2_APPROX (lines
329-480, Follows et al. 2006), iterated to convergence instead of one step per time step.
Units: dic, alk, po4, sit in mol m-3 (converted with permil = 1/1024.5, carbon_chem.F); t degC;
s g/kg; returns pCO2 in atm."""
from __future__ import annotations

import numpy as np

PERMIL = 1.0 / 1024.5


def coeffs(t, s) -> dict:
    t, s = np.asarray(t, float), np.asarray(s, float)
    tk = 273.15 + t
    tk100 = tk / 100.0
    tk1002 = tk100 * tk100
    invtk = 1.0 / tk
    dlogtk = np.log(tk)
    is_ = 19.924 * s / (1000.0 - 1.005 * s)
    is2, sqrtis, s2, sqrts, s15, scl = is_ * is_, np.sqrt(is_), s * s, np.sqrt(s), s ** 1.5, s / 1.80655
    P1atm, Rgas = 1.01325, 83.1451
    RT = Rgas * tk
    delta = 57.7 - 0.118 * tk
    B = -1636.75 + 12.0408 * tk - 0.0327957 * tk * tk + 3.16528 * tk * tk * tk * 0.00001
    c = {}
    c["fugf"] = np.exp((B + 2.0 * delta) * P1atm / RT)
    c["ff"] = np.exp(-162.8301 + 218.2968 / tk100 + 90.9241 * np.log(tk100) - 1.47696 * tk1002
                     + s * (.025695 - .025225 * tk100 + 0.0049867 * tk1002))
    c["k0"] = np.exp(93.4517 / tk100 - 60.2409 + 23.3585 * np.log(tk100)
                     + s * (0.023517 - 0.023656 * tk100 + 0.0047036 * tk1002))
    c["k1"] = 10.0 ** (-1.0 * (3670.7 * invtk - 62.008 + 9.7944 * dlogtk - 0.0118 * s + 0.000116 * s2))
    c["k2"] = 10.0 ** (-1.0 * (1394.7 * invtk + 4.777 - 0.0184 * s + 0.000118 * s2))
    c["kb"] = np.exp((-8966.90 - 2890.53 * sqrts - 77.942 * s + 1.728 * s15 - 0.0996 * s2) * invtk
                     + (148.0248 + 137.1942 * sqrts + 1.62142 * s)
                     + (-24.4344 - 25.085 * sqrts - 0.2474 * s) * dlogtk + 0.053105 * sqrts * tk)
    c["k1p"] = np.exp(-4576.752 * invtk + 115.525 - 18.453 * dlogtk + (-106.736 * invtk + 0.69171) * sqrts
                      + (-0.65643 * invtk - 0.01844) * s)
    c["k2p"] = np.exp(-8814.715 * invtk + 172.0883 - 27.927 * dlogtk + (-160.340 * invtk + 1.3566) * sqrts
                      + (0.37335 * invtk - 0.05778) * s)
    c["k3p"] = np.exp(-3070.75 * invtk - 18.141 + (17.27039 * invtk + 2.81197) * sqrts
                      + (-44.99486 * invtk - 0.09984) * s)
    c["ksi"] = np.exp(-8904.2 * invtk + 117.385 - 19.334 * dlogtk + (-458.79 * invtk + 3.5913) * sqrtis
                      + (188.74 * invtk - 1.5998) * is_ + (-12.1652 * invtk + 0.07871) * is2
                      + np.log(1.0 - 0.001005 * s))
    c["kw"] = np.exp(-13847.26 * invtk + 148.9652 - 23.6521 * dlogtk
                     + (118.67 * invtk - 5.977 + 1.0495 * dlogtk) * sqrts - 0.01615 * s)
    c["bt"] = 0.000232 * scl / 10.811
    return c


def pco2(dic, alk, t, s=34.5, po4=0.0, sit=0.0, ph0=8.0, n_iter=100):
    """(pCO2 atm, pH) from CALC_PCO2_APPROX iterated to convergence."""
    c = coeffs(t, s)
    dic, ta, pt, si = (np.asarray(x, float) * PERMIL for x in (dic, alk, po4, sit))
    ph = np.full(np.broadcast(dic, ta).shape, ph0, float)
    for _ in range(n_iter):
        h = 10.0 ** (-ph)
        bohg = c["bt"] * c["kb"] / (h + c["kb"])
        stuff = h ** 3 + c["k1p"] * h * h + c["k1p"] * c["k2p"] * h + c["k1p"] * c["k2p"] * c["k3p"]
        h3po4g = pt * h ** 3 / stuff
        hpo4g = pt * c["k1p"] * c["k2p"] * h / stuff
        po4g = pt * c["k1p"] * c["k2p"] * c["k3p"] / stuff
        siooh3g = si * c["ksi"] / (c["ksi"] + h)
        cag = ta - bohg - c["kw"] / h + h - hpo4g - 2.0 * po4g + h3po4g - siooh3g
        gamm = dic / cag
        disc = (1.0 - gamm) ** 2 * c["k1"] ** 2 - 4.0 * c["k1"] * c["k2"] * (1.0 - 2.0 * gamm)
        hnew = 0.5 * ((gamm - 1.0) * c["k1"] + np.sqrt(disc))
        ph_new = -np.log10(hnew)
        if np.allclose(ph_new, ph, atol=1e-10, rtol=0):
            ph = ph_new
            break
        ph = ph_new
    h = 10.0 ** (-ph)
    co2s = dic / (1.0 + c["k1"] / h + c["k1"] * c["k2"] / (h * h))
    return co2s / c["k0"] / c["fugf"], ph


def dic_at_pco2(pco2_atm, alk, t, s=34.5, po4=0.0, sit=0.0, lo=1.5, hi=2.6, n=60) -> float:
    """DIC (mol m-3) in equilibrium with pco2_atm at (alk, t, s): bisection on pco2()."""
    for _ in range(n):
        mid = 0.5 * (lo + hi)
        p, _ = pco2(mid, alk, t, s, po4, sit)
        lo, hi = (mid, hi) if p < pco2_atm else (lo, mid)
    return 0.5 * (lo + hi)
