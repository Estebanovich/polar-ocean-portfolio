"""Offline linear decomposition of Jarnikova et al. (2025, Sci. Adv. 11, eadr3589, Methods,
eqs. 1-5), as annual rates over the analysis region R (spec Sec. 7, H6b; Esteban: no iron, PO4
in place of Fe; the rate form is the design's proposal).

Inputs per state (parent = reference, branch = perturbed), units:
  mld_inv m (winter mixed-layer depth over R), area_R m2, g200 umol kg-1 m-1 (eq. 1 gradient),
  ddic1000 umol kg-1 (eq. 4), dpo4dz200 umol kg-1 m-1, dpo4_1000 umol kg-1, moc_sv Sv (upper cell
  upwelling in R), r_cp mol C/mol P (pkg/dic R_CP), e = export / PO4 supply in the parent.
Terms (mol C yr-1, + = DIC added to the surface layer of R):
  mld = dMLD g200 rho0 A                        (eq. 1, one winter entrainment per year)
  moc = dMOC dDIC1000 rho0 yr                   (eq. 4)
  bio = -R_CP e (dMLD dPO4/dz rho0 A + dMOC dPO4_1000 rho0 yr)   (eqs. 2, 3, 5 with PO4 for Fe)
Gradients, MLD, MOC and e come from the parent; dMLD and dMOC are branch - parent."""
from __future__ import annotations

import numpy as np

from . import carbonate

RHO0 = 1035.0
YEAR_S = 31536000.0
UMOL = 1e-6
SV = 1e6


def terms(parent: dict, branch: dict) -> dict:
    p = parent
    d_mld = branch["mld_inv"] - p["mld_inv"]
    d_moc = (branch["moc_sv"] - p["moc_sv"]) * SV
    mld = d_mld * p["g200"] * RHO0 * p["area_R"] * UMOL
    moc = d_moc * p["ddic1000"] * RHO0 * UMOL * YEAR_S
    supply = d_mld * p["dpo4dz200"] * RHO0 * p["area_R"] * UMOL + d_moc * p["dpo4_1000"] * RHO0 * UMOL * YEAR_S
    bio = -p["r_cp"] * p["e"] * supply
    return {"mld": float(mld), "moc": float(moc), "bio": float(bio), "sum_circ_bio": float(mld + moc + bio)}


# ---- temperature term (spec Sec. 7: acts through gas exchange; reported, not in sum_circ_bio) ----

def conductance(flux, dpco2, area) -> float:
    """Effective gas-exchange conductance of R (mol m-2 yr-1 atm-1): sum(A F) / sum(A dpCO2), with
    F the air-sea flux (+ into the ocean) and dpCO2 = pCO2atm - pCO2ocean (atm), both summed over
    the parent's window records. It stands for Kwexch ff / permil of pkg/dic; the change of ff with
    temperature on the atmospheric side is neglected (it multiplies the disequilibrium, not pCO2)."""
    a = np.asarray(area, float)
    return float((a * np.asarray(flux, float)).sum() / (a * np.asarray(dpco2, float)).sum())


def temperature_term(dic, alk, po4, t_parent, t_branch, area, cond) -> float:
    """Change of the air-sea flux into R (mol C yr-1, + into the ocean) from the branch - parent
    surface temperature alone: -cond sum A [pCO2(DIC, Alk, PO4, T_branch) - pCO2(..., T_parent)],
    with the parent's surface DIC, Alk and PO4 (mol m-3) held (analysis/carbonate.py)."""
    p_b, _ = carbonate.pco2(dic, alk, t_branch, po4=po4)
    p_p, _ = carbonate.pco2(dic, alk, t_parent, po4=po4)
    return float(-cond * (np.asarray(area, float) * (p_b - p_p)).sum())
