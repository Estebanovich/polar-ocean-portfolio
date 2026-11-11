#!/opt/anaconda3/envs/anemoi/bin/python
r"""POST HOC diagnostic of the main ``profile_test`` result of M7 (Exp. 3).

Usage::

    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_NUM_THREADS=1 \
    profile_test_diag.py [--out exp3_melt_emulator/figures/profile_test_diag.json]
                         [--diag-dir /Volumes/LaCie/portfolio_runs/exp3/diag]

Question: is the +413 % integrated error of M7 on ocean0_4km (main profile_test fold,
cv_results.json) an error of ours (inputs built differently for ocean0_4km, TF, scaling,
masks, units) or out-of-distribution behaviour of the frozen method?

Nothing of the protocol is re-run or modified. Steps:

1. Input consistency of ocean0_4km against the generator of Exp. 2 (gendata/profiles.py):
   the dataset T/S profile must equal both ``ocean0`` (AD16 Eq. 21-22, Table 6) and
   ``ramp(z_c = 360, d = 720)``; draft and shelf_mask must equal those of the s = 1.0
   training runs; TF of ``cell_features`` must equal an independent recomputation
   T - (lambda1 S + lambda2 + lambda3 p) (AD16 Table 4, Eq. 25 as data.freezing_temperature).
2. Marginal ranges of FEATURES on ocean0_4km against the 25 training runs.
3. Joint (conditional) range: TF against draft_minus_zc. g_theta sees x_g = (draft,
   draft_minus_zc, ...), never TF, so what matters is whether TF at a given position
   relative to the thermocline centre was seen in training.
4. One refit of M7 on the 25 training runs with the frozen defaults and tf_kind fixed to the
   value the unit selected ("local"); the member fit after the tf_kind selection does not
   depend on the selection fits (each fit forks the RNG), so the refit must reproduce the
   unit's per-member integrated melt exactly (checked). The prediction is decomposed per
   member into P = K TF+^n (Burgard 2022 Eq. 14 form with free n) and g_theta(x_g), and
   the excess integrated melt is attributed to bins of draft_minus_zc.
5. Counterfactual "required g": g_req = melt / P on ocean0_4km, against the learned g.

Integrated melt is the plain sum over the run's shelf_mask (uniform 4 km grid, as
metrics.py; Burgard 2022 Eq. 33-34 up to the constant cell area).
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from exp2_isomip_cavity.gendata.profiles import ocean0, ramp
from exp3_melt_emulator.data import FEATURES, cell_features, freezing_temperature, load_runs
from exp3_melt_emulator.methods.physics import C, QuadLocal

REPO = Path(__file__).resolve().parents[2]
UNIT = Path("/Volumes/LaCie/portfolio_runs/exp3/cv_units/main__profile_test__M7.json")
REF = "ocean0_4km"
# bins of draft_minus_zc (m, positive = deeper than the thermocline centre)
DZ_BINS = [-400, -250, -150, -50, 50, 150, 250, 400]


def _q(a):
    a = np.asarray(a, float)
    return {"min": float(a.min()), "p05": float(np.percentile(a, 5)),
            "median": float(np.median(a)), "p95": float(np.percentile(a, 95)),
            "max": float(a.max())}


def consistency(ref, train):
    """Step 1: ocean0_4km inputs built like the ensemble's."""
    T0, S0 = ocean0(ref.z)                       # AD16 Eq. 21-22 (Table 6)
    Tr, Sr = ramp(ref.z, 1.0, 360.0, 720.0)      # profiles.ramp docstring: z_c=360, d=720
    s1 = [r for r in train if r.s == 1.0]
    f = cell_features(ref)
    sh = ref.shelf_mask
    S = f["S"][sh]
    # independent TF: AD16 Table 4 liquidus via data.freezing_temperature (R2 pressure)
    tf_indep = f["T"][sh] - freezing_temperature(S, ref.draft[sh])
    # T, S at draft by direct evaluation of the analytic profile (no interpolation)
    Ta, Sa = ocean0(ref.draft[sh])
    return {
        "max_abs_T_profile_minus_ocean0": float(np.max(np.abs(ref.T_profile - T0))),
        "max_abs_S_profile_minus_ocean0": float(np.max(np.abs(ref.S_profile - S0))),
        "max_abs_T_profile_minus_ramp_zc360_d720": float(np.max(np.abs(ref.T_profile - Tr))),
        "max_abs_S_profile_minus_ramp_zc360_d720": float(np.max(np.abs(ref.S_profile - Sr))),
        "dataset_params": {"Tbot": ref.Tbot, "z_c": ref.z_c, "d": ref.d, "s": ref.s},
        "draft_equal_to_s1_runs": bool(all(np.array_equal(ref.draft, r.draft, equal_nan=True)
                                           for r in s1)),
        "shelf_mask_equal_to_s1_runs": bool(all(np.array_equal(ref.shelf_mask, r.shelf_mask)
                                                for r in s1)),
        "n_shelf_cells": int(sh.sum()),
        "max_abs_TF_minus_independent": float(np.max(np.abs(f["TF"][sh] - tf_indep))),
        "max_abs_T_interp_minus_analytic": float(np.max(np.abs(f["T"][sh] - Ta))),
        "draft_range_m": [float(ref.draft[sh].min()), float(ref.draft[sh].max())],
        "z_levels_range_m": [float(ref.z.min()), float(ref.z.max())],
        "melt_units_check_mean_m_per_a": float(ref.melt[sh].mean()),
    }


def marginal_ranges(ref, train):
    """Step 2: per-feature range of ocean0_4km vs the pooled training cells."""
    fr = cell_features(ref)
    pooled = {k: np.concatenate([cell_features(r)[k][r.shelf_mask] for r in train])
              for k in FEATURES}
    out = {}
    for k in FEATURES:
        v = fr[k][ref.shelf_mask]
        lo, hi = float(pooled[k].min()), float(pooled[k].max())
        out[k] = {"train_min": lo, "train_max": hi, "ref_min": float(v.min()),
                  "ref_max": float(v.max()),
                  "ref_frac_outside": float(np.mean((v < lo) | (v > hi)))}
    return out


def joint_tf(ref, train):
    """Step 3: TF in bins of draft_minus_zc, training cells vs ocean0_4km."""
    fr = cell_features(ref)
    sh = ref.shelf_mask
    dz_r, tf_r = fr["draft_minus_zc"][sh], fr["TF"][sh]
    dz_t = np.concatenate([cell_features(r)["draft_minus_zc"][r.shelf_mask] for r in train])
    tf_t = np.concatenate([cell_features(r)["TF"][r.shelf_mask] for r in train])
    bins = []
    above = np.zeros(sh.sum(), bool)
    for lo, hi in zip(DZ_BINS[:-1], DZ_BINS[1:]):
        mr = (dz_r >= lo) & (dz_r < hi)
        mt = (dz_t >= lo) & (dz_t < hi)
        rec = {"bin_m": [lo, hi], "n_ref": int(mr.sum()), "n_train": int(mt.sum())}
        if mt.any():
            rec["train_TF"] = _q(tf_t[mt])
        if mr.any():
            rec["ref_TF"] = _q(tf_r[mr])
        if mr.any() and mt.any():
            a = tf_r[mr] > tf_t[mt].max()
            above[np.flatnonzero(mr)[a]] = True
            rec["ref_frac_TF_above_train_max"] = float(a.mean())
        bins.append(rec)
    return {"bins": bins, "ref_frac_cells_TF_above_train_max_in_bin": float(above.mean()),
            "ref_cells_TF_above_train_max_in_bin": int(above.sum())}, above


def refit_and_decompose(ref, train_by_case, unit, diag_dir):
    """Step 4-5: refit M7 (tf_kind = local, frozen defaults) and decompose on ocean0_4km."""
    import torch
    from exp3_melt_emulator.methods.nets import Hybrid, _deterministic

    runs = [train_by_case[c] for c in unit["train"]]      # the unit's order
    kind = unit["info"]["tf_kind"]
    t0 = time.time()
    m = Hybrid(tf_kind=kind)
    m.fit(runs)
    fit_s = time.time() - t0
    cells = m._cells(ref)
    mus, _ = m.member_predictions(cells)
    with _deterministic(m.threads), torch.no_grad():
        g = np.stack([st.net.heads(m._tensors(cells, st)[0])[0].numpy() for st in m.members_])
    Ks = np.array([r["K"] for r in m.info["members"]])
    ns = np.array([r["n"] for r in m.info["members"]])
    b = m.basis(cells, kind)                               # TF+ (local)
    P = Ks[:, None] * b[None, :] ** ns[:, None]           # K TF+^n, m a-1
    recon = float(np.max(np.abs(P * g - mus)))
    y = ref.melt[ref.shelf_mask]

    # g on the training cells (same members), to check whether g on ocean0 is in range
    g_train = []
    for r in runs:
        c = m._cells(r)
        with _deterministic(m.threads), torch.no_grad():
            g_train.append(np.stack([st.net.heads(m._tensors(c, st)[0])[0].numpy()
                                     for st in m.members_]))
    g_train = np.concatenate(g_train, axis=1)

    integ = mus.sum(1)
    exp_integ = np.array(unit["results"][0]["members"]["integrated"])
    mu = mus.mean(0)
    dz = cells["draft_minus_zc"]
    bins = []
    for lo, hi in zip(DZ_BINS[:-1], DZ_BINS[1:]):
        mk = (dz >= lo) & (dz < hi)
        if not mk.any():
            continue
        Pm, gm = P[:, mk].mean(0), g[:, mk].mean(0)
        req = np.where(Pm > 1e-9, y[mk] / np.maximum(Pm, 1e-9), np.nan)
        bins.append({"bin_m": [lo, hi], "n_cells": int(mk.sum()),
                     "sum_true": float(y[mk].sum()), "sum_pred": float(mu[mk].sum()),
                     "excess": float(mu[mk].sum() - y[mk].sum()),
                     "median_TF": float(np.median(cells["TF"][mk])),
                     "median_P_KTFn": float(np.median(Pm)), "median_g": float(np.median(gm)),
                     "median_g_required": float(np.nanmedian(req))})
    excess_total = float(mu.sum() - y.sum())
    np.savez(Path(diag_dir) / "profile_test_M7_decomposition.npz",
             mus=mus, g=g, P=P, b=b, melt=y, draft_minus_zc=dz, TF=cells["TF"],
             draft=cells["draft"], g_train=g_train)
    # counterfactual: the same P with g = 1 (pure K TF+^n of each member)
    return {
        "tf_kind": kind, "fit_s": fit_s,
        "members_K": Ks.tolist(), "members_n": ns.tolist(),
        "reproduces_unit_integrated": bool(np.allclose(integ, exp_integ, rtol=1e-10, atol=0)),
        "max_rel_diff_vs_unit_integrated": float(np.max(np.abs(integ / exp_integ - 1))),
        "max_abs_reconstruction_P_times_g_minus_mu": recon,
        "integrated_true": float(y.sum()), "integrated_pred_mixture": float(mu.sum()),
        "integrated_rel_error": float(mu.sum() / y.sum() - 1),
        "integrated_P_only_g1_members": P.sum(1).tolist(),
        "integrated_rel_error_P_only_g1_mean": float(P.mean(0).sum() / y.sum() - 1),
        "g_ref": _q(g.mean(0)), "g_train": _q(g_train.mean(0)),
        "g_ref_max_member": float(g.max()), "g_train_max_member": float(g_train.max()),
        "bins": bins, "excess_total": excess_total,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--out", default=str(REPO / "exp3_melt_emulator/figures/profile_test_diag.json"))
    ap.add_argument("--diag-dir", default="/Volumes/LaCie/portfolio_runs/exp3/diag")
    ap.add_argument("--no-refit", action="store_true")
    a = ap.parse_args()
    Path(a.diag_dir).mkdir(parents=True, exist_ok=True)
    runs = load_runs()
    by = {r.case: r for r in runs}
    unit = json.loads(UNIT.read_text())
    ref = by[REF]
    train = [by[c] for c in unit["train"]]
    assert REF not in unit["train"] and len(train) == 25

    out = {"label": "POST HOC diagnostic (not pre-registered); protocol results unchanged",
           "unit_record": str(UNIT),
           "unit_metrics": unit["results"][0]["metrics"],
           "unit_spread_rel": unit["results"][0]["members"]["spread_rel"],
           "unit_members_K": [r["K"] for r in unit["info"]["members"]],
           "unit_members_n": [r["n"] for r in unit["info"]["members"]],
           "unit_tf_kind": unit["info"]["tf_kind"],
           "consistency": consistency(ref, train),
           "marginal_ranges": marginal_ranges(ref, train)}
    jt, above = joint_tf(ref, train)
    out["joint_TF_vs_draft_minus_zc"] = jt
    # M2 on the same fold (closed form, deterministic): reproduces the unit's +5 %
    m2 = QuadLocal()
    m2.fit(train)
    p2 = m2.predict(ref)[ref.shelf_mask]
    out["M2_check"] = {"coef": m2.coef, "integrated_rel_error":
                       float(p2.sum() / ref.melt[ref.shelf_mask].sum() - 1),
                       "K_equiv_at_n2": float(m2.coef * C**2)}
    dz = cell_features(ref)["draft_minus_zc"][ref.shelf_mask]
    y = ref.melt[ref.shelf_mask]
    out["M2_check"]["bins"] = [
        {"bin_m": [lo, hi], "sum_true": float(y[k].sum()), "sum_pred": float(p2[k].sum())}
        for lo, hi in zip(DZ_BINS[:-1], DZ_BINS[1:])
        for k in [(dz >= lo) & (dz < hi)] if k.any()]
    if not a.no_refit:
        d = refit_and_decompose(ref, by, unit, a.diag_dir)
        y = ref.melt[ref.shelf_mask]
        dec = np.load(Path(a.diag_dir) / "profile_test_M7_decomposition.npz")
        mu = dec["mus"].mean(0)
        d["excess_in_cells_TF_above_train_max_in_bin"] = float((mu[above] - y[above]).sum())
        d["excess_share_in_cells_TF_above_train_max_in_bin"] = float(
            (mu[above] - y[above]).sum() / (mu.sum() - y.sum()))
        out["M7_refit_decomposition"] = d
    Path(a.out).write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
