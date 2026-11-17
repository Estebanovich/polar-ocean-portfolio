"""analysis.heat_transport: J, the onshore heat transport across the J line (plan Task 3),
computed exactly as code_tap/cost_accumulate_mean.F does."""
import numpy as np
import pytest

from exp5_adjoint_slope.analysis import heat_transport as ht

NT, NR, NY, NX = 3, 4, 6, 5
RHO0, CP = 1027.5, 3994.0
DRF = np.array([10.0, 20.0, 30.0, 40.0])
DXG = 10e3
J = 2


def _hfac():
    return np.ones((NR, NY, NX))


def test_uniform_southward_flow_is_positive_and_exact():
    v = np.full((NT, NR, NY, NX), -0.01)
    th = np.zeros((NT, NR, NY, NX))
    got = ht.J(v, th, _hfac(), DRF, DXG, J, theta_f=-1.9, rho0=RHO0, cp=CP)
    area = NX * DXG * DRF.sum()
    want = RHO0 * CP * 0.01 * 1.9 * area * 1e-12
    assert got > 0
    assert got == pytest.approx(want, rel=1e-12)


def test_zero_flow_and_no_gm_flux():
    v = np.zeros((NT, NR, NY, NX))
    th = np.random.default_rng(0).normal(size=v.shape)
    assert ht.J(v, th, _hfac(), DRF, DXG, J, rho0=RHO0, cp=CP) == 0.0
    v = np.random.default_rng(1).normal(size=v.shape) * 0.01
    a = ht.J(v, th, _hfac(), DRF, DXG, J, rho0=RHO0, cp=CP)
    b = ht.J(v, th, _hfac(), DRF, DXG, J, gm_flux=np.zeros((NT, NR, NX)), rho0=RHO0, cp=CP)
    assert a == b


def test_theta_is_averaged_across_the_face_and_land_masked():
    v = np.full((1, NR, NY, NX), -0.01)
    th = np.zeros((1, NR, NY, NX))
    th[:, :, J] = 2.0          # cell north of the face
    th[:, :, J - 1] = 0.0      # cell south of the face
    hf = _hfac()
    hf[3, J] = 0.0             # bottom level of the face is land
    got = ht.J(v, th, hf, DRF, DXG, J, theta_f=-1.9, rho0=RHO0, cp=CP)
    want = RHO0 * CP * 0.01 * (1.0 + 1.9) * NX * DXG * DRF[:3].sum() * 1e-12
    assert got == pytest.approx(want, rel=1e-12)


DYC = 10e3


def test_gm_yflux_is_downgradient_kvy_term():
    # pkg/gmredi/gmredi_ytransport.F:128-134: df = -yA Kvy (T_j - T_j-1) / dyC, with
    # yA = dxG drF hFacS (pkg/generic_advdiff/gad_calc_rhs.F); units degC m3 s-1.
    th = np.zeros((1, NR, NY, NX))
    th[:, :, J] = 1.0                       # warmer north of the face
    kvy = np.full_like(th, 500.0)
    f = ht.gm_yflux(th, kvy, _hfac(), DRF, DXG, DYC, J)
    assert f.shape == (1, NR, NX)
    want = -DXG * DRF[:, None] * 500.0 * 1.0 / DYC
    assert np.allclose(f[0], want, rtol=1e-12)
    assert (f < 0).all()                    # southward: onshore


def test_gm_yflux_kvz_term():
    # gmredi_ytransport.F:154-186: df -= yA Kvz dTdz, dTdz the 4-point mean of the two
    # vertical differences at levels k and k+1 of cells j-1 and j (z positive up: T(k-1) - T(k)).
    th = np.zeros((1, NR, NY, NX))
    th[0] = (-np.arange(NR) * 1.0)[:, None, None]     # 1 degC colder per level downward
    drc = np.array([5.0, 15.0, 25.0, 35.0])
    kvz = np.full_like(th, 2.0)
    mask_c = np.ones((NR, NY, NX))
    f = ht.gm_yflux(th, np.zeros_like(th), _hfac(), DRF, DXG, DYC, J, kvz=kvz, drc=drc,
                    mask_c=mask_c)
    k = 1   # interior level: dTdz = 0.5 (1/drC(k) + 1/drC(k+1)) * 1 degC
    dtdz = 0.5 * (1.0 / drc[k] + 1.0 / drc[k + 1])
    assert f[0, k, 0] == pytest.approx(-DXG * DRF[k] * 2.0 * dtdz, rel=1e-12)
    # Top level: T(k-1) - T(k) with k-1 = k gives 0; only the lower half counts.
    assert f[0, 0, 0] == pytest.approx(-DXG * DRF[0] * 2.0 * 0.5 / drc[1], rel=1e-12)
    # Bottom level: no level below (maskp1 = 0), only the upper half.
    assert f[0, NR - 1, 0] == pytest.approx(-DXG * DRF[-1] * 2.0 * 0.5 / drc[-1], rel=1e-12)


def test_gm_flux_enters_J_onshore_positive():
    v = np.zeros((1, NR, NY, NX))
    th = np.zeros_like(v)
    gm = np.full((1, NR, NX), -3.0)         # degC m3 s-1, southward
    got = ht.J(v, th, _hfac(), DRF, DXG, J, gm_flux=gm, rho0=RHO0, cp=CP)
    assert got == pytest.approx(RHO0 * CP * 3.0 * NR * NX * 1e-12, rel=1e-12)


def test_gm_flux_land_masked():
    v = np.zeros((1, NR, NY, NX))
    hf = _hfac()
    hf[3, J] = 0.0
    th = np.zeros_like(v)
    th[:, :, J] = 1.0
    f = ht.gm_yflux(th, np.full_like(th, 500.0), hf, DRF, DXG, DYC, J)
    assert (f[0, 3] == 0).all()


def test_bolus_transport_from_psi_y():
    # pkg/gmredi/gmredi_residual_flow.F:76-81 (z-coordinates, flipSign4LHCoord = +1):
    # v* hFacS drF = PsiY(k+1) - PsiY(k), PsiY(Nr+1) = 0.
    psi = np.zeros((1, NR, NY, NX))
    psi[0, 1:, J] = 100.0
    vb = ht.bolus_transport(psi)
    assert vb[0, 0, J, 0] == pytest.approx(100.0)
    assert vb[0, 1, J, 0] == pytest.approx(0.0)
    assert vb[0, NR - 1, J, 0] == pytest.approx(-100.0)
    assert vb[0].sum(axis=0)[J, 0] == pytest.approx(0.0)


def test_bolus_heat_transport_sign():
    # Northward bolus flow in the warm deep cell, southward in the cold top: onshore heat < 0.
    v = np.zeros((1, NR, NY, NX))
    th = np.zeros_like(v)
    th[:, 0] = -1.9
    th[:, -1] = 1.0
    psi = np.zeros_like(v)
    psi[0, 1:, J] = -50.0
    got = ht.J(v, th, _hfac(), DRF, DXG, J, psi_y=psi, rho0=RHO0, cp=CP)
    want = -RHO0 * CP * (50.0 * 2.9) * NX * DXG * 1e-12
    assert got == pytest.approx(want, rel=1e-12)


def test_series_mean_equals_J():
    rng = np.random.default_rng(2)
    v = rng.normal(size=(NT, NR, NY, NX)) * 0.01
    th = rng.normal(size=v.shape)
    s = ht.J_series(v, th, _hfac(), DRF, DXG, J, rho0=RHO0, cp=CP)
    assert s.shape == (NT,)
    assert s.mean() == pytest.approx(ht.J(v, th, _hfac(), DRF, DXG, J, rho0=RHO0, cp=CP))


def test_window_iterations():
    # pkg/cost/cost_tile.F accumulates when myTime > endTime - lastinterval; at that call
    # myTime = n*deltaT is the end of step n (model/src/forward_step.F:801-802).
    from exp5_adjoint_slope.analysis.anchor import window_iters
    assert window_iters(n_total=144, last_steps=72) == list(range(73, 145))
    assert window_iters(n_total=10, last_steps=10) == list(range(1, 11))


def test_J_from_means_equals_J_of_products():
    # VVELTH is the time mean of vVel*(theta_j + theta_j-1)/2, GM_vbT of dxG*(PsiY(k+1)-PsiY(k))
    # *(theta_j + theta_j-1)/2 (gmredi_ytransport.F:246-255) and DFyE_TH of the Redi flux, so J
    # from the means equals the mean of the per-step J.
    rng = np.random.default_rng(3)
    v = rng.normal(size=(NT, NR, NY, NX)) * 0.01
    th = rng.normal(size=v.shape)
    psi = rng.normal(size=v.shape)
    psi[:, 0] = 0.0   # PsiY(k=1) = 0 (gmredi_calc_tensor.F:212-216, calc_psi_b.F:90): no net bolus
    gm = rng.normal(size=(NT, NR, NX))
    dfy = np.zeros((NT, NR, NY, NX))
    dfy[:, :, J, :] = gm
    tb = 0.5 * (th + np.roll(th, 1, axis=2))
    vth = (v * tb).mean(axis=0)
    vbt = (DXG * ht.bolus_transport(psi) * tb).mean(axis=0)
    got = ht.J_from_means(vth, v.mean(axis=0), _hfac(), DRF, DXG, J, dfy=dfy.mean(axis=0),
                          vbt=vbt, rho0=RHO0, cp=CP)
    want = ht.J(v, th, _hfac(), DRF, DXG, J, gm_flux=gm, psi_y=psi, rho0=RHO0, cp=CP)
    assert got == pytest.approx(want, rel=1e-12)


def test_rel_diff_with_extrapolated_last_record():
    # The model writes no snapshot after the last step, so the anchor has n-1 of the n records the
    # Fortran cost averages; the missing one is extrapolated linearly from the last two.
    from exp5_adjoint_slope.analysis import anchor as an
    full = np.linspace(1.0, 1.1, 72)
    fc = full.mean()
    assert an.rel_diff_extrapolated(full[:-1], fc, 72) == pytest.approx(0.0, abs=1e-14)
    assert an.rel_diff_extrapolated(full, fc, 72) == pytest.approx(0.0, abs=1e-14)
