"""Beam-width (FWHM) tests: the same checks as notebooks/09_beam_fwhm_tests.

The archive-cut tests need limTOD and the TRIS archive (``paths.tris_archive_dir``)
and are skipped without them.  The smoothing tests need only healpy.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

hp = pytest.importorskip("healpy")

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from moonrunes.stage2_beam_matching import extra_fwhm_deg, smooth_to_target  # noqa: E402

FWHM_E, FWHM_H = 19.155, 23.366
K = 2.0 * np.sqrt(2.0 * np.log(2.0))
HALF_DB = -10.0 * np.log10(2.0)


def _fwhm_from_profile(theta_deg, response):
    r = np.asarray(response) / response[0]
    i = int(np.argmax(r < 0.5))
    return 2 * (theta_deg[i - 1] + (0.5 - r[i - 1]) / (r[i] - r[i - 1]) * (theta_deg[i] - theta_deg[i - 1]))


@pytest.fixture(scope="module")
def tris_and_cuts():
    try:
        from moonrunes.config import load_config
        from moonrunes.stage1_tris_maps import import_limtod_tris

        config = load_config(REPO_ROOT / "configs" / "run_config.yaml")
        tris = import_limtod_tris(config)
        path = config.resolve_path("paths.tris_archive_dir") / "TRIS_Beam_Profile.txt"
        if not path.exists():
            pytest.skip("TRIS beam profile not downloaded")
        return tris, tris.read_tris_beam_cuts(path), config
    except ImportError as err:  # limTOD not installed
        pytest.skip(str(err))


def test_sigma_fwhm_constant():
    assert K == pytest.approx(2.3548200450309493, abs=1e-12)
    for fwhm in (FWHM_E, FWHM_H):
        assert np.exp(-0.5 * (fwhm / 2 / (fwhm / K)) ** 2) == pytest.approx(0.5, abs=1e-12)


def test_gaussian_cuts_give_back_their_fwhm(tris_and_cuts):
    tris, _, _ = tris_and_cuts
    angles = np.arange(0.0, 60.0, 0.01)
    db = lambda s: 10 * np.log10(np.exp(-0.5 * (angles / s) ** 2))
    cuts = tris.TRISPrincipalPlaneCuts(angle_deg=angles, h_plane_db=db(FWHM_H / K), e_plane_db=db(FWHM_E / K))
    e, h = cuts.half_power_full_width_deg()
    assert e == pytest.approx(FWHM_E, abs=1e-5)
    assert h == pytest.approx(FWHM_H, abs=1e-5)


def test_archive_widths_and_config_target(tris_and_cuts):
    _, cuts, config = tris_and_cuts
    e, h = cuts.half_power_full_width_deg()
    assert e == pytest.approx(FWHM_E, abs=5e-4)
    assert h == pytest.approx(FWHM_H, abs=5e-4)
    assert float(config.require("beam_matching.target_fwhm_deg")) == pytest.approx(h, abs=5e-4)


def test_gaussian_beam_map_widths_in_degrees(tris_and_cuts):
    # Catches deg/rad slips (x57) and sigma-for-FWHM mix-ups (x2.355) in the map builder.
    tris, _, _ = tris_and_cuts
    theta = np.arange(0.0, 40.0, 0.02)
    m = tris.approximate_tris_gaussian_beam_map(nside=256, fwhm_e_deg=FWHM_E, fwhm_h_deg=FWHM_H)
    prof = lambda phi: hp.get_interp_val(m, np.radians(theta), np.full(theta.size, np.radians(phi)))
    assert _fwhm_from_profile(theta, prof(0.0)) == pytest.approx(FWHM_E, abs=0.05)
    assert _fwhm_from_profile(theta, prof(90.0)) == pytest.approx(FWHM_H, abs=0.05)


@pytest.mark.parametrize("blend", ["db", "power"])
def test_cut_beam_principal_planes(tris_and_cuts, blend):
    tris, cuts, _ = tris_and_cuts
    a = cuts.angle_deg
    for phi, cut in ((0.0, cuts.e_plane_db), (180.0, cuts.e_plane_db),
                     (90.0, cuts.h_plane_db), (270.0, cuts.h_plane_db)):
        got = tris.tris_cut_beam_response(cuts, a, np.full(a.size, phi), blend=blend)
        np.testing.assert_allclose(got, 10 ** (cut / 10), rtol=1e-12)


def test_cut_beam_measured_fwhm(tris_and_cuts):
    tris, cuts, _ = tris_and_cuts
    ref_e, ref_h = cuts.half_power_full_width_deg()
    theta = np.arange(0.0, 40.0, 0.01)
    resp = lambda phi: tris.tris_cut_beam_response(cuts, theta, np.full(theta.size, phi))
    assert _fwhm_from_profile(theta, resp(0.0)) == pytest.approx(ref_e, abs=0.01)
    assert _fwhm_from_profile(theta, resp(90.0)) == pytest.approx(ref_h, abs=0.01)


def test_extra_fwhm_quadrature():
    assert extra_fwhm_deg(3.0, 5.0) == pytest.approx(4.0, abs=1e-12)
    with pytest.raises(ValueError):
        extra_fwhm_deg(5.0, 3.0)


def test_smoothing_width_and_quadrature_chain():
    # A point source smoothed to the target must have the exact spherical profile of a
    # harmonic Gaussian of that FWHM (23.425 deg real-space for 23.366), and going via 5 deg
    # must give the same answer.
    ns, target = 64, FWHM_H
    src = np.zeros(hp.nside2npix(ns))
    pix = hp.ang2pix(ns, 0.0, 0.0, lonlat=True)
    src[pix] = 1.0
    lon0, lat0 = hp.pix2ang(ns, pix, lonlat=True)
    half = np.arange(0.0, 40.0, 0.02)

    def width(m):
        peak = float(hp.get_interp_val(m, lon0, lat0, lonlat=True))
        up = hp.get_interp_val(m, lon0 + 0 * half, lat0 + half, lonlat=True)
        down = hp.get_interp_val(m, lon0 + 0 * half, lat0 - half, lonlat=True)
        return (_fwhm_from_profile(half, np.r_[peak, up[1:]]) + _fwhm_from_profile(half, np.r_[peak, down[1:]])) / 2

    t = np.arange(0.0, 40.0, 0.01)
    exact = _fwhm_from_profile(t, hp.bl2beam(hp.gauss_beam(np.radians(target), lmax=3 * ns - 1), np.radians(t)))
    direct = width(smooth_to_target(src, 0.0, target))
    chained = width(smooth_to_target(hp.smoothing(src, fwhm=np.radians(5.0), iter=3), 5.0, target))
    assert direct == pytest.approx(exact, abs=0.05)
    assert chained == pytest.approx(direct, abs=0.01)


@pytest.mark.parametrize("normalization", ["none", "peak"])
def test_round_trip_gaussian_map_through_fwhm_function(tris_and_cuts, normalization):
    # Test 2: Gaussian beam map (19.155/23.366) -> E/H cuts read off the map -> the FWHM
    # function -> the same two numbers to 0.1%.  Holds from nside 128 up (notebook 09, sec. 6).
    tris, _, _ = tris_and_cuts
    angles = np.arange(0.0, 60.0, 0.05)
    m = tris.approximate_tris_gaussian_beam_map(nside=128, fwhm_e_deg=FWHM_E, fwhm_h_deg=FWHM_H,
                                                normalization=normalization)
    th = np.radians(angles)
    e = hp.get_interp_val(m, th, np.zeros_like(th))
    h = hp.get_interp_val(m, th, np.full_like(th, np.pi / 2))
    cuts = tris.TRISPrincipalPlaneCuts(angle_deg=angles, e_plane_db=10 * np.log10(e),
                                       h_plane_db=10 * np.log10(h))
    e_fwhm, h_fwhm = cuts.half_power_full_width_deg()
    assert e_fwhm == pytest.approx(FWHM_E, rel=1e-3)
    assert h_fwhm == pytest.approx(FWHM_H, rel=1e-3)
