"""Tests for the elliptical TRIS-beam kernel (``moonrunes.elliptical_beam``).

The geometry and Gaussian tests are synthetic and need only healpy.  The cut-beam
tests need limTOD and the TRIS archive and are skipped without them.  The plots
that go with these checks are in ``notebooks/03_diagnostics.ipynb``, section
"The elliptical beam for the map route".
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

hp = pytest.importorskip("healpy")

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from moonrunes import elliptical_beam as eb  # noqa: E402
from moonrunes.config import load_config  # noqa: E402
from moonrunes import stage2_beam_matching as stage2  # noqa: E402

FWHM_E, FWHM_H = 19.155, 23.366


def _spec(kernel="elliptical_gaussian", **overrides):
    """A Gaussian-only spec that needs no archive (E plane due north unless overridden)."""
    kwargs = dict(kernel=kernel, fwhm_e_deg=FWHM_E, fwhm_h_deg=FWHM_H,
                  e_plane_east_of_meridian_deg=0.0, horizon_deg=90.0, nside_work=32,
                  eval_subdivision=1, max_broadening_fraction=0.002,
                  broad_native_kernel="elliptical_gaussian", chunk=64)
    kwargs.update(overrides)
    return eb.BeamSpec(**kwargs)


def _vec(ra_deg, dec_deg):
    return np.array(hp.ang2vec(ra_deg, dec_deg, lonlat=True), dtype=float)


def _fwhm(theta, response):
    r = np.asarray(response) / response[0]
    i = int(np.argmax(r < 0.5))
    return 2 * (theta[i - 1] + (0.5 - r[i - 1]) / (r[i] - r[i - 1]) * (theta[i] - theta[i - 1]))


@pytest.fixture(scope="module")
def cut_spec():
    try:
        config = load_config()
        return eb.beam_spec_from_config(config, kernel="tris_cuts")
    except Exception as err:  # limTOD or the archive missing
        pytest.skip("archive cuts unavailable: {}".format(err))


# ---------------------------------------------------------------------------
# geometry
# ---------------------------------------------------------------------------
def test_local_offsets_separation_and_position_angle():
    centre = _vec(0.0, 0.0)[None, :]
    points = np.array([_vec(10.0, 0.0), _vec(0.0, 10.0), _vec(350.0, 0.0), _vec(0.0, -10.0)])
    theta, psi = eb.local_offsets(centre, points)
    np.testing.assert_allclose(theta[:, 0], 10.0, atol=1e-9)
    # east (increasing RA) = 90, north = 0, west = -90, south = 180
    np.testing.assert_allclose(psi[:, 0], [90.0, 0.0, -90.0, 180.0], atol=1e-9)


def test_local_offsets_north_is_north_away_from_the_equator():
    # At Dec +42, a point due north (same RA, higher Dec) has position angle 0,
    # and a point at the same Dec but larger RA is to the east (psi near 90).
    centre = _vec(120.0, 42.0)[None, :]
    theta, psi = eb.local_offsets(centre, np.array([_vec(120.0, 52.0), _vec(125.0, 42.0)]))
    assert theta[0, 0] == pytest.approx(10.0, abs=1e-9)
    assert psi[0, 0] == pytest.approx(0.0, abs=1e-9)
    assert 80.0 < psi[1, 0] < 92.0


def test_eval_positions_lie_inside_their_output_pixel():
    pixels = np.array([0, 100, 3071])
    for sub in (1, 2, 4):
        vecs = eb.eval_positions(16, pixels, sub)
        assert vecs.shape == (pixels.size * sub * sub, 3)
        owners = hp.vec2pix(16, vecs[:, 0], vecs[:, 1], vecs[:, 2])
        np.testing.assert_array_equal(owners, np.repeat(pixels, sub * sub))


# ---------------------------------------------------------------------------
# the kernel shapes
# ---------------------------------------------------------------------------
def test_gaussian_kernel_is_half_power_at_half_the_fwhm_on_each_plane():
    spec = _spec(e_plane_east_of_meridian_deg=7.0)
    k = eb.kernel_for(spec, 0.0)
    # E plane at position angle 7, H plane at 97.
    r = eb.kernel_response(k, spec, np.array([FWHM_E / 2, FWHM_H / 2]), np.array([7.0, 97.0]))
    np.testing.assert_allclose(r, 0.5, rtol=1e-12)


def test_kernel_is_zero_beyond_the_horizon():
    spec = _spec(horizon_deg=90.0)
    k = eb.kernel_for(spec, 0.0)
    r = eb.kernel_response(k, spec, np.array([89.0, 91.0]), np.array([0.0, 0.0]))
    assert r[0] > 0.0 and r[1] == 0.0


def test_quadrature_kernel_for_a_gaussian_native_beam():
    k = eb.kernel_for(_spec(), 11.6)
    assert k.fwhm_e_deg == pytest.approx(np.sqrt(FWHM_E ** 2 - 11.6 ** 2))
    assert k.fwhm_h_deg == pytest.approx(np.sqrt(FWHM_H ** 2 - 11.6 ** 2))
    with pytest.raises(ValueError, match="narrow axis"):
        eb.kernel_for(_spec(), 20.0)            # wider than E: cannot be matched


def test_cut_response_matches_limtod_exactly(cut_spec):
    from moonrunes.stage1_tris_maps import import_limtod_tris

    tris = import_limtod_tris(load_config())
    cuts = tris.read_tris_beam_cuts(
        load_config().resolve_path("paths.tris_archive_dir") / eb.BEAM_PROFILE_FILE)
    rng = np.random.default_rng(3)
    theta, phi = rng.uniform(0, 180, 500), rng.uniform(-360, 360, 500)
    for blend in ("db", "power"):
        ours = eb.cut_response(cut_spec.cut_angle_deg, cut_spec.cut_e_db, cut_spec.cut_h_db,
                               theta, phi, blend=blend)
        theirs = tris.tris_cut_beam_response(cuts, theta, phi, blend=blend)
        np.testing.assert_allclose(ours, theirs, rtol=1e-12)


def test_cut_kernel_widths_on_its_principal_planes(cut_spec):
    k = eb.kernel_for(cut_spec, 0.0)
    theta = np.arange(0.0, 40.0, 0.01)
    pa_e = cut_spec.e_plane_east_of_meridian_deg
    e = eb.kernel_response(k, cut_spec, theta, np.full(theta.size, pa_e))
    h = eb.kernel_response(k, cut_spec, theta, np.full(theta.size, pa_e + 90.0))
    assert _fwhm(theta, e) == pytest.approx(FWHM_E, abs=0.01)
    assert _fwhm(theta, h) == pytest.approx(FWHM_H, abs=0.01)


def test_scaled_cuts_have_the_quadrature_widths(cut_spec):
    k = eb.kernel_for(cut_spec, 11.6)            # ARCADE 2: too broad to neglect
    assert k.shape == "tris_cuts" and k.scale_e > 1.0 and k.scale_h > 1.0
    theta = np.arange(0.0, 40.0, 0.01)
    pa_e = cut_spec.e_plane_east_of_meridian_deg
    e = eb.kernel_response(k, cut_spec, theta, np.full(theta.size, pa_e))
    h = eb.kernel_response(k, cut_spec, theta, np.full(theta.size, pa_e + 90.0))
    assert _fwhm(theta, e) == pytest.approx(np.sqrt(FWHM_E ** 2 - 11.6 ** 2), abs=0.02)
    assert _fwhm(theta, h) == pytest.approx(np.sqrt(FWHM_H ** 2 - 11.6 ** 2), abs=0.02)


def test_a_narrow_native_beam_gets_the_cut_beam_itself(cut_spec):
    k = eb.kernel_for(cut_spec, 0.93333)          # Haslam
    assert k.shape == "tris_cuts" and k.scale_e == 1.0 and k.scale_h == 1.0
    assert "neglected" in k.note


# ---------------------------------------------------------------------------
# the convolution
# ---------------------------------------------------------------------------
def test_a_constant_sky_smooths_to_the_same_constant():
    sky = np.full(hp.nside2npix(32), 7.5)
    out = eb.smooth_elliptical(sky, np.ones(sky.size, bool), 0.0, _spec(), 16)
    np.testing.assert_allclose(out["values"], 7.5, rtol=1e-10)
    np.testing.assert_allclose(out["weight"], 1.0, rtol=1e-10)


def test_masked_constant_does_not_sag_at_the_edge():
    # Normalised convolution: a constant over half the sky comes back constant
    # wherever it comes back at all, and the weight is ~0.5 on the boundary.
    nside = 32
    _, dec = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), lonlat=True)
    observed = dec > 0
    sky = np.where(observed, 5.0, 0.0)
    out = eb.smooth_elliptical(sky, observed, 0.0, _spec(), 16)
    good = out["weight"] > 0.05
    np.testing.assert_allclose(out["values"][good], 5.0, rtol=1e-9)
    _, dec16 = hp.pix2ang(16, np.arange(hp.nside2npix(16)), lonlat=True)
    near_equator = np.abs(dec16) < 2.0
    assert np.all(np.abs(out["weight"][near_equator] - 0.5) < 0.15)


def test_point_source_is_narrow_north_south_and_wide_east_west():
    # The orientation test: with the E plane due north (PA 0), a point source
    # smoothed with the elliptical Gaussian is 19.155 deg wide along Dec and
    # 23.366 deg along RA (the H plane is east-west).
    spec = _spec(nside_work=64, e_plane_east_of_meridian_deg=0.0)
    nside = 64
    sky = np.zeros(hp.nside2npix(nside))
    src = hp.ang2pix(nside, 180.0, 0.0, lonlat=True)
    sky[src] = 1.0
    ra0, dec0 = hp.pix2ang(nside, src, lonlat=True)
    offsets = np.arange(0.0, 30.0, 0.1)
    north = np.array([_vec(ra0, dec0 + d) for d in offsets])
    east = np.array([_vec(ra0 + d, dec0) for d in offsets])      # dec0 ~ 0: RA offset = angle
    v_n = eb.smooth_at_positions(sky, np.ones(sky.size, bool), 0.0, spec, north)["values"]
    v_e = eb.smooth_at_positions(sky, np.ones(sky.size, bool), 0.0, spec, east)["values"]
    # A single 0.92 deg working pixel adds ~0.6 deg in quadrature: < 0.02 deg here.
    assert _fwhm(offsets, v_n) == pytest.approx(FWHM_E, abs=0.1)
    assert _fwhm(offsets, v_e) == pytest.approx(FWHM_H, abs=0.1)


def test_match_to_tris_beam_keeps_only_observed_pixels_that_pass_the_floor():
    config = load_config()
    config._data["beam_matching"]["kernel"] = "elliptical_gaussian"
    nside = 32
    _, dec = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), lonlat=True)
    observed = (dec > 20) & (dec < 70)          # a 50 deg band, wider than the kernel
    sky = np.where(observed, 3.0, hp.UNSEEN)
    spec = _spec(nside_work=32)
    values, kept, info = stage2.match_to_tris_beam(sky, observed, 0.0, config, nside_out=16,
                                                   rule="all_children", spec=spec)
    assert kept.any()
    np.testing.assert_allclose(values[kept], 3.0, rtol=1e-9)
    assert np.all(values[~kept] == hp.UNSEEN)
    # every kept nside 16 pixel lies inside the observed band
    _, dec16 = hp.pix2ang(16, np.flatnonzero(kept), lonlat=True)
    assert dec16.min() > 20 and dec16.max() < 70
    assert info["kernel"] == "elliptical_gaussian" and info["extra_fwhm_deg"] is None


def test_spec_rejects_nonsense():
    with pytest.raises(ValueError):
        _spec(kernel="square")
    with pytest.raises(ValueError):
        _spec(eval_subdivision=3)
    with pytest.raises(ValueError, match="archive cuts"):
        _spec(kernel="tris_cuts")
