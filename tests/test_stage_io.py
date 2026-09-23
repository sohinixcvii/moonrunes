"""Stage 1 and stage 2 tests.

The expensive things (the archive, the FITS downloads) are skipped when they
are absent, so the algebra tests still run anywhere.  The algebra is the part
worth pinning: the Krylov solve, its preconditioner convention, and the
Woodbury posterior all have to agree with a dense reference, because a solve
that quietly returns the wrong answer with ``info == 0`` is the exact failure
mode this pipeline was written to avoid.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from moonrunes.config import ConfigError, load_config  # noqa: E402
from moonrunes import stage1_tris_maps as stage1  # noqa: E402
from moonrunes import stage2_beam_matching as stage2  # noqa: E402


@pytest.fixture(scope="module")
def config():
    return load_config()


def _random_problem(n_samples=12, n_parameters=40, seed=0):
    rng = np.random.default_rng(seed)
    operator = rng.normal(size=(n_samples, n_parameters))
    operator[:, -1] = 1.0  # the zero-level nuisance column
    noise_variance = rng.uniform(0.01, 0.05, size=n_samples)
    prior_mean = rng.normal(size=n_parameters)
    prior_variance = rng.uniform(0.5, 2.0, size=n_parameters)
    data = operator @ prior_mean + rng.normal(0, np.sqrt(noise_variance))
    return operator, data, noise_variance, prior_mean, prior_variance


def _dense_reference(operator, data, noise_variance, prior_mean, prior_variance):
    lhs = operator.T @ (operator / noise_variance[:, None]) + np.diag(
        1.0 / prior_variance
    )
    rhs = operator.T @ (data / noise_variance) + prior_mean / prior_variance
    covariance = np.linalg.inv(lhs)
    return covariance @ rhs, np.diag(covariance)


# ---------------------------------------------------------------------------
# config
# ---------------------------------------------------------------------------
def test_require_rejects_open_decisions(config):
    with pytest.raises(ConfigError, match="open decision"):
        config.require("tris.uncertainty_floor_k")
    with pytest.raises(ConfigError, match="not in"):
        config.require("tris.no_such_key")


def test_solver_decisions_are_the_recorded_ones(config):
    assert config.require("tris.mapmaker") == "prior_regularized"
    assert config.require("tris.solver.method") in stage1._SOLVERS
    assert config.require("gibbs.use_woodbury") is False  # Blocker 3


# ---------------------------------------------------------------------------
# Blocker 1: stage 1's grid has to survive the trip to bayesian_skymap
# ---------------------------------------------------------------------------
def test_pinned_nside_new_matches_bayesian_skymap(config):
    check = stage1._check_beam_against_bayesian_skymap(config)
    assert check["nside_new_pinned"] == check["nside_new_from_bayesian_func"] == 8


# ---------------------------------------------------------------------------
# the solve
# ---------------------------------------------------------------------------
def test_lhs_diagonal_matches_bayesian_skymap_probe(config):
    """Our analytic Jacobi diagonal is bayesian_skymap's probe, done cheaply."""
    bayesian_func = stage1.import_bayesian_func(config)
    operator, _data, noise_variance, _mean, prior_variance = _random_problem()
    inv_noise = 1.0 / noise_variance
    inv_prior = 1.0 / prior_variance

    def lhs_op(vector):
        return operator.T @ (inv_noise * (operator @ vector)) + inv_prior * vector

    probed = bayesian_func.estimate_diag_precond(lhs_op, operator.shape[1])
    analytic = stage1.lhs_diagonal(operator, inv_noise, inv_prior)
    np.testing.assert_allclose(analytic, probed, rtol=1e-12)


@pytest.mark.parametrize("method", ["gmres", "cg", "bicgstab"])
def test_krylov_solve_matches_the_dense_solve(method):
    problem = _random_problem()
    solution = stage1.krylov_map_solve(*problem, method=method, rtol=1e-12)
    expected, _variance = _dense_reference(*problem)
    assert solution.converged
    assert solution.relative_residual < 1e-9
    np.testing.assert_allclose(solution.parameters, expected, rtol=1e-7, atol=1e-10)


def test_woodbury_matches_the_dense_posterior():
    problem = _random_problem()
    mean, variance = stage1.woodbury_reference(*problem)
    expected_mean, expected_variance = _dense_reference(*problem)
    np.testing.assert_allclose(mean, expected_mean, rtol=1e-9, atol=1e-12)
    np.testing.assert_allclose(variance, expected_variance, rtol=1e-9, atol=1e-12)


def test_woodbury_stays_exact_when_pixels_outnumber_samples():
    """The real regime: 120 samples, tens of thousands of pixels."""
    problem = _random_problem(n_samples=8, n_parameters=300, seed=3)
    mean, variance = stage1.woodbury_reference(*problem)
    expected_mean, expected_variance = _dense_reference(*problem)
    np.testing.assert_allclose(mean, expected_mean, rtol=1e-8, atol=1e-11)
    np.testing.assert_allclose(variance, expected_variance, rtol=1e-8, atol=1e-11)


def test_krylov_solve_rejects_an_unknown_method():
    with pytest.raises(ValueError, match="solver must be one of"):
        stage1.krylov_map_solve(*_random_problem(), method="cgs")


# ---------------------------------------------------------------------------
# end to end, when the data is there
# ---------------------------------------------------------------------------
def test_beam_response_selection_is_per_observation_peak():
    """A pixel is kept if it clears the threshold at ANY one observation.

    Row 0 peaks at 10.0, so at threshold 0.1 its pixels must clear 1.0; row 1
    peaks at 1.0 and clears at 0.1.  Pixel 3 fails row 0 but passes row 1 --
    "at least one observation" is what keeps it.
    """
    operator = np.array([
        [10.0, 2.0, 0.5, 0.05],
        [ 0.1, 0.0, 0.0, 0.50],
    ])
    keep, info = stage1.beam_response_pixels(operator, 0.1)

    assert list(keep) == [0, 1, 3]
    assert info["n_pixels"] == 3
    assert info["threshold"] == 0.1
    assert info["sky_fraction"] == pytest.approx(3 / 4)


def test_beam_response_selection_tightens_monotonically():
    """A larger threshold can only ever keep fewer pixels."""
    rng = np.random.default_rng(0)
    operator = rng.uniform(0.0, 1.0, size=(20, 400))

    sizes = [stage1.beam_response_pixels(operator, t)[0].size
             for t in (0.003, 0.01, 0.03, 0.3)]

    assert sizes == sorted(sizes, reverse=True)


def test_beam_response_selection_reports_the_power_it_dropped():
    """The diagnostic that says whether the cut threw away real signal."""
    operator = np.array([[1.0, 0.5, 0.001]])

    keep, info = stage1.beam_response_pixels(operator, 0.01)

    assert list(keep) == [0, 1]
    assert info["beam_power_retained"] == pytest.approx(1.5 / 1.501)


def test_beam_response_selection_rejects_a_nonsense_threshold():
    operator = np.ones((2, 8))
    for bad in (0.0, 1.0, -0.1, 3.0):
        with pytest.raises(ValueError, match="beam_response_threshold"):
            stage1.beam_response_pixels(operator, bad)


def test_beam_response_selection_always_keeps_each_observation_peak():
    """The selection can never come back empty.

    Each row's own maximum clears any threshold below 1 by definition, so
    every observation contributes at least its peak pixel however severe the
    cut. That is why there is no empty-selection error to raise.
    """
    operator = np.array([[1.0, 1e-9], [1e-9, 2.0]])

    for threshold in (0.003, 0.01, 0.03, 0.9, 0.999):
        keep, info = stage1.beam_response_pixels(operator, threshold)
        assert list(keep) == [0, 1], threshold
        assert info["n_pixels"] == 2


def test_beam_response_selection_needs_a_positive_peak():
    """An all-zero observation means the geometry or mask is wrong."""
    operator = np.array([[1.0, 0.5], [0.0, 0.0]])
    with pytest.raises(RuntimeError, match="no positive beam response"):
        stage1.beam_response_pixels(operator, 0.01)


def _archive_or_skip(config):
    try:
        archive = config.resolve_path("paths.tris_archive_dir")
    except ConfigError:
        pytest.skip("paths.tris_archive_dir is not set")
    if not (archive / "TRIS_absolute_600.txt").exists():
        pytest.skip("no TRIS archive at {}".format(archive))
    return archive


def test_2500_mhz_is_refused_as_a_point_set(config):
    archive = _archive_or_skip(config)
    with pytest.raises(ValueError, match="point set"):
        stage1.solve_frequency(2500, config=config, archive_dir=archive, verbose=False)


def test_stage1_end_to_end_writes_the_expected_product(config, tmp_path):
    archive = _archive_or_skip(config)
    pytest.importorskip("pygdsm")

    # Same code path, a coarse grid so the test is seconds not minutes.
    fast = load_config(config.path)
    fast._data["tris"]["nside"] = 8
    fast._data["tris"]["nside_hires"] = 8
    fast._data["run"]["name"] = "pytest"

    manifest = stage1.run_stage1(
        config=fast, archive_dir=archive, output_dir=tmp_path, verbose=False
    )
    maps_path = tmp_path / "tris_maps_pytest.npz"
    assert maps_path.exists()
    assert (tmp_path / "stage1_manifest.json").exists()

    with np.load(maps_path) as bundle:
        assert int(bundle["nside"]) == 8
        n_freq = bundle["freq_mhz"].size
        n_pix = bundle["pixel_indices"].size
        assert n_freq == 2
        for key in ("sky_k", "sky_sigma_k", "prior_mean_k", "prior_sigma_k", "truth_k"):
            assert bundle[key].shape == (n_freq, n_pix), key
        assert np.all(bundle["sky_sigma_k"] > 0)
        # The posterior may not beat the prior by much -- 120 samples cannot
        # constrain hundreds of pixels -- but it must never be looser.
        assert np.all(bundle["sky_sigma_k"] <= bundle["prior_sigma_k"] * (1 + 1e-9))
        assert np.isfinite(bundle["sky_k"]).all()
        # Off-band pixels are NaN by design; stage 3 chooses how to fill them.
        assert np.isnan(bundle["sky_full_k"][:, ~bundle["band_mask"]]).all()

    for entry in manifest["frequencies"].values():
        # limTOD's dense solve has no iteration count to check; the Woodbury
        # cross-check is the whole convergence statement now.
        assert entry["solver"]["backend"] == "limtod"
        assert entry["solver"]["solver_vs_woodbury_sigma"] < 1e-4
    assert manifest["blocker_2_substitute"]["prior_equals_truth"] is True

    # The loader the notebooks and stages 3/5 read products through.
    products = stage1.load_stage1(
        fast, maps_path=maps_path, manifest_path=tmp_path / "stage1_manifest.json"
    )
    assert products.nside == 8
    assert list(products.frequencies_mhz) == [600.0, 820.0]
    assert products.index(820) == 1
    assert products.diagnostics(600)["solver"]["backend"] == "limtod"
    np.testing.assert_allclose(
        products.band_map("sky_k", 600)[products.pixel_indices],
        products["sky_k"][0],
    )
    with pytest.raises(KeyError, match="2500"):
        products.index(2500)


def test_load_stage1_says_what_to_run_when_nothing_is_there(tmp_path):
    config = load_config()
    with pytest.raises(FileNotFoundError, match="run_pipeline.py --stage 1"):
        stage1.load_stage1(config, maps_path=tmp_path / "absent.npz")


def test_stage1_refuses_to_clobber(config, tmp_path):
    archive = _archive_or_skip(config)
    pytest.importorskip("pygdsm")
    fast = load_config(config.path)
    fast._data["tris"]["nside"] = 4
    fast._data["tris"]["nside_hires"] = 4
    fast._data["run"]["name"] = "pytest_clobber"
    kwargs = dict(config=fast, archive_dir=archive, output_dir=tmp_path, verbose=False)

    stage1.run_stage1(**kwargs)
    with pytest.raises(FileExistsError):
        stage1.run_stage1(**kwargs)
    stage1.run_stage1(overwrite=True, **kwargs)  # explicit override is fine

# ===========================================================================
# stage 2 -- beam matching
# ===========================================================================
def test_extra_fwhm_subtracts_in_quadrature():
    """The kernel is sqrt(target^2 - native^2), not the target itself."""
    assert stage2.extra_fwhm_deg(3.0, 5.0) == pytest.approx(4.0)
    # Smoothing an already-matched map is a no-op, not an error.
    assert stage2.extra_fwhm_deg(5.0, 5.0) == pytest.approx(0.0)


def test_extra_fwhm_refuses_to_sharpen():
    """A native beam coarser than the target is a bookkeeping error.

    np.sqrt of a negative returns NaN, which would propagate through
    hp.smoothing into an all-NaN map with no exception anywhere.
    """
    with pytest.raises(ValueError, match="coarsest resolution"):
        stage2.extra_fwhm_deg(30.0, stage2.TARGET_FWHM_DEG)


def test_target_is_the_larger_tris_axis():
    """Step 0.1: the conservative choice, so nothing stays under-smoothed."""
    assert stage2.TARGET_FWHM_DEG == stage2.TRIS_HPBW_H_DEG
    assert stage2.TARGET_FWHM_DEG > stage2.TRIS_HPBW_E_DEG


def test_haslam_beam_is_configured_in_degrees_not_arcminutes(config):
    """BEAMSIZE = 56.0 is in arcmin; as degrees it would be 60x too coarse."""
    haslam_fwhm = float(config.require("beam_matching.haslam_native_fwhm_deg"))
    assert haslam_fwhm == pytest.approx(56.0 / 60.0, rel=1e-4)
    assert haslam_fwhm < 1.0


def test_beam_matching_decisions_are_the_recorded_ones(config):
    """The config is the source of truth; the module may not re-derive it."""
    target = float(config.require("beam_matching.target_fwhm_deg"))
    assert target == pytest.approx(stage2.TRIS_HPBW_H_DEG)
    assert float(config.require("beam_matching.weight_floor")) == 0.5
    # Every dataset's beam comes from a config key, not a literal.
    for dataset in stage2.DATASETS:
        assert float(config.require(dataset.fwhm_key)) > 0.0


def test_observed_mask_knows_each_datasets_convention():
    """Three datasets, three conventions, none of them hp.UNSEEN on disk."""
    hp = pytest.importorskip("healpy")
    sky = np.array([1.0, 0.0, 2.0, np.nan, hp.UNSEEN])

    zeros = stage2.observed_mask(sky, "zeros")      # ARCADE 2
    assert list(zeros) == [True, False, True, False, False]

    nans = stage2.observed_mask(sky, "nan")         # stage 1 TRIS maps
    assert list(nans) == [True, True, True, False, False]

    # A NaN is unobserved under every convention, including ARCADE 2's.
    assert not stage2.observed_mask(np.array([np.nan]), "zeros")[0]

    full = stage2.observed_mask(np.array([1.0, 0.0, 2.0]), "full")   # Haslam
    assert full.all()


def test_observed_mask_rejects_an_unknown_convention():
    with pytest.raises(ValueError, match="unknown mask convention"):
        stage2.observed_mask(np.zeros(12), "sentinel")


def test_smoothing_conserves_the_mean_of_a_full_sky_map():
    """A convolution redistributes power, it does not create or destroy it."""
    hp = pytest.importorskip("healpy")
    rng = np.random.default_rng(0)
    sky = rng.normal(10.0, 1.0, size=hp.nside2npix(32))

    smoothed = stage2.smooth_to_target(sky, 1.0, 20.0)

    # Not exact: map2alm's quadrature leaves a few parts in 1e6 at nside 32.
    assert smoothed.mean() == pytest.approx(sky.mean(), rel=1e-4)
    assert smoothed.std() < sky.std()      # and it really did smooth


def test_smoothing_does_not_modify_its_input():
    """The old script mutated the arrays hp.read_map handed it."""
    hp = pytest.importorskip("healpy")
    sky = np.linspace(1.0, 2.0, hp.nside2npix(16))
    mask = np.ones_like(sky, dtype=bool)
    mask[:100] = False
    before = sky.copy()

    stage2.smooth_to_target(sky, 5.0, 20.0, weight_mask=mask)

    assert np.array_equal(sky, before)


def test_masked_smoothing_does_not_drag_the_mask_into_the_data():
    """Step 0.3: the weight-map division is what makes this correct.

    A constant sky over a partial mask must smooth back to that constant.
    Without the division the unobserved zeros dilute every pixel near the
    boundary, which is exactly the bleed the procedure exists to prevent.
    """
    hp = pytest.importorskip("healpy")
    npix = hp.nside2npix(32)
    theta = hp.pix2ang(32, np.arange(npix))[0]
    mask = theta < np.radians(90)          # one hemisphere observed
    sky = np.where(mask, 5.0, 0.0)

    corrected = stage2.smooth_to_target(sky, 5.0, 20.0, weight_mask=mask)
    naive = stage2.smooth_to_target(sky, 5.0, 20.0)

    survived = corrected != hp.UNSEEN
    assert survived.sum() > 0
    assert corrected[survived] == pytest.approx(5.0, abs=1e-3)
    # The uncorrected version sags towards zero across the boundary.
    assert naive[mask].min() < 4.0


def test_a_feature_narrower_than_the_kernel_cannot_survive():
    """The weight floor is what stops a thin sliver being smoothed into data.

    A band 10 deg wide, smoothed with a 20 deg kernel, never accumulates half
    its own weight anywhere -- so every pixel of it is dropped rather than
    handed on as a diluted temperature.  This is the protection Step 0.3's
    threshold exists to give.
    """
    hp = pytest.importorskip("healpy")
    theta = hp.pix2ang(32, np.arange(hp.nside2npix(32)))[0]
    mask = np.abs(np.degrees(theta) - 90) < 5.0
    sky = np.where(mask, 3.0, 0.0)

    corrected = stage2.smooth_to_target(sky, 5.0, 20.0, weight_mask=mask)

    assert np.all(corrected == hp.UNSEEN)


def test_a_mask_wider_than_the_kernel_is_not_eroded():
    """The flip side: no spurious loss where the weight is genuinely high.

    The cut sits at weight 0.5, which for a symmetric kernel is the mask
    boundary itself, so a region comfortably wider than the kernel comes back
    with its footprint intact.  Erosion is a statement about the mask's shape,
    not an unavoidable tax on every smoothing.
    """
    hp = pytest.importorskip("healpy")
    theta = hp.pix2ang(32, np.arange(hp.nside2npix(32)))[0]
    mask = np.abs(np.degrees(theta) - 90) < 30.0
    sky = np.where(mask, 3.0, 0.0)

    corrected = stage2.smooth_to_target(sky, 5.0, 20.0, weight_mask=mask)

    assert np.array_equal(corrected != hp.UNSEEN, mask)
    assert corrected[mask] == pytest.approx(3.0, abs=1e-3)


def test_weight_mask_shape_is_checked():
    hp = pytest.importorskip("healpy")
    sky = np.zeros(hp.nside2npix(16))
    with pytest.raises(ValueError, match="weight_mask has shape"):
        stage2.smooth_to_target(sky, 1.0, 20.0, weight_mask=np.ones(5, dtype=bool))


# ---------------------------------------------------------------------------
# stage 2 end to end -- needs the FITS downloads
# ---------------------------------------------------------------------------
def _maps_or_skip(config):
    missing = [
        dataset.path_key
        for dataset in stage2.DATASETS
        if not config.resolve_path(dataset.path_key).exists()
    ]
    if missing:
        pytest.skip("no local maps for {} -- see DATA_SOURCES.md".format(missing))


def test_beam_match_runs_over_every_dataset(config):
    hp = pytest.importorskip("healpy")
    _maps_or_skip(config)

    results = stage2.beam_match(config, verbose=False)

    assert set(results) == {d.name for d in stage2.DATASETS}
    for name, entry in results.items():
        smoothed = entry["map"]
        good = smoothed != hp.UNSEEN
        assert good.any(), name
        assert np.isfinite(smoothed[good]).all(), name
        # Every kernel is the quadrature difference, never the raw target.
        assert entry["extra_fwhm_deg"] < stage2.TARGET_FWHM_DEG
        assert entry["extra_fwhm_deg"] == pytest.approx(
            stage2.extra_fwhm_deg(entry["native_fwhm_deg"], stage2.TARGET_FWHM_DEG)
        )


def test_beam_match_erodes_only_the_masked_datasets(config):
    """Haslam is full sky and must keep every pixel; ARCADE 2 loses its edge."""
    pytest.importorskip("healpy")
    _maps_or_skip(config)

    results = stage2.beam_match(config, verbose=False)

    haslam = results["haslam_408"]
    assert haslam["observed_after"] == haslam["observed_before"]
    for name in ("arcade2_3150", "arcade2_3410"):
        entry = results[name]
        assert entry["observed_after"] < entry["observed_before"]


def test_beam_match_says_when_a_map_is_missing(config, tmp_path):
    """The error names the config key to set, not just the path."""
    pytest.importorskip("healpy")
    config = load_config()
    config._data["paths"]["arcade2_map_3150"] = str(tmp_path / "absent.fits")
    with pytest.raises(FileNotFoundError, match="paths.arcade2_map_3150"):
        stage2.beam_match(config, verbose=False)
