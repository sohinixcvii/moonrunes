"""Stage 1 tests.

The expensive things (the archive, pyGDSM's data cache) are skipped when they
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
        assert entry["solver"]["info"] == 0
        assert entry["solver"]["krylov_vs_woodbury_sigma"] < 1e-4
    assert manifest["blocker_2_substitute"]["prior_equals_truth"] is True

    # The loader the notebooks and stages 3/5 read products through.
    products = stage1.load_stage1(
        fast, maps_path=maps_path, manifest_path=tmp_path / "stage1_manifest.json"
    )
    assert products.nside == 8
    assert list(products.frequencies_mhz) == [600.0, 820.0]
    assert products.index(820) == 1
    assert products.diagnostics(600)["solver"]["info"] == 0
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
