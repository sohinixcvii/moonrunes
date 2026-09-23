"""Stage 1 -- calibrated TRIS sky maps from the public archive.

What this stage produces
------------------------
One prior-regularized map per TRIS frequency (600 and 820 MHz; 2.5 GHz is
dropped, it is a point set and is flagged poor quality in the literature
review), on a declination band around the ring, together with the fitted
zero-level offset and a full posterior standard deviation.  These maps are what
stage 3 writes into ``bayesian_skymap``'s ``sky_lowres``/``data1`` key: TRIS is
an absolute radiometer, so it is the absolutely-calibrated leg of the
recalibration.

How it is built
---------------
The forward model is the one verified in ``limTOD``'s TRIS walkthrough
(``examples/TRIS/tris_limtod_walkthrough.ipynb``) and nothing here re-derives
it: the beam comes from the archive's own E/H principal-plane cuts (HPBW
19.155 / 23.366 deg), the pointing is the parked-zenith geometry
(``RA = LST``, ``dec = lat = 42 deg 26'``) with the E plane rolled 7 deg east
via ``selfrot = -7``, and the below-horizon response is masked.
``limTOD.tris.build_tris_mapmaking_inputs`` assembles the operator, the data
and the noise; the zero level is carried as an explicit nuisance parameter with
a Gaussian prior, which is exactly equivalent to the rank-1 noise term and also
hands back the fitted offset.

The solve
---------
**The solve is limTOD's own.**  ``TRISMapMakingInputs.solve`` wraps
``limTOD.wiener_filter_map``, which forms the dense normal matrix
``A^T N^-1 A + S^-1 + eps I`` and factorises it directly.  limTOD builds the
prior from ``prior_map``/``prior_sigma_k`` and takes the zero-level column's
width from the ``zero_level_sigma_k`` handed to
``build_tris_mapmaking_inputs``, so nothing is stacked by hand here.

That path is dense in the number of parameters, and the cost is not academic:
at ``tris.nside = 64`` the band carries 25705 parameters, so each of the
several ``n x n`` float64 arrays it allocates is ~4.9 GB and the peak is
~25 GB.  It runs comfortably at nside 32 and below.  See
``TRIS_MAP_PROVENANCE.md``.

:func:`krylov_map_solve` is retained as the matrix-free reference
implementation the test suite pins against a dense solve; it is no longer what
makes the maps.  ``bayesian_skymap.bayesian_func`` is still imported for real:
its ``nside_for_beam`` validates the pinned ``beam.nside_new`` against the TRIS
beam width (Blocker 1), and its ``estimate_diag_precond`` defines the
preconditioner convention that :func:`lhs_diagonal` reproduces analytically.

Because "it returned an array" is not "it solved the system" -- a dense
factorisation with a regularizer added to the diagonal can come back biased
without raising -- every solve is cross-checked against the exact Woodbury
identity

    (A^T N^-1 A + P^-1)^-1 = P - P A^T (N + A P A^T)^-1 A P,

which is cheap here (the ring has 120 samples, so the inner matrix is 120x120)
and is *not* the approximation Blocker 3 rules out -- that one drops
cross-frequency coupling in the multi-frequency Gibbs solver; this one is an
algebraic identity at a single frequency.  The same identity supplies the
posterior standard deviation, which a Krylov solve alone cannot give.

Blocker 2, honestly
-------------------
There is no true sky for real data, and Berkhuijsen (1972) is not in hand yet,
so GSM2008 stands in as both the prior template and the reference "truth".
That is circular and the manifest says so (``prior_equals_truth``).  The
walkthrough measured GSM2008 as ~1.9 K cold at 600 MHz; with a template that
wrong, the fitted zero level absorbs the model offset and the reduced
chi-square climbs.  Both numbers are reported per frequency rather than
smoothed over -- they are the honest product of one ring plus a template, and
stage 5 reads them.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Optional, Sequence, Tuple

import numpy as np

from .config import REPO_ROOT, Config, ConfigError, load_config

STAGE = "stage1"

# The archive ships one file per frequency; 2.5 GHz is a point set, not a ring,
# and this stage only makes maps from rings.
_RING_FILENAMES = {600: "TRIS_absolute_600.txt", 820: "TRIS_absolute_820.txt"}
_POINT_SET_FILENAMES = {2500: "TRIS_absolute_2500MHz.txt"}
_BEAM_FILENAME = "TRIS_Beam_Profile.txt"

_SOLVERS = ("gmres", "cg", "bicgstab")


# ---------------------------------------------------------------------------
# imports that live outside this package
# ---------------------------------------------------------------------------
def import_limtod_tris(config: Optional[Config] = None):
    """Import ``limTOD.tris``, falling back to the sibling source checkout.

    limTOD is a normal pip dependency, but the walkthrough runs against the
    checkout next door and environments here often have only that.
    """
    try:
        from limTOD import tris  # noqa: F401
    except ImportError:
        candidate = None
        if config is not None:
            try:
                candidate = config.resolve_path("paths.limtod")
            except ConfigError:
                candidate = None
        if candidate is None:
            candidate = (REPO_ROOT.parent / "limTOD").resolve()
        if not (candidate / "limTOD" / "__init__.py").exists():
            raise ImportError(
                "limTOD is neither installed nor at {} -- pip install limTOD, "
                "or point paths.limtod at the checkout".format(candidate)
            )
        sys.path.insert(0, str(candidate))
        from limTOD import tris  # noqa: F811
    return tris


def import_bayesian_func(config: Optional[Config] = None):
    """Import ``bayesian_func`` from the ``bayesian_skymap`` submodule."""
    checkout = (
        config.resolve_path("paths.bayesian_skymap")
        if config is not None
        else (REPO_ROOT / "external" / "bayesian_skymap")
    )
    if not (checkout / "bayesian_func.py").exists():
        raise ImportError(
            "no bayesian_func.py under {} -- is the submodule checked "
            "out? (git submodule update --init)".format(checkout)
        )
    if str(checkout) not in sys.path:
        sys.path.insert(0, str(checkout))
    import bayesian_func  # noqa: F401

    return bayesian_func


# ---------------------------------------------------------------------------
# the Krylov solve (the bayesian_skymap pattern), and its exact reference
# ---------------------------------------------------------------------------
def lhs_diagonal(
    operator: np.ndarray, inv_noise: np.ndarray, inv_prior: np.ndarray
) -> np.ndarray:
    """Diagonal of ``A^T N^-1 A + P^-1``, for Jacobi preconditioning.

    This is the analytic form of what ``bayesian_func.estimate_diag_precond``
    probes column by column.  The probe costs one matvec per parameter, which
    is ~25k dense matvecs at the working nside; this is one einsum.  The test
    suite pins the two against each other.
    """
    return np.einsum("ij,ij->j", operator, inv_noise[:, None] * operator) + inv_prior


@dataclass(frozen=True)
class KrylovSolution:
    """Result of one MAP solve, with everything stage 5 needs to judge it."""

    parameters: np.ndarray
    info: int
    iterations: int
    relative_residual: float
    method: str
    restart: Optional[int]
    seconds: float

    @property
    def converged(self) -> bool:
        return self.info == 0


def krylov_map_solve(
    operator: np.ndarray,
    data: np.ndarray,
    noise_variance: np.ndarray,
    prior_mean: np.ndarray,
    prior_variance: np.ndarray,
    *,
    method: str = "gmres",
    rtol: float = 1e-9,
    maxiter: int = 1000,
    restart: Optional[int] = None,
) -> KrylovSolution:
    """MAP solve of ``(A^T N^-1 A + P^-1) x = A^T N^-1 d + P^-1 mu``.

    Deliberately matrix-free, mirroring ``bayesian_skymap``'s amplitude step:
    a ``LinearOperator`` around the matvec, a Jacobi preconditioner, and a
    Krylov solve.  No dense normal matrix is ever formed -- at the working
    nside it would be 25k x 25k.
    """
    import time

    from scipy.sparse import linalg as splinalg

    if method not in _SOLVERS:
        raise ValueError(
            "solver must be one of {}, got {!r}".format(_SOLVERS, method)
        )
    operator = np.asarray(operator, dtype=float)
    n_samples, n_parameters = operator.shape
    inv_noise = 1.0 / np.asarray(noise_variance, dtype=float)
    inv_prior = 1.0 / np.asarray(prior_variance, dtype=float)
    if inv_prior.shape != (n_parameters,):
        raise ValueError("prior_variance must have one entry per parameter")
    if prior_mean.shape != (n_parameters,):
        raise ValueError("prior_mean must have one entry per parameter")

    def matvec(vector: np.ndarray) -> np.ndarray:
        return operator.T @ (inv_noise * (operator @ vector)) + inv_prior * vector

    lhs = splinalg.LinearOperator(
        matvec=matvec, shape=(n_parameters, n_parameters), dtype=float
    )
    rhs = operator.T @ (inv_noise * data) + inv_prior * prior_mean

    diagonal = lhs_diagonal(operator, inv_noise, inv_prior)
    preconditioner = splinalg.LinearOperator(
        matvec=lambda v: v / diagonal,
        shape=(n_parameters, n_parameters),
        dtype=float,
    )

    iterations = 0

    def callback(_state) -> None:
        nonlocal iterations
        iterations += 1

    kwargs: Dict[str, Any] = dict(
        M=preconditioner, rtol=rtol, atol=0.0, maxiter=maxiter, callback=callback
    )
    if method == "gmres":
        kwargs["callback_type"] = "pr_norm"
        # The data term has rank n_samples, so a preconditioned GMRES cycle
        # longer than that terminates inside one restart instead of stalling
        # at scipy's default restart of 20.
        restart = int(restart) if restart else min(n_parameters, n_samples + 40)
        kwargs["restart"] = restart
    else:
        restart = None

    solver = getattr(splinalg, method)
    start = time.perf_counter()
    parameters, info = solver(lhs, rhs, **kwargs)
    seconds = time.perf_counter() - start

    residual = rhs - matvec(parameters)
    denominator = np.linalg.norm(rhs)
    relative = float(np.linalg.norm(residual) / denominator) if denominator else 0.0
    return KrylovSolution(
        parameters=parameters,
        info=int(info),
        iterations=int(iterations),
        relative_residual=relative,
        method=method,
        restart=restart,
        seconds=float(seconds),
    )


def woodbury_reference(
    operator: np.ndarray,
    data: np.ndarray,
    noise_variance: np.ndarray,
    prior_mean: np.ndarray,
    prior_variance: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Exact MAP mean and posterior variance via the Woodbury identity.

    ``S = P - P A^T (N + A P A^T)^-1 A P`` and
    ``x = mu + P A^T (N + A P A^T)^-1 (d - A mu)``.  The inner matrix is
    ``n_samples x n_samples`` (120 for a TRIS ring), so this is exact and cheap
    even with 25k sky pixels -- it is an algebraic identity at one frequency,
    not the multi-frequency Woodbury acceleration Blocker 3 rules out.

    Returns ``(mean, variance)``, both over all parameters.
    """
    from scipy.linalg import cho_factor, cho_solve

    operator = np.asarray(operator, dtype=float)
    prior_variance = np.asarray(prior_variance, dtype=float)
    scaled = operator * prior_variance  # A P
    inner = scaled @ operator.T
    inner[np.diag_indices_from(inner)] += np.asarray(noise_variance, dtype=float)

    factor = cho_factor(inner, lower=True)
    mean = prior_mean + prior_variance * (
        operator.T @ cho_solve(factor, data - operator @ prior_mean)
    )
    quadratic = np.einsum("ij,ij->j", operator, cho_solve(factor, operator))
    variance = prior_variance - prior_variance**2 * quadratic
    return mean, np.maximum(variance, 0.0)


# ---------------------------------------------------------------------------
# the sky template (Blocker 2 substitute)
# ---------------------------------------------------------------------------
def gsm2008_template(
    frequency_mhz: float, nside: int, *, tris_convention: bool = True
) -> np.ndarray:
    """GSM2008 at ``frequency_mhz``, equatorial RING, in the TRIS convention.

    Pass the ring's *effective* frequency.  The order here (degrade, then
    rotate, then add the CMB monopole back) is the walkthrough's, so the
    template matches the numbers that notebook reports.
    """
    import healpy as hp
    from pygdsm import GlobalSkyModel

    tris = import_limtod_tris()
    galactic = hp.ud_grade(GlobalSkyModel().generate(float(frequency_mhz)), nside)
    equatorial = hp.Rotator(coord=["G", "C"]).rotate_map_pixel(galactic)
    if not tris_convention:
        return equatorial
    return tris.to_tris_temperature_convention(equatorial, float(frequency_mhz))


# ---------------------------------------------------------------------------
# per-frequency product
# ---------------------------------------------------------------------------
@dataclass
class TRISFrequencyMap:
    """One solved TRIS ring: the map, its uncertainty, and its diagnostics."""

    frequency_mhz: float
    effective_frequency_mhz: float
    nside: int
    pixel_indices: np.ndarray
    sky_k: np.ndarray
    sky_sigma_k: np.ndarray
    prior_mean_k: np.ndarray
    prior_sigma_k: np.ndarray
    truth_k: np.ndarray
    zero_level_k: float
    zero_level_sigma_k: float
    zero_level_prior_sigma_k: float
    ra_deg: np.ndarray
    data_k: np.ndarray
    residual_k: np.ndarray
    beam_coverage: np.ndarray
    diagnostics: Dict[str, Any] = field(default_factory=dict)

    def healpix_map(self, fill: float = np.nan) -> np.ndarray:
        import healpy as hp

        full = np.full(hp.nside2npix(self.nside), fill, dtype=float)
        full[self.pixel_indices] = self.sky_k
        return full

    def healpix_sigma(self, fill: float = np.nan) -> np.ndarray:
        import healpy as hp

        full = np.full(hp.nside2npix(self.nside), fill, dtype=float)
        full[self.pixel_indices] = self.sky_sigma_k
        return full


def beam_response_pixels(
    operator: np.ndarray, threshold: float
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """Pixels the beam actually saw, per the collaborator's criterion.

    A pixel is retained when the beam response at it exceeds ``threshold`` of
    the peak response **for at least one observation**.  ``operator`` is the
    full-sky sky-to-sample projection: row ``t`` holds the (normalized) beam
    weight of every pixel at pointing ``t``, so its own maximum is that
    observation's peak response and the comparison is per observation, not
    against a global maximum.

    This replaces the declination-band cut, which selected on geometry alone
    and so kept pixels the beam never illuminated and cut pixels in the wings
    that it did.  The returned diagnostics record what the threshold cost.

    Returns ``(pixel_indices, diagnostics)``, indices sorted ascending.
    """
    if not 0.0 < float(threshold) < 1.0:
        raise ValueError(
            "tris.beam_response_threshold must be in (0, 1), got {!r}".format(
                threshold
            )
        )
    operator = np.asarray(operator, dtype=float)
    peak_per_sample = operator.max(axis=1, keepdims=True)
    if not np.all(peak_per_sample > 0):
        raise RuntimeError(
            "an observation has no positive beam response anywhere on the sky"
        )
    # Every observation's own peak pixel clears any threshold below 1, so the
    # selection is never empty once the peak check above has passed.
    above = operator > (float(threshold) * peak_per_sample)
    keep = np.flatnonzero(above.any(axis=0))
    retained_power = float(
        operator[:, keep].sum() / operator.sum()
    )
    return keep, {
        "threshold": float(threshold),
        "n_pixels": int(keep.size),
        "sky_fraction": float(keep.size / operator.shape[1]),
        "samples_per_pixel": float(operator.shape[0] / keep.size),
        "beam_power_retained": retained_power,
    }


def solve_frequency(
    frequency_mhz: int,
    *,
    config: Config,
    archive_dir: Path,
    verbose: bool = True,
) -> TRISFrequencyMap:
    """Read one ring, build its operator and prior, and solve for the map."""
    import healpy as hp

    tris = import_limtod_tris(config)

    if int(frequency_mhz) in _POINT_SET_FILENAMES:
        raise ValueError(
            "{} MHz is a point set in the archive, not a drift ring, and is "
            "flagged poor quality in the TRIS literature review -- stage 1 "
            "makes maps from rings only".format(frequency_mhz)
        )
    try:
        filename = _RING_FILENAMES[int(frequency_mhz)]
    except KeyError:
        raise ValueError(
            "no TRIS ring for {} MHz; known rings: {}".format(
                frequency_mhz, sorted(_RING_FILENAMES)
            )
        )

    ring = tris.read_tris_ring(archive_dir / filename)
    cuts = tris.read_tris_beam_cuts(archive_dir / _BEAM_FILENAME)
    effective_mhz = float(ring.effective_frequency_mhz)

    nside = int(config.require("tris.nside"))
    nside_hires = int(config.require("tris.nside_hires"))
    if nside_hires < nside:
        raise ConfigError("tris.nside_hires must be >= tris.nside")

    # One row per ring has a 0.000 K statistical error; a floor is not optional.
    floor = config.get_path("tris.uncertainty_floor_k")
    positive = ring.statistical_uncertainty_k[ring.statistical_uncertainty_k > 0]
    floor = float(positive.min()) if floor is None else float(floor)

    zero_prior = float(config.require("tris.zero_level_sigma_k"))

    # ---- pixel selection ---------------------------------------------------
    # Two passes.  The first projects the beam onto the WHOLE sky so the
    # response at every pixel is known; the second rebuilds the operator on
    # just the pixels that pass the beam-response threshold.  Selecting on the
    # beam rather than on a declination band is the collaborator's criterion
    # and is the physical statement: keep what the horn actually saw.
    #
    # The full-sky pass is affordable only because this is a deliberately
    # coarse map -- 768 pixels at nside 8, 3072 at nside 16.
    common = dict(
        cuts=cuts,
        uncertainty_floor_k=floor,
        zero_level_sigma_k=zero_prior,
        apply_horizon_mask=bool(config.require("beam.apply_horizon_mask")),
        nside_hires=nside_hires,
    )
    threshold = float(config.require("tris.beam_response_threshold"))
    full_sky = tris.build_tris_mapmaking_inputs(
        ring,
        nside=nside,
        pixel_indices=np.arange(hp.nside2npix(nside)),
        **common,
    )
    # Drop the trailing zero-level column before thresholding: it is a
    # nuisance parameter, not a pixel.
    sky_columns = np.asarray(full_sky.operator)[:, : hp.nside2npix(nside)]
    pixel_indices, selection = beam_response_pixels(sky_columns, threshold)

    inputs = tris.build_tris_mapmaking_inputs(
        ring,
        nside=nside,
        pixel_indices=pixel_indices,
        **common,
    )
    n_sky = int(inputs.sky_parameter_count)
    if not inputs.has_zero_level:  # pragma: no cover - zero_level_sigma_k is required
        raise RuntimeError("expected a zero-level nuisance column")

    # ---- prior / truth (Blocker 2: GSM2008 stands in for both) -------------
    template_name = str(config.require("tris.prior.template")).lower()
    if template_name != "gsm2008":
        raise ConfigError(
            "tris.prior.template only implements 'gsm2008' (the Blocker 2 "
            "substitute); got {!r}".format(template_name)
        )
    template_full = gsm2008_template(effective_mhz, nside)
    prior_mean_sky, prior_sigma_sky = tris.tris_prior_from_template(
        template_full,
        inputs.pixel_indices,
        relative_sigma=float(config.require("tris.prior.relative_sigma")),
        absolute_sigma_k=float(config.require("tris.prior.absolute_sigma_k")),
        floor_sigma_k=float(config.require("tris.prior.floor_sigma_k")),
    )

    # ---- solve -------------------------------------------------------------
    # limTOD owns the solve: TRISMapMakingInputs.solve wraps
    # limTOD.wiener_filter_map, which forms the dense normal matrix and
    # factorises it directly.  It builds the prior itself from prior_map /
    # prior_sigma_k and takes the zero-level column's width from the
    # zero_level_sigma_k passed to build_tris_mapmaking_inputs, so nothing is
    # stacked by hand here.
    regularization = float(config.require("tris.solver.regularization"))
    import time as _time

    _start = _time.perf_counter()
    solved = inputs.solve(
        prior_map=prior_mean_sky,
        prior_sigma_k=prior_sigma_sky,
        regularization=regularization,
    )
    solver_seconds = _time.perf_counter() - _start

    sky_k = np.asarray(solved.sky_k, dtype=float)
    sky_sigma_k = np.asarray(solved.sky_uncertainty_k, dtype=float)
    zero_level_k = float(solved.zero_level_k)
    zero_level_sigma_k = float(solved.zero_level_uncertainty_k)
    residual = np.asarray(solved.residual_k, dtype=float)
    chi_square = float(solved.chi_square)
    dof = int(solved.degrees_of_freedom)

    # ---- verification ------------------------------------------------------
    # "It returned an array" is not "it solved the system".  The dense path
    # adds `regularization * I` to the normal matrix and inverts a
    # 25k x 25k operator, so it can come back biased or ill-conditioned
    # without raising.  The Woodbury identity is exact at a single frequency
    # (the inner matrix is n_samples x n_samples) and is what says whether it
    # did: same check the Krylov path carried, now pointed at limTOD's answer.
    prior_mean = np.concatenate([prior_mean_sky, [0.0]])
    prior_variance = np.concatenate([prior_sigma_sky**2, [zero_prior**2]])
    reference_mean, posterior_variance = woodbury_reference(
        inputs.operator,
        inputs.data_k,
        inputs.noise.variance_k2,
        prior_mean,
        prior_variance,
    )
    # The meaningful unit is the posterior width, not the residual norm: the
    # prior width spans an order of magnitude between the Galactic plane and
    # the cold sky, so an absolute difference means nothing on its own.
    posterior_sigma = np.sqrt(np.maximum(posterior_variance, 0.0))
    parameters = np.concatenate([sky_k, [zero_level_k]])
    difference = np.abs(parameters - reference_mean)
    solver_vs_exact_sigma = float(
        np.max(difference / np.maximum(posterior_sigma, 1e-300))
    )
    solver_vs_exact_rel = float(
        difference.max() / max(float(np.max(np.abs(reference_mean))), 1e-300)
    )
    cross_check_max_sigma = float(config.require("tris.solver.cross_check_max_sigma"))
    if bool(config.require("tris.solver.cross_check")) and (
        solver_vs_exact_sigma > cross_check_max_sigma
    ):
        raise RuntimeError(
            "limTOD's solve at {} MHz sits {:.3e} posterior sigma from the "
            "exact Woodbury reference (tolerance {:.1e}); it raised nothing, "
            "which is exactly the silently-wrong case stage 5 exists to catch "
            "-- lower tris.solver.regularization or check the "
            "conditioning".format(
                frequency_mhz, solver_vs_exact_sigma, cross_check_max_sigma,
            )
        )

    implied = float(inputs.implied_monopole_prior_sigma_k(prior_sigma_sky))
    shrinkage = sky_sigma_k / prior_sigma_sky
    truth_k = template_full[inputs.pixel_indices]

    diagnostics: Dict[str, Any] = {
        "n_samples": dof,
        "n_sky_pixels": n_sky,
        "pixel_selection": selection,
        "n_parameters": int(inputs.parameter_count),
        "uncertainty_floor_k": floor,
        "archive_zero_level_uncertainty_k": _describe_zero_level(ring),
        "solver": {
            "backend": "limtod",
            "method": "wiener_filter_map (dense, direct)",
            "regularization": regularization,
            "seconds": solver_seconds,
            "solver_vs_woodbury_sigma": solver_vs_exact_sigma,
            "solver_vs_woodbury_rel": solver_vs_exact_rel,
            "cross_check_max_sigma": cross_check_max_sigma,
        },
        "chi_square": chi_square,
        "degrees_of_freedom": dof,
        "reduced_chi_square": chi_square / max(dof, 1),
        "residual_rms_k": float(residual.std()),
        "zero_level_k": zero_level_k,
        "zero_level_sigma_k": zero_level_sigma_k,
        "zero_level_prior_sigma_k": zero_prior,
        "zero_level_pull_prior_sigma": zero_level_k / zero_prior,
        "monopole_degeneracy": float(inputs.monopole_degeneracy),
        "implied_monopole_prior_sigma_k": implied,
        "monopole_split_is_prior_driven": bool(implied >= 0.3 * zero_prior),
        "beam_coverage_min": float(inputs.beam_coverage.min()),
        "prior_sigma_median_k": float(np.median(prior_sigma_sky)),
        "posterior_sigma_median_k": float(np.median(sky_sigma_k)),
        "sigma_shrinkage_median": float(np.median(shrinkage)),
        "pixels_shrunk_gt_5pct": int((shrinkage < 0.95).sum()),
        "mean_offset_from_template_k": float(np.mean(sky_k - truth_k)),
    }

    if verbose:
        _report(frequency_mhz, effective_mhz, diagnostics)

    return TRISFrequencyMap(
        frequency_mhz=float(frequency_mhz),
        effective_frequency_mhz=effective_mhz,
        nside=nside,
        pixel_indices=np.asarray(inputs.pixel_indices, dtype=np.int64),
        sky_k=sky_k,
        sky_sigma_k=sky_sigma_k,
        prior_mean_k=prior_mean_sky,
        prior_sigma_k=prior_sigma_sky,
        truth_k=truth_k,
        zero_level_k=zero_level_k,
        zero_level_sigma_k=zero_level_sigma_k,
        zero_level_prior_sigma_k=zero_prior,
        ra_deg=np.asarray(ring.ra_deg, dtype=float),
        data_k=np.asarray(inputs.data_k, dtype=float),
        residual_k=residual,
        beam_coverage=np.asarray(inputs.beam_coverage, dtype=float),
        diagnostics=diagnostics,
    )


def _describe_zero_level(ring) -> Any:
    """The archive's own zero level, symmetric at 600 MHz, asymmetric at 820."""
    value = ring.zero_level_uncertainty_k
    if np.isscalar(value):
        return float(value)
    return {"positive_k": float(value.positive_k), "negative_k": float(value.negative_k)}


def _report(frequency_mhz: float, effective_mhz: float, diagnostics: Dict) -> None:
    solver = diagnostics["solver"]
    selection = diagnostics["pixel_selection"]
    print(
        "  {:>4.0f} MHz (nu_eff {:.1f}): {} solve, vs exact {:.2e} sigma "
        "in {:.2f}s".format(
            frequency_mhz,
            effective_mhz,
            solver["backend"],
            solver["solver_vs_woodbury_sigma"],
            solver["seconds"],
        )
    )
    print(
        "            beam response > {:g} of peak: {} pixels ({:.1%} of sky), "
        "{:.4%} of beam power kept".format(
            selection["threshold"],
            selection["n_pixels"],
            selection["sky_fraction"],
            selection["beam_power_retained"],
        )
    )
    print(
        "            reduced chi2 {:.2f}   zero level {:+.4f} +- {:.4f} K "
        "({:+.1f} prior sigma)".format(
            diagnostics["reduced_chi_square"],
            diagnostics["zero_level_k"],
            diagnostics["zero_level_sigma_k"],
            diagnostics["zero_level_pull_prior_sigma"],
        )
    )
    print(
        "            sigma shrinkage median {:.3f}  ({} of {} pixels tightened "
        ">5%)  implied monopole prior {:.4f} K".format(
            diagnostics["sigma_shrinkage_median"],
            diagnostics["pixels_shrunk_gt_5pct"],
            diagnostics["n_sky_pixels"],
            diagnostics["implied_monopole_prior_sigma_k"],
        )
    )


# ---------------------------------------------------------------------------
# stage entry point
# ---------------------------------------------------------------------------
def _check_beam_against_bayesian_skymap(config: Config) -> Dict[str, Any]:
    """Blocker 1: the pinned ``nside_new`` must be what the fork's heuristic gives.

    Stage 1's own working grid is free, but its product has to land on the grid
    stage 3 hands to ``bayesian_skymap``, so the reconciliation belongs here --
    where it can still be fixed -- rather than three stages downstream.
    """
    bayesian_func = import_bayesian_func(config)

    e_plane = float(config.require("beam.e_plane_hpbw_deg"))
    h_plane = float(config.require("beam.h_plane_hpbw_deg"))
    beam_deg = float(config.require("beam.beam_deg"))
    pinned = int(config.require("beam.nside_new"))
    pixels_per_fwhm = float(config.require("beam.pixels_per_fwhm"))

    expected_mean = 0.5 * (e_plane + h_plane)
    if not np.isclose(beam_deg, expected_mean, rtol=0, atol=1e-3):
        raise ConfigError(
            "beam.beam_deg ({}) is not the mean of the E/H HPBWs ({:.4f}); the "
            "config documents it as that mean".format(beam_deg, expected_mean)
        )
    derived = int(bayesian_func.nside_for_beam(beam_deg, pixels_per_fwhm))
    if derived != pinned:
        raise ConfigError(
            "beam.nside_new is pinned at {} but bayesian_func.nside_for_beam("
            "{}) gives {} -- Blocker 1's grid decision has moved under the "
            "pipeline".format(pinned, beam_deg, derived)
        )
    return {
        "beam_deg": beam_deg,
        "nside_new_pinned": pinned,
        "nside_new_from_bayesian_func": derived,
        "pixels_per_fwhm": pixels_per_fwhm,
        "bayesian_func": str(Path(bayesian_func.__file__).resolve()),
    }


def _output_dir(config: Config) -> Path:
    root = Path(str(config.require("run.outputs_dir")))
    if not root.is_absolute():
        root = REPO_ROOT / root
    return root / STAGE


def run_stage1(
    config: Optional[Config] = None,
    *,
    config_path: Optional[Path] = None,
    archive_dir: Optional[Path] = None,
    output_dir: Optional[Path] = None,
    overwrite: Optional[bool] = None,
    verbose: bool = True,
) -> Dict[str, Any]:
    """Run stage 1 end to end and write its products.

    Returns the manifest dict.  Products land in ``outputs/stage1/``:
    ``tris_maps_<run name>.npz`` (the maps) and ``stage1_manifest.json`` (the
    decisions and diagnostics that produced them).
    """
    import healpy as hp

    config = config if config is not None else load_config(config_path)
    archive = Path(archive_dir) if archive_dir else config.resolve_path(
        "paths.tris_archive_dir"
    )
    if not archive.exists():
        raise FileNotFoundError(
            "no TRIS archive at {} -- set paths.tris_archive_dir to the local "
            "copy of the four LAMBDA text files".format(archive)
        )

    destination = Path(output_dir) if output_dir else _output_dir(config)
    destination.mkdir(parents=True, exist_ok=True)
    run_name = str(config.require("run.name"))
    maps_path = destination / "tris_maps_{}.npz".format(run_name)
    manifest_path = destination / "{}_manifest.json".format(STAGE)
    allow_overwrite = (
        bool(config.require("run.overwrite")) if overwrite is None else bool(overwrite)
    )
    for path in (maps_path, manifest_path):
        if path.exists() and not allow_overwrite:
            raise FileExistsError(
                "{} exists; set run.overwrite or pass --overwrite".format(path)
            )

    beam_check = _check_beam_against_bayesian_skymap(config)

    frequencies = [int(f) for f in config.require("tris.frequencies_mhz")]
    mapmaker = str(config.require("tris.mapmaker"))
    if mapmaker != "prior_regularized":
        raise ConfigError(
            "stage 1 implements tris.mapmaker='prior_regularized' (the "
            "gameplan's stage 1.4 decision); 'beamwidth_binning' is "
            "limTOD.tris.tris_beamwidth_binning_map and is not wired up here"
        )
    if verbose:
        print(
            "stage 1: {} at nside {} for {} MHz".format(
                mapmaker, config.require("tris.nside"), frequencies
            )
        )

    results = [
        solve_frequency(
            frequency, config=config, archive_dir=archive, verbose=verbose
        )
        for frequency in frequencies
    ]

    # Both rings share an RA grid and a boresight declination, so they share a
    # pixel band.  Stage 3 stacks them into one (Nfreq, npix) array, so a
    # mismatch here would be a silent misalignment downstream.
    reference = results[0]
    for result in results[1:]:
        if not np.array_equal(result.pixel_indices, reference.pixel_indices):
            raise RuntimeError(
                "the {} MHz and {} MHz rings solved different pixel sets; "
                "stage 3 cannot stack them".format(
                    reference.frequency_mhz, result.frequency_mhz
                )
            )

    nside = reference.nside
    npix = hp.nside2npix(nside)
    band_mask = np.zeros(npix, dtype=bool)
    band_mask[reference.pixel_indices] = True

    payload = {
        "nside": np.int64(nside),
        "pixel_indices": reference.pixel_indices,
        "band_mask": band_mask,
        "freq_mhz": np.array([r.frequency_mhz for r in results]),
        "effective_freq_mhz": np.array([r.effective_frequency_mhz for r in results]),
        "ref_freq_mhz": np.float64(config.require("tris.ref_freq_mhz")),
        "sky_k": np.vstack([r.sky_k for r in results]),
        "sky_sigma_k": np.vstack([r.sky_sigma_k for r in results]),
        "prior_mean_k": np.vstack([r.prior_mean_k for r in results]),
        "prior_sigma_k": np.vstack([r.prior_sigma_k for r in results]),
        "truth_k": np.vstack([r.truth_k for r in results]),
        "zero_level_k": np.array([r.zero_level_k for r in results]),
        "zero_level_sigma_k": np.array([r.zero_level_sigma_k for r in results]),
        "zero_level_prior_sigma_k": np.array(
            [r.zero_level_prior_sigma_k for r in results]
        ),
        "ra_deg": np.vstack([r.ra_deg for r in results]),
        "data_k": np.vstack([r.data_k for r in results]),
        "residual_k": np.vstack([r.residual_k for r in results]),
        "beam_coverage": np.vstack([r.beam_coverage for r in results]),
        # Full-sky maps with NaN outside the solved band.  Stage 3 decides how
        # to fill before degrading to nside_new -- that is its call, not ours.
        "sky_full_k": np.vstack([r.healpix_map() for r in results]),
        "sky_sigma_full_k": np.vstack([r.healpix_sigma() for r in results]),
    }
    np.savez_compressed(maps_path, **payload)

    manifest = {
        "stage": STAGE,
        "run_name": run_name,
        "written_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "config_path": str(config.path) if config.path else None,
        "archive_dir": str(archive),
        "products": {
            "maps": str(maps_path),
            "manifest": str(manifest_path),
        },
        "config_used": {
            "run": config.section("run"),
            "beam": config.section("beam"),
            "tris": config.section("tris"),
        },
        "blocker_1_beam_reconciliation": beam_check,
        "blocker_2_substitute": {
            "prior_template": "gsm2008",
            "truth_template": "gsm2008",
            "prior_equals_truth": True,
            "note": (
                "No true sky exists for real TRIS data and Berkhuijsen (1972) "
                "is not in hand, so GSM2008 is both the prior template and the "
                "reference truth. Agreement between the stage 1 map and "
                "truth_k is therefore circular and is NOT evidence the "
                "pipeline recovered anything. The walkthrough measured GSM2008 "
                "as ~1.9 K cold at 600 MHz, so a large fitted zero level and a "
                "large reduced chi-square are the expected, honest signature "
                "of that template offset -- read the offset as a "
                "model-mismatch amplitude, not as the archive's zero point."
            ),
        },
        "blocker_3_note": (
            "The map solve is the direct Krylov solve (LinearOperator + "
            "{}), never the multi-frequency Woodbury acceleration. The "
            "Woodbury identity used for the posterior variance is exact at a "
            "single frequency and is cross-checked against the Krylov "
            "mean.".format(config.require("tris.solver.method"))
        ),
        "frequencies": {
            "{:.0f}".format(r.frequency_mhz): {
                "effective_frequency_mhz": r.effective_frequency_mhz,
                **r.diagnostics,
            }
            for r in results
        },
        "environment": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "healpy": hp.__version__,
            "limtod": str(Path(import_limtod_tris(config).__file__).resolve()),
        },
    }
    with open(manifest_path, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=False)

    if verbose:
        print("stage 1: wrote {}".format(maps_path))
        print("stage 1: wrote {}".format(manifest_path))
    return manifest


# ---------------------------------------------------------------------------
# reading stage 1 back (notebooks, stage 3, stage 5)
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Stage1Products:
    """Everything stage 1 wrote, keyed the way the notebooks want it.

    ``arrays`` holds the npz contents (stacked over frequency, in
    ``freq_mhz`` order); ``manifest`` holds the decisions and diagnostics.
    Use :meth:`index` rather than assuming an order.
    """

    arrays: Dict[str, np.ndarray]
    manifest: Dict[str, Any]
    maps_path: Path
    manifest_path: Path

    def __getitem__(self, key: str) -> np.ndarray:
        return self.arrays[key]

    @property
    def nside(self) -> int:
        return int(self.arrays["nside"])

    @property
    def frequencies_mhz(self) -> np.ndarray:
        return self.arrays["freq_mhz"]

    @property
    def pixel_indices(self) -> np.ndarray:
        return self.arrays["pixel_indices"]

    def index(self, frequency_mhz: float) -> int:
        """Row of the stacked arrays for a nominal frequency."""
        matches = np.flatnonzero(
            np.isclose(self.arrays["freq_mhz"], float(frequency_mhz))
        )
        if matches.size != 1:
            raise KeyError(
                "{} MHz is not one of {}".format(
                    frequency_mhz, list(self.arrays["freq_mhz"])
                )
            )
        return int(matches[0])

    def diagnostics(self, frequency_mhz: float) -> Dict[str, Any]:
        return self.manifest["frequencies"]["{:.0f}".format(float(frequency_mhz))]

    def band_map(self, key: str, frequency_mhz: float, fill: float = np.nan):
        """Scatter one band array back onto a full-sky RING map."""
        import healpy as hp

        full = np.full(hp.nside2npix(self.nside), fill, dtype=float)
        full[self.pixel_indices] = self.arrays[key][self.index(frequency_mhz)]
        return full


def stage1_product_paths(config: Config) -> Tuple[Path, Path]:
    """``(maps npz, manifest json)`` for this run, whether or not they exist."""
    destination = _output_dir(config)
    run_name = str(config.require("run.name"))
    return (
        destination / "tris_maps_{}.npz".format(run_name),
        destination / "{}_manifest.json".format(STAGE),
    )


def load_stage1(
    config: Optional[Config] = None,
    *,
    config_path: Optional[Path] = None,
    maps_path: Optional[Path] = None,
    manifest_path: Optional[Path] = None,
) -> Stage1Products:
    """Load stage 1's products, defaulting to this run's own output paths."""
    config = config if config is not None else load_config(config_path)
    default_maps, default_manifest = stage1_product_paths(config)
    maps_path = Path(maps_path) if maps_path else default_maps
    manifest_path = Path(manifest_path) if manifest_path else default_manifest
    if not maps_path.exists():
        raise FileNotFoundError(
            "no stage 1 maps at {} -- run `python run_pipeline.py --stage 1` "
            "first".format(maps_path)
        )
    with np.load(maps_path) as bundle:
        arrays = {key: bundle[key] for key in bundle.files}
    with open(manifest_path, "r", encoding="utf-8") as handle:
        manifest = json.load(handle)
    return Stage1Products(
        arrays=arrays,
        manifest=manifest,
        maps_path=maps_path,
        manifest_path=manifest_path,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="stage1_tris_maps",
        description="Stage 1 -- prior-regularized TRIS maps from the public archive.",
    )
    parser.add_argument("--config", type=Path, default=None, help="run config path")
    parser.add_argument(
        "--archive-dir",
        type=Path,
        default=None,
        help="override paths.tris_archive_dir",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=None, help="override outputs/stage1"
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="replace existing stage 1 products"
    )
    parser.add_argument("--quiet", action="store_true", help="suppress progress output")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    run_stage1(
        config_path=args.config,
        archive_dir=args.archive_dir,
        output_dir=args.output_dir,
        overwrite=True if args.overwrite else None,
        verbose=not args.quiet,
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
