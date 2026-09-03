"""Stage 2 -- prepare Haslam, the gain operator, and the prior covariance.

What this stage produces
------------------------
The three things ``bayesian_skymap`` needs on the miscalibrated side of the
problem, all on its fixed ``nside = 128`` grid and all in **galactic**
coordinates (the decision recorded in ``haslam.frame``):

``sky_gain``
    The reprocessed Haslam 408 MHz map -- destriped and desourced, Remazeilles
    et al. (2015) -- degraded to nside 128 with the CMB monopole removed.
    ``bayesian_skymap`` models this as ``gain * s * (freq/ref)**beta``, a pure
    power law with no monopole, so the monopole has to come off first or the
    forward model is wrong by 2.7 K everywhere.
``operator`` / ``operator_spect``
    The gain-parameter -> pixel map and the spectral-index-region -> pixel map,
    both ``(npix, nregions)`` 0/1 indicator matrices built the way
    ``bayesian_skymap``'s own generator builds them.  12 gain regions and 6 beta
    regions: the gain deliberately has more freedom than beta.
``covA``
    The per-pixel prior *variance* on the sky amplitudes -- the substitute for
    the truth-informed ``covA = (0.1 * true_s)**2`` the paper's scripts
    hard-code (Blocker 2).

Blocker 2, and why the default prior is circular
------------------------------------------------
The prior the pipeline document wants is Berkhuijsen (1972) at 820 MHz, because
it is the one basis independent of Haslam itself.  It is not in hand, and this
repo never touches the network.  Three sources are implemented and selected by
``haslam.prior.source``:

``gsm2008`` (default)
    **Not Haslam-independent.**  GSM2008 is built from a set of surveys that
    includes Haslam 408 MHz, and at 408 MHz its reconstruction is dominated by
    it -- ``notebooks/02`` measures a log-space correlation of 0.993 and a
    median ratio of 1.005 -- so a GSM2008 prior on the true sky is a Haslam
    prior in all but name.  This stage sets
    ``prior_is_haslam_derived`` in its manifest so that nothing downstream can
    read a recovered gain built on it as an independent measurement.  It exists
    to let stages 3-5 be exercised end to end.
``tris_stage1``
    Stage 1's TRIS maps extrapolated to 408 MHz.  Genuinely independent -- TRIS
    is an absolute radiometer -- but only across the declination band it
    observed; outside it the width falls back to something deliberately loose.
``berkhuijsen_1972``
    The real answer.  Raises until ``paths.berkhuijsen_map`` points at a file.

One thing worth knowing before reading too much into the prior *mean*:
``gibbs_*.py`` sets ``sprior_mean = 0`` and only ever consumes ``covA``.  The
template therefore sets the prior's **scale**, not its centre.  Stage 2 emits
the mean anyway (stage 4 may want to wire it in), but ``covA`` is the operative
product.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional, Sequence, Tuple

import numpy as np

from .config import REPO_ROOT, Config, ConfigError, load_config
from .stage1_tris_maps import (
    import_bayesian_func,
    import_limtod_tris,
    load_stage1,
)

STAGE = "stage2"

#: The frequency the Haslam map is at, and the frequency every prior source is
#: extrapolated to.  Not configurable: it is a property of the map.
HASLAM_FREQ_MHZ = 408.0

#: pygdsm's cached copy of the destriped/desourced Remazeilles (2015) map.
HASLAM_URL = (
    "https://lambda.gsfc.nasa.gov/data/foregrounds/haslam_2014/"
    "haslam408_dsds_Remazeilles2014.fits"
)

#: Why the default prior is not the independent basis Blocker 2 asks for.
#: Measured in ``notebooks/02``, not asserted.
_GSM2008_CIRCULARITY = (
    "GSM2008 is built from a set of surveys that includes Haslam 408 MHz, and "
    "at 408 MHz its reconstruction is dominated by it (log-space r = 0.993, "
    "median ratio 1.005), so this prior is NOT Haslam-independent"
)

_PRIOR_SOURCES = ("gsm2008", "tris_stage1", "berkhuijsen_1972")
_REGION_GEOMETRIES = ("theta_bands", "equal_area")
_S0_INITS = ("haslam", "tris_extrap")


# ---------------------------------------------------------------------------
# region operators
# ---------------------------------------------------------------------------
def region_operator(
    nside: int, nregions: int, geometry: str = "theta_bands"
) -> np.ndarray:
    """``(npix, nregions)`` 0/1 indicator matrix of latitude bands.

    ``theta_bands`` reproduces ``bayesian_skymap``'s own
    ``tests/make_test_data.region_operator`` exactly -- equal-spaced in
    colatitude, contiguous in galactic latitude -- and the test suite pins the
    two against each other.  ``equal_area`` splits at equal pixel counts
    instead, so no gain parameter ends up carried by a tiny polar band.
    """
    import healpy as hp

    if geometry not in _REGION_GEOMETRIES:
        raise ValueError(
            "region geometry must be one of {}, got {!r}".format(
                _REGION_GEOMETRIES, geometry
            )
        )
    nregions = int(nregions)
    if nregions < 1:
        raise ValueError("nregions must be >= 1")

    npix = hp.nside2npix(nside)
    theta, _phi = hp.pix2ang(nside, np.arange(npix))
    if geometry == "theta_bands":
        edges = np.linspace(0, np.pi, nregions + 1)
    else:
        # Equal area on the sphere is equal spacing in cos(theta).
        edges = np.arccos(np.linspace(1.0, -1.0, nregions + 1))
    index = np.clip(np.digitize(theta, edges[1:-1]), 0, nregions - 1)
    operator = np.zeros((npix, nregions))
    operator[np.arange(npix), index] = 1.0
    return operator


def region_index(operator: np.ndarray) -> np.ndarray:
    """Which region each pixel belongs to, for plotting and diagnostics."""
    return np.argmax(operator, axis=1).astype(np.int16)


# ---------------------------------------------------------------------------
# the maps
# ---------------------------------------------------------------------------
def load_haslam(config: Config) -> Tuple[np.ndarray, Dict[str, Any]]:
    """The reprocessed Haslam 408 MHz map at ``haslam.nside``, galactic, in K.

    ``paths.haslam_map`` overrides; ``null`` means pygdsm's cached copy of the
    destriped/desourced Remazeilles (2015) map, which is the reprocessed map
    the pipeline document asks for.  Nothing here reaches the network: the
    astropy cache is consulted, and a cache miss is reported as such.
    """
    import healpy as hp

    frame = str(config.require("haslam.frame")).lower()
    if frame != "galactic":
        raise ConfigError(
            "stage 2 implements haslam.frame = 'galactic' (the recorded "
            "decision -- Haslam is native there and the regions are galactic "
            "latitude bands); got {!r}".format(frame)
        )

    override = config.get_path("paths.haslam_map")
    if override is not None:
        path = Path(str(override))
        if not path.is_absolute():
            path = REPO_ROOT / path
        source = "paths.haslam_map"
    else:
        from astropy.utils.data import download_file

        try:
            path = Path(download_file(HASLAM_URL, cache=True, show_progress=False))
        except Exception as error:  # pragma: no cover - needs an empty cache
            raise FileNotFoundError(
                "the reprocessed Haslam map is not in the astropy cache and "
                "could not be fetched ({}). Download {} once, or point "
                "paths.haslam_map at a local copy.".format(error, HASLAM_URL)
            )
        source = "pygdsm astropy cache (Remazeilles 2015 dsds)"

    native = hp.read_map(path, dtype=np.float64)
    native_nside = hp.get_nside(native)
    nside = int(config.require("haslam.nside"))
    sky = hp.ud_grade(native, nside) if native_nside != nside else native.copy()

    monopole_k = 0.0
    if bool(config.require("haslam.remove_cmb_monopole")):
        tris = import_limtod_tris(config)
        monopole_k = float(tris.cmb_monopole_rj_k(HASLAM_FREQ_MHZ))
        sky = sky - monopole_k

    if np.any(sky <= 0):
        raise RuntimeError(
            "{} Haslam pixels are non-positive after monopole removal; covA is "
            "inverted per pixel, so a zero there would blow up the "
            "solve".format(int((sky <= 0).sum()))
        )

    meta = {
        "path": str(path),
        "source": source,
        "native_nside": int(native_nside),
        "nside": nside,
        "frame": frame,
        "frequency_mhz": HASLAM_FREQ_MHZ,
        "cmb_monopole_removed_k": monopole_k,
        "min_k": float(sky.min()),
        "median_k": float(np.median(sky)),
        "max_k": float(sky.max()),
    }
    return sky, meta


def _gsm2008_galactic(frequency_mhz: float, nside: int) -> np.ndarray:
    """GSM2008 at one frequency, galactic RING, Galactic-only (no monopole)."""
    import healpy as hp
    from pygdsm import GlobalSkyModel

    return hp.ud_grade(GlobalSkyModel().generate(float(frequency_mhz)), nside)


def _tris_at_408(config: Config, nside: int) -> Tuple[np.ndarray, np.ndarray]:
    """Stage 1's TRIS map extrapolated to 408 MHz, resampled into galactic.

    Returns ``(template, covered)``: the template is NaN where the TRIS band
    never reached, and ``covered`` is the boolean galactic footprint of the
    band.  Resampling is nearest-neighbour by construction -- rotate each
    galactic pixel centre into equatorial and look up the stage 1 pixel it
    lands in -- so no interpolation smears the band edge.
    """
    import healpy as hp

    tris = import_limtod_tris(config)
    products = load_stage1(config)
    ref_freq = float(config.require("tris.ref_freq_mhz"))
    row = products.index(ref_freq)
    effective = float(products["effective_freq_mhz"][row])

    # TRIS temperatures include the CMB monopole; the Haslam side is
    # Galactic-only, so it has to come off before the power-law extrapolation.
    band = products["sky_full_k"][row] - tris.cmb_monopole_rj_k(effective)
    beta = float(config.require("haslam.prior.extrapolation_beta"))
    band = band * (HASLAM_FREQ_MHZ / effective) ** beta

    theta, phi = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
    theta_eq, phi_eq = hp.Rotator(coord=["G", "C"])(theta, phi)
    source_pixels = hp.ang2pix(products.nside, theta_eq, phi_eq)
    template = band[source_pixels]
    covered = np.isfinite(template) & (template > 0)
    return template, covered


def build_prior(
    config: Config, haslam_k: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, Dict[str, Any]]:
    """``(prior_mean_k, prior_sigma_k, meta)`` for the sky amplitudes.

    ``covA = prior_sigma_k**2`` is what ``bayesian_skymap`` consumes.
    """
    import healpy as hp

    source = str(config.require("haslam.prior.source")).lower()
    if source not in _PRIOR_SOURCES:
        raise ConfigError(
            "haslam.prior.source must be one of {}, got {!r}".format(
                _PRIOR_SOURCES, source
            )
        )
    nside = int(config.require("haslam.nside"))
    frac = float(config.require("haslam.prior.frac_sigma"))
    floor = float(config.require("haslam.prior.floor_k"))

    meta: Dict[str, Any] = {
        "source": source,
        "frac_sigma": frac,
        "floor_k": floor,
        "extrapolated_to_mhz": HASLAM_FREQ_MHZ,
    }

    if source == "gsm2008":
        template = _gsm2008_galactic(HASLAM_FREQ_MHZ, nside)
        sigma = np.maximum(frac * np.abs(template), floor)
        meta.update(
            template_freq_mhz=HASLAM_FREQ_MHZ,
            extrapolation_beta=None,
            haslam_derived=True,
            note=(
                "{}. A gain recovered against it is circular. Berkhuijsen "
                "(1972) is the real answer.".format(_GSM2008_CIRCULARITY)
            ),
        )

    elif source == "berkhuijsen_1972":
        raw = config.get_path("paths.berkhuijsen_map")
        if raw is None:
            raise ConfigError(
                "haslam.prior.source is 'berkhuijsen_1972' but "
                "paths.berkhuijsen_map is null -- supply the 820 MHz survey "
                "map, or pick another source (see the config's notes on why "
                "gsm2008 is not a substitute for it)"
            )
        path = Path(str(raw))
        if not path.is_absolute():
            path = REPO_ROOT / path
        if not path.exists():
            raise FileNotFoundError("no Berkhuijsen map at {}".format(path))
        survey_freq = float(config.require("haslam.prior.survey_freq_mhz"))
        beta = float(config.require("haslam.prior.extrapolation_beta"))
        native = hp.read_map(path, dtype=np.float64)
        template = hp.ud_grade(native, nside) * (
            HASLAM_FREQ_MHZ / survey_freq
        ) ** beta
        sigma = np.maximum(frac * np.abs(template), floor)
        meta.update(
            path=str(path),
            template_freq_mhz=survey_freq,
            extrapolation_beta=beta,
            haslam_derived=False,
            note="Haslam-independent: an external survey, extrapolated in "
            "frequency only.",
        )

    else:  # tris_stage1
        wide = float(config.require("haslam.prior.out_of_band_frac_sigma"))
        beta = float(config.require("haslam.prior.extrapolation_beta"))
        template, covered = _tris_at_408(config, nside)
        # Outside the band TRIS says nothing, so the only scale available is
        # Haslam's own -- taken deliberately loose so it barely constrains.
        template = np.where(covered, template, haslam_k)
        sigma = np.where(
            covered,
            np.maximum(frac * np.abs(template), floor),
            np.maximum(wide * np.abs(haslam_k), floor),
        )
        meta.update(
            template_freq_mhz=float(config.require("tris.ref_freq_mhz")),
            extrapolation_beta=beta,
            out_of_band_frac_sigma=wide,
            band_pixels=int(covered.sum()),
            band_sky_fraction=float(covered.mean()),
            haslam_derived=False,
            note=(
                "Haslam-independent inside the TRIS declination band "
                "({:.0%} of the sky); outside it the width falls back to "
                "{:g}x the local Haslam value, which is deliberately "
                "uninformative.".format(covered.mean(), wide)
            ),
        )
        meta["covered"] = covered

    covered = meta.pop("covered", None)
    prior_mean = np.asarray(template, dtype=float)
    prior_sigma = np.asarray(sigma, dtype=float)
    if not np.all(np.isfinite(prior_sigma)) or np.any(prior_sigma <= 0):
        raise RuntimeError("the prior sigma is not strictly positive everywhere")
    meta.update(
        prior_sigma_median_k=float(np.median(prior_sigma)),
        prior_sigma_min_k=float(prior_sigma.min()),
        prior_sigma_max_k=float(prior_sigma.max()),
    )
    if covered is not None:
        return prior_mean, prior_sigma, {**meta, "_covered": covered}
    return prior_mean, prior_sigma, meta


def build_s0_init(
    config: Config, haslam_k: np.ndarray, prior_mean_k: np.ndarray
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """Blocker 2's substitute for ``nstart=0 -> initialise s0 at the true sky``.

    There is no true sky for real data.  ``haslam`` initialises at Haslam
    itself, which is crude but not circular *for an initialisation* -- only for
    a prior.  ``tris_extrap`` uses the TRIS extrapolation instead.
    """
    choice = str(config.require("haslam.s0_init")).lower()
    if choice not in _S0_INITS:
        raise ConfigError(
            "haslam.s0_init must be one of {}, got {!r}".format(_S0_INITS, choice)
        )
    if choice == "haslam":
        return haslam_k.copy(), {
            "s0_init": choice,
            "note": "Haslam itself. Crude, and circular only if read as a "
            "prior rather than a starting point.",
        }
    nside = int(config.require("haslam.nside"))
    template, covered = _tris_at_408(config, nside)
    initial = np.where(covered, template, haslam_k)
    return initial, {
        "s0_init": choice,
        "band_pixels": int(covered.sum()),
        "note": "TRIS extrapolated to 408 MHz inside its band, Haslam "
        "elsewhere.",
    }


# ---------------------------------------------------------------------------
# stage entry point
# ---------------------------------------------------------------------------
def _output_dir(config: Config) -> Path:
    root = Path(str(config.require("run.outputs_dir")))
    if not root.is_absolute():
        root = REPO_ROOT / root
    return root / STAGE


def stage2_product_paths(config: Config) -> Tuple[Path, Path]:
    """``(products npz, manifest json)`` for this run, whether or not they exist."""
    destination = _output_dir(config)
    run_name = str(config.require("run.name"))
    return (
        destination / "haslam_prep_{}.npz".format(run_name),
        destination / "{}_manifest.json".format(STAGE),
    )


def run_stage2(
    config: Optional[Config] = None,
    *,
    config_path: Optional[Path] = None,
    output_dir: Optional[Path] = None,
    overwrite: Optional[bool] = None,
    verbose: bool = True,
) -> Dict[str, Any]:
    """Run stage 2 end to end and write its products.  Returns the manifest."""
    import healpy as hp

    config = config if config is not None else load_config(config_path)

    if int(config.require("haslam.op_alm")) != 0:
        raise ConfigError(
            "stage 2 implements the region-based gain parametrisation "
            "(haslam.op_alm = 0, the recorded decision). op_alm = 1 builds no "
            "operator matrix at all -- bayesian_skymap synthesises the gain "
            "from alms internally -- so there is nothing for this stage to "
            "produce."
        )

    destination = Path(output_dir) if output_dir else _output_dir(config)
    destination.mkdir(parents=True, exist_ok=True)
    products_path, manifest_path = stage2_product_paths(config)
    if output_dir:
        products_path = destination / products_path.name
        manifest_path = destination / manifest_path.name
    allow_overwrite = (
        bool(config.require("run.overwrite")) if overwrite is None else bool(overwrite)
    )
    for path in (products_path, manifest_path):
        if path.exists() and not allow_overwrite:
            raise FileExistsError(
                "{} exists; set run.overwrite or pass --overwrite".format(path)
            )

    nside = int(config.require("haslam.nside"))
    geometry = str(config.require("haslam.region_geometry"))
    n_gain = int(config.require("haslam.n_gain_regions"))
    n_beta = int(config.require("haslam.n_beta_regions"))
    beta_init = np.asarray(config.require("assemble.beta_init"), dtype=float)
    if beta_init.size != n_beta:
        raise ConfigError(
            "assemble.beta_init has {} entries but haslam.n_beta_regions is "
            "{} -- one initial index per spectral region".format(
                beta_init.size, n_beta
            )
        )

    if verbose:
        print(
            "stage 2: Haslam at nside {} ({}), {} gain / {} beta {} regions".format(
                nside, config.require("haslam.frame"), n_gain, n_beta, geometry
            )
        )

    haslam_k, haslam_meta = load_haslam(config)
    if verbose:
        print(
            "  Haslam: {} -> nside {}, CMB monopole {:.4f} K removed, "
            "{:.2f} .. {:.1f} K".format(
                haslam_meta["native_nside"],
                nside,
                haslam_meta["cmb_monopole_removed_k"],
                haslam_meta["min_k"],
                haslam_meta["max_k"],
            )
        )

    operator = region_operator(nside, n_gain, geometry)
    operator_spect = region_operator(nside, n_beta, geometry)
    beta_map = operator_spect @ beta_init

    prior_mean_k, prior_sigma_k, prior_meta = build_prior(config, haslam_k)
    prior_meta.pop("_covered", None)
    cov_a = prior_sigma_k**2

    gain_frac = float(config.require("haslam.prior.gain_frac_sigma"))
    cov_g_diag = np.full(n_gain, gain_frac**2)

    s0_init_k, s0_meta = build_s0_init(config, haslam_k, prior_mean_k)

    if verbose:
        print(
            "  prior '{}': sigma median {:.3f} K  ({})".format(
                prior_meta["source"],
                prior_meta["prior_sigma_median_k"],
                "HASLAM-DERIVED, not independent"
                if prior_meta.get("haslam_derived")
                else "Haslam-independent",
            )
        )

    # A region that holds no pixels would be an unconstrained gain parameter.
    gain_counts = operator.sum(axis=0).astype(int)
    beta_counts = operator_spect.sum(axis=0).astype(int)
    if gain_counts.min() == 0 or beta_counts.min() == 0:
        raise RuntimeError("a region holds no pixels; reduce the region count")

    payload = {
        "nside": np.int64(nside),
        "frame": str(config.require("haslam.frame")),
        "haslam_freq_mhz": np.float64(HASLAM_FREQ_MHZ),
        "sky_gain": haslam_k,
        "operator": operator,
        "operator_spect": operator_spect,
        "gain_region_index": region_index(operator),
        "beta_region_index": region_index(operator_spect),
        "gain_region_pixels": gain_counts,
        "beta_region_pixels": beta_counts,
        "covA": cov_a,
        "prior_mean_k": prior_mean_k,
        "prior_sigma_k": prior_sigma_k,
        "covG_diag": cov_g_diag,
        "s0_init_k": s0_init_k,
        "beta_init": beta_init,
        "beta_init_map": beta_map,
    }
    np.savez_compressed(products_path, **payload)

    bayesian_func = import_bayesian_func(config)
    manifest = {
        "stage": STAGE,
        "run_name": str(config.require("run.name")),
        "written_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "config_path": str(config.path) if config.path else None,
        "products": {"maps": str(products_path), "manifest": str(manifest_path)},
        "config_used": {"haslam": config.section("haslam")},
        "haslam": haslam_meta,
        "regions": {
            "geometry": geometry,
            "n_gain_regions": n_gain,
            "n_beta_regions": n_beta,
            "gain_region_pixels": gain_counts.tolist(),
            "beta_region_pixels": beta_counts.tolist(),
            "note": (
                "12 gain / 6 beta regions follow bayesian_skymap's own "
                "generator (tests/make_test_data.py). The six regions that "
                "repository initialises are the SPECTRAL INDEX regions, not "
                "the gain regions."
            ),
        },
        "prior": prior_meta,
        "prior_is_haslam_derived": bool(prior_meta.get("haslam_derived", False)),
        "gain_prior": {
            "covG_diag": cov_g_diag.tolist(),
            "gain_frac_sigma": gain_frac,
            "note": (
                "gibbs_*.py builds covG = (0.15 * (1 + g_true))**2 from the "
                "TRUE gain, which does not exist for real data. This is the "
                "same width evaluated at g = 0. Stage 4 has to wire it in; "
                "the script will otherwise reach for a truth it does not have."
            ),
        },
        "s0_init": s0_meta,
        "blocker_2_status": (
            "UNRESOLVED. The prior source is '{}'. {}".format(
                prior_meta["source"], prior_meta.get("note", "")
            )
        ),
        "environment": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "healpy": hp.__version__,
            "bayesian_func": str(Path(bayesian_func.__file__).resolve()),
        },
    }
    with open(manifest_path, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=False)

    if verbose:
        print("stage 2: wrote {}".format(products_path))
        print("stage 2: wrote {}".format(manifest_path))
        if manifest["prior_is_haslam_derived"]:
            print(
                "stage 2: WARNING -- the prior is Haslam-derived, so a gain "
                "recovered from it is circular. See blocker_2_status."
            )
    return manifest


# ---------------------------------------------------------------------------
# reading stage 2 back
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Stage2Products:
    """Everything stage 2 wrote, plus its manifest."""

    arrays: Dict[str, np.ndarray]
    manifest: Dict[str, Any]
    products_path: Path
    manifest_path: Path

    def __getitem__(self, key: str) -> np.ndarray:
        return self.arrays[key]

    @property
    def nside(self) -> int:
        return int(self.arrays["nside"])

    @property
    def frame(self) -> str:
        return str(self.arrays["frame"])

    @property
    def prior_is_haslam_derived(self) -> bool:
        return bool(self.manifest["prior_is_haslam_derived"])


def load_stage2(
    config: Optional[Config] = None,
    *,
    config_path: Optional[Path] = None,
    products_path: Optional[Path] = None,
    manifest_path: Optional[Path] = None,
) -> Stage2Products:
    """Load stage 2's products, defaulting to this run's own output paths."""
    config = config if config is not None else load_config(config_path)
    default_products, default_manifest = stage2_product_paths(config)
    products_path = Path(products_path) if products_path else default_products
    manifest_path = Path(manifest_path) if manifest_path else default_manifest
    if not products_path.exists():
        raise FileNotFoundError(
            "no stage 2 products at {} -- run `python run_pipeline.py "
            "--stage 2` first".format(products_path)
        )
    with np.load(products_path, allow_pickle=False) as bundle:
        arrays = {key: bundle[key] for key in bundle.files}
    with open(manifest_path, "r", encoding="utf-8") as handle:
        manifest = json.load(handle)
    return Stage2Products(
        arrays=arrays,
        manifest=manifest,
        products_path=products_path,
        manifest_path=manifest_path,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="stage2_haslam_prep",
        description="Stage 2 -- Haslam, the gain operator, and the prior covariance.",
    )
    parser.add_argument("--config", type=Path, default=None, help="run config path")
    parser.add_argument(
        "--output-dir", type=Path, default=None, help="override outputs/stage2"
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="replace existing stage 2 products"
    )
    parser.add_argument("--quiet", action="store_true", help="suppress progress output")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    run_stage2(
        config_path=args.config,
        output_dir=args.output_dir,
        overwrite=True if args.overwrite else None,
        verbose=not args.quiet,
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
