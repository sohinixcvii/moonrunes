"""Stage 2 -- beam-match the three datasets onto one common resolution.

What this stage produces
------------------------
Haslam and both ARCADE 2 bands, smoothed so that every map carries the same
effective beam as TRIS, which is the coarsest of the three.  Stage 3 assembles
per-pixel SEDs from them, and an SED is only meaningful if each point in it was
measured through the same beam -- otherwise the spectral index picks up the
resolution difference between surveys as if it were sky.

This is Step 0 of ``tris_haslam_pipeline.md``, applied.

The target FWHM (Step 0.1)
--------------------------
``hp.smoothing`` takes one symmetric FWHM, so TRIS's elliptical beam
(19.155 deg E / 23.366 deg H) cannot be carried through exactly.  The target is
the **larger** axis: that guarantees the other two surveys end up at least as
coarse as TRIS along both of its axes, so nothing survives under-smoothed in
the narrow direction.  A proper elliptical convolution is more correct and
substantially more work; the pipeline document flags this as a deliberate
simplification, not an oversight.

The quadrature correction (Step 0.2)
------------------------------------
Gaussian beams add in quadrature, so a map already at ``native`` needs an
*extra* kernel of ``sqrt(target^2 - native^2)`` -- not ``target``.  Smoothing
straight to the target over-smooths, by 3.5% for ARCADE 2 here and by 0.1% for
Haslam.

Masks before smoothing (Step 0.3)
---------------------------------
``hp.smoothing`` works through a global spherical harmonic transform and knows
nothing about masks, so unobserved pixels bleed across the mask boundary into
real data.  The fix is the standard one: zero the unobserved pixels, smooth a
binary weight map with the same kernel, divide, and drop any pixel whose
smoothed weight falls below a threshold -- those are edge pixels the division
cannot rescue.

Each dataset marks unobserved sky differently and none of them uses
``hp.UNSEEN`` on disk (measured in ``notebooks/01_explore_data.ipynb``):
ARCADE 2 uses exact zeros, Haslam is full-sky, and stage 1's TRIS maps use NaN
alongside an explicit ``band_mask``.  :func:`observed_mask` is the one place
that knows this.

What this stage does NOT do yet
-------------------------------
Step 0.4-0.5 -- regridding everything onto the common ``nside`` grid, and the
coordinate-frame reconciliation that has to happen with it.  Haslam is
galactic, stage 1's TRIS maps are equatorial, and ARCADE 2's headers do not say
(no ``COORDSYS``), which is the open TODO in the pipeline document.  Smoothing
is frame-independent, so this stage is correct as far as it goes; assembling
SEDs across these maps is not, until that is resolved.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional

import numpy as np

from .config import Config, load_config

STAGE = "stage2"

# TRIS's real, asymmetric half-power beam widths, from the archive's own E/H
# principal-plane cuts (``TRIS_Beam_Profile.txt``); the ring header's "18 deg"
# is 6% off and is not what limTOD measures.  Section 01 of
# ``notebooks/01_explore_data.ipynb`` re-derives both from the cuts.
TRIS_HPBW_E_DEG = 19.155
TRIS_HPBW_H_DEG = 23.366

#: Step 0.1 -- the common target, the larger of TRIS's two axes.  The live
#: value is ``beam_matching.target_fwhm_deg``; this is what that key is set to
#: and what the pure-algebra helpers default to.
TARGET_FWHM_DEG = max(TRIS_HPBW_E_DEG, TRIS_HPBW_H_DEG)

#: Below this smoothed weight a pixel is mask boundary, not data (Step 0.3).
#: The live value is ``beam_matching.weight_floor``.
DEFAULT_WEIGHT_FLOOR = 0.5


# ---------------------------------------------------------------------------
# the beam algebra
# ---------------------------------------------------------------------------
def extra_fwhm_deg(native_fwhm_deg: float, target_fwhm_deg: float) -> float:
    """The kernel that takes ``native`` to ``target``, in degrees.

    Gaussian beams add in quadrature.  Raises rather than returning NaN when
    the map is already coarser than the target: that is a real error in the
    caller's resolution bookkeeping, and ``np.sqrt`` of a negative would
    otherwise propagate silently through ``hp.smoothing`` into an all-NaN map.
    """
    if target_fwhm_deg < native_fwhm_deg:
        raise ValueError(
            "cannot smooth a {:.4f} deg beam to {:.4f} deg -- the target must "
            "be the coarsest resolution in play".format(
                native_fwhm_deg, target_fwhm_deg
            )
        )
    return float(np.sqrt(target_fwhm_deg**2 - native_fwhm_deg**2))


def smooth_to_target(
    sky: np.ndarray,
    native_fwhm_deg: float,
    target_fwhm_deg: float = TARGET_FWHM_DEG,
    *,
    weight_mask: Optional[np.ndarray] = None,
    weight_floor: float = DEFAULT_WEIGHT_FLOOR,
    smoothing_iter: int = 3,
) -> np.ndarray:
    """Smooth ``sky`` from its native resolution down to ``target_fwhm_deg``.

    ``weight_mask`` is True where the map was actually observed.  Pass it for
    any partial-sky map: without it the unobserved pixels are convolved into
    the real ones and the result is wrong near every mask edge (Step 0.3).
    Pixels whose smoothed weight lands below ``weight_floor`` come back as
    ``hp.UNSEEN``.

    ``smoothing_iter`` is healpy's own default of 3 map2alm iterations, which
    is the accurate choice.  The ``iter=0`` this repo pins elsewhere
    (``gibbs.smoothing_iter``) is a Gibbs-solver requirement -- it keeps the
    forward operator symmetric so the Krylov solve does not diverge -- and does
    not transfer to a one-shot smoothing like this one.

    The input is never modified.
    """
    import healpy as hp

    fwhm_rad = np.radians(extra_fwhm_deg(native_fwhm_deg, target_fwhm_deg))

    if weight_mask is None:
        return hp.smoothing(sky, fwhm=fwhm_rad, iter=smoothing_iter)

    weight_mask = np.asarray(weight_mask, dtype=bool)
    if weight_mask.shape != sky.shape:
        raise ValueError(
            "weight_mask has shape {} but the map has shape {}".format(
                weight_mask.shape, sky.shape
            )
        )

    # Zero rather than UNSEEN: the sentinel is -1.6e30 and would dominate the
    # transform.  The division by the smoothed weight is what undoes the
    # dilution this introduces.
    zeroed = np.where(weight_mask, sky, 0.0)
    smoothed = hp.smoothing(zeroed, fwhm=fwhm_rad, iter=smoothing_iter)
    weight = hp.smoothing(weight_mask.astype(float), fwhm=fwhm_rad,
                          iter=smoothing_iter)

    with np.errstate(invalid="ignore", divide="ignore"):
        corrected = smoothed / weight
    corrected[weight < weight_floor] = hp.UNSEEN
    return corrected


def observed_mask(sky: np.ndarray, convention: str) -> np.ndarray:
    """True where a map was observed, for each dataset's own convention.

    ``zeros``   ARCADE 2 -- unobserved sky is exact 0.0, with no sentinel.
    ``nan``     stage 1's TRIS maps -- NaN off the declination band.
    ``full``    Haslam -- all sky observed, so everything is True.

    ``hp.UNSEEN`` and NaN are unobserved under every convention -- only the
    extra rule differs -- so a map that has already been through
    :func:`smooth_to_target` masks correctly whatever it started as.
    """
    import healpy as hp

    sky = np.asarray(sky)
    # A NaN or a sentinel is unobserved whatever the dataset's own convention
    # is; only the extra rule differs between them.
    mask = np.isfinite(sky) & (sky != hp.UNSEEN)
    if convention == "zeros":
        return mask & (sky != 0.0)
    if convention in ("nan", "full"):
        return mask
    raise ValueError(
        "unknown mask convention {!r} -- expected 'zeros', 'nan' or "
        "'full'".format(convention)
    )


# ---------------------------------------------------------------------------
# the inputs
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Dataset:
    """One map to beam-match, and everything needed to do it correctly."""

    name: str
    path_key: str          # dotted config key holding the file path
    fwhm_key: str          # dotted config key holding the native beam FWHM
    mask_convention: str   # see observed_mask
    frame: str             # as the file states it, or "unstated"
    fwhm_source: str       # where that FWHM came from


DATASETS = (
    Dataset("haslam_408", "paths.haslam_map_fits",
            "beam_matching.haslam_native_fwhm_deg", "full", "galactic",
            "BEAMSIZE = 56.0 arcmin, HDU 1 of the FITS file"),
    Dataset("arcade2_3150", "paths.arcade2_map_3150",
            "beam_matching.arcade2_native_fwhm_deg", "zeros", "unstated",
            "external -- no beam keyword in the FITS header"),
    Dataset("arcade2_3410", "paths.arcade2_map_3410",
            "beam_matching.arcade2_native_fwhm_deg", "zeros", "unstated",
            "external -- no beam keyword in the FITS header"),
)


def load_input_map(dataset: Dataset, config: Config) -> np.ndarray:
    """Read one dataset's map, as RING-ordered float."""
    import healpy as hp

    path = config.resolve_path(dataset.path_key)
    if not path.exists():
        raise FileNotFoundError(
            "no {} map at {} -- set {} to the local download "
            "(see DATA_SOURCES.md)".format(dataset.name, path, dataset.path_key)
        )
    return hp.read_map(path)


# ---------------------------------------------------------------------------
# the stage
# ---------------------------------------------------------------------------
def beam_match(
    config: Optional[Config] = None,
    *,
    config_path: Optional[Path] = None,
    target_fwhm_deg: Optional[float] = None,
    verbose: bool = True,
) -> Dict[str, Dict[str, object]]:
    """Beam-match every dataset in :data:`DATASETS` to ``target_fwhm_deg``.

    Every number comes from the ``beam_matching`` config block; nothing here
    re-derives a decision the config records.  ``target_fwhm_deg`` overrides it
    for a one-off comparison.

    Returns one entry per dataset: the smoothed map, the mask it was smoothed
    under, and the numbers a manifest would want.  Writes nothing.
    """
    import healpy as hp

    config = config if config is not None else load_config(config_path)
    if target_fwhm_deg is None:
        target_fwhm_deg = float(config.require("beam_matching.target_fwhm_deg"))
    weight_floor = float(config.require("beam_matching.weight_floor"))

    if verbose:
        print("TRIS beam: {}° (E) x {}° (H)".format(
            TRIS_HPBW_E_DEG, TRIS_HPBW_H_DEG))
        print("common target FWHM: {:.4f}° (the larger axis)\n".format(
            target_fwhm_deg))

    results: Dict[str, Dict[str, object]] = {}
    for dataset in DATASETS:
        native_fwhm_deg = float(config.require(dataset.fwhm_key))
        sky = load_input_map(dataset, config)
        mask = observed_mask(sky, dataset.mask_convention)
        kernel = extra_fwhm_deg(native_fwhm_deg, target_fwhm_deg)

        smoothed = smooth_to_target(
            sky,
            native_fwhm_deg,
            target_fwhm_deg,
            weight_mask=None if dataset.mask_convention == "full" else mask,
            weight_floor=weight_floor,
        )

        observed_before = int(mask.sum())
        observed_after = int(np.sum(smoothed != hp.UNSEEN))
        results[dataset.name] = {
            "map": smoothed,
            "mask": mask,
            "nside": int(hp.get_nside(sky)),
            "native_fwhm_deg": native_fwhm_deg,
            "fwhm_source": dataset.fwhm_source,
            "extra_fwhm_deg": kernel,
            "frame": dataset.frame,
            "observed_before": observed_before,
            "observed_after": observed_after,
        }

        if verbose:
            print("{:<14} nside {:<4} native {:7.4f}° -> extra kernel "
                  "{:7.4f}°".format(dataset.name, results[dataset.name]["nside"],
                                    native_fwhm_deg, kernel))
            print("{:<14} observed {} -> {} pixels ({:+d} at the mask edge), "
                  "frame {}".format("", observed_before, observed_after,
                                    observed_after - observed_before,
                                    dataset.frame))

    if verbose:
        print("\nNot done here: the regrid onto the common nside and the frame "
              "reconciliation\n(Step 0.4-0.5).  Haslam is galactic, TRIS is "
              "equatorial, ARCADE 2 does not say.")
    return results


def plot_beam_matched(results: Dict[str, Dict[str, object]]) -> None:
    """Mollview of each beam-matched map.  Import-safe: call it explicitly."""
    import healpy as hp

    for name, entry in results.items():
        hp.mollview(entry["map"], unit="K",
                    title="{}, beam-matched".format(name))


def main(argv=None) -> int:
    import argparse

    parser = argparse.ArgumentParser(prog="stage2_beam_matching")
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)

    beam_match(config_path=args.config, verbose=not args.quiet)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
