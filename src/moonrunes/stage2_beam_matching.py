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

The frame and the grid (Step 0.4-0.5)
-------------------------------------
Everything is matched **to TRIS**: TRIS's beam and TRIS's grid, and the TRIS
maps themselves are not touched.  Each map is rotated into equatorial (the
frame of stage 1's TRIS maps) with :mod:`moonrunes.frames`, smoothed, cut to
the pixels the survey actually observed whose smoothed weight reached the
floor, and regridded onto ``beam.nside_new`` (= ``tris.nside``).
:func:`run_stage2` writes the result to ``outputs/stage2/``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np

from .config import Config, ConfigError, load_config

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
    fwhm_source: str       # where that FWHM came from
    units: str             # what the temperatures are, for the manifest


DATASETS = (
    Dataset("haslam_408", "paths.haslam_map_fits",
            "beam_matching.haslam_native_fwhm_deg", "full",
            "header BEAMSIZE = 56.0 arcmin (HDU 1)",
            "K, Rayleigh-Jeans, CMB monopole included"),
    Dataset("arcade2_3150", "paths.arcade2_map_3150",
            "beam_matching.arcade2_native_fwhm_deg", "zeros",
            "header BEAMSZ = 11.6 deg (HDU 0); LAMBDA; Singal et al. 2011",
            "K, thermodynamic, CMB monopole included (dipole removed)"),
    Dataset("arcade2_3410", "paths.arcade2_map_3410",
            "beam_matching.arcade2_native_fwhm_deg", "zeros",
            "header BEAMSZ = 11.6 deg (HDU 0); LAMBDA; Singal et al. 2011",
            "K, thermodynamic, CMB monopole included (dipole removed)"),
)

#: The frame every stage 2 product is written in: the frame of stage 1's TRIS
#: maps, which stage 2 never touches (everything is matched to TRIS).
OUTPUT_FRAME = "C"


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
# onto the common grid
# ---------------------------------------------------------------------------
def regrid_to_nside(
    sky: np.ndarray, observed: np.ndarray, nside_out: int, rule: str = "all_children"
) -> Tuple[np.ndarray, np.ndarray]:
    """Average a RING map onto a coarser (or the same) grid.

    A coarse pixel gets a value only if **all** of its children were observed
    (``rule = "all_children"``), so every value is a mean of observed sky and
    nothing is extrapolated into a pixel the survey only partly saw.  For a
    full-sky map this is plain ``hp.ud_grade`` averaging.  Going to a *finer*
    grid would invent resolution, so it raises.

    Returns ``(values, observed)`` with ``hp.UNSEEN`` where unobserved.
    """
    import healpy as hp

    if rule != "all_children":
        raise ValueError(
            "beam_matching.coarse_pixel_rule {!r} is not implemented; only "
            "'all_children'".format(rule)
        )
    sky = np.asarray(sky, dtype=float)
    observed = np.asarray(observed, dtype=bool)
    nside_in = hp.get_nside(sky)
    if nside_out > nside_in:
        raise ValueError(
            "cannot regrid nside {} up to {}: that would invent resolution the "
            "map does not have".format(nside_in, nside_out)
        )
    if nside_out == nside_in:
        return np.where(observed, sky, hp.UNSEEN), observed.copy()

    filled = np.where(observed, sky, 0.0)
    total = hp.ud_grade(filled, nside_out, order_in="RING", order_out="RING")
    fraction = hp.ud_grade(observed.astype(float), nside_out,
                           order_in="RING", order_out="RING")
    keep = fraction >= 1.0 - 1e-9
    values = np.full(total.size, hp.UNSEEN)
    values[keep] = total[keep] / fraction[keep]
    return values, keep


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
    """Match every dataset in :data:`DATASETS` to the TRIS beam, on the TRIS grid.

    Per dataset, in this order:

    1. read it, and check the frame the file itself declares
       (:func:`moonrunes.frames.declared_frame`);
    2. rotate it into equatorial, the frame of stage 1's TRIS maps -- harmonic
       space for a full-sky map, mask-aware pixel space for a partial one;
    3. smooth it to ``beam_matching.target_fwhm_deg`` (TRIS's beam) with the
       quadrature kernel, mask-aware for a partial map;
    4. keep a pixel only if it was observed before smoothing **and** its
       smoothed weight reached ``beam_matching.weight_floor`` -- never a pixel
       the survey did not see;
    5. regrid onto ``beam.nside_new`` with :func:`regrid_to_nside`.

    TRIS itself is not touched: everything else is matched to it.  Returns one
    entry per dataset with the final map and the numbers the manifest records.
    Writes nothing; :func:`run_stage2` does.
    """
    import healpy as hp
    from . import frames

    config = config if config is not None else load_config(config_path)
    if target_fwhm_deg is None:
        target_fwhm_deg = float(config.require("beam_matching.target_fwhm_deg"))
    weight_floor = float(config.require("beam_matching.weight_floor"))
    nside_out = int(config.require("beam.nside_new"))
    rule = str(config.require("beam_matching.coarse_pixel_rule"))
    tris_nside = int(config.require("tris.nside"))
    if nside_out != tris_nside:
        raise ConfigError(
            "beam.nside_new ({}) must equal tris.nside ({}): stage 2 puts every "
            "map on the TRIS grid and leaves TRIS itself alone".format(
                nside_out, tris_nside)
        )

    if verbose:
        print("stage 2: match to the TRIS beam, {}° (E) x {}° (H); target FWHM "
              "{:.4f}° (the larger axis)".format(
                  TRIS_HPBW_E_DEG, TRIS_HPBW_H_DEG, target_fwhm_deg))
        print("         output: equatorial, nside {} (the TRIS map grid)\n".format(
            nside_out))

    results: Dict[str, Dict[str, object]] = {}
    for dataset in DATASETS:
        path = config.resolve_path(dataset.path_key)
        native_fwhm_deg = float(config.require(dataset.fwhm_key))
        sky = load_input_map(dataset, config)
        declared = frames.declared_frame(path)
        frame = declared["frame"]
        if frame not in ("G", "C"):
            raise ValueError(
                "{} declares frame {!r}; stage 2 handles Galactic or equatorial "
                "maps only".format(path.name, frame)
            )
        native_nside = int(hp.get_nside(sky))
        observed_native = observed_mask(sky, dataset.mask_convention)

        # 2. into equatorial
        if frame == OUTPUT_FRAME:
            rotated, observed = np.where(observed_native, sky, hp.UNSEEN), observed_native
            rotation = "none (already equatorial)"
        elif dataset.mask_convention == "full":
            rotated = frames.rotate_full_sky(sky, (frame, OUTPUT_FRAME))
            observed = np.ones(rotated.size, dtype=bool)
            rotation = "harmonic (rotate_map_alms), full sky"
        else:
            rotated, observed, _ = frames.rotate_masked(
                sky, observed_native, (frame, OUTPUT_FRAME))
            rotation = "pixel space, mask-aware (frames.rotate_masked)"

        # 3. to the TRIS beam
        kernel = extra_fwhm_deg(native_fwhm_deg, target_fwhm_deg)
        smoothed = smooth_to_target(
            rotated, native_fwhm_deg, target_fwhm_deg,
            weight_mask=None if dataset.mask_convention == "full" else observed,
            weight_floor=weight_floor,
        )

        # 4. only pixels the survey saw, and the floor kept
        passed_floor = smoothed != hp.UNSEEN
        final = observed & passed_floor
        extrapolated = int((passed_floor & ~observed).sum())

        # 5. onto the TRIS grid
        values, observed_out = regrid_to_nside(
            np.where(final, smoothed, hp.UNSEEN), final, nside_out, rule)

        results[dataset.name] = {
            "map": values,
            "observed": observed_out,
            "input_path": str(path),
            "declared_frame": frame,
            "frame_statement": declared["statements"],
            "rotation": rotation,
            "native_nside": native_nside,
            "native_fwhm_deg": native_fwhm_deg,
            "fwhm_source": dataset.fwhm_source,
            "extra_fwhm_deg": kernel,
            "units": dataset.units,
            "observed_before": int(observed.sum()),
            "observed_after": int(final.sum()),
            "extrapolated_removed": extrapolated,
            "observed_on_output_grid": int(observed_out.sum()),
        }

        if verbose:
            r = results[dataset.name]
            print("{:<14} {} nside {:<4} native {:7.4f}° -> extra kernel {:7.4f}°".format(
                dataset.name, frame, native_nside, native_fwhm_deg, kernel))
            print("{:<14} observed {} -> {} pixels after smoothing ({} extrapolated "
                  "removed); {} on the nside {} grid".format(
                      "", r["observed_before"], r["observed_after"],
                      extrapolated, r["observed_on_output_grid"], nside_out))
    return results


def _output_dir(config: Config) -> Path:
    from .config import REPO_ROOT

    root = Path(str(config.require("run.outputs_dir")))
    if not root.is_absolute():
        root = REPO_ROOT / root
    return root / STAGE


def stage2_product_paths(config: Config, output_dir: Optional[Path] = None) -> Tuple[Path, Path]:
    """``(maps npz, manifest json)`` for this run, whether or not they exist."""
    destination = Path(output_dir) if output_dir else _output_dir(config)
    run_name = str(config.require("run.name"))
    return (destination / "beam_matched_{}.npz".format(run_name),
            destination / "{}_manifest.json".format(STAGE))


def run_stage2(
    config: Optional[Config] = None,
    *,
    config_path: Optional[Path] = None,
    output_dir: Optional[Path] = None,
    overwrite: Optional[bool] = None,
    verbose: bool = True,
) -> Dict[str, object]:
    """Run stage 2 and write its product.  Returns the manifest dict.

    ``outputs/stage2/beam_matched_<run name>.npz`` holds, stacked in
    :data:`DATASETS` order: ``maps`` (``hp.UNSEEN`` where unobserved) and
    ``observed`` on the ``beam.nside_new`` RING grid, in equatorial, plus
    ``names``, ``nside`` and ``target_fwhm_deg``.  ``stage2_manifest.json``
    records the config used and every number :func:`beam_match` reports.
    """
    import datetime as _dt
    import json
    import sys

    import healpy as hp

    config = config if config is not None else load_config(config_path)
    maps_path, manifest_path = stage2_product_paths(config, output_dir)
    maps_path.parent.mkdir(parents=True, exist_ok=True)
    allow_overwrite = (bool(config.require("run.overwrite"))
                       if overwrite is None else bool(overwrite))
    for path in (maps_path, manifest_path):
        if path.exists() and not allow_overwrite:
            raise FileExistsError(
                "{} exists; set run.overwrite or pass --overwrite".format(path))

    results = beam_match(config, verbose=verbose)
    names = [d.name for d in DATASETS]
    nside_out = int(config.require("beam.nside_new"))
    np.savez_compressed(
        maps_path,
        names=np.array(names),
        nside=np.int64(nside_out),
        frame=np.array(OUTPUT_FRAME),
        target_fwhm_deg=np.float64(config.require("beam_matching.target_fwhm_deg")),
        maps=np.vstack([results[n]["map"] for n in names]),
        observed=np.vstack([results[n]["observed"] for n in names]),
    )

    manifest = {
        "stage": STAGE,
        "run_name": str(config.require("run.name")),
        "written_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "config_path": str(config.path) if config.path else None,
        "products": {"maps": str(maps_path), "manifest": str(manifest_path)},
        "output_grid": {"nside": nside_out, "ordering": "RING",
                        "frame": "equatorial (C), the frame of the stage 1 TRIS maps",
                        "unobserved": "hp.UNSEEN"},
        "tris": "not touched: every other map is matched to the TRIS beam and grid",
        "config_used": {
            "run": config.section("run"),
            "beam_matching": config.section("beam_matching"),
            "beam_nside_new": nside_out,
            "tris_nside": int(config.require("tris.nside")),
        },
        "datasets": {n: {k: v for k, v in results[n].items()
                         if k not in ("map", "observed")} for n in names},
        "environment": {"python": sys.version.split()[0], "numpy": np.__version__,
                        "healpy": hp.__version__},
    }
    with open(manifest_path, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
    if verbose:
        print("\nstage 2: wrote {}".format(maps_path))
        print("stage 2: wrote {}".format(manifest_path))
    return manifest


@dataclass(frozen=True)
class Stage2Products:
    """Stage 2's product: beam-matched maps on the TRIS grid, in equatorial."""

    arrays: Dict[str, np.ndarray]
    manifest: Dict[str, object]
    maps_path: Path
    manifest_path: Path

    @property
    def names(self) -> list:
        return [str(n) for n in self.arrays["names"]]

    @property
    def nside(self) -> int:
        return int(self.arrays["nside"])

    def map(self, name: str) -> np.ndarray:
        """One dataset's beam-matched map, ``hp.UNSEEN`` where unobserved."""
        return self.arrays["maps"][self.names.index(name)]

    def observed(self, name: str) -> np.ndarray:
        return self.arrays["observed"][self.names.index(name)]


def load_stage2(
    config: Optional[Config] = None,
    *,
    config_path: Optional[Path] = None,
    maps_path: Optional[Path] = None,
    manifest_path: Optional[Path] = None,
) -> Stage2Products:
    """Load stage 2's product, defaulting to this run's own output paths."""
    import json

    config = config if config is not None else load_config(config_path)
    default_maps, default_manifest = stage2_product_paths(config)
    maps_path = Path(maps_path) if maps_path else default_maps
    manifest_path = Path(manifest_path) if manifest_path else default_manifest
    if not maps_path.exists():
        raise FileNotFoundError(
            "no stage 2 product at {} -- run `python run_pipeline.py --stage 2` "
            "first".format(maps_path))
    with np.load(maps_path) as bundle:
        arrays = {key: bundle[key] for key in bundle.files}
    with open(manifest_path, "r", encoding="utf-8") as handle:
        manifest = json.load(handle)
    return Stage2Products(arrays=arrays, manifest=manifest,
                          maps_path=maps_path, manifest_path=manifest_path)


def plot_beam_matched(results: Dict[str, Dict[str, object]]) -> None:
    """Mollview of each beam-matched map.  Import-safe: call it explicitly."""
    import healpy as hp

    for name, entry in results.items():
        hp.mollview(entry["map"], unit="K", coord="C",
                    title="{}, beam-matched".format(name))


def main(argv=None) -> int:
    import argparse

    parser = argparse.ArgumentParser(prog="stage2_beam_matching")
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)

    run_stage2(config_path=args.config, output_dir=args.output_dir,
               overwrite=True if args.overwrite else None, verbose=not args.quiet)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
