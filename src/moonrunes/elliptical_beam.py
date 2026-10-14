"""The elliptical TRIS beam as a smoothing kernel, for the map route.

Why this module exists
----------------------
Stage 2 and the map route (:mod:`moonrunes.sed_prototype`) smooth every survey
to the TRIS beam so that each point of a per-pixel SED was seen through the
same beam.  ``hp.smoothing`` can only apply a kernel that is the same in every
direction around its centre (it works through spherical harmonics, where such a
kernel is one number per multipole).  TRIS's beam is not: it is 19.155 deg wide
in its E plane and 23.366 deg in its H plane, and it has non-Gaussian
shoulders.  The original stage 2 therefore used a *circular* Gaussian at the
larger axis, 23.366 deg, which leaves every map over-smoothed by 22% along the
narrow axis.

This module applies the real shape instead, directly on the sphere, as a
weighted sum over pixels.  It is selected with ``beam_matching.kernel``:

``circular_gaussian``
    The original behaviour (``hp.smoothing`` at ``target_fwhm_deg``).  Nothing
    in this module is used for it, except for side-by-side diagnostics.
``elliptical_gaussian``
    A Gaussian with FWHM 19.155 deg (E) x 23.366 deg (H), correctly oriented.
``tris_cuts``
    The beam limTOD uses for the TOD route and stage 1: the archive's measured
    E and H cuts, blended with ``cos^2 phi`` / ``sin^2 phi`` weights
    (:func:`limTOD.tris.tris_cut_beam_response`).  Elliptical *and* with the
    measured shoulders.  The recommended choice.

The kernel, mathematically
--------------------------
For an output position :math:`p` (a unit vector) and an input pixel :math:`q`,
let :math:`\\theta` be their angular separation and :math:`\\psi` the position
angle of :math:`q` seen from :math:`p`, measured from north through east (the
usual astronomical convention).  The kernel is

.. math::

    K(p, q) = B(\\theta, \\phi), \\qquad \\phi = \\psi - \\psi_E,

where :math:`B(\\theta, \\phi)` is the beam in its own frame (E plane at
:math:`\\phi = 0, 180`, H plane at :math:`\\phi = 90, 270`) and
:math:`\\psi_E` is the position angle of the E plane on the sky.

**Orientation.**  TRIS is a zenith-pointing drift scan.  limTOD's
:func:`limTOD.tris.tris_zenith_geometry` puts the E plane along the meridian,
rolled ``TRIS_E_PLANE_EAST_OF_MERIDIAN_DEG`` = 7 deg east of it, as
``TRIS_Beam_Profile.txt`` and both ring headers state.  At the zenith the local
meridian is the hour circle, so on the equatorial sky the E plane points
north-south (rolled 7 deg towards east) and the H plane east-west.  That is
:math:`\\psi_E` = ``beam_matching.elliptical.e_plane_east_of_meridian_deg``.
The cut beam is symmetric under :math:`\\phi \\to -\\phi`, so the handedness of
the beam's own frame cannot matter -- only :math:`\\psi_E` does.

TRIS only ever pointed at Dec = +42.4 deg.  At every other declination the
kernel uses the same orientation relative to local north: the beam as TRIS
would have seen it had it been pointing there.

**Horizon.**  limTOD blocks beam response below the horizon (the horn sees
ground there, not sky).  At the zenith that is :math:`\\theta > 90` deg, so the
kernel is set to zero beyond ``horizon_deg``.  It only matters for the cut
beam, whose archive floor is about -60 dB all the way to 176 deg; the Gaussian
is negligible there anyway.

**Normalised convolution.**  Exactly as for the circular kernel (``PIPELINE.md``
Appendix A.6): with :math:`M(q)` the observed fraction of input pixel
:math:`q`,

.. math::

    T_s(p) = \\frac{\\sum_q K(p,q) M(q) T(q)}{\\sum_q K(p,q) M(q)}, \\qquad
    w(p) = \\frac{\\sum_q K(p,q) M(q)}{\\sum_q K(p,q)} .

:math:`w` is the fraction of the kernel's weight on observed sky, compared with
``beam_matching.weight_floor`` exactly as before.

The native beam of each survey
------------------------------
The effective beam of a smoothed map is the survey's own beam convolved with
the kernel.  For Gaussian beams this is exact per axis: the kernel needs
FWHM :math:`\\sqrt{\\mathrm{target}^2 - \\mathrm{native}^2}` along E and along H
(the quadrature rule, ``PIPELINE.md`` A.4).  The cut beam is not a Gaussian and
has no exact "subtract the native beam" operation, so :func:`kernel_for`
decides per map:

* **Native beam negligible** -- if convolving the cut beam with it would widen
  the narrow (E) axis by at most ``max_broadening_fraction`` (0.2%: natives up
  to about 1.2 deg; Haslam 0.93 deg widens it by 0.12%), the kernel is the cut
  beam itself.
* **Native beam too broad** (LWA1 35/80 MHz, Maipu 45 MHz, ARCADE 2) -- set by
  ``broad_native_kernel``:

  - ``scaled_cuts`` (default): the cut beam, compressed in angle separately
    along E and H so that its half-power widths are the quadrature widths.
    Convolved with the native beam, the core comes back to 19.155 x 23.366 deg
    and the shoulders keep the measured shape, approximately.  Exact only if
    both were Gaussian; notebook 03 measures how close it is.
  - ``elliptical_gaussian``: the per-axis quadrature Gaussian.  The core is
    exact (for a Gaussian native beam) but the shoulders are Gaussian, not
    TRIS's.

How it is computed
------------------
A direct sum over a *working grid* (``nside_work``, default 64: 0.92 deg
pixels, 20x finer than the 19 deg beam), evaluated at ``eval_subdivision**2``
points inside each output pixel and averaged, so that each output value is a
pixel average like the circular path's ``all_children`` regrid.  Each block of
evaluation points is one matrix product, so the cost is
``(output pixels x eval points) x (working pixels)`` kernel evaluations:
3072 x 4 x 49152 = 6e8 for a full-sky nside 16 map at the defaults.  The
runtime has not been measured yet.

Nothing here writes files; :func:`moonrunes.stage2_beam_matching.match_to_tris_beam`
is the one entry point the pipeline uses.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Dict, Optional

import numpy as np

from .config import Config, ConfigError

#: Every value ``beam_matching.kernel`` may take.
KERNELS = ("circular_gaussian", "elliptical_gaussian", "tris_cuts")
#: Every value ``beam_matching.elliptical.broad_native_kernel`` may take.
BROAD_NATIVE_KERNELS = ("scaled_cuts", "elliptical_gaussian")
#: FWHM = 2 sqrt(2 ln 2) sigma (notebook 09, section 1).
FWHM_PER_SIGMA = 2.0 * np.sqrt(2.0 * np.log(2.0))
#: The archive file holding the E/H principal-plane cuts.
BEAM_PROFILE_FILE = "TRIS_Beam_Profile.txt"


# ---------------------------------------------------------------------------
# configuration
# ---------------------------------------------------------------------------
@dataclass(frozen=True, eq=False)
class BeamSpec:
    """Everything that defines the elliptical smoothing, from the config.

    The three ``cut_*`` arrays are the archive's E/H cuts (sorted by angle), and
    are only needed for ``kernel = "tris_cuts"``.
    """

    kernel: str
    fwhm_e_deg: float
    fwhm_h_deg: float
    e_plane_east_of_meridian_deg: float
    horizon_deg: float
    nside_work: int
    eval_subdivision: int
    max_broadening_fraction: float
    broad_native_kernel: str
    blend: str = "db"
    chunk: int = 128
    cut_angle_deg: Optional[np.ndarray] = field(default=None, repr=False)
    cut_e_db: Optional[np.ndarray] = field(default=None, repr=False)
    cut_h_db: Optional[np.ndarray] = field(default=None, repr=False)

    def __post_init__(self) -> None:
        import healpy as hp

        if self.kernel not in KERNELS:
            raise ValueError("kernel {!r} is not one of {}".format(self.kernel, KERNELS))
        if self.broad_native_kernel not in BROAD_NATIVE_KERNELS:
            raise ValueError("broad_native_kernel {!r} is not one of {}".format(
                self.broad_native_kernel, BROAD_NATIVE_KERNELS))
        if self.blend not in ("db", "power"):
            raise ValueError('blend must be "db" or "power"')
        if not (self.fwhm_e_deg > 0 and self.fwhm_h_deg > 0):
            raise ValueError("the target FWHMs must be positive")
        if not 0.0 < self.horizon_deg <= 180.0:
            raise ValueError("horizon_deg must be in (0, 180]")
        if not hp.isnsideok(self.nside_work):
            raise ValueError("nside_work {} is not a valid nside".format(self.nside_work))
        s = int(self.eval_subdivision)
        if s < 1 or (s & (s - 1)) != 0:
            raise ValueError("eval_subdivision must be a power of 2 (1, 2, 4, ...)")
        if self.max_broadening_fraction < 0:
            raise ValueError("max_broadening_fraction must be >= 0")
        if self.chunk < 1:
            raise ValueError("chunk must be >= 1")
        if self.kernel == "tris_cuts" and self.cut_angle_deg is None:
            raise ValueError('kernel "tris_cuts" needs the archive cuts (cut_angle_deg, '
                             'cut_e_db, cut_h_db); build the spec with beam_spec_from_config')


def beam_spec_from_config(config: Config, kernel: Optional[str] = None) -> BeamSpec:
    """Build the :class:`BeamSpec` from ``beam_matching.*``.

    ``kernel`` overrides ``beam_matching.kernel`` (the diagnostics notebook uses
    it to put the three kernels side by side).  For ``tris_cuts`` the archive
    cuts are read through limTOD, and the configured FWHMs must agree with what
    the cuts themselves give to 0.001 deg -- otherwise the config and the
    kernel would silently describe two different beams.
    """
    kernel = str(config.require("beam_matching.kernel")) if kernel is None else str(kernel)
    if kernel not in KERNELS:
        raise ConfigError("beam_matching.kernel {!r} is not one of {}".format(kernel, KERNELS))

    def get(key):
        return config.require("beam_matching.elliptical." + key)

    kwargs = dict(
        kernel=kernel,
        fwhm_e_deg=float(get("fwhm_e_deg")),
        fwhm_h_deg=float(get("fwhm_h_deg")),
        e_plane_east_of_meridian_deg=float(get("e_plane_east_of_meridian_deg")),
        horizon_deg=float(get("horizon_deg")),
        nside_work=int(get("nside_work")),
        eval_subdivision=int(get("eval_subdivision")),
        max_broadening_fraction=float(get("max_broadening_fraction")),
        broad_native_kernel=str(get("broad_native_kernel")),
        blend=str(get("blend")),
        chunk=int(get("chunk")),
    )
    if kernel == "tris_cuts":
        from .stage1_tris_maps import import_limtod_tris

        tris = import_limtod_tris(config)
        path = config.resolve_path("paths.tris_archive_dir") / BEAM_PROFILE_FILE
        cuts = tris.read_tris_beam_cuts(path)
        measured_e, measured_h = cuts.half_power_full_width_deg()
        if (abs(measured_e - kwargs["fwhm_e_deg"]) > 1e-3
                or abs(measured_h - kwargs["fwhm_h_deg"]) > 1e-3):
            raise ConfigError(
                "beam_matching.elliptical.fwhm_e_deg/fwhm_h_deg ({}, {}) do not match "
                "the archive cuts ({:.4f}, {:.4f}) in {}".format(
                    kwargs["fwhm_e_deg"], kwargs["fwhm_h_deg"], measured_e, measured_h, path))
        order = np.argsort(cuts.angle_deg)
        kwargs.update(cut_angle_deg=np.asarray(cuts.angle_deg, dtype=float)[order],
                      cut_e_db=np.asarray(cuts.e_plane_db, dtype=float)[order],
                      cut_h_db=np.asarray(cuts.h_plane_db, dtype=float)[order])
    return BeamSpec(**kwargs)


def with_kernel(spec: BeamSpec, kernel: str) -> BeamSpec:
    """The same spec with another ``kernel`` (``tris_cuts`` needs the cuts loaded)."""
    return replace(spec, kernel=kernel)


# ---------------------------------------------------------------------------
# the beam, in its own frame
# ---------------------------------------------------------------------------
def gaussian_response(theta_deg, phi_deg, fwhm_e_deg: float, fwhm_h_deg: float) -> np.ndarray:
    """Elliptical Gaussian: ``exp(-((theta cos phi / sigma_E)^2 + (theta sin phi / sigma_H)^2) / 2)``.

    The same expression limTOD's ``approximate_tris_gaussian_beam_map`` uses, with
    ``theta`` the angle from the beam centre and ``phi`` the beam's own azimuth
    (0 on the E plane).  Peak 1.  Degrees in, degrees throughout.
    """
    theta = np.asarray(theta_deg, dtype=float)
    phi = np.deg2rad(np.asarray(phi_deg, dtype=float))
    sigma_e = float(fwhm_e_deg) / FWHM_PER_SIGMA
    sigma_h = float(fwhm_h_deg) / FWHM_PER_SIGMA
    return np.exp(-0.5 * ((theta * np.cos(phi) / sigma_e) ** 2
                          + (theta * np.sin(phi) / sigma_h) ** 2))


def cut_response(angle_deg, e_db, h_db, theta_deg, phi_deg, *, scale_e: float = 1.0,
                 scale_h: float = 1.0, blend: str = "db") -> np.ndarray:
    """The archive-cut beam: limTOD's ``tris_cut_beam_response``, plus an optional scale.

    Each cut is interpolated in angle (linearly in dB, held constant beyond the
    last tabulated angle, as limTOD does), and the two are blended with
    ``cos^2 phi`` (E) and ``sin^2 phi`` (H) weights -- in dB for ``blend="db"``,
    in linear power for ``"power"``.  With ``scale_e = scale_h = 1`` this is
    limTOD's function exactly (``tests/test_elliptical_beam.py`` checks it).

    ``scale_e`` / ``scale_h`` compress the E / H cut in angle: the cut is read at
    ``theta * scale``, so its half-power width shrinks by ``1 / scale``.  That is
    how :func:`kernel_for` makes a cut-shaped kernel narrower than the beam, for
    surveys whose own beam is too broad to ignore.
    """
    theta = np.asarray(theta_deg, dtype=float)
    phi = np.deg2rad(np.asarray(phi_deg, dtype=float))
    e_at = np.interp(theta * scale_e, angle_deg, e_db, left=e_db[0], right=e_db[-1])
    h_at = np.interp(theta * scale_h, angle_deg, h_db, left=h_db[0], right=h_db[-1])
    weight_e = np.cos(phi) ** 2
    weight_h = np.sin(phi) ** 2
    if blend == "db":
        return 10.0 ** ((e_at * weight_e + h_at * weight_h) / 10.0)
    return 10.0 ** (e_at / 10.0) * weight_e + 10.0 ** (h_at / 10.0) * weight_h


# ---------------------------------------------------------------------------
# the kernel for one survey
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Kernel:
    """The kernel applied to one map, and what effective beam it should give.

    ``shape`` is ``"tris_cuts"`` or ``"elliptical_gaussian"`` (a circular
    Gaussian is an elliptical one with equal axes).  ``fwhm_e_deg`` /
    ``fwhm_h_deg`` are the kernel's own half-power widths; ``scale_e`` /
    ``scale_h`` are the cut compression factors (1 = the measured cut).
    ``expected_e_deg`` / ``expected_h_deg`` are the widths of native beam (x)
    kernel -- exact for Gaussians, a core-width estimate for the cut beam.
    """

    shape: str
    fwhm_e_deg: float
    fwhm_h_deg: float
    native_fwhm_deg: float
    expected_e_deg: float
    expected_h_deg: float
    note: str
    scale_e: float = 1.0
    scale_h: float = 1.0

    def as_dict(self) -> Dict[str, object]:
        return {
            "kernel_shape": self.shape,
            "kernel_fwhm_e_deg": float(self.fwhm_e_deg),
            "kernel_fwhm_h_deg": float(self.fwhm_h_deg),
            "cut_scale_e": float(self.scale_e),
            "cut_scale_h": float(self.scale_h),
            "expected_beam_e_deg": float(self.expected_e_deg),
            "expected_beam_h_deg": float(self.expected_h_deg),
            "kernel_note": self.note,
        }


def kernel_for(spec: BeamSpec, native_fwhm_deg: float) -> Kernel:
    """Choose the kernel for a survey with a (circular, Gaussian) native beam.

    See the module docstring, "The native beam of each survey".  Raises if the
    native beam is at least as wide as the narrow target axis: such a map is
    already coarser than TRIS along E and cannot be matched by smoothing.
    """
    n = float(native_fwhm_deg)
    e, h = float(spec.fwhm_e_deg), float(spec.fwhm_h_deg)
    if not np.isfinite(n) or n < 0:
        raise ValueError("native_fwhm_deg must be a finite, non-negative number")
    if spec.kernel == "circular_gaussian":
        # Diagnostics only: stage 2's circular path is hp.smoothing, not this.
        c = max(e, h)
        if n >= c:
            raise ValueError("cannot smooth a {:.3f} deg beam to {:.3f} deg".format(n, c))
        k = float(np.sqrt(c ** 2 - n ** 2))
        return Kernel("elliptical_gaussian", k, k, n, c, c,
                      "circular Gaussian at the larger axis, quadrature kernel")
    if n >= min(e, h):
        raise ValueError(
            "cannot smooth a {:.3f} deg beam to {:.3f} x {:.3f} deg: the native beam is "
            "already wider than TRIS's narrow axis".format(n, e, h))
    k_e, k_h = float(np.sqrt(e ** 2 - n ** 2)), float(np.sqrt(h ** 2 - n ** 2))

    if spec.kernel == "elliptical_gaussian":
        return Kernel("elliptical_gaussian", k_e, k_h, n, e, h,
                      "elliptical Gaussian, per-axis quadrature kernel (exact for a "
                      "Gaussian native beam)")

    # tris_cuts
    broadening = float(np.hypot(e, n) / e - 1.0)      # largest on the narrow axis
    if broadening <= spec.max_broadening_fraction:
        return Kernel("tris_cuts", e, h, n, float(np.hypot(e, n)), float(np.hypot(h, n)),
                      "measured cut beam; native beam neglected (widens the E core by "
                      "{:.2f}%)".format(100.0 * broadening))
    if spec.broad_native_kernel == "elliptical_gaussian":
        return Kernel("elliptical_gaussian", k_e, k_h, n, e, h,
                      "native beam too broad to neglect ({:.2f}% on E): elliptical-Gaussian "
                      "quadrature kernel, so the core matches TRIS but the shoulders are "
                      "Gaussian".format(100.0 * broadening))
    return Kernel("tris_cuts", k_e, k_h, n, e, h,
                  "native beam too broad to neglect ({:.2f}% on E): cut beam compressed to "
                  "the quadrature core widths (approximate; exact only for Gaussian "
                  "cores)".format(100.0 * broadening),
                  scale_e=e / k_e, scale_h=h / k_h)


def kernel_response(kernel: Kernel, spec: BeamSpec, theta_deg, position_angle_deg) -> np.ndarray:
    """The kernel's weight at separation ``theta`` and sky position angle ``psi``.

    ``psi`` is measured from north through east; the beam's own azimuth is
    ``phi = psi - e_plane_east_of_meridian_deg``.  Zero beyond ``horizon_deg``.
    """
    theta = np.asarray(theta_deg, dtype=float)
    phi = np.asarray(position_angle_deg, dtype=float) - spec.e_plane_east_of_meridian_deg
    if kernel.shape == "tris_cuts":
        response = cut_response(spec.cut_angle_deg, spec.cut_e_db, spec.cut_h_db, theta, phi,
                                scale_e=kernel.scale_e, scale_h=kernel.scale_h, blend=spec.blend)
    elif kernel.shape == "elliptical_gaussian":
        response = gaussian_response(theta, phi, kernel.fwhm_e_deg, kernel.fwhm_h_deg)
    else:
        raise ValueError("unknown kernel shape {!r}".format(kernel.shape))
    return np.where(theta <= spec.horizon_deg, response, 0.0)


# ---------------------------------------------------------------------------
# geometry
# ---------------------------------------------------------------------------
def local_offsets(centres: np.ndarray, points: np.ndarray):
    """Separation and position angle of every point, seen from every centre.

    ``centres`` is ``(k, 3)`` and ``points`` ``(n, 3)``, unit vectors in the
    map's frame (equatorial here: x to RA 0, y to RA 90, z to the north pole).
    Returns ``(theta_deg, position_angle_deg)``, each ``(n, k)``.

    At each centre, ``east = (-sin a, cos a, 0)`` and
    ``north = (-sin d cos a, -sin d sin a, cos d)``.  With ``e``, ``n``, ``r``
    the components of a point along east, north and the centre itself,
    ``theta = atan2(hypot(e, n), r)`` (accurate at small and large separations,
    unlike ``arccos``) and ``psi = atan2(e, n)`` (north = 0, east = 90).
    """
    c = np.atleast_2d(np.asarray(centres, dtype=float))
    x0, y0, z0 = c[:, 0], c[:, 1], c[:, 2]
    rho = np.maximum(np.hypot(x0, y0), 1e-12)          # cos(dec); no centre sits on a pole
    east = np.stack([-y0 / rho, x0 / rho, np.zeros_like(x0)], axis=1)
    north = np.stack([-z0 * x0 / rho, -z0 * y0 / rho, rho], axis=1)
    pts = np.asarray(points, dtype=float)
    e = pts @ east.T
    n = pts @ north.T
    r = pts @ c.T
    theta = np.degrees(np.arctan2(np.hypot(e, n), r))
    position_angle = np.degrees(np.arctan2(e, n))
    return theta, position_angle


def to_work_grid(sky: np.ndarray, observed: np.ndarray, nside_work: int):
    """Put a RING map on the working grid, keeping track of what was observed.

    Returns ``(filled, fraction)``: ``fraction`` is the observed fraction of each
    working pixel, and ``filled`` = fraction x (mean observed temperature in it).
    So ``filled`` is already the ``M * T`` of the normalised convolution.
    Degrading averages the children; upgrading (ARCADE 2 is native nside 16)
    copies each parent value into its children, which adds no information and
    loses none.
    """
    import healpy as hp

    sky = np.asarray(sky, dtype=float)
    observed = np.asarray(observed, dtype=bool)
    if sky.shape != observed.shape:
        raise ValueError("sky {} and observed {} differ in shape".format(sky.shape, observed.shape))
    filled = np.where(observed, sky, 0.0)
    fraction = observed.astype(float)
    if hp.get_nside(sky) != nside_work:
        filled = hp.ud_grade(filled, nside_work, order_in="RING", order_out="RING")
        fraction = hp.ud_grade(fraction, nside_work, order_in="RING", order_out="RING")
    return filled, fraction


def eval_positions(nside_out: int, pixels: np.ndarray, subdivision: int) -> np.ndarray:
    """Unit vectors of ``subdivision**2`` points inside each output pixel.

    They are the centres of the pixel's children at ``nside_out * subdivision``
    (in NESTED numbering the children of pixel ``p`` are ``p*s^2 ... p*s^2 + s^2 - 1``),
    returned as ``(len(pixels) * subdivision**2, 3)``, grouped by output pixel.
    """
    import healpy as hp

    pixels = np.asarray(pixels, dtype=np.int64)
    n_sub = int(subdivision) ** 2
    nested = hp.ring2nest(nside_out, pixels)
    children = nested[:, None] * n_sub + np.arange(n_sub)[None, :]
    ring = hp.nest2ring(nside_out * int(subdivision), children.ravel())
    return np.array(hp.pix2vec(nside_out * int(subdivision), ring)).T


# ---------------------------------------------------------------------------
# the convolution
# ---------------------------------------------------------------------------
def _sums_at(centres, filled, fraction, points, kernel, spec, progress=False):
    """``(sum K M T, sum K M, sum K)`` at each centre, in blocks of ``spec.chunk``."""
    n = len(centres)
    num, den, tot = np.empty(n), np.empty(n), np.empty(n)
    for start in range(0, n, spec.chunk):
        block = slice(start, min(start + spec.chunk, n))
        theta, psi = local_offsets(centres[block], points)
        weights = kernel_response(kernel, spec, theta, psi)        # (n_points, k)
        num[block] = filled @ weights
        den[block] = fraction @ weights
        tot[block] = weights.sum(axis=0)
        if progress and (start // spec.chunk) % 20 == 0:
            print("  elliptical smoothing: {}/{} evaluation points".format(block.stop, n))
    return num, den, tot


def smooth_elliptical(sky: np.ndarray, observed: np.ndarray, native_fwhm_deg: float,
                      spec: BeamSpec, nside_out: int, *, pixels: Optional[np.ndarray] = None,
                      kernel: Optional[Kernel] = None, progress: bool = False) -> Dict[str, object]:
    """Smooth a RING map (equatorial) with the elliptical kernel, onto ``nside_out``.

    ``sky`` may be at any nside; it is put on ``spec.nside_work`` first
    (:func:`to_work_grid`).  ``pixels`` restricts the output to some nside_out
    pixels (default: all).  Returns a dict with, as full nside_out RING maps
    (NaN where not computed):

    * ``values`` -- the smoothed temperature, the normalised convolution averaged
      over the evaluation points of each pixel;
    * ``weight`` -- the mean, and ``weight_min`` the smallest, observed-weight
      fraction ``w`` over a pixel's evaluation points;

    plus ``kernel`` (the :class:`Kernel` used) and ``pixels``.  No masking is
    applied here; :func:`moonrunes.stage2_beam_matching.match_to_tris_beam` does it.
    """
    import healpy as hp

    kernel = kernel if kernel is not None else kernel_for(spec, native_fwhm_deg)
    filled, fraction = to_work_grid(sky, observed, spec.nside_work)
    points = np.array(hp.pix2vec(spec.nside_work, np.arange(filled.size))).T
    npix_out = hp.nside2npix(nside_out)
    pixels = np.arange(npix_out) if pixels is None else np.unique(np.asarray(pixels, dtype=np.int64))
    centres = eval_positions(nside_out, pixels, spec.eval_subdivision)

    num, den, tot = _sums_at(centres, filled, fraction, points, kernel, spec, progress)
    n_sub = spec.eval_subdivision ** 2
    num, den, tot = (a.reshape(len(pixels), n_sub) for a in (num, den, tot))
    with np.errstate(invalid="ignore", divide="ignore"):
        values_p = num.sum(axis=1) / den.sum(axis=1)
        weight_points = den / tot

    out = {name: np.full(npix_out, np.nan) for name in ("values", "weight", "weight_min")}
    out["values"][pixels] = values_p
    out["weight"][pixels] = weight_points.mean(axis=1)
    out["weight_min"][pixels] = weight_points.min(axis=1)
    out.update(kernel=kernel, pixels=pixels)
    return out


def smooth_at_positions(sky: np.ndarray, observed: np.ndarray, native_fwhm_deg: float,
                        spec: BeamSpec, vectors: np.ndarray, *,
                        kernel: Optional[Kernel] = None) -> Dict[str, object]:
    """The smoothed map evaluated at arbitrary unit vectors (no pixel averaging).

    For diagnostics: sampling a smoothed point source along a line through it
    gives the effective beam profile directly, without a coarse output grid in
    the way.  Returns ``values``, ``weight`` (arrays over ``vectors``) and ``kernel``.
    """
    import healpy as hp

    kernel = kernel if kernel is not None else kernel_for(spec, native_fwhm_deg)
    filled, fraction = to_work_grid(sky, observed, spec.nside_work)
    points = np.array(hp.pix2vec(spec.nside_work, np.arange(filled.size))).T
    vectors = np.atleast_2d(np.asarray(vectors, dtype=float))
    num, den, tot = _sums_at(vectors, filled, fraction, points, kernel, spec)
    with np.errstate(invalid="ignore", divide="ignore"):
        return {"values": num / den, "weight": den / tot, "kernel": kernel}


def kernel_image(kernel: Kernel, spec: BeamSpec, half_width_deg: float = 40.0,
                 step_deg: float = 0.25):
    """The kernel on a flat (east, north) offset grid, for plotting.

    The grid is the azimuthal-equidistant projection about the kernel centre:
    a point at offset ``(x_east, y_north)`` is at separation ``hypot(x, y)`` and
    position angle ``atan2(x, y)``.  The kernel depends on exactly those two
    numbers, so the picture is exact, not a small-angle approximation.
    Returns ``(x, y, response)``; ``x`` increases to the east.
    """
    axis = np.arange(-half_width_deg, half_width_deg + step_deg / 2, step_deg)
    x, y = np.meshgrid(axis, axis)
    theta = np.hypot(x, y)
    position_angle = np.degrees(np.arctan2(x, y))
    return x, y, kernel_response(kernel, spec, theta, position_angle)
