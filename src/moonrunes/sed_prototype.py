"""Stage 3 prototype: the TOD route and the map route of the SED calibration.

The same calculation ``notebooks/07_sed_calibration`` develops step by step,
packaged so that plotting notebooks (``notebooks/08_paper_plots``) can call it.
It is a prototype, not a pipeline stage: it writes nothing, and the choices it
makes are the module constants below.

* **Map route** (:func:`map_route`): every usable map matched to the TRIS beam
  on the nside grid (stage 2 procedure), TRIS from the stage 1 maps + the fitted
  zero level, and a power law fitted at every TRIS-stripe pixel.
* **TOD route** (:func:`tod_route`): limTOD observes every map at its native
  resolution through the TRIS beam, a power law is fitted at every ring sample
  against the real TRIS data, and the correction TOD (predicted - target) is
  mapped with limTOD's Wiener filter and a zero-mean prior.

Temperatures in the fits are Rayleigh-Jeans with the CMB removed
(:func:`moonrunes.sky_maps.rj_excess`).  Predictions are quoted in the target's
own convention (fit + T_CMB,RJ(408 MHz)).
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from . import frames, sky_maps
from .config import Config
from .stage1_tris_maps import import_limtod_tris, load_stage1
from .stage2_beam_matching import regrid_to_nside, smooth_to_target

NU0_MHZ = 408.0
#: Error assumed for maps with no published per-pixel error (fraction of T).
ASSUMED_FRACTION = 0.10
#: ARCADE 2 calibration uncertainty, Fixsen et al. 2011 Table 3 (K).
ARCADE_CAL_K = {3150.0: 7.5e-3, 3410.0: 6.0e-3}
MIN_COVERAGE = 0.5
#: Zero-mean prior on the TOD-route correction: (fraction of target, constant K).
CORRECTION_PRIOR = (0.10, 0.91)


def power_law(nu, amplitude, beta):
    return amplitude * (np.asarray(nu, dtype=float) / NU0_MHZ) ** beta


def fit_power_law(nu, t, s) -> Optional[Dict[str, object]]:
    """Weighted power-law fit to the points with ``t > 0``.  ``None`` if < 2 points."""
    from scipy.optimize import curve_fit

    nu, t, s = (np.asarray(x, dtype=float) for x in (nu, t, s))
    ok = np.isfinite(t) & np.isfinite(s) & (t > 0) & (s > 0)
    if ok.sum() < 2:
        return None
    slope, intercept = np.polyfit(np.log(nu[ok] / NU0_MHZ), np.log(t[ok]), 1, w=t[ok] / s[ok])
    try:
        popt, pcov = curve_fit(power_law, nu[ok], t[ok], p0=[np.exp(intercept), slope],
                               sigma=s[ok], absolute_sigma=True, maxfev=20000)
    except RuntimeError:
        return None
    model = power_law(nu, *popt)
    chi2 = float(np.sum(((t[ok] - model[ok]) / s[ok]) ** 2))
    n = int(ok.sum())
    err = np.sqrt(np.diag(pcov)) if np.all(np.isfinite(pcov)) else np.array([np.nan, np.nan])
    return dict(A=float(popt[0]), beta=float(popt[1]), sigma_A=float(err[0]),
                sigma_beta=float(err[1]), chi2=chi2, n=n,
                chi2_dof=chi2 / (n - 2) if n > 2 else np.nan, used=ok, model=model)


def sigma_for(m: Dict[str, object], t_excess):
    """Error bar for a non-TRIS point: ARCADE 2 calibration, else the assumed fraction."""
    t = np.asarray(t_excess, dtype=float)
    if m["freq_mhz"] in ARCADE_CAL_K:
        return np.full_like(t, ARCADE_CAL_K[m["freq_mhz"]])
    return ASSUMED_FRACTION * np.abs(t)


def ring_sigma(ring) -> np.ndarray:
    """TRIS ring error: statistical (floored) ⊕ published zero level (larger side)."""
    stat = np.asarray(ring.statistical_uncertainty_k, dtype=float)
    stat = np.where(stat > 0, stat, stat[stat > 0].min())
    z = ring.zero_level_uncertainty_k
    zero = float(z) if np.isscalar(z) else max(z.positive_k, z.negative_k)
    return np.hypot(stat, zero)


class Inputs:
    """Everything both routes need, built once."""

    def __init__(self, config: Config, res_dir: Path,
                 target_file: str = "HASLAM_DESTRIPED_ONLY.fits", verbose: bool = True):
        import healpy as hp

        self.config = config
        self.tris = import_limtod_tris(config)
        self.nside = int(config.require("tris.nside"))
        self.target_fwhm = float(config.require("beam_matching.target_fwhm_deg"))
        self.weight_floor = float(config.require("beam_matching.weight_floor"))
        self.rule = str(config.require("beam_matching.coarse_pixel_rule"))
        self.stage1 = load_stage1(config)
        self.band = self.stage1["band_mask"].astype(bool)
        self.all_maps = sky_maps.discover(res_dir)
        self.target = next(m for m in self.all_maps if m["file"] == target_file)
        self.fit_maps = [m for m in self.all_maps if m["usable"] and m["role"] == "fit"]

        if verbose:
            print("building the TRIS-grid (beam-matched) and raw nside {} maps...".format(self.nside))
        self.grid: Dict[str, np.ndarray] = {}
        self.raw16: Dict[str, np.ndarray] = {}
        for m in self.fit_maps + [self.target]:
            self.grid[m["name"]], self.raw16[m["name"]] = self._to_grid(m)

        archive = config.resolve_path("paths.tris_archive_dir")
        self.rings = {f: self.tris.read_tris_ring(archive / "TRIS_absolute_{}.txt".format(f))
                      for f in (600, 820)}
        self.ra = np.asarray(self.rings[600].ra_deg)
        self.cuts = self.tris.read_tris_beam_cuts(archive / "TRIS_Beam_Profile.txt")
        self.geom = self.tris.tris_zenith_geometry(self.ra)
        self.hires = int(config.require("tris.nside_hires"))
        self.beam = self.tris.tris_cut_beam_map(self.cuts, nside=self.nside, normalization="peak")
        self.horizon = self.tris.tris_horizon_mask(self.nside)
        if verbose:
            print("simulating TRIS TODs of every map with limTOD...")
        self.tods: Dict[str, Dict[str, np.ndarray]] = {}
        for m in self.fit_maps + [self.target]:
            tod, cov = self._native_tod(m)
            self.tods[m["name"]] = dict(tod=tod, coverage=cov, good=cov >= MIN_COVERAGE)

        self.tris_map: Dict[float, tuple] = {}
        for f in (600, 820):
            i = self.stage1.index(f)
            nu = float(self.stage1["effective_freq_mhz"][i])
            t = self.stage1.band_map("sky_k", f) + float(self.stage1["zero_level_k"][i])
            s = np.hypot(self.stage1.band_map("sky_sigma_k", f),
                         float(self.stage1["zero_level_sigma_k"][i]))
            self.tris_map[nu] = (t, s, f)
        self.npix = hp.nside2npix(self.nside)

    # -- maps onto the TRIS grid ------------------------------------------------
    def _rotated(self, m):
        import healpy as hp

        sky = m["sky"]
        observed = np.isfinite(sky)
        if m["frame"] == "C":
            return np.where(observed, sky, hp.UNSEEN), observed
        if observed.all():
            return frames.rotate_full_sky(sky, (m["frame"], "C")), observed
        rotated, obs_rot, _ = frames.rotate_masked(sky, observed, (m["frame"], "C"))
        return rotated, obs_rot

    def _to_grid(self, m):
        import healpy as hp

        rotated, observed = self._rotated(m)
        raw_values, raw_kept = regrid_to_nside(np.where(observed, rotated, hp.UNSEEN),
                                               observed, self.nside, self.rule)
        smoothed = smooth_to_target(rotated, m["fwhm_deg"], self.target_fwhm,
                                    weight_mask=None if observed.all() else observed,
                                    weight_floor=self.weight_floor)
        final = observed & (smoothed != hp.UNSEEN)
        values, kept = regrid_to_nside(np.where(final, smoothed, hp.UNSEEN), final,
                                       self.nside, self.rule)
        return np.where(kept, values, np.nan), np.where(raw_kept, raw_values, np.nan)

    # -- TODs -----------------------------------------------------------------------
    def observe(self, sky_eq):
        from limTOD import generate_TOD_sky

        g = self.geom
        return np.asarray(generate_TOD_sky(
            self.beam, sky_eq, g.lst_deg, g.latitude_deg, g.azimuth_deg, g.elevation_deg,
            g.selfrot_deg, nside_hires=self.hires, normalize_beam=True,
            horizontal_mask=self.horizon))

    def _native_tod(self, m):
        import healpy as hp

        observed = np.isfinite(m["sky"])
        filled, weight = np.where(observed, m["sky"], 0.0), observed.astype(float)
        if m["nside"] > 64:
            filled, weight = hp.ud_grade(filled, 64), hp.ud_grade(weight, 64)
        if m["frame"] != "C":
            rot = hp.Rotator(coord=[m["frame"], "C"])
            filled, weight = rot.rotate_map_pixel(filled), rot.rotate_map_pixel(weight)
        summed = self.observe(hp.ud_grade(filled, self.nside))
        coverage = self.observe(hp.ud_grade(weight, self.nside))
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(coverage > 0, summed / coverage, np.nan), coverage


def _points_at_pixel(inp: Inputs, p: int):
    nus, ts, ss, names = [], [], [], []
    for nu, (t, s, nominal) in inp.tris_map.items():
        nus.append(nu); ts.append(t[p] - sky_maps.cmb_rj_k(nu)); ss.append(s[p])
        names.append("TRIS {} MHz (stage 1 map)".format(nominal))
    for m in inp.fit_maps:
        v = inp.grid[m["name"]][p]
        if np.isfinite(v):
            t_ex = float(sky_maps.rj_excess(v, m))
            nus.append(m["freq_mhz"]); ts.append(t_ex); ss.append(float(sigma_for(m, t_ex)))
            names.append(m["name"])
    return np.array(nus), np.array(ts), np.array(ss), names


def map_route(inp: Inputs) -> Dict[str, object]:
    """Per-pixel fit over the TRIS stripe on the TRIS grid."""
    cmb408 = sky_maps.cmb_rj_k(NU0_MHZ)
    npix = inp.npix
    out = {k: np.full(npix, np.nan) for k in ("pred", "beta", "sigma_beta", "sigma_pred",
                                               "chi2_dof", "n")}
    resid_frac: Dict[str, np.ndarray] = {}
    resid_norm: Dict[str, np.ndarray] = {}
    fit_minus_map: Dict[str, np.ndarray] = {}
    points: Dict[int, tuple] = {}
    for p in np.flatnonzero(inp.band):
        nu, t, s, names = _points_at_pixel(inp, p)
        r = fit_power_law(nu, t, s)
        if r is None:
            continue
        out["pred"][p] = r["A"] + cmb408
        out["beta"][p], out["sigma_beta"][p] = r["beta"], r["sigma_beta"]
        out["sigma_pred"][p], out["chi2_dof"][p], out["n"][p] = r["sigma_A"], r["chi2_dof"], r["n"]
        for name, ti, si, mi in zip(names, t, s, r["model"]):
            resid_frac.setdefault(name, np.full(npix, np.nan))[p] = (ti - mi) / mi
            resid_norm.setdefault(name, np.full(npix, np.nan))[p] = (ti - mi) / si
            fit_minus_map.setdefault(name, np.full(npix, np.nan))[p] = mi - ti
        points[p] = (nu, t, s, names, r)
    target = inp.grid[inp.target["name"]]
    out.update(old=target, delta=out["pred"] - target, resid_frac=resid_frac,
               resid_norm=resid_norm, fit_minus_map=fit_minus_map, points=points)
    return out


def tod_route(inp: Inputs) -> Dict[str, object]:
    """Per-sample fit along the ring, and the correction mapped by the Wiener filter."""
    from limTOD import wiener_filter_map

    cmb408 = sky_maps.cmb_rj_k(NU0_MHZ)
    n_t = inp.ra.size
    pred, sig, beta, chi2, n = (np.full(n_t, np.nan) for _ in range(5))
    points: Dict[int, tuple] = {}
    for k in range(n_t):
        nus, ts, ss, names = [], [], [], []
        for f, ring in inp.rings.items():
            nu = ring.effective_frequency_mhz
            nus.append(nu); ts.append(ring.temperature_k[k] - sky_maps.cmb_rj_k(nu))
            ss.append(ring_sigma(ring)[k]); names.append("TRIS {} MHz (ring data)".format(f))
        for m in inp.fit_maps:
            d = inp.tods[m["name"]]
            if d["good"][k]:
                t_ex = float(sky_maps.rj_excess(d["tod"][k], m))
                nus.append(m["freq_mhz"]); ts.append(t_ex); ss.append(float(sigma_for(m, t_ex)))
                names.append(m["name"])
        r = fit_power_law(nus, ts, ss)
        if r:
            points[k] = (np.array(nus), np.array(ts), np.array(ss), names, r)
            pred[k], sig[k], beta[k], chi2[k], n[k] = (r["A"] + cmb408, r["sigma_A"], r["beta"],
                                                       r["chi2_dof"], r["n"])
    target_tod = inp.tods[inp.target["name"]]["tod"]
    delta_tod = pred - target_tod

    tris = inp.tris
    inputs = tris.build_tris_mapmaking_inputs(
        inp.rings[600], nside=inp.nside, pixel_indices=inp.stage1.pixel_indices, cuts=inp.cuts,
        uncertainty_floor_k=0.004,
        zero_level_sigma_k=float(inp.config.require("tris.zero_level_sigma_k")),
        apply_horizon_mask=bool(inp.config.require("beam.apply_horizon_mask")),
        nside_hires=inp.hires)
    B = np.asarray(inputs.operator)[:, :-1]
    pix = inp.stage1.pixel_indices
    old = inp.grid[inp.target["name"]]
    prior_sigma = np.hypot(CORRECTION_PRIOR[0] * old[pix], CORRECTION_PRIOR[1])
    ok = np.isfinite(delta_tod) & np.isfinite(sig) & (sig > 0)
    delta_band, delta_sigma = wiener_filter_map(
        delta_tod[ok], B[ok], noise_variance=sig[ok] ** 2,
        prior_inv_cov=prior_sigma ** -2.0, guess=np.zeros(pix.size))
    delta_map = np.full(inp.npix, np.nan); delta_map[pix] = delta_band
    sigma_map = np.full(inp.npix, np.nan); sigma_map[pix] = delta_sigma
    return dict(pred=pred, sigma_pred=sig, beta=beta, chi2_dof=chi2, n=n, points=points,
                target_tod=target_tod, delta_tod=delta_tod, delta=delta_map,
                sigma_delta=sigma_map, old=old, pred_map=old + delta_map)


def apply_full_resolution(inp: Inputs, delta_map_eq: np.ndarray) -> Dict[str, np.ndarray]:
    """Corrected full-resolution target: target + the nside correction, interpolated.

    The target is rotated into equatorial at its native nside (harmonic), the
    nside-``inp.nside`` correction is interpolated bilinearly onto the native grid
    (so it is smooth, not blocky), and added inside the TRIS stripe; outside it,
    the target is left as it was.
    """
    import healpy as hp

    native = frames.rotate_full_sky(inp.target["sky"], (inp.target["frame"], "C"))
    nside = hp.get_nside(native)
    theta, phi = hp.pix2ang(nside, np.arange(native.size))
    stripe = inp.band[hp.ang2pix(inp.nside, theta, phi)]
    filled = np.where(np.isfinite(delta_map_eq), delta_map_eq, 0.0)
    correction = hp.get_interp_val(filled, theta, phi)
    corrected = np.where(stripe, native + correction, native)
    return dict(old=native, corrected=corrected, correction=np.where(stripe, correction, np.nan),
                stripe=stripe)
