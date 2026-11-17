"""The end-to-end pipeline, run with ``python run_pipeline.py --step <name>`` (or ``--step all``).

Steps (``STEPS``), each reading the products of the steps before it:

0. ``tris_maps``   TRIS maps from the two archive rings (stage 1, :mod:`moonrunes.stage1_tris_maps`).  Its band
                   mask is used by everything that builds :class:`moonrunes.sed_prototype.Inputs`, and its maps by
                   the paper's Figure 3.  Written to ``outputs/stage1/`` (``run.outputs_dir``).
1. ``beam_match``  every allowed map (``res/`` minus ``sed.exclude_maps``, plus the Haslam target) beam-matched to
                   the TRIS beam on the TRIS grid (``pipeline.beam_match_kernel``, the measured ``tris_cuts`` beam).
2. ``ring_tods``   ring TODs (the 120 TRIS ring samples) of every allowed map, the TRIS rings and Haslam, with
                   their error bars: the SED points of every later step.  Also the sky models used for comparison
                   (GSM2008, GSM2016, LFSM and the Carter / Irfan / GINAN corrections) at 408 MHz and at the TRIS
                   frequencies, through the same beam.  The TODs use the same kernel at the ring pointings
                   (``sed.tod_beam: beam_matched``); they are not sampled from the step 1 maps.
3. ``compare``     every spectral model of ``configs/model_comp_config.yaml`` fitted at every ring sample to the
                   allowed maps, **TRIS excluded** (``comparison.exclude_tris``); UltraNest evidence by default
                   (``comparison.statistic``).
4. ``select``      the best model from step 3 (``pipeline.select``).
5. ``fit``         the best model fitted at every ring sample to every allowed data set, **TRIS included**,
                   with UltraNest (posterior kept).  These fits are our spectral model.
6. ``predict_408`` our model's 408 MHz TOD through the TRIS beam, against Haslam and the 408 MHz sky models:
                   residuals and a scale-offset (York) fit.
7. ``rhino``       our model at the RHINO frequencies: through the TRIS beam, deconvolved (Wiener, prior = a sky
                   model) and observed with the RHINO beam at Jodrell Bank, with the sky models through both beams.

``report`` collects the summary every step writes into its ``manifest.json`` and writes ``PIPELINE_REPORT.md``.

Products go to ``pipeline.output_dir`` (``paper_outputs/pipeline/``), one folder per step (``STEP_DIRS``).  The
file formats are the ones the notebooks already read (``notebooks/10_paper_plots.ipynb``,
``notebooks/test_results.ipynb``), so a product can be swapped for its notebook-made predecessor.
Everything is configured in ``configs/run_pipeline.yaml`` (the data, beam, ring, corrections and ``pipeline:``) and,
for the models, priors, sampler, error bars and the comparison's statistic, ``configs/model_comp_config.yaml``
(``pipeline.model_comparison_config``).  No setting is in both.
"""

from __future__ import annotations

import copy
import csv
import datetime as _dt
import hashlib
import json
import logging
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

from .config import (REPO_ROOT, Config, ConfigError, error_bars_path, load_config, load_model_comp_config,
                     model_comp_config_path)

os.environ.setdefault("TQDM_DISABLE", "1")

STEPS = ("tris_maps", "beam_match", "ring_tods", "compare", "select", "fit", "predict_408", "rhino")
STEP_DIRS = {
    "beam_match": "step1_beam_matched_maps",
    "ring_tods": "step2_ring_tods",
    "compare": "step3_model_comparison",
    "select": "step4_model_selection",
    "fit": "step5_fit",
    "predict_408": "step6_predict_408",
    "rhino": "step7_rhino",
}
REPORT_NAME = "PIPELINE_REPORT.md"

# The 408 MHz sky models, in the order every reader of sky_tods.npz expects
# (test_results.ipynb Parts B, E-H: 0 Haslam, 1 Carter, 2 Irfan, 3 GSM2008, 4 GSM2016, 5 GINAN, 6 LFSM).
SKY_408 = ["Haslam (uncorrected)", "Haslam, Carter correction", "Haslam, Irfan correction",
           "GSM2008", "GSM2016", "GSM2016, GINAN correction", "LFSM"]
# The sky models carried to the TRIS frequencies (the paper's Figure 4).
SKY_TRIS = ["GSM2008", "GSM2016", "GSM2008, Carter correction", "GSM2008, Irfan correction",
            "GSM2016, Carter correction", "GSM2016, Irfan correction", "GSM2016, GINAN correction",
            "LFSM, stretched"]


# ---------------------------------------------------------------------------
# configuration and bookkeeping
# ---------------------------------------------------------------------------
def pcfg(config: Config, dotted: str):
    """``pipeline.<dotted>`` from run_pipeline.yaml (raises if missing: an open decision)."""
    return config.require("pipeline." + dotted)


def output_root(config: Config) -> Path:
    return config.resolve_path("pipeline.output_dir")


def step_dir(config: Config, step: str, create: bool = True) -> Path:
    d = output_root(config) / STEP_DIRS[step]
    if create:
        d.mkdir(parents=True, exist_ok=True)
    return d


def mc_config(config: Config) -> dict:
    """``configs/model_comp_config.yaml`` (``pipeline.model_comparison_config``)."""
    return load_model_comp_config(config=config)


def mc_errors(config: Config, rc) -> dict:
    """The error-bar rows (``data.error_bars_csv``, relative to the model-comparison config's folder)."""
    return rc.read_error_bars(error_bars_path(mc_config(config), model_comp_config_path(config)))


def _mc_modules():
    """``scripts/model_comparison/models.py`` and ``run_comparison.py`` (data helpers; their main() is not run)."""
    mc_dir = str(REPO_ROOT / "scripts" / "model_comparison")
    if mc_dir not in sys.path:
        sys.path.insert(0, mc_dir)
    import models  # noqa: E402
    import run_comparison as rc  # noqa: E402

    return models, rc


def _jsonable(x):
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, (np.floating, np.integer, np.bool_)):
        return x.item()
    if isinstance(x, Path):
        return str(x)
    return x


def write_manifest(d: Path, step: str, config: Config, summary: dict, products: List[str],
                   extra: Optional[dict] = None) -> dict:
    """``manifest.json`` of a step: what made it, with which settings, and its summary (read by ``report``)."""
    m = {
        "step": step,
        "step_number": STEPS.index(step),
        "written_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "produced_by": "run_pipeline.py --step {}".format(step),
        "config_path": str(config.path) if config.path else None,
        "pipeline_config": config.section("pipeline"),
        "products": products,
        "summary": summary,
    }
    if extra:
        m.update(extra)
    (d / "manifest.json").write_text(json.dumps(_jsonable(m), indent=2) + "\n")
    return m


def read_manifest(d: Path) -> Optional[dict]:
    p = d / "manifest.json"
    return json.loads(p.read_text()) if p.exists() else None


def _require(path: Path, step: str) -> Path:
    if not path.exists():
        raise FileNotFoundError("{} is missing: run `python run_pipeline.py --step {}` first".format(path, step))
    return path


def with_kernel(config: Config, kernel: str) -> Config:
    """A copy of ``config`` with ``beam_matching.kernel`` set (run_pipeline.yaml is not edited)."""
    data = copy.deepcopy(dict(config))
    data["beam_matching"]["kernel"] = kernel
    return Config(data, config.path)


class Context:
    """One pipeline run: the config, and the TOD-route inputs built once and shared by steps 1, 2 and 7."""

    def __init__(self, config: Config, verbose: bool = True, overwrite: bool = False):
        self.config = config
        self.verbose = verbose
        self.overwrite = overwrite
        self._inputs = None

    def say(self, *a):
        if self.verbose:
            print(*a, flush=True)

    def inputs(self):
        if self._inputs is None:
            from .sed_prototype import Inputs

            kernel = str(pcfg(self.config, "beam_match_kernel"))
            self.say("building the TOD-route inputs (every allowed map beam-matched; a few minutes)...")
            self._inputs = Inputs(with_kernel(self.config, kernel), REPO_ROOT / "res",
                                  target_file=str(pcfg(self.config, "target_file")), verbose=self.verbose)
        return self._inputs


# ---------------------------------------------------------------------------
# step 0: TRIS maps (stage 1)
# ---------------------------------------------------------------------------
def step_tris_maps(ctx: Context) -> dict:
    from .stage1_tris_maps import run_stage1

    manifest = run_stage1(config=ctx.config, overwrite=True if ctx.overwrite else None, verbose=ctx.verbose)
    return {"note": "stage 1 manifest: outputs/stage1/stage1_manifest.json",
            "frequencies": {f: {k: d.get(k) for k in ("effective_frequency_mhz", "zero_level_k", "reduced_chi_square")}
                            for f, d in manifest["frequencies"].items()}}


# ---------------------------------------------------------------------------
# step 1: beam-matched maps
# ---------------------------------------------------------------------------
def step_beam_match(ctx: Context) -> dict:
    inp = ctx.inputs()
    d = step_dir(ctx.config, "beam_match")
    maps = inp.fit_maps + [inp.target]
    names = [m["name"] for m in maps]
    grid = np.vstack([inp.grid[n] for n in names])
    np.savez_compressed(
        d / "beam_matched_maps.npz", names=np.array(names), files=np.array([m["file"] for m in maps]),
        freq_mhz=np.array([float(m["freq_mhz"]) for m in maps]),
        native_fwhm_deg=np.array([float(m["fwhm_deg"]) for m in maps]),
        nside=np.int64(inp.nside), frame=np.array("C"), kernel=np.array(inp.kernel),
        maps=grid, band_mask=inp.band)
    rows = []
    for m, row in zip(maps, grid):
        rows.append({"map": m["name"], "file": m["file"], "freq_mhz": float(m["freq_mhz"]),
                     "native_fwhm_deg": float(m["fwhm_deg"]), "role": "target" if m is inp.target else "fit",
                     "pixels_kept": int(np.isfinite(row).sum()),
                     "pixels_kept_in_tris_band": int(np.isfinite(row[inp.band]).sum()),
                     "band_pixels": int(inp.band.sum())})
    _write_csv(d / "beam_matched_maps_summary.csv", rows)
    summary = {"kernel": inp.kernel, "nside": inp.nside, "maps": len(maps),
               "excluded_by_config": inp.excluded,
               "per_map": {r["file"]: {"pixels_kept_in_tris_band": r["pixels_kept_in_tris_band"],
                                       "of": r["band_pixels"]} for r in rows}}
    write_manifest(d, "beam_match", ctx.config, summary,
                   ["beam_matched_maps.npz", "beam_matched_maps_summary.csv"],
                   {"smoothing": _jsonable(inp.smoothing)})
    return summary


# ---------------------------------------------------------------------------
# step 2: ring TODs
# ---------------------------------------------------------------------------
_IRFAN_MAPS: Dict[tuple, tuple] = {}


def irfan_key(sky_cfg) -> dict:
    """What defines the Irfan correction (for cache keys)."""
    return {"zero_level_k": float(sky_cfg["irfan_monopole_k"]), "nu0_mhz": float(sky_cfg["irfan_nu0_mhz"]),
            "beta_map": str(sky_cfg["irfan_beta_map"]), "cs_map": str(sky_cfg["irfan_cs_map"])}


def irfan_scaling(sky_cfg, nside: int, nu_mhz: float, from_mhz: float = 408.0) -> np.ndarray:
    """Per-pixel T(nu) / T(from_mhz) of Irfan & Puglisi's synchrotron model (their eqs. 1, 2 and 11):
    T(p, nu) = T(p, nu0) (nu/nu0)^(beta_s(p) + c_s(p) ln(nu/nu0)), nu0 = 45 MHz, with their beta_s and c_s maps.
    The maps are HEALPix nside 256, RING, Galactic (the frame checked against Haslam's small-scale structure);
    they are brought to ``nside`` with hp.ud_grade (pixel copy upward, mean downward).  Their eq. 6 applies this
    to (T_tot - 8.9 K) at 408 MHz, which is what the corrected 408 MHz maps here are."""
    import healpy as hp

    key = (str(sky_cfg["irfan_beta_map"]), str(sky_cfg["irfan_cs_map"]), int(nside))
    if key not in _IRFAN_MAPS:
        beta = np.load(REPO_ROOT / sky_cfg["irfan_beta_map"]).astype(float)
        cs = np.load(REPO_ROOT / sky_cfg["irfan_cs_map"]).astype(float)
        _IRFAN_MAPS[key] = (hp.ud_grade(beta, nside), hp.ud_grade(cs, nside))
    beta, cs = _IRFAN_MAPS[key]
    nu0 = float(sky_cfg["irfan_nu0_mhz"])
    x1, x0 = np.log(float(nu_mhz) / nu0), np.log(float(from_mhz) / nu0)
    return np.exp(beta * (x1 - x0) + cs * (x1 ** 2 - x0 ** 2))


def _sky_model_maps_408(inp, sky_cfg) -> Dict[str, tuple]:
    """name -> (CMB-removed RJ map, frame, native FWHM) at 408 MHz, in SKY_408 order (as test_results Part A/B)."""
    from pygdsm import GlobalSkyModel, GlobalSkyModel16, LowFrequencySkyModel
    from pygdsm.mckay25 import generate_mckay25_correction, T_CMB

    from . import sky_maps

    h = inp.target
    e_h = np.asarray(sky_maps.rj_excess(h["sky"], h), float)
    cmb408 = float(sky_maps.cmb_rj_k(408.0))
    carter, irfan = sky_cfg["carter"], float(sky_cfg["irfan_monopole_k"])
    ginan_f = float(sky_cfg["ginan_freq_mhz"])
    native = sky_cfg["native_fwhm_deg"]
    g16 = GlobalSkyModel16(freq_unit="MHz", data_unit="TRJ", include_cmb=False)
    g16_408 = np.array(g16.generate(408.0), float)
    g16_lo = np.array(g16.generate(ginan_f), float)
    off, scl = generate_mckay25_correction(np.array([ginan_f / 1e3]))
    g16_cmb = np.array(GlobalSkyModel16(freq_unit="MHz", data_unit="TCMB", include_cmb=True).generate(ginan_f), float)
    hf = float(h["fwhm_deg"])
    return {
        "Haslam (uncorrected)": (e_h, h["frame"], hf),
        "Haslam, Carter correction": (float(carter["scale"]) * e_h + float(carter["offset_k"]), h["frame"], hf),
        "Haslam, Irfan correction": (e_h - (irfan - cmb408), h["frame"], hf),
        "GSM2008": (np.array(GlobalSkyModel(freq_unit="MHz", include_cmb=False).generate(408.0), float), "G",
                    float(native["GSM2008"])),
        "GSM2016": (g16_408, "G", float(native["GSM2016"])),
        "GSM2016, GINAN correction": (((g16_cmb - off[0]) * scl[0] - T_CMB) * g16_408 / g16_lo, "G",
                                      float(native["GSM2016"])),
        "LFSM": (np.array(LowFrequencySkyModel(freq_unit="MHz", include_cmb=False).generate(408.0), float), "G",
                 float(native["LFSM"])),
    }


def _sky_model_maps_tris(inp, sky_cfg, maps408, nu_tris: Dict[str, float]) -> Dict[str, Dict[str, tuple]]:
    """key_t -> name -> (map, frame, native FWHM): GSM2008/16 at the TRIS frequencies; the Carter and GINAN
    corrections carried there from 408 MHz with the base model's per-pixel spectrum, the Irfan correction with
    Irfan & Puglisi's own beta_s and c_s maps (``irfan_scaling``); LFSM extended with its own 300-408 MHz index."""
    import healpy as hp
    from pygdsm import GlobalSkyModel, GlobalSkyModel16, LowFrequencySkyModel

    from . import sky_maps

    native = sky_cfg["native_fwhm_deg"]
    carter, irfan = sky_cfg["carter"], float(sky_cfg["irfan_monopole_k"])
    cmb408 = float(sky_maps.cmb_rj_k(408.0))
    g08_408 = maps408["GSM2008"][0]
    g16_408 = maps408["GSM2016"][0]
    base408 = {"GSM2008": g08_408, "GSM2016": g16_408}
    corrected408 = {
        ("GSM2008", "Carter correction"): float(carter["scale"]) * g08_408 + float(carter["offset_k"]),
        ("GSM2008", "Irfan correction"): g08_408 - (irfan - cmb408),
        ("GSM2016", "Carter correction"): float(carter["scale"]) * g16_408 + float(carter["offset_k"]),
        ("GSM2016", "Irfan correction"): g16_408 - (irfan - cmb408),
        ("GSM2016", "GINAN correction"): maps408["GSM2016, GINAN correction"][0],
    }
    lo = float(sky_cfg["lfsm_beta_from_mhz"])
    lfsm = LowFrequencySkyModel(freq_unit="MHz", include_cmb=False)
    lfsm_408 = np.array(lfsm.generate(408.0), float)
    lfsm_beta = np.log(lfsm_408 / np.array(lfsm.generate(lo), float)) / np.log(408.0 / lo)
    g08 = GlobalSkyModel(freq_unit="MHz", include_cmb=False)
    g16 = GlobalSkyModel16(freq_unit="MHz", data_unit="TRJ", include_cmb=False)
    out = {}
    for key_t, nu in nu_tris.items():
        at_nu = {"GSM2008": np.array(g08.generate(nu), float), "GSM2016": np.array(g16.generate(nu), float)}
        skies = {n: (v, "G", float(native[n])) for n, v in at_nu.items()}
        for (base, corr), m408 in corrected408.items():
            if corr == "Irfan correction":
                scale = irfan_scaling(sky_cfg, hp.npix2nside(m408.size), nu)
            else:
                scale = at_nu[base] / base408[base]
            skies["{}, {}".format(base, corr)] = (m408 * scale, "G", float(native[base]))
        skies["LFSM, stretched"] = (lfsm_408 * (nu / 408.0) ** lfsm_beta, "G", float(native["LFSM"]))
        out[key_t] = {n: skies[n] for n in SKY_TRIS}
    return out


def step_ring_tods(ctx: Context) -> dict:
    import healpy as hp

    from . import elliptical_beam, sky_maps
    from .sed_prototype import MIN_COVERAGE, sigma_for, tod_method_key

    config = ctx.config
    _, rc = _mc_modules()
    mc = mc_config(config)
    inp = ctx.inputs()
    d = step_dir(config, "ring_tods")
    errors = mc_errors(config, rc)
    missing = [m["file"] for m in inp.fit_maps if m["file"] not in errors]
    if missing:
        raise KeyError("no error-bar row in {} for: {}".format(
            error_bars_path(mc, model_comp_config_path(config)), missing))

    # ---- the data: TRIS rings + every allowed map, as (nu, T, sigma) per ring sample (run_comparison.tod_units)
    ra = np.asarray(inp.ra, dtype=float)
    n = ra.size
    keys = [rc.tris_key(f) for f in rc.TRIS_BANDS] + [m["file"] for m in inp.fit_maps]
    nu = {}
    T = {k: np.full(n, np.nan) for k in keys}
    S = {k: np.full(n, np.nan) for k in keys}
    for uid, _, pts in rc.tod_units(inp, errors):
        for key, f, t, s in pts:
            nu[key], T[key][uid], S[key][uid] = float(f), float(t), float(s)
    for m in inp.fit_maps:                                   # a map with no valid sample still needs its frequency
        nu.setdefault(m["file"], float(m["freq_mhz"]))
    tgt = inp.target
    h_ex = np.asarray(sky_maps.rj_excess(np.asarray(inp.tods[tgt["name"]]["tod"], float), tgt), dtype=float)
    h_good = np.asarray(inp.tods[tgt["name"]]["good"], dtype=bool)
    h_T = np.where(h_good, h_ex, np.nan)
    h_S = np.where(h_good, np.asarray(sigma_for(tgt, h_ex), dtype=float), np.nan)
    np.savez(d / "ring_tods.npz", ra=ra, keys=np.array(keys), nu=np.array([nu[k] for k in keys]),
             T=np.array([T[k] for k in keys]), S=np.array([S[k] for k in keys]),
             target_file=np.array(tgt["file"]), haslam_nu=np.float64(tgt["freq_mhz"]), haslam_T=h_T, haslam_S=h_S,
             tod_key=np.array(tod_method_key(config)))
    long = {"sample": np.arange(n), "ra_deg": ra}
    for k in keys:
        long["T_K " + k], long["sigma_K " + k] = T[k], S[k]
    long["T_K " + tgt["file"]], long["sigma_K " + tgt["file"]] = h_T, h_S
    _write_csv(d / "ring_tods.csv", [dict(zip(long, vals)) for vals in zip(*long.values())])

    # ---- every map and every sky model through the same beam, in notebook 18's format (beam_matched_tods.npz)
    spec = elliptical_beam.beam_spec_from_config(config, kernel="tris_cuts")
    floor = float(config.require("beam_matching.weight_floor"))
    out, meta = {"ra_deg": ra}, []

    def record(kind, name, freq, native, tod, w, key):
        out[key] = tod
        k = elliptical_beam.kernel_for(spec, float(native))
        meta.append(dict(kind=kind, name=name, freq_mhz=float(freq), native_fwhm_deg=float(native), kernel=k.shape,
                         kernel_e_deg=k.fwhm_e_deg, kernel_h_deg=k.fwhm_h_deg, n_valid=int(np.isfinite(tod).sum())))

    for m in inp.fit_maps + [tgt]:
        tod = np.asarray(inp.tods[m["name"]]["tod"], float)
        w = np.asarray(inp.tods[m["name"]]["coverage"], float)
        t = np.where(w >= floor, np.asarray(sky_maps.rj_excess(tod, m), float), np.nan)
        out["weight::survey::" + m["file"]] = w
        record("survey", m["file"], m["freq_mhz"], m["fwhm_deg"], t, w, "survey::" + m["file"])

    def observe(sky, frame, native):
        m = dict(sky=np.asarray(sky, float), nside=hp.npix2nside(len(sky)), frame=frame, fwhm_deg=float(native))
        return inp.ring_tod(m)

    sky_cfg = config.section("pipeline")["sky_models"]
    ctx.say("sky models at 408 MHz through the TRIS beam...")
    maps408 = _sky_model_maps_408(inp, sky_cfg)
    raw408 = {}
    for name in SKY_408:
        sky, fr, native = maps408[name]
        tod, w = observe(sky, fr, native)
        raw408[name] = (np.asarray(tod, float), np.asarray(w, float))
        out["weight::model408::" + name] = raw408[name][1]
        record("model 408", name, 408.0, native, np.where(raw408[name][1] >= floor, raw408[name][0], np.nan),
               raw408[name][1], "model408::" + name)
    nu_tris = {rc.tris_key(f): float(inp.rings[f].effective_frequency_mhz) for f in rc.TRIS_BANDS}
    ctx.say("sky models at the TRIS frequencies {} through the TRIS beam...".format(list(nu_tris.values())))
    for key_t, skies in _sky_model_maps_tris(inp, sky_cfg, maps408, nu_tris).items():
        for name, (sky, fr, native) in skies.items():
            tod, w = observe(sky, fr, native)
            record("model " + key_t, name, nu_tris[key_t], native, np.where(np.asarray(w) >= floor, tod, np.nan),
                   w, "modeltris::{}::{}".format(name, key_t))
    np.savez(d / "beam_matched_tods.npz", **out)
    _write_csv(d / "beam_matched_tods_meta.csv", meta)
    _write_csv(d / "beam_matched_tods.csv",
               [dict(zip(["ra_deg"] + [k for k in out if not k.startswith("weight::") and k != "ra_deg"], vals))
                for vals in zip(ra, *[out[k] for k in out if not k.startswith("weight::") and k != "ra_deg"])])
    # the 408 MHz sky models in the format test_results.ipynb Parts B and E-I read (raw TOD + coverage)
    np.savez(d / "sky_tods.npz", names=np.array(SKY_408),
             **{"tod_{}".format(i): raw408[n][0] for i, n in enumerate(SKY_408)},
             **{"cov_{}".format(i): raw408[n][1] for i, n in enumerate(SKY_408)})

    summary = {
        "ring_samples": int(n), "ra_deg": [float(ra.min()), float(ra.max())],
        "data_sets": {k: {"freq_mhz": nu[k], "samples_with_a_point": int(np.isfinite(T[k]).sum()),
                          "median_T_K": _nanmedian(T[k]), "median_sigma_K": _nanmedian(S[k])} for k in keys},
        "haslam": {"file": tgt["file"], "samples": int(np.isfinite(h_T).sum()), "median_T_K": _nanmedian(h_T),
                   "median_sigma_K": _nanmedian(h_S)},
        "sky_models_408": {n: _nanmedian(np.where(raw408[n][1] >= MIN_COVERAGE, raw408[n][0], np.nan))
                           for n in SKY_408},
    }
    write_manifest(d, "ring_tods", config, summary,
                   ["ring_tods.npz", "ring_tods.csv", "beam_matched_tods.npz", "beam_matched_tods.csv",
                    "beam_matched_tods_meta.csv", "sky_tods.npz"])
    return summary


def load_ring_tods(config: Config) -> dict:
    z = np.load(_require(step_dir(config, "ring_tods", False) / "ring_tods.npz", "ring_tods"), allow_pickle=False)
    keys = [str(k) for k in z["keys"]]
    return dict(ra=z["ra"], keys=keys, nu={k: float(v) for k, v in zip(keys, z["nu"])},
                T={k: z["T"][i] for i, k in enumerate(keys)}, S={k: z["S"][i] for i, k in enumerate(keys)},
                target_file=str(z["target_file"]), haslam_T=z["haslam_T"], haslam_S=z["haslam_S"],
                haslam_nu=float(z["haslam_nu"]))


def points(data: dict, unit: int, keys: List[str]):
    """The usable (key, nu, T, sigma) points of one ring sample, in data-set order (run_comparison.usable)."""
    _, rc = _mc_modules()
    return rc.usable([(k, data["nu"][k], float(data["T"][k][unit]), float(data["S"][k][unit]))
                      for k in keys if np.isfinite(data["T"][k][unit])])


# ---------------------------------------------------------------------------
# step 3: model comparison (TRIS excluded)
# ---------------------------------------------------------------------------
def step_compare(ctx: Context) -> dict:
    config = ctx.config
    models_mod, rc = _mc_modules()
    mc = mc_config(config)
    data = load_ring_tods(config)
    d = step_dir(config, "compare")
    ccfg = mc.get("comparison") or {}
    if "statistic" not in ccfg or "exclude_tris" not in ccfg:
        raise ConfigError("comparison.statistic and comparison.exclude_tris must be set in {}".format(
            model_comp_config_path(config)))
    stat = str(ccfg["statistic"])
    if stat not in ("evidence", "max_likelihood"):
        raise ConfigError("comparison.statistic must be evidence or max_likelihood, got {!r}".format(stat))
    exclude_tris = bool(ccfg["exclude_tris"])
    keys = [k for k in data["keys"] if not (exclude_tris and k.startswith("TRIS_"))]
    models = models_mod.build_models(mc)
    min_points = int(mc["data"]["min_points"])
    seed = int(pcfg(config, "compare.seed"))
    rng = np.random.default_rng(seed)
    logging.getLogger("ultranest").setLevel(logging.WARNING)
    ctx.say("step 3: {} for {} on {} (TRIS {})".format(stat, ", ".join(models), ", ".join(keys),
                                                        "excluded" if exclude_tris else "included"))
    rows, skipped, n_points, used = [], 0, 0, set()
    t0 = time.time()
    for u in range(data["ra"].size):
        pts = points(data, u, keys)
        if len(pts) < min_points:
            skipped += 1
            continue
        nu, t, s = (np.array([q[i] for q in pts], dtype=float) for i in (1, 2, 3))
        used.update(q[0] for q in pts)
        n_points += len(pts)
        row = dict(unit=u, ra_deg=float(data["ra"][u]), n_points=len(pts), inputs=";".join(q[0] for q in pts))
        for name, (pnames, fn) in models.items():
            tr = models_mod.prior_transform(pnames, mc["priors"])
            if stat == "max_likelihood":
                lnl, c2, p = rc.max_likelihood(pnames, fn, tr, nu, t, s,
                                               int((mc.get("comparison") or {}).get("ml_starts", 24)), rng)
                row["chi2 " + name], row["lnLmax " + name], row["dof " + name] = c2, lnl, len(pts) - len(pnames)
                for pn, pv in zip(pnames, p):
                    row["{} {}".format(pn, name)] = float(pv)
            else:
                np.random.seed(seed + u)                      # UltraNest draws from numpy's global RNG
                lz, lze = rc.log_evidence(pnames, fn, tr, nu, t, s, mc["sampler"])
                row["logZ " + name], row["logZerr " + name] = lz, lze
        rows.append(row)
        ctx.say("  unit {:>4}: {} points -- {:.0f} s".format(u, len(pts), time.time() - t0))
    if not rows:
        raise RuntimeError("no ring sample had >= {} usable points".format(min_points))
    used_keys = [k for k in keys if k in used]
    meta = dict(n_points=n_points, add_zhat_to_tris=bool(mc["zero_level"]["add_zhat_to_tris"]), placeholder_points=0,
                moment_beta_fixed=mc["models"].get("moment_beta_fixed") or "free",
                tod_beam=str(config.require("sed.tod_beam")), datasets=";".join(used_keys))
    if stat == "max_likelihood":
        rc.write_maxlike(rows, models, d, "tod", meta)
        totals = None
    else:
        _write_csv(d / "evidence_per_unit_tod.csv", rows)
        totals = {m: float(sum(r["logZ " + m] for r in rows)) for m in models}
        errs = {m: float(np.sqrt(sum(r["logZerr " + m] ** 2 for r in rows))) for m in models}
        wins = {m: sum(1 for r in rows if max(models, key=lambda q: r["logZ " + q]) == m) for m in models}
        _write_csv(d / "evidence_summary_tod.csv", [
            dict(model=m, sum_logZ=totals[m], sum_logZ_err=errs[m], logZ_minus_power_law=totals[m] - totals["power law"],
                 units_won=wins[m], n_units=len(rows), route="tod", add_zhat_to_tris=meta["add_zhat_to_tris"],
                 placeholder_points=0, n_points=n_points, moment_beta_fixed=meta["moment_beta_fixed"],
                 tod_beam=meta["tod_beam"], datasets=meta["datasets"]) for m in models])
    summary = {"statistic": stat, "tris_excluded": exclude_tris, "data_sets": used_keys, "units": len(rows),
               "skipped_units": skipped, "points_per_unit": sorted({r["n_points"] for r in rows}),
               "sum_logZ": totals}
    write_manifest(d, "compare", config, summary,
                   ["evidence_per_unit_tod.csv", "evidence_summary_tod.csv"] if stat == "evidence"
                   else ["maxlike_per_unit_tod.csv", "maxlike_summary_tod.csv"])
    return summary


# ---------------------------------------------------------------------------
# step 4: model selection
# ---------------------------------------------------------------------------
def step_select(ctx: Context) -> dict:
    """The best model of step 3.  Evidence: the highest ln Z summed over the ring samples; a model within
    ``pipeline.select.tie_lnz`` of it with fewer parameters is preferred.  Maximum likelihood (no Occam
    penalty of its own): the lowest summed BIC = chi2 + k ln n, with the same tie rule."""
    import pandas as pd

    config = ctx.config
    models_mod, _ = _mc_modules()
    models = models_mod.build_models(mc_config(config))
    d3 = step_dir(config, "compare", False)
    tie = float(pcfg(config, "select.tie_lnz"))
    npar = {m: len(p) for m, (p, _) in models.items()}
    if (d3 / "evidence_per_unit_tod.csv").exists():
        ev = pd.read_csv(d3 / "evidence_per_unit_tod.csv")
        score = {m: float(ev["logZ " + m].sum()) for m in models}            # higher is better
        wins = {m: int((ev[["logZ " + q for q in models]].idxmax(axis=1) == "logZ " + m).sum()) for m in models}
        statistic = "evidence (sum of ln Z over {} ring samples)".format(len(ev))
    else:
        ml = pd.read_csv(_require(d3 / "maxlike_per_unit_tod.csv", "compare"))
        score = {m: float(-(ml["chi2 " + m] + npar[m] * np.log(ml["n_points"])).sum()) for m in models}   # -BIC
        wins = {}
        statistic = "maximum likelihood (minus the summed BIC over {} ring samples)".format(len(ml))
    best = max(score, key=score.get)
    close = [m for m in models if score[best] - score[m] <= tie]
    chosen = min(close, key=lambda m: (npar[m], -score[m]))
    sel = {"model": chosen, "parameters": models[chosen][0], "statistic": statistic, "rule":
           "highest score; any model within tie_lnz = {:g} of it with fewer parameters wins".format(tie),
           "scores": score, "score_minus_best": {m: score[m] - score[best] for m in models},
           "units_won": wins, "candidates_within_tie": close}
    d = step_dir(config, "select")
    (d / "selection.json").write_text(json.dumps(_jsonable(sel), indent=2) + "\n")
    write_manifest(d, "select", config, sel, ["selection.json"])
    ctx.say("step 4: best model = {}".format(chosen))
    return sel


def selected_model(config: Config) -> str:
    return json.loads(_require(step_dir(config, "select", False) / "selection.json", "select").read_text())["model"]


# ---------------------------------------------------------------------------
# step 5: the best model with TRIS included (our spectral model)
# ---------------------------------------------------------------------------
def fit_path(config: Config, unit: int, model: str) -> Path:
    return step_dir(config, "fit", False) / "fits" / "unit{:03d}__{}.npz".format(unit, model.replace(" ", "_"))


def prior_key(mc: dict) -> str:
    """The key test_results.ipynb Part B checks a fit file against (its three base parameters)."""
    return json.dumps({n: mc["priors"][n] for n in ("A_K", "beta", "C")}, sort_keys=True)


def step_fit(ctx: Context) -> dict:
    from .sed_prototype import tod_method_key

    config = ctx.config
    models_mod, _ = _mc_modules()
    mc = mc_config(config)
    model = selected_model(config)
    names, fn = models_mod.build_models(mc)[model]
    transform = models_mod.prior_transform(names, mc["priors"])
    data = load_ring_tods(config)
    d = step_dir(config, "fit")
    (d / "fits").mkdir(exist_ok=True)
    np.save(d / "ring_ra.npy", np.asarray(data["ra"], dtype=float))
    pkey, tkey = prior_key(mc), tod_method_key(config)
    seed = int(pcfg(config, "fit.seed"))
    min_points = int(mc["data"]["min_points"])
    scfg = mc["sampler"]
    logging.getLogger("ultranest").setLevel(logging.WARNING)
    rows, reused, t0 = [], 0, time.time()
    for u in range(data["ra"].size):
        pts = points(data, u, data["keys"])                       # TRIS included
        if len(pts) < min_points:
            continue
        nu, t, s = (np.array([q[i] for q in pts], dtype=float) for i in (1, 2, 3))
        keys = np.array([q[0] for q in pts])
        p = fit_path(config, u, model)
        r = None
        if p.exists() and not ctx.overwrite:
            z = np.load(p, allow_pickle=False)
            if (str(z["prior_key"]) == pkey and str(z["tod_key"]) == tkey and list(z["keys"]) == list(keys)
                    and np.allclose(z["t"], t) and np.allclose(z["s"], s)):
                r, reused = {k: z[k] for k in z.files}, reused + 1
        if r is None:
            import ultranest

            norm = -np.sum(np.log(s * np.sqrt(2 * np.pi)))

            def loglike(params):
                res = (t - fn(nu, params)) / s
                return -0.5 * np.sum(res * res, axis=-1) + norm

            np.random.seed(seed + u)
            sampler = ultranest.ReactiveNestedSampler(names, loglike, transform, vectorized=True,
                                                      log_dir=scfg.get("log_dir"))
            logging.getLogger("ultranest").setLevel(logging.WARNING)
            res = sampler.run(min_num_live_points=int(scfg["min_num_live_points"]), dlogz=float(scfg["dlogz"]),
                              show_status=False, viz_callback=False)
            r = dict(logz=res["logz"], logzerr=res["logzerr"], samples=res["samples"],
                     maxpoint=np.asarray(res["maximum_likelihood"]["point"]), names=np.array(names))
            np.savez(p, **r, nu=nu, t=t, s=s, keys=keys, prior_key=pkey, tod_key=tkey)
            ctx.say("  unit {:>4}: fitted ({} points) -- {:.0f} s".format(u, len(pts), time.time() - t0))
        best = np.asarray(r["maxpoint"], float)
        chi2 = float(np.sum(((t - fn(nu, best[None, :])[0]) / s) ** 2))
        row = dict(unit=u, ra_deg=float(data["ra"][u]), n_points=len(pts), inputs=";".join(keys),
                   logZ=float(r["logz"]), logZerr=float(r["logzerr"]), chi2_maxlike=chi2, dof=len(pts) - len(names))
        for i, pn in enumerate(names):
            smp = np.asarray(r["samples"], float)[:, i]
            row[pn + " maxlike"], row[pn + " median"], row[pn + " std"] = (float(best[i]), float(np.median(smp)),
                                                                           float(np.std(smp)))
        rows.append(row)
    _write_csv(d / "fit_summary.csv", rows)
    summary = {"model": model, "parameters": names, "units_fitted": len(rows), "reused_from_cache": reused,
               "points_per_unit": sorted({r["n_points"] for r in rows}),
               "data_sets": sorted({k for r in rows for k in r["inputs"].split(";")}),
               "median_parameters_maxlike": {pn: _nanmedian(np.array([r[pn + " maxlike"] for r in rows])) for pn in names},
               "median_parameter_std": {pn: _nanmedian(np.array([r[pn + " std"] for r in rows])) for pn in names},
               "median_chi2_per_dof": _nanmedian(np.array([r["chi2_maxlike"] / r["dof"] if r["dof"] > 0 else np.nan
                                                            for r in rows]))}
    write_manifest(d, "fit", config, summary, ["fits/unit###__{}.npz".format(model.replace(" ", "_")), "ring_ra.npy",
                                                "fit_summary.csv"])
    return summary


DATA_LABELS = {"TRIS_600": "TRIS 600", "TRIS_820": "TRIS 820", "LWA_35MHz.fits": "LWA1 35",
               "LWA1_38.fits": "LWA1 38", "LWA1_40.fits": "LWA1 40", "LWA1_45.fits": "LWA1 45",
               "MAIPU_45mhZ.fits": "Maipu 45", "LWA1_50.fits": "LWA1 50", "LWA1_60.fits": "LWA1 60",
               "LWA1_70.fits": "LWA1 70", "LWA1_74.fits": "LWA1 74", "LWA_80MHz.fits": "LWA1 80",
               "CHIPASS_1400MHz.fits": "CHIPASS 1400", "STOCKERT_VILLA_1400MHZ.fits": "Stockert 1420",
               "RHODES_HARTRAO_2326MHz.fits": "Rhodes 2326"}


def data_label(key: str) -> str:
    """Short plot / table label of a data-set key (``TRIS_600`` or a map file name)."""
    return DATA_LABELS.get(key, Path(key).stem)


def fit_data_sets(config: Optional[Config] = None, model: str = "power law") -> Dict[str, str]:
    """The data sets of step 5's fits (the pipeline's data sets with data on the ring, TRIS included), low frequency
    first: key -> label.  Read from the fit files themselves, so it is exactly what step 5 fitted."""
    config = config if config is not None else load_config()
    nu = {}
    for u in range(int(np.load(_require(step_dir(config, "fit", False) / "ring_ra.npy", "fit")).size)):
        p = fit_path(config, u, model)
        if p.exists():
            z = np.load(p, allow_pickle=False)
            nu.update({str(k): float(f) for k, f in zip(z["keys"], z["nu"])})
    if not nu:
        raise ConfigError("no step 5 fits found: run `python run_pipeline.py --step fit`")
    return {k: data_label(k) for k in sorted(nu, key=lambda q: (nu[q], q))}


def our_model(config: Config):
    """(model name, model function, {unit: fit}) from step 5: our spectral model."""
    models_mod, _ = _mc_modules()
    model = selected_model(config)
    _, fn = models_mod.build_models(mc_config(config))[model]
    fits = {}
    n = load_ring_tods(config)["ra"].size
    for u in range(n):
        p = fit_path(config, u, model)
        if p.exists():
            fits[u] = dict(np.load(p, allow_pickle=False))
    if not fits:
        raise FileNotFoundError("no step 5 fits: run `python run_pipeline.py --step fit` first")
    return model, fn, fits


def predict(config: Config, freqs_mhz, n_units: int):
    """Our model's ring TOD at ``freqs_mhz``: (T, sigma), each (nfreq, n_units), K, RJ, CMB removed.
    T = the model at UltraNest's maximum-likelihood point; sigma = the standard deviation of the model over the
    posterior samples (at 408 MHz, for the power law, exactly the posterior std of A)."""
    _, fn, fits = our_model(config)
    nu = np.atleast_1d(np.asarray(freqs_mhz, dtype=float))
    T = np.full((nu.size, n_units), np.nan)
    S = np.full((nu.size, n_units), np.nan)
    for u, r in fits.items():
        T[:, u] = fn(nu, np.asarray(r["maxpoint"], float)[None, :])[0]
        S[:, u] = np.std(fn(nu, np.asarray(r["samples"], float)), axis=0)
    return T, S


# ---------------------------------------------------------------------------
# step 6: the 408 MHz prediction against Haslam and the 408 MHz sky models
# ---------------------------------------------------------------------------
def york(x, y, sx, sy, inflate=True, n_iter=200, tol=1e-12):
    """y = a x + b with errors on both axes (York et al. 2004, Am. J. Phys. 72, 367; no x-y correlation);
    sigma inflated by sqrt(chi2/dof) when > 1.  The same code as test_results.ipynb Part F."""
    vx, vy = np.asarray(sx, float) ** 2, np.asarray(sy, float) ** 2
    a = np.polyfit(x, y, 1)[0]
    for _ in range(n_iter):
        W = 1.0 / (vy + a * a * vx)
        xb, yb = np.sum(W * x) / W.sum(), np.sum(W * y) / W.sum()
        U, V = x - xb, y - yb
        bet = W * (U * vy + a * V * vx)
        a_new = np.sum(W * bet * V) / np.sum(W * bet * U)
        done = abs(a_new - a) <= tol * abs(a_new)
        a = a_new
        if done:
            break
    W = 1.0 / (vy + a * a * vx)
    xb, yb = np.sum(W * x) / W.sum(), np.sum(W * y) / W.sum()
    b = yb - a * xb
    xa = xb + W * ((x - xb) * vy + a * (y - yb) * vx)
    xab = np.sum(W * xa) / W.sum()
    sa = np.sqrt(1.0 / np.sum(W * (xa - xab) ** 2))
    sb = np.sqrt(1.0 / W.sum() + xab ** 2 * sa ** 2)
    chi2 = float(np.sum(W * (y - a * x - b) ** 2))
    dof = x.size - 2
    f = np.sqrt(max(1.0, chi2 / dof)) if inflate else 1.0
    return dict(a=a, sa=sa * f, b=b, sb=sb * f, chi2_dof=chi2 / dof)


def step_predict_408(ctx: Context) -> dict:
    from .sed_prototype import MIN_COVERAGE, sigma_for

    config = ctx.config
    data = load_ring_tods(config)
    ra, n = data["ra"], data["ra"].size
    model = selected_model(config)
    P, SP = (v[0] for v in predict(config, [408.0], n))
    z = np.load(_require(step_dir(config, "ring_tods", False) / "sky_tods.npz", "ring_tods"), allow_pickle=False)
    names = [str(x) for x in z["names"]]
    target = data["target_file"]
    T_MOD = {nm: np.where(z["cov_{}".format(i)] >= MIN_COVERAGE, z["tod_{}".format(i)], np.nan)
             for i, nm in enumerate(names)}
    S_MOD = {nm: np.asarray(sigma_for(dict(freq_mhz=408.0, file=target), T_MOD[nm]), dtype=float) for nm in names}
    d = step_dir(config, "predict_408")

    pred = {"sample": np.arange(n), "ra_deg": ra, "prediction_K": P, "prediction_sigma_K": SP,
            "haslam_K": data["haslam_T"], "haslam_sigma_K": data["haslam_S"],
            "haslam_minus_prediction_K": data["haslam_T"] - P}
    for nm in names:
        pred["T_K " + nm] = T_MOD[nm]
    _write_csv(d / "prediction_408.csv", [dict(zip(pred, vals)) for vals in zip(*pred.values())])

    # model - prediction, chi2 with the prediction's sigma (test_results Part B's table, for our model)
    rows = []
    for nm in names:
        r = T_MOD[nm] - P
        ok = np.isfinite(r) & np.isfinite(SP) & (SP > 0)
        chi2 = float(np.sum((r[ok] / SP[ok]) ** 2))
        rows.append({"fit": model, "sky model": nm, "N samples": int(ok.sum()), "chi2": chi2,
                     "chi2 / N": chi2 / ok.sum(), "mean of residual [K]": float(np.mean(r[ok])),
                     "median of residual [K]": float(np.median(r[ok])), "std of residual [K]": float(np.std(r[ok])),
                     "rms of residual [K]": float(np.sqrt(np.mean(r[ok] ** 2))),
                     "samples within fit 1σ": int(np.sum(np.abs(r[ok]) <= SP[ok]))})
    rows.sort(key=lambda q: q["chi2"])
    for i, q in enumerate(rows, 1):
        q["rank"] = i
    _write_csv(d / "residuals_408.csv", rows)

    # model = a P + b, York, errors on both axes (test_results Part F), all samples and Cygnus A excluded
    cyg_ra, cyg_half = float(pcfg(config, "predict_408.cygnus_ra_deg")), float(pcfg(config, "predict_408.cygnus_half_width_deg"))
    cyg = np.abs((ra - cyg_ra + 180.0) % 360.0 - 180.0) <= cyg_half
    sp0 = np.nan_to_num(SP, nan=0.0)
    so = []
    for ex in (False, True):
        for nm in names:
            ok = np.isfinite(P) & np.isfinite(T_MOD[nm]) & np.isfinite(S_MOD[nm])
            if ex:
                ok &= ~cyg
            f = york(P[ok], T_MOD[nm][ok], sp0[ok], S_MOD[nm][ok])
            ratio = T_MOD[nm][ok] / P[ok]
            so.append({"samples": "excluding Cygnus A ± {:g}°".format(cyg_half) if ex else "all", "model": nm,
                       "slope a": f["a"], "σ_a": f["sa"], "intercept b [K]": f["b"], "σ_b [K]": f["sb"],
                       "χ²/ν": f["chi2_dof"], "N samples": int(ok.sum()), "mean ratio": float(np.mean(ratio)),
                       "median ratio": float(np.median(ratio))})
    _write_csv(d / "scale_offset_408.csv", so)
    h = data["haslam_T"] - P
    summary = {"model": model, "median_prediction_K": _nanmedian(P), "median_prediction_sigma_K": _nanmedian(SP),
               "prediction_sigma_range_K": [float(np.nanmin(SP)), float(np.nanmax(SP))],
               "haslam_minus_prediction": {"median_K": _nanmedian(h), "rms_K": float(np.sqrt(np.nanmean(h ** 2)))},
               "residuals_408": {q["sky model"]: {"median_K": q["median of residual [K]"], "rms_K": q["rms of residual [K]"],
                                                  "chi2_per_N": q["chi2 / N"], "rank": q["rank"]} for q in rows},
               "scale_offset_haslam": {q["samples"]: {"a": q["slope a"], "sigma_a": q["σ_a"], "b_K": q["intercept b [K]"],
                                                      "sigma_b_K": q["σ_b [K]"], "N": q["N samples"]}
                                       for q in so if q["model"] == "Haslam (uncorrected)"},
               "cygnus_window_samples": int(cyg.sum())}
    write_manifest(d, "predict_408", config, summary, ["prediction_408.csv", "residuals_408.csv", "scale_offset_408.csv"])
    return summary


# ---------------------------------------------------------------------------
# step 7: RHINO
# ---------------------------------------------------------------------------
def step_rhino(ctx: Context) -> dict:
    """test_results.ipynb Part J Test 6 (formerly notebook 19), fed by our model (step 5)."""
    import healpy as hp
    from scipy.interpolate import RegularGridInterpolator
    from limTOD.simulator import pointing_beam_in_eq_sys

    from . import elliptical_beam, frames, sky_maps
    from .sed_prototype import tod_method_key
    from .stage1_tris_maps import import_limtod_tris

    config = ctx.config
    rh = config.section("pipeline")["rhino"]
    sky_cfg = config.section("pipeline")["sky_models"]
    freqs = [float(f) for f in rh["freqs_mhz"]]
    data = load_ring_tods(config)
    RA = np.asarray(data["ra"], dtype=float)
    n = RA.size
    d = step_dir(config, "rhino")
    model = selected_model(config)
    T_TRIS, S_TRIS = predict(config, freqs, n)
    if not (np.all(np.isfinite(T_TRIS)) and np.all(S_TRIS > 0)):
        raise RuntimeError("our model is missing at some ring samples: rerun step 5 (fit)")

    NSIDE = int(config.require("tris.nside"))
    NPIX = hp.nside2npix(NSIDE)
    PRIOR = {k: float(config.require("tris.prior." + k)) for k in ("relative_sigma", "absolute_sigma_k", "floor_sigma_k")}
    geom = import_limtod_tris(config).tris_zenith_geometry(RA)
    ring_vec = hp.ang2vec(np.radians(90.0 - float(geom.latitude_deg)) * np.ones_like(RA),
                          np.radians(np.asarray(geom.lst_deg, dtype=float)))
    spec = elliptical_beam.beam_spec_from_config(config, kernel="tris_cuts")
    kernel = elliptical_beam.kernel_for(spec, 0.0)
    nw = spec.nside_work
    pts = np.array(hp.pix2vec(nw, np.arange(hp.nside2npix(nw)))).T
    th, psi = elliptical_beam.local_offsets(ring_vec, pts)
    A64 = elliptical_beam.kernel_response(kernel, spec, th, psi).T
    A64 = A64 / A64.sum(axis=1, keepdims=True)
    parent = hp.ang2pix(NSIDE, *hp.pix2ang(nw, np.arange(hp.nside2npix(nw))))
    A16 = np.zeros((n, NPIX))
    np.add.at(A16.T, parent, A64.T)
    del pts, th, psi

    sky_models = list(rh["sky_models"])
    native = {q: float(sky_cfg["native_fwhm_deg"][q.split(",")[0]]) for q in sky_models}
    carter, irfan = sky_cfg["carter"], float(sky_cfg["irfan_monopole_k"])
    # the cache key: every setting the sky-model TODs depend on (Irfan: irfan_key)
    sky_key = json.dumps({"freqs": freqs, "models": sky_models, "native": native,
                          "carter": [float(carter["scale"]), float(carter["offset_k"])], "irfan": irfan,
                          "irfan_method": irfan_key(sky_cfg),
                          "tod": tod_method_key(config), "ra": hashlib.md5(RA.tobytes()).hexdigest(), "nside": NSIDE},
                         sort_keys=True)
    (d / "cache").mkdir(exist_ok=True)
    sky_cache = d / "cache" / "rhino_sky_models.npz"

    def tris_tod(sky_g, native_fwhm):
        eq = frames.rotate_full_sky(np.asarray(sky_g, dtype=float), ("G", "C"))
        r = elliptical_beam.smooth_at_positions(eq, np.ones(eq.size, bool), native_fwhm, spec, ring_vec)
        return np.asarray(r["values"], dtype=float), eq

    if sky_cache.exists() and str(np.load(sky_cache)["key"]) == sky_key:
        zc = np.load(sky_cache)
        SKY = {q: (zc["tod_" + q], zc["map_" + q]) for q in sky_models}
        ctx.say("step 7: sky models read from the cache")
    else:
        from pygdsm import GlobalSkyModel, GlobalSkyModel16, LowFrequencySkyModel
        from pygdsm.mckay25 import generate_mckay25_correction, T_CMB

        g08 = GlobalSkyModel(freq_unit="MHz", include_cmb=False)
        g08_408 = np.asarray(g08.generate(408.0), dtype=float)
        irfan_off = -(irfan - float(sky_maps.cmb_rj_k(408.0)))
        g16_cmb = GlobalSkyModel16(freq_unit="MHz", data_unit="TCMB", include_cmb=True)

        def ginan(f):
            off, scl = generate_mckay25_correction(np.array([f / 1e3]))
            return (np.asarray(g16_cmb.generate(f), dtype=float) - off[0]) * scl[0] - T_CMB

        gen = {"GSM2008": g08.generate,
               "GSM2016": GlobalSkyModel16(freq_unit="MHz", data_unit="TRJ", include_cmb=False).generate,
               "LFSM": LowFrequencySkyModel(freq_unit="MHz", include_cmb=False).generate,
               "GSM2008, Carter correction": lambda f: ((float(carter["scale"]) * g08_408 + float(carter["offset_k"]))
                                                        * np.asarray(g08.generate(f), dtype=float) / g08_408),
               "GSM2008, Irfan correction": lambda f: ((g08_408 + irfan_off)
                                                       * irfan_scaling(sky_cfg, hp.npix2nside(g08_408.size), f)),
               "GSM2016, GINAN correction": ginan}
        SKY = {}
        for q in sky_models:
            tods, maps = [], []
            for f in freqs:
                tod, eq = tris_tod(gen[q](f), native[q])
                tods.append(tod)
                maps.append(hp.ud_grade(eq, NSIDE))
            SKY[q] = (np.array(tods), np.array(maps))
        np.savez(sky_cache, key=sky_key, **{"tod_" + q: v[0] for q, v in SKY.items()},
                 **{"map_" + q: v[1] for q, v in SKY.items()})

    def prior_sigma(mu):
        return np.maximum(np.hypot(PRIOR["relative_sigma"] * mu, PRIOR["absolute_sigma_k"]), PRIOR["floor_sigma_k"])

    def deconvolve(dd, s, mu):
        P = prior_sigma(mu) ** 2
        PAt = P[:, None] * A16.T
        Minv = np.linalg.inv(A16 @ PAt + np.diag(s ** 2))
        return mu + PAt @ (Minv @ (dd - A16 @ mu)), PAt, Minv, P

    beam_file, horizon_file = config.resolve_path("pipeline.rhino.beam_file"), config.resolve_path("pipeline.rhino.horizon_file")
    rb = np.load(_require(beam_file, "rhino (the RHINO beam is outside the repo)"))
    za, az = rb["zenith_angle"], rb["azimuth_angle"]
    az_ext = np.append(az, 2 * np.pi)
    horizon = hp.ud_grade(hp.read_map(_require(horizon_file, "rhino (the horizon mask is outside the repo)")), 64)
    horizon = np.where(horizon >= 0.5, 1.0, 0.0)
    lat = float(rh["latitude_deg"])

    def rhino_beam_map(f, nside=64):
        i = int(np.flatnonzero(np.isclose(rb["freqs_mhz"], f))[0])
        g = 10 ** (rb["abs_dir_dbi"][i] / 10)
        g = np.concatenate([g, g[:, :1]], axis=1)
        tt, pp = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
        return RegularGridInterpolator((za, az_ext), g)(np.stack([tt, pp], axis=1))

    def rhino_operator(f):
        alm = hp.map2alm(rhino_beam_map(f))
        return np.array([pointing_beam_in_eq_sys(alm, lst, lat, 0.0, 90.0, 0.0, NSIDE, normalize=True,
                                                 horizontal_mask=horizon) for lst in RA])

    prior_model = str(rh["deconvolution_prior"])
    support = (A16 / A16.max(axis=1, keepdims=True) > float(config.require("tris.beam_response_threshold"))).any(axis=0)
    RH = {k: [] for k in ("ours", "ours_sigma", "closure", "support", "below_horizon")}
    RH.update({q: [] for q in sky_models})
    for j, f in enumerate(freqs):
        R = rhino_operator(f)
        if not np.allclose(R.sum(axis=1), 1.0):
            raise RuntimeError("the RHINO operator is not normalised at {:g} MHz".format(f))
        bm = rhino_beam_map(f)
        RH["below_horizon"].append(bm[horizon == 0].sum() / bm.sum())
        RH["support"].append(R @ support)
        x, PAt, Minv, P = deconvolve(T_TRIS[j], S_TRIS[j], SKY[prior_model][1][j])
        RPAt = R @ PAt
        var = np.einsum("ip,ip->i", R * P, R) - np.einsum("ik,kl,il->i", RPAt, Minv, RPAt)
        RH["ours"].append(R @ x)
        RH["ours_sigma"].append(np.sqrt(var))
        RH["closure"].append(A16 @ x - T_TRIS[j])
        for q in sky_models:
            RH[q].append(R @ SKY[q][1][j])
        ctx.say("  {:g} MHz done".format(f))
    RH = {k: np.array(v) for k, v in RH.items()}

    for beam, ours, sig, mods in (("tris_beam", T_TRIS, S_TRIS, {q: SKY[q][0] for q in sky_models}),
                                  ("rhino_beam", RH["ours"], RH["ours_sigma"], {q: RH[q] for q in sky_models})):
        cols = {"sample": np.arange(n), "lst_deg": RA}
        for j, f in enumerate(freqs):
            cols["T_K_{:g}MHz_this_work".format(f)], cols["sigma_K_{:g}MHz_this_work".format(f)] = ours[j], sig[j]
            for q in sky_models:
                cols["T_K_{:g}MHz_{}".format(f, q)] = mods[q][j]
        _write_csv(d / "rhino_tods_{}.csv".format(beam), [dict(zip(cols, v)) for v in zip(*cols.values())])
    rows, by_beam = [], {}
    for beam, ours, sig, mods in (("TRIS", T_TRIS, S_TRIS, {q: SKY[q][0] for q in sky_models}),
                                  ("RHINO", RH["ours"], RH["ours_sigma"], {q: RH[q] for q in sky_models})):
        for j, f in enumerate(freqs):
            row = {"beam": beam, "frequency (MHz)": f, "median this work (K)": float(np.median(ours[j])),
                   "median σ this work (K)": float(np.median(sig[j]))}
            for q in sky_models:
                dq = mods[q][j] - ours[j]
                row["median {} − this work (K)".format(q)] = float(np.median(dq))
                row["median ({} − this work) / this work".format(q)] = float(np.median(dq / ours[j]))
            rows.append(row)
            by_beam.setdefault(beam, {})["{:g}".format(f)] = {"median_K": row["median this work (K)"],
                                                               "median_sigma_K": row["median σ this work (K)"]}
    _write_csv(d / "rhino_summary.csv", rows)
    checks = [{"frequency (MHz)": f, "max |A x - d| (K)": float(np.abs(RH["closure"][j]).max()),
               "median σ_d (K)": float(np.median(S_TRIS[j])),
               "RHINO beam on TRIS-covered pixels, min over LST": float(RH["support"][j].min()),
               "RHINO beam below the horizon": float(RH["below_horizon"][j])} for j, f in enumerate(freqs)]
    _write_csv(d / "rhino_checks.csv", checks)
    summary = {"model": model, "frequencies_mhz": freqs, "deconvolution_prior": prior_model,
               "this_work": by_beam,
               "median_model_minus_this_work_fraction_rhino_beam": {
                   q: [r["median ({} − this work) / this work".format(q)] for r in rows if r["beam"] == "RHINO"]
                   for q in sky_models},
               "max_closure_K": float(np.abs(RH["closure"]).max())}
    write_manifest(d, "rhino", config, summary,
                   ["rhino_tods_tris_beam.csv", "rhino_tods_rhino_beam.csv", "rhino_summary.csv", "rhino_checks.csv",
                    "cache/rhino_sky_models.npz"])
    return summary


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
TITLES = {"tris_maps": "Step 0 — TRIS maps", "beam_match": "Step 1 — Beam matching",
          "ring_tods": "Step 2 — Ring TODs", "compare": "Step 3 — Model comparison (TRIS excluded)",
          "select": "Step 4 — Model selection", "fit": "Step 5 — Our spectral model (TRIS included)",
          "predict_408": "Step 6 — Haslam at 408 MHz", "rhino": "Step 7 — RHINO prediction"}


def _md_table(header: List[str], rows: List[list]) -> List[str]:
    def fmt(v):
        if isinstance(v, float):
            v = "{:.2f}".format(v) if abs(v) >= 100 else "{:.4g}".format(v)
        return str(v).replace("|", "\\|")
    return (["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
            + ["| " + " | ".join(fmt(v) for v in r) + " |" for r in rows] + [""])


def _read_csv_rows(path: Path) -> List[dict]:
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def _report_compare(config: Config) -> List[str]:
    """Step 3: the evidence (or maximum-likelihood) table of every model, from step 3's summary CSV."""
    d = step_dir(config, "compare", False)
    models_mod, _ = _mc_modules()
    npar = {m: len(p) for m, (p, _) in models_mod.build_models(mc_config(config)).items()}
    if (d / "evidence_summary_tod.csv").exists():
        rows = _read_csv_rows(d / "evidence_summary_tod.csv")
        rows.sort(key=lambda r: -float(r["sum_logZ"]))
        best = float(rows[0]["sum_logZ"])
        n = rows[0]["n_units"]
        out = ["**Evidence of every model** (UltraNest ln Z summed over {} ring samples; {} points in all; TRIS "
               "excluded):".format(n, rows[0].get("n_points", "?")), ""]
        return out + _md_table(
            ["model", "parameters", "Σ ln Z", "± (error)", "Δ ln Z vs best", "samples where it has the highest ln Z"],
            [[r["model"], npar.get(r["model"], "?"), float(r["sum_logZ"]), float(r["sum_logZ_err"]),
              float(r["sum_logZ"]) - best, "{} of {}".format(r["units_won"], n)] for r in rows])
    if (d / "maxlike_summary_tod.csv").exists():
        rows = _read_csv_rows(d / "maxlike_summary_tod.csv")
        return (["**Maximum-likelihood summary of every model** (no evidence was computed; columns as in "
                 "`maxlike_summary_tod.csv`):", ""]
                + _md_table(list(rows[0]), [list(r.values()) for r in rows]))
    return ["No step 3 summary CSV found.", ""]


def _report_select(config: Config) -> List[str]:
    """Step 4: the criteria the model was chosen on, and how each model fared against them."""
    p = step_dir(config, "select", False) / "selection.json"
    if not p.exists():
        return []
    sel = json.loads(p.read_text())
    models_mod, _ = _mc_modules()
    npar = {m: len(q) for m, (q, _) in models_mod.build_models(mc_config(config)).items()}
    tie = pcfg(config, "select.tie_lnz")
    best = max(sel["scores"], key=sel["scores"].get)
    crit = _md_table(["criterion", "setting / result"], [
        ["statistic compared", sel["statistic"]],
        ["1. highest score", "{} (score {:.2f})".format(best, sel["scores"][best])],
        ["2. models within `select.tie_lnz` of the best", "tie_lnz = {}: {}".format(
            tie, ", ".join(sel["candidates_within_tie"]))],
        ["3. tie-break among those", "fewest parameters, then highest score"],
        ["chosen model", "{} (parameters: {})".format(sel["model"], ", ".join(sel["parameters"]))]])
    per = _md_table(["model", "parameters", "score", "score − best", "within tie_lnz", "samples won", "chosen"],
                    [[m, npar.get(m, "?"), float(s), float(sel["score_minus_best"][m]),
                      "yes" if m in sel["candidates_within_tie"] else "no",
                      sel.get("units_won", {}).get(m, "—"), "**yes**" if m == sel["model"] else ""]
                     for m, s in sorted(sel["scores"].items(), key=lambda kv: -kv[1])])
    return ["**Selection criteria** (rule: {}):".format(sel["rule"]), ""] + crit + per


def _report_predict_408(config: Config) -> List[str]:
    """Step 6: the three sky models closest to our 408 MHz prediction (ranked by χ², as residuals_408.csv)."""
    p = step_dir(config, "predict_408", False) / "residuals_408.csv"
    if not p.exists():
        return []
    rows = sorted(_read_csv_rows(p), key=lambda r: int(r["rank"]))
    return (["**Top 3 sky models against our prediction at 408 MHz** (model − prediction along the ring, TRIS beam; "
             "ranked by χ² with the prediction's σ; all {} in `residuals_408.csv`):".format(len(rows)), ""]
            + _md_table(["rank", "sky model", "χ²/N", "median (K)", "mean (K)", "rms (K)", "samples within 1σ", "N"],
                        [[int(r["rank"]), r["sky model"], float(r["chi2 / N"]), float(r["median of residual [K]"]),
                          float(r["mean of residual [K]"]), float(r["rms of residual [K]"]),
                          int(r["samples within fit 1σ"]), int(r["N samples"])] for r in rows[:3]]))


def _report_rhino(config: Config) -> List[str]:
    """Step 7: per beam, the three sky models closest to our prediction over the RHINO frequencies, ranked by χ²
    (model − this work with this work's σ, summed over every sample and frequency), as step 6 ranks at 408 MHz."""
    d = step_dir(config, "rhino", False)
    freqs = [float(f) for f in pcfg(config, "rhino.freqs_mhz")]
    out = []
    for beam, fname in (("TRIS beam (ring at Dec +42.4°)", "rhino_tods_tris_beam.csv"),
                        ("RHINO beam (Jodrell Bank zenith)", "rhino_tods_rhino_beam.csv")):
        if not (d / fname).exists():
            continue
        rows = _read_csv_rows(d / fname)
        col = lambda name: np.array([float(r[name]) if r[name] not in ("", "nan") else np.nan for r in rows])
        models = [k[len("T_K_{:g}MHz_".format(freqs[0])):] for k in rows[0]
                  if k.startswith("T_K_{:g}MHz_".format(freqs[0])) and not k.endswith("_this_work")]
        res = []
        for q in models:
            chi2, n, med, all_d = 0.0, 0, [], []
            for f in freqs:
                ours, sig = col("T_K_{:g}MHz_this_work".format(f)), col("sigma_K_{:g}MHz_this_work".format(f))
                dq = col("T_K_{:g}MHz_{}".format(f, q)) - ours
                ok = np.isfinite(dq) & np.isfinite(sig) & (sig > 0)
                chi2 += float(np.sum((dq[ok] / sig[ok]) ** 2)); n += int(ok.sum())
                med.append(float(np.median(dq[ok]))); all_d.append(dq[ok])
            all_d = np.concatenate(all_d)
            res.append([q, chi2 / n, *med, float(np.sqrt(np.mean(all_d ** 2))), n])
        res.sort(key=lambda r: r[1])
        out += ["**Top 3 sky models against our prediction, {}** (model − this work; ranked by χ²/N with this work's "
                "σ over every sample and frequency; {} models in `{}`):".format(beam, len(res), fname), ""]
        out += _md_table(["rank", "sky model", "χ²/N"] + ["median {:g} MHz (K)".format(f) for f in freqs]
                         + ["rms, all (K)", "N"], [[i] + r for i, r in enumerate(res[:3], 1)])
    return out


REPORT_TABLES = {"compare": _report_compare, "select": _report_select, "predict_408": _report_predict_408,
                 "rhino": _report_rhino}


def write_report(config: Config, path: Optional[Path] = None) -> Path:
    """PIPELINE_REPORT.md from every step's manifest.json (steps never run are listed as such), with plain-text
    tables for steps 3, 4, 6 and 7 read from their products.  ``path`` overrides the output file."""
    root = output_root(config)
    root.mkdir(parents=True, exist_ok=True)
    lines = ["# Pipeline report", "",
             "Written {} by `python run_pipeline.py --step report` from each step's `manifest.json`.".format(
                 _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")), ""]
    for step in STEPS[1:]:
        m = read_manifest(step_dir(config, step, False))
        lines += ["## " + TITLES[step], ""]
        if m is None:
            lines += ["Not run yet.", ""]
            continue
        lines += ["Written {} ({}).".format(m.get("written_utc"), m.get("produced_by")), ""]
        if step in REPORT_TABLES:
            try:
                lines += REPORT_TABLES[step](config)
            except Exception as e:                       # a table must never stop the report
                lines += ["(table not written: {}: {})".format(type(e).__name__, e), ""]
        lines += ["```json", json.dumps(m.get("summary"), indent=2), "```", ""]
    p = Path(path) if path is not None else root / REPORT_NAME
    p.write_text("\n".join(lines))
    return p


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
RUNNERS = {"tris_maps": step_tris_maps, "beam_match": step_beam_match, "ring_tods": step_ring_tods,
           "compare": step_compare, "select": step_select, "fit": step_fit, "predict_408": step_predict_408,
           "rhino": step_rhino}


def run(steps: List[str], config: Optional[Config] = None, verbose: bool = True, overwrite: bool = False,
        report: bool = True) -> Dict[str, Any]:
    config = config if config is not None else load_config()
    ctx = Context(config, verbose=verbose, overwrite=overwrite)
    out = {}
    for step in steps:
        t0 = time.time()
        ctx.say("\n=== step {} — {} ===".format(STEPS.index(step), step))
        out[step] = RUNNERS[step](ctx)
        ctx.say("=== {} done in {:.0f} s ===".format(step, time.time() - t0))
    if report:
        ctx.say("report: {}".format(write_report(config)))
    return out


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _nanmedian(a) -> Optional[float]:
    a = np.asarray(a, dtype=float)
    return float(np.nanmedian(a)) if np.isfinite(a).any() else None


def _write_csv(path: Path, rows: List[dict]) -> None:
    if not rows:
        path.write_text("")
        return
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        for r in rows:
            w.writerow({k: (v.item() if isinstance(v, np.generic) else v) for k, v in r.items()})
