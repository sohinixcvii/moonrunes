"""Does beam averaging change the SHAPE of T vs LST across frequencies?  LWA1 35 MHz and Stockert 1420 MHz.

The same test as ``beam_averaging_blobs/beam_averaging_blobs.py`` (run there on Haslam 408 MHz), with the
map at its own frequency nu_map as the template instead of Haslam:

    T_nu = T_map * (nu / nu_map) ** beta(pixel),   nu = nu_map, 600, 820 MHz,

beta = the same background beta0 = -2.7 plus the same four 4 deg blobs on the TRIS stripe, imported unchanged
from ``beam_averaging_blobs.py`` (not tuned).  Each map is observed along the TRIS ring with the same observer
(``sed_prototype.Inputs._native_tod`` / ``Inputs.observe``: limTOD's TRIS cut beam, zenith pointing, horizon
mask, nside 16), and compared with a control at uniform beta0.  Shape change = 100 (S_blob / S_control - 1),
S = T / mean over the ring.

Two variants:

* ``as given``     -- the power law applied to the map as it is (CMB included), as in the Haslam test;
* ``CMB removed``  -- the CMB's RJ temperature at nu_map subtracted first (and not added back), so only the
                      Galactic emission is scaled.  For Stockert the CMB is most of T away from the plane.

The control stays exactly shape-invariant in both (one constant factor per frequency).  LWA1 35 MHz does not
cover the far south; the observer averages observed pixels only (its ``coverage``), as the TOD route does.

Outputs, ``test_outputs/lwa35_stockert/beam_averaging_blobs/``:

* ``<map stem>/summary.csv``, ``<map stem>/t_vs_lst.csv``, ``<map stem>/t_vs_lst_comparison_<variant>.png``
* ``comparison_with_haslam.csv`` / ``.png`` (Haslam read from ``test_outputs/beam_averaging_blobs/HASLAM_DESTRIPED_ONLY/``,
  not rerun)

Run (``limtod`` conda env; under a minute):

    python beam_tests_lwa35_stockert1420/beam_averaging_blobs_maps.py
"""

from __future__ import annotations

import csv
import os
import sys
from pathlib import Path
from types import SimpleNamespace

os.environ.setdefault("TQDM_DISABLE", "1")
sys.dont_write_bytecode = True

import numpy as np
import healpy as hp
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]                 # scripts/beam_tests_lwa35_stockert1420 -> the repo
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts" / "beam_averaging_blobs"))

import beam_averaging_blobs as bab  # noqa: E402  (read only: blob parameters and beta_map)
from moonrunes import sky_maps  # noqa: E402
from moonrunes.config import load_config  # noqa: E402
from moonrunes.sed_prototype import Inputs  # noqa: E402
from moonrunes.stage1_tris_maps import import_limtod_tris  # noqa: E402

MAPS = ["LWA_35MHz.fits", "STOCKERT_VILLA_1400MHZ.fits"]
VARIANTS = ["as given", "CMB removed"]
VSLUG = {"as given": "as_given", "CMB removed": "cmb_removed"}
OUT_ROOT = REPO / "test_outputs" / "lwa35_stockert" / "beam_averaging_blobs"
HASLAM_SUMMARY = REPO / "test_outputs" / "beam_averaging_blobs" / "HASLAM_DESTRIPED_ONLY" / "summary.csv"
BLOB_COL, CTRL_COL = "#b2432a", "#256abf"


def make_observer(config, tris):
    """The attributes sed_prototype.Inputs builds for the TOD route (as beam_averaging_blobs.py)."""
    archive = config.resolve_path("paths.tris_archive_dir")
    ring = tris.read_tris_ring(archive / "TRIS_absolute_600.txt")
    nside = int(config.require("tris.nside"))
    obs = SimpleNamespace(
        nside=nside,
        geom=tris.tris_zenith_geometry(np.asarray(ring.ra_deg)),
        beam=tris.tris_cut_beam_map(tris.read_tris_beam_cuts(archive / "TRIS_Beam_Profile.txt"),
                                    nside=nside, normalization="peak"),
        horizon=tris.tris_horizon_mask(nside),
        hires=int(config.require("tris.nside_hires")))
    obs.observe = lambda sky: Inputs.observe(obs, sky)
    return obs


def run_map(fname, config, obs, tris_freqs):
    d = sky_maps.describe(config.resolve_path("paths.haslam_map_fits").parent / fname)
    assert d["in_kelvin"] and not d["thermodynamic"]
    nu0 = float(d["freq_mhz"])
    freqs = [nu0] + tris_freqs
    stem = Path(fname).stem
    out_dir = OUT_ROOT / stem
    out_dir.mkdir(parents=True, exist_ok=True)
    lst = np.asarray(obs.geom.lst_deg, dtype=float)
    stripe_dec = float(obs.geom.latitude_deg)
    dec_c = [stripe_dec] * bab.N_BLOBS if bab.BLOB_DEC_DEG is None else list(bab.BLOB_DEC_DEG)
    cmb = sky_maps.cmb_rj_k(nu0)

    bar = "=" * 78
    print(bar)
    print(f"MAP:          {d['name']}  ({d['path']})")
    print(f"              nside {d['nside']}, frame {d['frame']}, native FWHM {d['fwhm_deg']:.3f}°, "
          f"{d['n_unobserved']} unobserved pixels, CMB {d['cmb']} (RJ {cmb:.3f} K at {nu0:g} MHz)")
    print(f"FREQUENCIES:  {', '.join(f'{f:g}' for f in freqs)} MHz  ({nu0:g} = reference = the map's own)")
    print(f"BLOBS:        beta0 = {bab.BETA0}; " + "; ".join(
        f"RA {bab.BLOB_RA_DEG[i]:g}°, FWHM {bab.BLOB_FWHM_DEG[i]:g}°, β {bab.BLOB_BETA[i]:+.2f}" for i in range(bab.N_BLOBS)))
    print(bar)

    beta = bab.beta_map(d["nside"], d["frame"], bab.BLOB_RA_DEG, dec_c, bab.BLOB_FWHM_DEG, bab.BLOB_BETA, bab.BETA0)
    ring_vec = np.array(hp.ang2vec(lst, np.full(lst.size, stripe_dec), lonlat=True)).T
    if d["frame"] != "C":
        ring_vec = hp.Rotator(coord=["C", d["frame"]])(ring_vec)
    ring_th, ring_ph = hp.vec2ang(ring_vec.T)

    rows, curves = [], {}
    for variant in VARIANTS:
        t0 = d["sky"] - (cmb if variant == "CMB removed" else 0.0)
        tod, pencil, cov_min = {}, {}, 1.0
        for case, b in (("blob", beta), ("control", np.full_like(beta, bab.BETA0))):
            for nu in freqs:
                t_nu = t0 * (nu / nu0) ** b
                tod[case, nu], coverage = Inputs._native_tod(obs, dict(sky=t_nu, nside=d["nside"], frame=d["frame"]))
                assert np.all(coverage > 0)
                cov_min = min(cov_min, float(np.min(coverage)))
                pencil[case, nu] = hp.get_interp_val(t_nu, ring_th, ring_ph)
        S = {k: v / np.mean(v) for k, v in tod.items()}
        P = {k: v / np.mean(v) for k, v in pencil.items()}
        ctrl_dev = {nu: float(np.max(np.abs(S["control", nu] / S["control", nu0] - 1))) for nu in freqs}
        print(f"\n[{variant}] min observed beam weight along the ring {cov_min:.4f}; control max |S_nu/S_ref − 1| "
              f"= {max(ctrl_dev.values()):.2e}")
        print(f"{'nu [MHz]':>9}{'ln(nu/nu0)':>11}{'max shape Δ':>13}{'at LST':>9}{'rms':>8}{'norm Δ':>10}{'no-beam':>10}")
        for nu in freqs:
            shape = 100 * (S["blob", nu] / S["control", nu] - 1)
            k = int(np.argmax(np.abs(shape)))
            pshape = 100 * (P["blob", nu] / P["control", nu] - 1)
            lnr = float(np.log(nu / nu0))
            row = dict(map=d["name"], variant=variant, ref_freq_mhz=nu0, freq_mhz=nu, ln_freq_ratio=lnr,
                       max_abs_shape_change_pct=float(np.max(np.abs(shape))), signed_at_max_pct=float(shape[k]),
                       at_lst_deg=float(lst[k]), rms_shape_change_pct=float(np.sqrt(np.mean(shape ** 2))),
                       max_abs_shape_change_per_unit_ln_ratio_pct=(float(np.max(np.abs(shape)) / abs(lnr)) if lnr else 0.0),
                       normalisation_change_pct=float(100 * (np.mean(tod["blob", nu]) / np.mean(tod["control", nu]) - 1)),
                       control_shape_vs_ref_max_abs=ctrl_dev[nu],
                       no_beam_max_abs_shape_change_pct=float(np.max(np.abs(pshape))),
                       min_beam_coverage=cov_min)
            rows.append(row)
            curves[variant, nu] = shape
            print(f"{nu:9.1f}{lnr:+11.3f}{row['max_abs_shape_change_pct']:12.3f}%{row['at_lst_deg']:8.1f}°"
                  f"{row['rms_shape_change_pct']:7.3f}%{row['normalisation_change_pct']:+9.3f}%"
                  f"{row['no_beam_max_abs_shape_change_pct']:9.2f}%")

        with open(out_dir / f"t_vs_lst_{VSLUG[variant]}.csv", "w") as f:
            cols = ["lst_deg"]
            for nu in freqs:
                cols += [f"T_blob_{nu:g}_K", f"T_control_{nu:g}_K", f"S_blob_{nu:g}", f"S_control_{nu:g}",
                         f"shape_change_pct_{nu:g}", f"T_blob_minus_control_{nu:g}_K"]
            f.write(f"# map = {d['file']}; variant = {variant}; T in K (RJ), "
                    f"{'CMB included' if variant == 'as given' else f'CMB ({cmb:.4f} K at {nu0:g} MHz) removed before scaling'}\n")
            f.write(",".join(cols) + "\n")
            for j in range(lst.size):
                vals = [lst[j]]
                for nu in freqs:
                    vals += [tod["blob", nu][j], tod["control", nu][j], S["blob", nu][j], S["control", nu][j],
                             100 * (S["blob", nu][j] / S["control", nu][j] - 1), tod["blob", nu][j] - tod["control", nu][j]]
                f.write(",".join(f"{v:.8g}" for v in vals) + "\n")

        # Plot: as beam_averaging_blobs.py; columns = frequencies, rows = S, T_nu / T_ref, shape change.
        fig, axes = plt.subplots(3, len(freqs), figsize=(4.6 * len(freqs), 9.5), sharex=True)
        for c in range(1, len(freqs)):
            axes[0, c].sharey(axes[0, 0]); axes[2, c].sharey(axes[2, 0])
        for c, nu in enumerate(freqs):
            axes[0, c].plot(lst, S["control", nu], color=CTRL_COL, lw=2, label=f"control (β = {bab.BETA0})")
            axes[0, c].plot(lst, S["blob", nu], color=BLOB_COL, lw=1.4, ls="--", label="blobs")
            axes[0, c].set_title(f"{nu:g} MHz", loc="left", fontsize=10)
            axes[1, c].plot(lst, tod["control", nu] / tod["control", nu0], color=CTRL_COL, lw=2)
            axes[1, c].plot(lst, tod["blob", nu] / tod["blob", nu0], color=BLOB_COL, lw=1.4, ls="--")
            axes[2, c].plot(lst, curves[variant, nu], color=BLOB_COL, lw=1.8)
            axes[2, c].axhline(0, color="0.4", lw=0.8)
            for a in axes[:, c]:
                for i in range(bab.N_BLOBS):
                    a.axvline(bab.BLOB_RA_DEG[i], color="0.82", lw=0.9, zorder=0)
                a.grid(axis="y", alpha=0.2)
            for i in range(bab.N_BLOBS):
                axes[2, c].annotate(f"β {bab.BLOB_BETA[i]:+.1f}", (bab.BLOB_RA_DEG[i], 1), xycoords=("data", "axes fraction"),
                                    xytext=(2, -10), textcoords="offset points", fontsize=7, color="0.35")
            axes[2, c].set_xlabel("LST [deg]")
            axes[2, c].set_xlim(0, 360); axes[2, c].set_xticks(np.arange(0, 361, 60))
        axes[0, 0].set_ylabel("T / mean(T)")
        axes[1, 0].set_ylabel(f"T_ν / T_{nu0:g}")
        axes[2, 0].set_ylabel("shape change: S_blob / S_control − 1 [%]")
        axes[0, 0].legend(fontsize=8, frameon=False)
        fig.suptitle(f"{d['name']} ({d['file']}, {variant}) extrapolated with the blob β map, observed with the TRIS "
                     "beam along the ring (grey lines: blob LSTs)", x=0.01, ha="left", fontsize=10)
        plt.tight_layout()
        fig.savefig(out_dir / f"t_vs_lst_comparison_{VSLUG[variant]}.png", dpi=150, facecolor="white")
        plt.close(fig)

    with open(out_dir / "summary.csv", "w") as f:
        f.write(f"# map = {d['file']}; reference = the map's own {nu0:g} MHz; shape change = 100*(S_blob/S_control - 1), "
                "S = T/mean_LST(T); beam = limTOD TRIS cut beam (sed_prototype.Inputs._native_tod)\n")
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.8g}" if isinstance(v, float) else v) for k, v in r.items()})
    print(f"\nwrote {out_dir}")
    return d, rows, curves, lst


def comparison(results, tris_freqs):
    """Haslam (earlier run, read from its CSV) next to the two maps: per-frequency numbers and the curves."""
    hrows = list(csv.DictReader(line for line in open(HASLAM_SUMMARY) if not line.startswith("#")))
    hcurve_rows = list(csv.DictReader(open(HASLAM_SUMMARY.parent / "t_vs_lst.csv")))
    hl = np.array([float(r["lst_deg"]) for r in hcurve_rows])
    out = []
    for r in hrows:
        nu = float(r["freq_mhz"])
        if nu in tris_freqs:
            lnr = float(np.log(nu / 408.0))
            out.append(dict(map="Haslam 408 MHz (destriped only)", variant="as given", ref=408.0, nu=nu, lnr=lnr,
                            mx=float(r["max_abs_shape_change_pct"]), at=float(r["at_lst_deg"]),
                            rms=float(r["rms_shape_change_pct"]), per=float(r["max_abs_shape_change_pct"]) / abs(lnr),
                            norm=float(r["normalisation_change_pct"]), nob=float(r["no_beam_max_abs_shape_change_pct"])))
    for d, rows, _, _ in results:
        for r in rows:
            if r["freq_mhz"] in tris_freqs:
                out.append(dict(map=r["map"], variant=r["variant"], ref=r["ref_freq_mhz"], nu=r["freq_mhz"],
                                lnr=r["ln_freq_ratio"], mx=r["max_abs_shape_change_pct"], at=r["at_lst_deg"],
                                rms=r["rms_shape_change_pct"], per=r["max_abs_shape_change_per_unit_ln_ratio_pct"],
                                norm=r["normalisation_change_pct"], nob=r["no_beam_max_abs_shape_change_pct"]))
    with open(OUT_ROOT / "comparison_with_haslam.csv", "w", newline="") as f:
        f.write("# Haslam rows read from test_outputs/beam_averaging_blobs/HASLAM_DESTRIPED_ONLY/summary.csv (not rerun). "
                "shape change = 100*(S_blob/S_control - 1); per unit ln ratio = max |shape change| / |ln(nu/nu_ref)|\n")
        w = csv.writer(f)
        w.writerow(["map", "variant", "reference frequency [MHz]", "frequency [MHz]", "ln(nu / nu_ref)",
                    "max |shape change| [%]", "at LST [deg]", "rms shape change [%]",
                    "max |shape change| per unit |ln(nu/nu_ref)| [%]", "normalisation change [%]",
                    "max |shape change| without the beam [%]"])
        for o in out:
            w.writerow([o["map"], o["variant"], f"{o['ref']:g}", f"{o['nu']:g}"] +
                       [f"{o[k]:.6g}" for k in ("lnr", "mx", "at", "rms", "per", "norm", "nob")])

    fig, axes = plt.subplots(1, len(tris_freqs), figsize=(7 * len(tris_freqs), 4.6), sharey=True)
    for ax, nu in zip(np.atleast_1d(axes), tris_freqs):
        ax.plot(hl, [float(r[f"shape_change_pct_{nu:g}"]) for r in hcurve_rows], color="0.45", lw=1.8,
                label="Haslam 408 MHz (as given; earlier run)")
        for (d, rows, curves, lst), col in zip(results, ("#7b3294", "#008837")):
            for variant, ls in zip(VARIANTS, ("-", "--")):
                ax.plot(lst, curves[variant, nu], color=col, lw=1.6, ls=ls, label=f"{d['name']} ({variant})")
        for i in range(bab.N_BLOBS):
            ax.axvline(bab.BLOB_RA_DEG[i], color="0.85", lw=0.9, zorder=0)
            ax.annotate(f"β {bab.BLOB_BETA[i]:+.1f}", (bab.BLOB_RA_DEG[i], 1), xycoords=("data", "axes fraction"),
                        xytext=(2, -10), textcoords="offset points", fontsize=7, color="0.35")
        ax.axhline(0, color="0.4", lw=0.8)
        ax.set_xlim(0, 360); ax.set_xticks(np.arange(0, 361, 60)); ax.set_xlabel("LST [deg]")
        ax.set_title(f"{nu:g} MHz (each map relative to its own frequency)", loc="left", fontsize=10)
        ax.grid(axis="y", alpha=0.2)
    np.atleast_1d(axes)[0].set_ylabel("shape change: S_blob / S_control − 1 [%]")
    np.atleast_1d(axes)[0].legend(fontsize=8, frameon=False, loc="lower left")
    fig.suptitle("Shape change from the same β blobs, three template maps (grey lines: blob LSTs)",
                 x=0.01, ha="left", fontsize=11)
    plt.tight_layout()
    fig.savefig(OUT_ROOT / "comparison_with_haslam.png", dpi=150, facecolor="white")
    print(f"wrote {OUT_ROOT / 'comparison_with_haslam.csv'} and .png")


def main():
    config = load_config(REPO / "configs" / "run_pipeline.yaml")
    tris = import_limtod_tris(config)
    obs = make_observer(config, tris)
    tris_freqs = [float(f) for f in config.require("tris.frequencies_mhz")]
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    results = [run_map(f, config, obs, tris_freqs) for f in MAPS]
    comparison(results, tris_freqs)
    # the beta map is the same as in test_outputs/beam_averaging_blobs/<run>/beta_map.png


if __name__ == "__main__":
    main()
