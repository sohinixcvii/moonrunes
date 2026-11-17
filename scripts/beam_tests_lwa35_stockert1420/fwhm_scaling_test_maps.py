"""Two-step vs one-step beam matching with four FWHM scalings, on LWA1 35 MHz and Stockert 1420 MHz.

The same test as ``scripts/fwhm_scaling_test.py`` (which runs it on Haslam 408 MHz; the test itself is
``notebooks/test_results.ipynb`` Part K2, "Two-step against one-step beam matching", ``PIPELINE.md`` App. E):

1. step 1: the map is smoothed with the made-up "complicated" survey beam (5 x 3 deg oval main lobe at
   PA 30 deg, a -15 dB ring at 7 deg, a -25 dB 20 deg halo);
2. two-step: that map is matched to the TRIS beam with ``match_to_tris_beam``, telling the pipeline one
   FWHM for the survey beam (major, minor, arithmetic or geometric mean of the main lobe);
3. one-step: the map is matched to the TRIS beam directly, telling the pipeline the map's own beam;
4. metric: 100 * (two-step / one-step - 1) [%] and two-step - one-step [K], over the sky and in the TRIS stripe.

Imported unchanged from ``scripts/fwhm_scaling_test.py``: the complicated beam, its measured widths, the
statistics, the kernel set-up and the CSV header texts.  Two things differ, because these maps are not
full-sky Haslam:

* **Masks.**  LWA1 35 MHz covers 82% of the sky and Stockert misses 47 pixels.  The map is degraded to
  nside 64 with the pipeline's ``all_children`` rule, step 1 averages only observed pixels and keeps a pixel
  only where at least ``beam_matching.weight_floor`` of the survey beam's weight was observed, and the
  pipeline is given that observed mask.  Statistics use the pixels kept in both the two-step and one-step maps.
* **Two variants of the FWHM given to the pipeline.**  ``survey only``, as in the Haslam test (Haslam's own
  0.93 deg beam was ignored there), and ``survey + native``: sqrt(width^2 + native^2), because the step-1
  map has really been smoothed by the map's own beam as well.  For LWA1 35 (native 4.65 deg) the two differ
  a lot; for Stockert (0.59 deg) hardly.

Outputs, ``test_outputs/lwa35_stockert/fwhm_scaling/<map stem>/``:

* ``fwhm_scaling_table.csv``            one row per variant x kernel x scaling
* ``fwhm_scaling_comparison.png``       median / std / max |delta| [%] against the scaling
* ``fwhm_scaling_residual_maps_<variant>.png`` and ``fwhm_scaling_residual_maps.npz``

and, after both maps, ``fwhm_scaling/comparison_with_haslam.csv/.png`` (Haslam from
``test_outputs/fwhm_scaling_test/fwhm_scaling_table.csv``, not rerun).

Run (``limtod`` conda env; about 15 minutes per map):

    python beam_tests_lwa35_stockert1420/fwhm_scaling_test_maps.py [--maps LWA_35MHz.fits STOCKERT_VILLA_1400MHZ.fits]
    python beam_tests_lwa35_stockert1420/fwhm_scaling_test_maps.py --replot-maps   # residual maps only, seconds
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import time
import warnings
from pathlib import Path

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
sys.path.insert(0, str(REPO / "scripts"))

import fwhm_scaling_test as fst  # noqa: E402  (read only: functions and constants)
from moonrunes import elliptical_beam as eb, frames, sky_maps  # noqa: E402
from moonrunes.config import load_config  # noqa: E402
from moonrunes.stage1_tris_maps import load_stage1  # noqa: E402
from moonrunes.stage2_beam_matching import match_to_tris_beam, regrid_to_nside  # noqa: E402

warnings.filterwarnings("ignore", category=UserWarning)

MAPS = ["LWA_35MHz.fits", "STOCKERT_VILLA_1400MHZ.fits"]
VARIANTS = ["survey only", "survey + native"]
VSLUG = {"survey only": "survey_only", "survey + native": "survey_plus_native"}
NS64, NSIDE_OUT = fst.NS64, fst.NSIDE_OUT
OUT_ROOT = REPO / "test_outputs" / "lwa35_stockert" / "fwhm_scaling"
HASLAM_CSV = REPO / "test_outputs" / "fwhm_scaling_test" / "fwhm_scaling_table.csv"

CSV_HEADERS = {"map": "map", "freq_mhz": "frequency [MHz]", "native_fwhm_deg": "map's own beam FWHM [deg]",
               "variant": "FWHM given to the pipeline: survey beam only, or survey beam + map's own beam in quadrature",
               **fst.CSV_HEADERS, "n_pixels": "pixels compared (nside 16, kept in both maps)",
               "n_stripe_pixels": "TRIS-stripe pixels compared"}


def to_nside64(d):
    """The map in equatorial on the nside 64 grid, with its observed mask (pipeline rotation and regrid)."""
    sky = np.asarray(d["sky"], dtype=float)
    observed = np.isfinite(sky)
    if d["frame"] == "C":
        rotated, obs_rot = np.where(observed, sky, hp.UNSEEN), observed
    elif observed.all():
        rotated, obs_rot = frames.rotate_full_sky(sky, (d["frame"], "C")), observed
    else:
        rotated, obs_rot, _ = frames.rotate_masked(sky, observed, (d["frame"], "C"))
    if not hp.isnpixok(rotated.size) or hp.npix2nside(rotated.size) < NS64:
        raise ValueError("map coarser than nside 64")
    v, kept = regrid_to_nside(np.where(obs_rot, rotated, hp.UNSEEN), obs_rot, NS64, "all_children")
    return np.where(kept, v, 0.0), np.asarray(kept, bool)


def convolve_masked(sky, observed, beam_fn, pts, floor, chunk=256):
    """Step 1 on observed pixels only: out(p) = sum_q B T obs / sum_q B obs; kept where that observed
    beam weight is >= floor of the full beam weight (fst.convolve_native, with a mask)."""
    out, frac = np.empty(sky.size), np.empty(sky.size)
    t, o = np.where(observed, sky, 0.0), observed.astype(float)
    for s in range(0, len(pts), chunk):
        th, psi = eb.local_offsets(pts[s:s + chunk], pts)
        w = beam_fn(th, psi)
        wo = o @ w
        with np.errstate(invalid="ignore", divide="ignore"):
            out[s:s + chunk] = (t @ w) / wo
        frac[s:s + chunk] = wo / w.sum(axis=0)
    keep = observed & (frac >= floor)
    return np.where(keep, out, 0.0), keep


def run_map(fname, config, specs, configs, band, scalings, labels):
    d = sky_maps.describe(config.resolve_path("paths.haslam_map_fits").parent / fname)
    native = float(d["fwhm_deg"])
    stem = Path(fname).stem
    out_dir = OUT_ROOT / stem
    out_dir.mkdir(parents=True, exist_ok=True)
    target_circ = float(config.require("beam_matching.target_fwhm_deg"))
    fwhm_e = float(config.require("beam_matching.elliptical.fwhm_e_deg"))
    fwhm_h = float(config.require("beam_matching.elliptical.fwhm_h_deg"))
    floor = float(config.require("beam_matching.weight_floor"))

    bar = "=" * 78
    print(bar)
    print(f"MAP:        {d['name']}  ({d['path']})")
    print(f"            nside {d['nside']}, frame {d['frame']}, {d['unit']}, native FWHM {native:.3f} deg, "
          f"{d['n_unobserved']} unobserved pixels")
    print(f"FREQUENCY:  {d['freq_mhz']:g} MHz")
    print(bar)

    sky64, obs64 = to_nside64(d)
    pts64 = np.array(hp.pix2vec(NS64, np.arange(sky64.size))).T
    print(f"nside 64: {obs64.mean():.1%} of pixels observed")

    def pipeline(sky, observed, given, label):
        spec = None if fst.KERNEL_KEY[label] == "circular_gaussian" else specs[label]
        v, kept, info = match_to_tris_beam(sky, observed, given, configs[label], nside_out=NSIDE_OUT,
                                           rule="all_children", spec=spec, full_sky=bool(observed.all()))
        return np.where(kept, v, np.nan), info

    t0 = time.time()
    pre, pre_obs = convolve_masked(sky64, obs64, fst.complicated_beam, pts64, floor)
    print(f"step 1 (complicated survey beam) done in {time.time() - t0:.0f} s; "
          f"{pre_obs.mean():.1%} of pixels kept")
    direct = {label: pipeline(sky64, obs64, native, label)[0] for label in labels}
    print(f"one-step reference done in {time.time() - t0:.0f} s")

    rows, resid = [], {}
    for variant in VARIANTS:
        for label in labels:
            for name in fst.SCALING_ORDER:
                given = scalings[name] if variant == "survey only" else float(np.hypot(scalings[name], native))
                two, info = pipeline(pre, pre_obs, given, label)
                r = two - direct[label]
                both = np.isfinite(r)
                resid[variant, label, name] = r
                pct = fst.stats(100 * (two / direct[label] - 1), band)
                kel = {"mean": float(np.nanmean(r)), **fst.stats(r, band)}
                axis = fst.fwhm_axis_text(name, scalings)
                if variant == "survey + native":
                    axis += f", plus the map's own {native:.2f} deg beam in quadrature = {given:.3f} deg"
                rows.append(dict(map=d["name"], freq_mhz=d["freq_mhz"], native_fwhm_deg=native, variant=variant,
                                 kernel=label, scaling=name, survey_beam=fst.survey_beam_text(),
                                 target_beam=fst.target_beam_text(fst.KERNEL_KEY[label], fwhm_e, fwhm_h, target_circ),
                                 fwhm_axis=axis, fwhm_given_deg=given,
                                 compared_against=(f"{d['name']} smoothed straight to the same target beam "
                                                   f"(pipeline told the map's own {native:.3f} deg beam)"),
                                 kernel_applied=fst.kernel_applied_text(info),
                                 **{f"pct_{k}": v for k, v in pct.items()},
                                 **{f"K_{k}": v for k, v in kel.items()},
                                 n_pixels=int(both.sum()), n_stripe_pixels=int((both & band).sum())))
            print(f"{variant}, {label} done in {time.time() - t0:.0f} s")

    cols = list(CSV_HEADERS)
    with open(out_dir / "fwhm_scaling_table.csv", "w", newline="") as f:
        f.write(f"# Each row: {d['name']} ({d['file']}) is smoothed with the survey beam, then the pipeline smooths "
                f"it to the target beam, told the survey beam has the width given (= two-step). That is compared with "
                f"the map smoothed straight to the target beam (= one-step). residual = two-step - one-step [K]; "
                f"fractional difference = two-step / one-step - 1 [%]; std = standard deviation; nside {NSIDE_OUT}; "
                f"pixels kept in both maps only. Test from notebooks/03, as scripts/fwhm_scaling_test.py.\n")
        w = csv.writer(f)
        w.writerow([CSV_HEADERS[c] for c in cols])
        for r in rows:
            w.writerow([f"{r[c]:.6g}" if isinstance(r[c], float) else r[c] for c in cols])
    np.savez(out_dir / "fwhm_scaling_residual_maps.npz",
             **{f"{VSLUG[v]}__{fst.KERNEL_KEY[l]}__{n}": m for (v, l, n), m in resid.items()})

    print(f"\n{'variant':<17}{'kernel':<26}{'scaling':<12}{'FWHM':>8}{'median':>9}{'std':>8}{'max|Δ|':>9}"
          f"{'stripe std':>12}{'mean [K]':>11}{'std [K]':>10}")
    for r in rows:
        print(f"{r['variant']:<17}{r['kernel']:<26}{r['scaling']:<12}{r['fwhm_given_deg']:7.3f}°"
              f"{r['pct_median']:+8.2f}%{r['pct_rms']:7.2f}%{r['pct_max_abs']:8.2f}%{r['pct_stripe_rms']:11.2f}%"
              f"{r['K_mean']:+11.4f}{r['K_rms']:10.4f}")

    # Plot: as the Haslam figure; solid = survey only, dashed = survey + native.
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4))
    panels = [("pct_median", "median"), ("pct_rms", "std over the sky"), ("pct_stripe_rms", "std in the TRIS stripe")]
    x = np.arange(len(fst.SCALING_ORDER))
    for ax, (key, title) in zip(axes, panels):
        for label in labels:
            for variant, ls in zip(VARIANTS, ("-", "--")):
                rs = {r["scaling"]: r for r in rows if r["kernel"] == label and r["variant"] == variant}
                ax.plot(x, [rs[n][key] for n in fst.SCALING_ORDER], ls, marker="o", color=fst.KCOL[label],
                        lw=1.6, ms=5, label=f"{label}, {variant}")
        if key == "pct_median":
            ax.axhline(0, color="0.4", lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{n}\n{scalings[n]:.2f}°" for n in fst.SCALING_ORDER], fontsize=8)
        ax.set_xlabel("survey-beam width (dashed: ⊕ the map's own beam)")
        ax.set_ylabel("two-step / one-step − 1 [%]")
        ax.set_title(title, loc="left", fontsize=10)
        ax.grid(axis="y", alpha=0.2)
    axes[0].legend(fontsize=7, frameon=False)
    fig.suptitle(f"{d['name']} ({native:.2f}° beam), complicated survey beam: FWHM scaling (nside {NS64} → {NSIDE_OUT}); "
                 f"solid: survey width only, dashed: survey ⊕ own beam", fontsize=10, x=0.01, ha="left")
    plt.tight_layout()
    fig.savefig(out_dir / "fwhm_scaling_comparison.png", dpi=150, facecolor="white")
    plt.close(fig)

    plot_maps(resid, scalings, native, labels, d["name"], out_dir)
    print(f"\nwrote {out_dir}")
    return rows


def plot_maps(resid, scalings, native, labels, name, out_dir):
    # One residual-map figure per variant (fst.plot_residual_maps), titled with the width actually given.
    for variant in VARIANTS:
        sub = {(l, n): resid[variant, l, n] for l in labels for n in fst.SCALING_ORDER}
        given = scalings if variant == "survey only" else {n: float(np.hypot(w, native)) for n, w in scalings.items()}
        path = out_dir / f"fwhm_scaling_residual_maps_{VSLUG[variant]}.png"
        fst.plot_residual_maps(sub, given, labels, path)
        fig = plt.gcf()                  # fst's title names Haslam; replace it with this map's
        fig._suptitle.set_text(f"Residual two-step − one-step [K], {name}, FWHM given: {variant}"
                               + ("" if variant == "survey only" else f" (survey width ⊕ {native:.2f}°)")
                               + f", nside {NSIDE_OUT} (grey: not compared; colour scale = 99th percentile of |residual|)")
        fig.savefig(path, dpi=150, facecolor="white")
        plt.close("all")


def replot_maps(fname, config, scalings):
    # Redraw the residual-map figures of one map from its saved .npz (no pipeline rerun).
    d = sky_maps.describe(config.resolve_path("paths.haslam_map_fits").parent / fname)
    out_dir = OUT_ROOT / Path(fname).stem
    z = np.load(out_dir / "fwhm_scaling_residual_maps.npz")
    key_label = {v: k for k, v in fst.KERNEL_KEY.items()}
    vlab = {v: k for k, v in VSLUG.items()}
    resid = {(vlab[k.split("__")[0]], key_label[k.split("__")[1]], k.split("__")[2]): z[k] for k in z.files}
    labels = [l for l in fst.KERNEL_KEY if any(k[1] == l for k in resid)]
    plot_maps(resid, scalings, float(d["fwhm_deg"]), labels, d["name"], out_dir)
    print(f"redrew the residual maps in {out_dir}")


def read_haslam():
    with open(HASLAM_CSV) as f:
        rows = list(csv.DictReader(line for line in f if not line.startswith("#")))
    H = fst.CSV_HEADERS
    return [dict(map="Haslam 408 MHz", variant="survey only", kernel=r[H["kernel"]], scaling=r[H["scaling"]],
                 pct_median=float(r[H["pct_median"]]), pct_rms=float(r[H["pct_rms"]]),
                 pct_stripe_rms=float(r[H["pct_stripe_rms"]]), pct_max_abs=float(r[H["pct_max_abs"]])) for r in rows]


def comparison(all_rows, labels):
    rows = read_haslam() + [r for rr in all_rows for r in rr]
    keys = ["map", "variant", "kernel", "scaling", "pct_median", "pct_rms", "pct_stripe_rms", "pct_max_abs"]
    with open(OUT_ROOT / "comparison_with_haslam.csv", "w", newline="") as f:
        f.write("# Haslam rows read from test_outputs/fwhm_scaling_test/fwhm_scaling_table.csv (not rerun; survey width "
                "only, Haslam's own beam ignored). pct = two-step / one-step - 1 [%].\n")
        w = csv.writer(f)
        w.writerow(["map", "FWHM given to the pipeline", "pipeline kernel", "FWHM scaling",
                    "median of fractional difference [%]", "std of fractional difference [%]",
                    "std of fractional difference in TRIS stripe [%]", "max |fractional difference| [%]"])
        for r in rows:
            w.writerow([f"{r[k]:.6g}" if isinstance(r[k], float) else r[k] for k in keys])

    maps = ["Haslam 408 MHz"] + list(dict.fromkeys(r["map"] for rr in all_rows for r in rr))
    series = [(m, v) for m in maps for v in VARIANTS if any(r["map"] == m and r["variant"] == v for r in rows)]
    styles = {}
    mk = ["o", "s", "^", "D", "v"]
    for i, s in enumerate(series):
        styles[s] = dict(marker=mk[maps.index(s[0]) % len(mk)], ls="-" if s[1] == "survey only" else "--")
    fig, axes = plt.subplots(2, len(labels), figsize=(5.2 * len(labels), 8), sharex=True)
    axes = np.atleast_2d(axes).reshape(2, len(labels))
    x = np.arange(len(fst.SCALING_ORDER))
    for c, label in enumerate(labels):
        for rI, key in enumerate(("pct_median", "pct_stripe_rms")):
            ax = axes[rI, c]
            for s in series:
                rs = {r["scaling"]: r for r in rows if r["map"] == s[0] and r["variant"] == s[1] and r["kernel"] == label}
                if len(rs) == len(fst.SCALING_ORDER):
                    ax.plot(x, [rs[n][key] for n in fst.SCALING_ORDER], color=fst.KCOL[label], lw=1.4, ms=5,
                            alpha=1.0 if s[0] != "Haslam 408 MHz" else 0.6, label=f"{s[0]}, {s[1]}", **styles[s])
            if key == "pct_median":
                ax.axhline(0, color="0.4", lw=0.8)
            ax.grid(axis="y", alpha=0.2)
            ax.set_title(f"{label}: {'median' if key == 'pct_median' else 'std in the TRIS stripe'}",
                         loc="left", fontsize=10)
            ax.set_ylabel("two-step / one-step − 1 [%]")
        axes[1, c].set_xticks(x)
        axes[1, c].set_xticklabels(fst.SCALING_ORDER)
        axes[1, c].set_xlabel("survey-beam width told to the pipeline")
    axes[0, 0].legend(fontsize=7.5, frameon=False)
    fig.suptitle("Two-step vs one-step beam matching on three maps (Haslam: earlier run, survey width only)",
                 x=0.01, ha="left", fontsize=11)
    plt.tight_layout()
    fig.savefig(OUT_ROOT / "comparison_with_haslam.png", dpi=150, facecolor="white")
    print(f"wrote {OUT_ROOT / 'comparison_with_haslam.csv'} and .png")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--maps", nargs="+", default=MAPS)
    ap.add_argument("--replot-maps", action="store_true", help="only redraw the residual-map figures from the .npz")
    args = ap.parse_args()
    config = load_config(REPO / "configs" / "run_pipeline.yaml")

    theta_line = np.arange(0.0, 40.0, 0.01)
    major = fst.fwhm_from_profile(theta_line, fst.complicated_beam(theta_line, np.full(theta_line.size, fst.CB["pa"])))
    minor = fst.fwhm_from_profile(theta_line, fst.complicated_beam(theta_line, np.full(theta_line.size, fst.CB["pa"] + 90)))
    scalings = {"major": major, "minor": minor, "arithmetic": 0.5 * (major + minor), "geometric": float(np.sqrt(major * minor))}
    print("survey-beam widths (measured main lobe): " + ", ".join(f"{n} {scalings[n]:.4f}°" for n in fst.SCALING_ORDER))

    if args.replot_maps:
        for f in args.maps:
            replot_maps(f, config, scalings)
        return
    spec_cuts = eb.beam_spec_from_config(config, kernel="tris_cuts")
    specs = {"circular Gaussian (old)": eb.with_kernel(spec_cuts, "circular_gaussian"),
             "elliptical Gaussian": eb.with_kernel(spec_cuts, "elliptical_gaussian"),
             "TRIS cuts (measured)": spec_cuts}
    labels = list(specs)

    def config_with_kernel(kernel):
        c = load_config(config.path)
        c._data["beam_matching"]["kernel"] = kernel
        return c

    configs = {label: config_with_kernel(fst.KERNEL_KEY[label]) for label in labels}
    band = load_stage1(config)["band_mask"].astype(bool)
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    all_rows = [run_map(f, config, specs, configs, band, scalings, labels) for f in args.maps]
    comparison(all_rows, labels)


if __name__ == "__main__":
    main()
