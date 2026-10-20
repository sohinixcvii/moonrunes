"""Two-step vs one-step beam matching, rerun with four FWHM scalings.

The existing test is ``notebooks/03_diagnostics.ipynb``, section "Two-step against
one-step beam matching" (``PIPELINE.md`` Appendix E).  It convolves Haslam with a
complicated, non-Gaussian survey beam (elliptical 5 x 3 deg main lobe at PA 30 deg,
a -15 dB sidelobe ring at 7 deg, a -25 dB 20 deg error beam), then matches it to the
TRIS target with ``match_to_tris_beam``, telling the pipeline one FWHM for the survey
beam.  That two-step map is compared with Haslam matched directly (one step).

The notebook tells the pipeline the main lobe's *major-axis* FWHM.  This script
reruns the same test with the FWHM given to the pipeline set to:

* major      -- the measured major-axis FWHM (what the notebook uses)
* minor      -- the measured minor-axis FWHM
* arithmetic -- (major + minor) / 2
* geometric  -- sqrt(major * minor)

and reports the notebook's discrepancy metric, the fractional difference
100 * (two-step / one-step - 1) [%]: its median, std (``np.nanstd``, which the notebook
calls "rms"), max |delta| over the sky, and std in the TRIS stripe.  The same
statistics of the residual two-step - one-step [K] are added alongside, plus its mean.
Nothing in ``src/`` is modified; the beam and the step-1 convolution are copied from
the notebook.

Outputs (``--out``, default ``outputs/fwhm_scaling_test/``):

* ``fwhm_scaling_table.csv``         -- one row per kernel x scaling; headers name each statistic
* ``fwhm_scaling_comparison.png``    -- median / std / max |delta| [%] against the scaling
* ``fwhm_scaling_residual_maps.png`` -- the residual map [K] of every test, one shared scale
* ``fwhm_scaling_residual_maps.npz`` -- those maps (nside 16), keyed ``<kernel>__<scaling>``

Run (``limtod`` conda env; ~7.5 min, most of it the elliptical and cut-beam kernels):

    python scripts/fwhm_scaling_test.py [--out DIR]
    python scripts/fwhm_scaling_test.py --replot-maps   # redraw the residual maps only, seconds
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

import numpy as np
import healpy as hp
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from moonrunes import elliptical_beam as eb, frames  # noqa: E402
from moonrunes.config import load_config  # noqa: E402
from moonrunes.stage1_tris_maps import load_stage1  # noqa: E402
from moonrunes.stage2_beam_matching import match_to_tris_beam  # noqa: E402

warnings.filterwarnings("ignore", category=UserWarning)

K_ = 2.0 * np.sqrt(2.0 * np.log(2.0))
NS64 = 64
NSIDE_OUT = 16
HASLAM_FREQ_MHZ = 408.0

# The complicated survey beam, exactly as in notebook 03.
CB = dict(major=5.0, minor=3.0, pa=30.0, ring_r=7.0, ring_w=1.5, ring_db=-15.0,
          shoulder_fwhm=20.0, shoulder_db=-25.0)

KERNEL_KEY = {"circular Gaussian (old)": "circular_gaussian",
              "elliptical Gaussian": "elliptical_gaussian",
              "TRIS cuts (measured)": "tris_cuts"}
KCOL = {"circular Gaussian (old)": "#6b7280", "elliptical Gaussian": "#eb6834",
        "TRIS cuts (measured)": "#256abf"}
SCALING_ORDER = ["minor", "geometric", "arithmetic", "major"]
KSHORT = {"circular Gaussian (old)": "circular", "elliptical Gaussian": "elliptical",
          "TRIS cuts (measured)": "TRIS cuts"}

# CSV headers.  "residual" = two-step - one-step [K]; "fractional difference" =
# two-step / one-step - 1 [%] (the notebook's metric).  "std" is np.nanstd, i.e. the
# rms about the mean -- what the notebook labels "rms".
CSV_HEADERS = {
    "kernel": "pipeline kernel",
    "scaling": "FWHM scaling",
    "fwhm_given_deg": "FWHM given to pipeline [deg]",
    "pct_median": "median of fractional difference (two-step/one-step - 1) [%]",
    "pct_rms": "std of fractional difference [%]",
    "pct_max_abs": "max |fractional difference| [%]",
    "pct_stripe_rms": "std of fractional difference in TRIS stripe [%]",
    "K_mean": "mean of residual (two-step - one-step) [K]",
    "K_median": "median of residual [K]",
    "K_rms": "std of residual [K]",
    "K_max_abs": "max |residual| [K]",
    "K_stripe_rms": "std of residual in TRIS stripe [K]",
    "kernel_applied_deg": "kernel applied by pipeline [deg]",
}


def complicated_beam(theta, psi):
    phi = np.radians(np.asarray(psi) - CB["pa"])
    sa, sb = CB["major"] / K_, CB["minor"] / K_
    main = np.exp(-0.5 * ((theta * np.cos(phi) / sa) ** 2 + (theta * np.sin(phi) / sb) ** 2))
    ring = 10 ** (CB["ring_db"] / 10) * np.exp(-0.5 * ((theta - CB["ring_r"]) / (CB["ring_w"] / K_)) ** 2)
    shoulder = 10 ** (CB["shoulder_db"] / 10) * np.exp(-0.5 * (theta / (CB["shoulder_fwhm"] / K_)) ** 2)
    return np.where(theta <= 90.0, main + ring + shoulder, 0.0)


def fwhm_from_profile(offsets, response):
    r = np.asarray(response, dtype=float) / response[0]
    i = int(np.argmax(r < 0.5))
    return 2 * (offsets[i - 1] + (0.5 - r[i - 1]) / (r[i] - r[i - 1]) * (offsets[i] - offsets[i - 1]))


def convolve_native(sky, beam_fn, pts, chunk=256):
    # Direct sum on the nside 64 grid: out(p) = sum_q B(p, q) T(q) / sum_q B(p, q).
    out = np.empty(sky.size)
    for s in range(0, len(pts), chunk):
        th, psi = eb.local_offsets(pts[s:s + chunk], pts)
        w = beam_fn(th, psi)
        out[s:s + chunk] = (sky @ w) / w.sum(axis=0)
    return out


def stats(d, band):
    return dict(median=float(np.nanmedian(d)), rms=float(np.nanstd(d)),
                max_abs=float(np.nanmax(np.abs(d))), stripe_rms=float(np.nanstd(d[band])))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=str(REPO / "outputs" / "fwhm_scaling_test"))
    ap.add_argument("--replot-maps", action="store_true",
                    help="only redraw the residual-map figure from the saved .npz and CSV")
    args = ap.parse_args()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    if args.replot_maps:
        return replot_maps(out_dir)

    config = load_config(REPO / "configs" / "run_config.yaml")
    haslam_path = config.resolve_path("paths.haslam_map_fits")
    haslam_native = float(config.require("beam_matching.haslam_native_fwhm_deg"))

    bar = "=" * 78
    print(bar)
    print(f"MAP:        Haslam all-sky  ({haslam_path})")
    print(f"FREQUENCY:  {HASLAM_FREQ_MHZ:.0f} MHz")
    print(f"Grid:       equatorial, nside {NS64} in, nside {NSIDE_OUT} out; "
          f"Haslam native FWHM {haslam_native}°")
    print(f"Test:       notebooks/03 'Two-step against one-step beam matching' "
          f"(PIPELINE.md App. E), complicated survey beam")
    print(f"Metric:     100 * (two-step / one-step - 1) [%]: median, rms, max |Δ|, rms in TRIS stripe")
    print(bar + "\n")

    theta_line = np.arange(0.0, 40.0, 0.01)
    major = fwhm_from_profile(theta_line, complicated_beam(theta_line, np.full(theta_line.size, CB["pa"])))
    minor = fwhm_from_profile(theta_line, complicated_beam(theta_line, np.full(theta_line.size, CB["pa"] + 90)))
    scalings = {"major": major, "minor": minor,
                "arithmetic": 0.5 * (major + minor), "geometric": float(np.sqrt(major * minor))}
    print("FWHM given to the pipeline for the complicated beam (measured main-lobe widths):")
    for name in SCALING_ORDER:
        print(f"  {name:<11} {scalings[name]:.4f}°")
    print()

    # Kernels: the three the notebook ran.  tris_cuts needs the archive cuts.
    specs = {}
    try:
        spec_cuts = eb.beam_spec_from_config(config, kernel="tris_cuts")
        specs = {"circular Gaussian (old)": eb.with_kernel(spec_cuts, "circular_gaussian"),
                 "elliptical Gaussian": eb.with_kernel(spec_cuts, "elliptical_gaussian"),
                 "TRIS cuts (measured)": spec_cuts}
    except Exception as err:  # archive or limTOD missing
        print(f"WARNING: TRIS archive cuts unavailable ({err}); running circular_gaussian only.\n")
        specs = {"circular Gaussian (old)": None}

    def config_with_kernel(kernel):
        c = load_config(config.path)
        c._data["beam_matching"]["kernel"] = kernel
        return c

    configs = {label: config_with_kernel(KERNEL_KEY[label]) for label in specs}

    sky64 = hp.ud_grade(frames.rotate_full_sky(hp.read_map(haslam_path), ("G", "C")), NS64)
    ones64 = np.ones(sky64.size, bool)
    pts64 = np.array(hp.pix2vec(NS64, np.arange(sky64.size))).T
    band = load_stage1(config)["band_mask"].astype(bool)

    def pipeline(sky, native_fwhm, label):
        spec = None if KERNEL_KEY[label] == "circular_gaussian" else specs[label]
        v, kept, info = match_to_tris_beam(sky, ones64, native_fwhm, configs[label], nside_out=NSIDE_OUT,
                                           rule="all_children", spec=spec, full_sky=True)
        return np.where(kept, v, np.nan), info

    t0 = time.time()
    pre = convolve_native(sky64, complicated_beam, pts64)
    print(f"step 1 (complicated survey beam) done in {time.time() - t0:.0f} s")
    direct = {label: pipeline(sky64, haslam_native, label)[0] for label in specs}
    print(f"one-step reference done in {time.time() - t0:.0f} s")

    rows, resid = [], {}
    for label in specs:
        for name in SCALING_ORDER:
            two, info = pipeline(pre, scalings[name], label)
            resid[label, name] = two - direct[label]                     # residual map [K]
            pct = stats(100 * (two / direct[label] - 1), band)
            kel = {"mean": float(np.nanmean(resid[label, name])), **stats(resid[label, name], band)}
            chose = (f"circular {info['extra_fwhm_deg']:.3f}"
                     if info["kernel"] == "circular_gaussian"
                     else f"{info['kernel_shape']} {info['kernel_fwhm_e_deg']:.3f} x {info['kernel_fwhm_h_deg']:.3f}")
            rows.append(dict(kernel=label, scaling=name, fwhm_given_deg=scalings[name],
                             **{f"pct_{k}": v for k, v in pct.items()},
                             **{f"K_{k}": v for k, v in kel.items()},
                             kernel_applied_deg=chose))
        print(f"{label} done in {time.time() - t0:.0f} s")
    print()

    # CSV (header names in CSV_HEADERS; the metadata line comes first, prefixed with #)
    csv_path = out_dir / "fwhm_scaling_table.csv"
    cols = list(CSV_HEADERS)
    with open(csv_path, "w", newline="") as f:
        f.write(f"# map=Haslam {HASLAM_FREQ_MHZ:.0f} MHz ({haslam_path.name}); test=notebook 03 two-step vs one-step"
                f" (complicated survey beam); residual = two-step - one-step [K]; nside {NSIDE_OUT} maps\n")
        w = csv.writer(f)
        w.writerow([CSV_HEADERS[c] for c in cols])
        for r in rows:
            w.writerow([f"{r[c]:.6g}" if isinstance(r[c], float) else r[c] for c in cols])

    # Residual maps, so the plot can be redone without the 7-minute rerun.
    np.savez(out_dir / "fwhm_scaling_residual_maps.npz",
             **{f"{KERNEL_KEY[label]}__{name}": m for (label, name), m in resid.items()})

    # Table to stdout
    print(f"{'kernel':<26}{'scaling':<12}{'FWHM':>8}{'median':>9}{'std':>8}{'max|Δ|':>9}{'stripe std':>12}"
          f"{'mean [K]':>10}{'std [K]':>9}{'max|Δ| [K]':>12}")
    for r in rows:
        print(f"{r['kernel']:<26}{r['scaling']:<12}{r['fwhm_given_deg']:7.3f}°{r['pct_median']:+8.2f}%"
              f"{r['pct_rms']:7.2f}%{r['pct_max_abs']:8.2f}%{r['pct_stripe_rms']:11.2f}%"
              f"{r['K_mean']:+10.4f}{r['K_rms']:9.4f}{r['K_max_abs']:12.4f}")

    # Plot: the metric for each scaling (ordered by the FWHM given), one line per kernel.
    # Evenly spaced categories: geometric and arithmetic are only 0.13 deg apart.
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    panels = [("pct_median", "median"), ("pct_rms", "rms (dotted: TRIS stripe)"), ("pct_max_abs", "max |Δ|")]
    x = np.arange(len(SCALING_ORDER))
    for ax, (key, title) in zip(axes, panels):
        for label in specs:
            rs = {r["scaling"]: r for r in rows if r["kernel"] == label}
            ax.plot(x, [rs[n][key] for n in SCALING_ORDER], "-o", color=KCOL[label], lw=2, ms=6, label=label)
            if key == "pct_rms":
                ax.plot(x, [rs[n]["pct_stripe_rms"] for n in SCALING_ORDER], ":o", color=KCOL[label], lw=1.3, ms=4)
        if key == "pct_median":
            ax.axhline(0, color="0.4", lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{n}\n{scalings[n]:.2f}°" for n in SCALING_ORDER], fontsize=7.5)
        ax.set_xlabel("FWHM given to the pipeline")
        ax.set_ylabel("two-step / one-step − 1 [%]")
        ax.set_title(title, loc="left", fontsize=9)
        ax.grid(axis="y", alpha=0.2)
    axes[0].legend(fontsize=7.5, frameon=False)
    fig.suptitle(f"Haslam {HASLAM_FREQ_MHZ:.0f} MHz, complicated survey beam: FWHM scaling "
                 f"(nside {NS64} → {NSIDE_OUT})", fontsize=10, x=0.01, ha="left")
    plt.tight_layout()
    png_path = out_dir / "fwhm_scaling_comparison.png"
    fig.savefig(png_path, dpi=150, facecolor="white")

    res_path = out_dir / "fwhm_scaling_residual_maps.png"
    plot_residual_maps(resid, scalings, list(specs), res_path)

    print(f"\nwrote {csv_path}\nwrote {png_path}\nwrote {res_path}"
          f"\nwrote {out_dir / 'fwhm_scaling_residual_maps.npz'}")


def plot_residual_maps(resid, scalings, labels, path):
    """One panel per test (rows = kernels, columns = scalings), one shared K scale."""
    divb = LinearSegmentedColormap.from_list("divb", ["#1c5cab", "#86b6ef", "#eeeeec", "#f2a582", "#b2432a"])
    divb.set_bad("#e8e8e6")
    span = float(np.nanpercentile(np.abs(np.concatenate([m for m in resid.values()])), 99))
    nr, nc = len(labels), len(SCALING_ORDER)
    top_in, bottom_in = 0.5, 0.9                         # room for the title and one shared colour bar
    fig = plt.figure(figsize=(4.4 * nc, 3.0 * nr + top_in + bottom_in))
    H = fig.get_figheight()
    w, h = 1.0 / nc, (1.0 - (top_in + bottom_in) / H) / nr
    y0 = bottom_in / H
    for r_, label in enumerate(labels):
        for c_, name in enumerate(SCALING_ORDER):
            m = resid[label, name]
            # Panels placed explicitly (sub=(nr, nc, i) leaves no room for a title or a shared colour bar).
            hp.mollview(m, fig=fig.number, cmap=divb, min=-span, max=span, cbar=False, notext=True,
                        title=f"{KSHORT[label]}, {name} ({scalings[name]:.2f}°)\n"
                              f"mean {np.nanmean(m):+.3f} K, std {np.nanstd(m):.3f} K")
            plt.gca().set_position([c_ * w + 0.005, y0 + (nr - 1 - r_) * h, w - 0.01, h - 0.07])
            hp.graticule(dpar=30, dmer=60, color="#9aa0a6", lw=0.3)
    cax = fig.add_axes([0.3, 0.5 / H, 0.4, 0.16 / H])
    fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(-span, span), cmap=divb), cax=cax,
                 orientation="horizontal", label="residual, two-step − one-step [K]")
    fig.suptitle(f"Residual two-step − one-step [K], Haslam {HASLAM_FREQ_MHZ:.0f} MHz, complicated survey beam, "
                 f"nside {NSIDE_OUT}  (colour scale ±{span:.2f} K = 99th percentile of |residual|)",
                 x=0.01, y=0.995, ha="left", va="top", fontsize=11)
    fig.savefig(path, dpi=150, facecolor="white")


def replot_maps(out_dir):
    """Redraw the residual-map figure from the saved .npz and CSV (no pipeline rerun)."""
    z = np.load(out_dir / "fwhm_scaling_residual_maps.npz")
    key_label = {v: k for k, v in KERNEL_KEY.items()}
    resid = {(key_label[k.split("__")[0]], k.split("__")[1]): z[k] for k in z.files}
    with open(out_dir / "fwhm_scaling_table.csv") as f:
        rows = list(csv.DictReader(line for line in f if not line.startswith("#")))
    scalings = {r[CSV_HEADERS["scaling"]]: float(r[CSV_HEADERS["fwhm_given_deg"]]) for r in rows}
    labels = [l for l in KERNEL_KEY if any(k[0] == l for k in resid)]
    path = out_dir / "fwhm_scaling_residual_maps.png"
    plot_residual_maps(resid, scalings, labels, path)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
