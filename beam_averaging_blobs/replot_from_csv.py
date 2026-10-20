"""Replot the beam-averaging blob results from the saved CSVs only (no analysis is rerun).

Reads <run>/t_vs_lst.csv and <run>/summary.csv; blob positions and betas are the
parameters at the top of beam_averaging_blobs.py (imported, not recomputed).

    python beam_averaging_blobs/replot_from_csv.py [RUN]      # RUN defaults to HASLAM_DESTRIPED_ONLY
"""

from __future__ import annotations

import sys
import textwrap
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from beam_averaging_blobs import BLOB_BETA, BLOB_RA_DEG  # noqa: E402  (parameters only)

RUN = sys.argv[1] if len(sys.argv) > 1 else "HASLAM_DESTRIPED_ONLY"
RUN_DIR = HERE / RUN
FREQS = [600, 820]
STYLE = {600: dict(color="#256abf", ls="-"), 820: dict(color="#eb6834", ls="--")}

plt.rcParams.update({"font.size": 14, "axes.titlesize": 15, "axes.labelsize": 14,
                     "xtick.labelsize": 13, "ytick.labelsize": 13, "legend.fontsize": 13})

d = np.genfromtxt(RUN_DIR / "t_vs_lst.csv", delimiter=",", names=True)
s = np.genfromtxt(RUN_DIR / "summary.csv", delimiter=",", names=True, skip_header=1)   # line 1 = metadata
summary = {int(r["freq_mhz"]): r for r in np.atleast_1d(s)}
lst = d["lst_deg"]


def caption(fig, text, y=0.01):
    fig.text(0.01, y, textwrap.fill(text, 125), ha="left", va="bottom", fontsize=12, color="0.25")


# ---------------------------------------------------------------------------
# Plot A: change in shape vs LST, 600 and 820 MHz.
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(13, 7.2))
for ra, beta in zip(BLOB_RA_DEG, BLOB_BETA):
    ax.axvline(ra, color="0.75", lw=1.2, zorder=0)
    ax.annotate(f"blob, β = {beta:.1f}", (ra, 1.0), xycoords=("data", "axes fraction"), xytext=(5, -18),
                textcoords="offset points", fontsize=12, color="0.35")
ax.axhline(0, color="0.4", lw=0.9)
for nu in FREQS:
    y = d[f"shape_change_pct_{nu}"]
    ax.plot(lst, y, lw=2.4, label=f"{nu} MHz", **STYLE[nu])
    k = int(np.argmax(np.abs(y)))
    assert abs(abs(y[k]) - summary[nu]["max_abs_shape_change_pct"]) < 1e-6   # same number as summary.csv
    ax.plot(lst[k], y[k], "o", ms=9, color=STYLE[nu]["color"], zorder=5)
    ax.annotate(f"largest at {nu} MHz:\n{y[k]:+.2f}% at LST {lst[k]:.0f}°", (lst[k], y[k]),
                xytext=(-200, -10 if y[k] < 0 else 10), textcoords="offset points", fontsize=12,
                color=STYLE[nu]["color"], arrowprops=dict(arrowstyle="-", color=STYLE[nu]["color"], lw=1))
ax.set_xlim(0, 360)
ax.set_xticks(np.arange(0, 361, 60))
ax.set_xlabel("LST (degrees)")
ax.set_ylabel("Change in shape of T vs LST (%)")
ax.set_title("Change in the shape of T vs LST caused by the spectral-index blobs, seen through the TRIS beam",
             loc="left", fontsize=14)
ax.legend(frameon=False, loc="lower left")
ax.grid(axis="y", alpha=0.25)
fig.subplots_adjust(left=0.09, right=0.98, top=0.92, bottom=0.27)
caption(fig, "Each curve compares the blob sky with a sky that has one spectral index (−2.7) everywhere. "
             "Both are normalised by their mean along the ring, so only differences in shape remain; "
             "0% means the shape is unchanged. Each blob leaves a bump (flatter spectrum) or a dip (steeper "
             "spectrum) about one beam wide, larger at 820 MHz because it is further from 408 MHz. "
             f"Map: {RUN}.fits.", y=0.015)
out_a = RUN_DIR / "plotA_shape_change_vs_lst.png"
fig.savefig(out_a, dpi=150, facecolor="white")

# ---------------------------------------------------------------------------
# Plot B: largest change without the beam vs with the beam (summary.csv only).
# ---------------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 6.6))
x = np.arange(len(FREQS))
for ax, key, ylab, title in (
        (ax1, "no_beam_max_abs_shape_change_pct", "Largest change in shape, no beam (%)",
         "Ring sampled without the beam"),
        (ax2, "max_abs_shape_change_pct", "Largest change in shape, with TRIS beam (%)",
         "Ring seen through the TRIS beam")):
    vals = [summary[nu][key] for nu in FREQS]
    bars = ax.bar(x, vals, width=0.55, color=[STYLE[nu]["color"] for nu in FREQS])
    for b, v in zip(bars, vals):
        ax.annotate(f"{v:.2f}%", (b.get_x() + b.get_width() / 2, v), xytext=(0, 4),
                    textcoords="offset points", ha="center", fontsize=13)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{nu} MHz" for nu in FREQS])
    ax.set_ylabel(ylab)
    ax.set_title(title, loc="left")
    ax.set_ylim(0, 1.18 * max(vals))
    ax.grid(axis="y", alpha=0.25)
nb = {nu: summary[nu]["no_beam_max_abs_shape_change_pct"] for nu in FREQS}
wb = {nu: summary[nu]["max_abs_shape_change_pct"] for nu in FREQS}
fig.suptitle("Largest change in the shape of T vs LST, without and with the TRIS beam (note the different y-axes)",
             x=0.01, ha="left", fontsize=14)
fig.subplots_adjust(left=0.08, right=0.98, top=0.86, bottom=0.25, wspace=0.3)
caption(fig, "Left: the largest change in shape if the ring is sampled at single points, with no beam. "
             "Right: the same quantity seen through the TRIS beam. The blobs are 4° wide and the beam is "
             "about 19° × 23°, so each blob fills only a small part of the beam and its effect is diluted: "
             f"the largest change drops from {nb[600]:.1f}% to {wb[600]:.2f}% at 600 MHz and from {nb[820]:.1f}% "
             f"to {wb[820]:.2f}% at 820 MHz. The two largest values need not come from the same blob. Only "
             f"these maximum values were saved, so the no-beam result is shown as bars, not as a curve "
             f"against LST. Map: {RUN}.fits.", y=0.015)
out_b = RUN_DIR / "plotB_beam_dilution.png"
fig.savefig(out_b, dpi=150, facecolor="white")

print(f"wrote {out_a}\nwrote {out_b}")
