"""Does beam averaging change the SHAPE of T vs LST across frequencies?

Haslam 408 MHz is extrapolated to each TRIS frequency with a position-dependent
spectral index beta(pixel) -- a uniform background plus a few small Gaussian blobs
on the TRIS stripe -- and observed through the TRIS beam along the TRIS ring.  A
control run uses the background beta everywhere.  For the control, T_nu / <T_nu>
must be the same curve at every frequency (the scaling is one constant and the
beam is linear).  For the blob run, any difference from the control in the
normalised curve is a change of shape, not of normalisation.

Repo code used, read only:
* ``moonrunes.sky_maps.describe``  -- loads Haslam (the repo's map loader);
* ``moonrunes.sed_prototype.Inputs._native_tod`` / ``Inputs.observe`` -- the TOD
  route's own "observe a native-resolution map with the TRIS beam along the ring"
  (limTOD's measured cut beam, zenith pointing, horizon mask), called unchanged on
  a stand-in object that carries only the attributes those two methods read;
* ``limTOD.tris`` (via ``import_limtod_tris``) -- ring LSTs, beam, horizon.

Outputs go to ``beam_averaging_blobs/<map stem>/`` (set by ``HASLAM_FILE``):
summary.csv, t_vs_lst.csv, t_vs_lst_comparison.png, beta_map.png.  See README.md in
this folder; ``replot_from_csv.py`` redraws the summary figures from the CSVs.

Run (``limtod`` conda env; a few seconds):

    python beam_averaging_blobs/beam_averaging_blobs.py
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from types import SimpleNamespace

os.environ.setdefault("TQDM_DISABLE", "1")

import numpy as np
import healpy as hp
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(REPO / "src"))

from moonrunes import sky_maps  # noqa: E402
from moonrunes.config import load_config  # noqa: E402
from moonrunes.sed_prototype import Inputs  # noqa: E402
from moonrunes.stage1_tris_maps import import_limtod_tris  # noqa: E402

# ---------------------------------------------------------------------------
# PARAMETERS -- fixed before the first run; not tuned to the result.
# ---------------------------------------------------------------------------
HASLAM_FILE = "HASLAM_DESTRIPED_ONLY.fits"   # in the repo's res/ folder; "HASLAM.fits" was the first run
NU_REF_MHZ = 408.0
BETA0 = -2.7                         # background spectral index
N_BLOBS = 4                           # 3-5 allowed; the lists below must have this length
BLOB_RA_DEG = [30.0, 120.0, 210.0, 300.0]     # = LST of the blob centre (zenith pointing)
BLOB_DEC_DEG = None                   # None -> on the stripe centre (the TRIS site latitude)
BLOB_FWHM_DEG = [4.0, 4.0, 4.0, 4.0]          # Gaussian FWHM of each blob; TRIS beam is 19.2 x 23.4
BLOB_BETA = [-2.5, -2.9, -2.6, -2.8]          # beta at each blob centre
NORMALISATION = "mean"                # T_nu(LST) / mean over the ring's LST samples
OUT_DIR = HERE / Path(HASLAM_FILE).stem       # one subfolder per input map

assert len(BLOB_RA_DEG) == len(BLOB_FWHM_DEG) == len(BLOB_BETA) == N_BLOBS and 3 <= N_BLOBS <= 5
assert all(-2.9 <= b <= -2.5 for b in BLOB_BETA)


def beta_map(nside, frame, ra_c, dec_c, fwhm, betas, beta0):
    """beta at every pixel of a map in ``frame``: beta0 + sum_i (beta_i - beta0) G_i(r_i)."""
    vec = np.array(hp.pix2vec(nside, np.arange(hp.nside2npix(nside))))
    if frame != "C":
        vec = hp.Rotator(coord=[frame, "C"])(vec)            # pixel centres in equatorial
    beta = np.full(vec.shape[1], beta0, dtype=float)
    for ra, dec, fw, b in zip(ra_c, dec_c, fwhm, betas):
        c = np.array(hp.ang2vec(ra, dec, lonlat=True))
        r = np.degrees(np.arccos(np.clip(c @ vec, -1.0, 1.0)))
        beta += (b - beta0) * np.exp(-0.5 * (r / (fw / (2 * np.sqrt(2 * np.log(2))))) ** 2)
    return beta


def main():
    config = load_config(REPO / "configs" / "run_config.yaml")
    tris = import_limtod_tris(config)
    haslam = sky_maps.describe(config.resolve_path("paths.haslam_map_fits").parent / HASLAM_FILE)
    assert abs(haslam["freq_mhz"] - NU_REF_MHZ) < 1 and haslam["in_kelvin"]
    assert haslam["n_unobserved"] == 0, "the β extrapolation and no-beam sampling assume a full-sky map"
    OUT_DIR.mkdir(exist_ok=True)
    tris_freqs = [float(f) for f in config.require("tris.frequencies_mhz")]
    freqs = [NU_REF_MHZ] + tris_freqs

    # The observer: the same attributes sed_prototype.Inputs builds for the TOD route.
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
    lst = np.asarray(obs.geom.lst_deg, dtype=float)
    stripe_dec = float(obs.geom.latitude_deg)
    dec_c = [stripe_dec] * N_BLOBS if BLOB_DEC_DEG is None else list(BLOB_DEC_DEG)

    bar = "=" * 78
    print(bar)
    print(f"MAP:          {haslam['name']}  ({haslam['path']})")
    print(f"              nside {haslam['nside']}, frame {haslam['frame']}, {haslam['unit']}, "
          f"CMB {haslam['cmb']}, native FWHM {haslam['fwhm_deg']:.3f}°")
    print(f"FREQUENCIES:  {', '.join(f'{f:g}' for f in freqs)} MHz  (408 = reference; TRIS from "
          f"config tris.frequencies_mhz)")
    print(f"BEAM:         limTOD TRIS cut beam (archive), zenith pointing at Dec {stripe_dec:.3f}°, "
          f"selfrot {obs.geom.selfrot_deg[0]:+.0f}°, observed at nside {nside} (nside_hires {obs.hires})")
    print(f"STRIPE:       {lst.size} LST samples, {lst.min():.2f}–{lst.max():.2f}° (LST = RA)")
    print(f"BLOBS:        beta0 = {BETA0}")
    for i in range(N_BLOBS):
        print(f"              #{i + 1}: RA {BLOB_RA_DEG[i]:6.1f}°, Dec {dec_c[i]:.2f}°, "
              f"FWHM {BLOB_FWHM_DEG[i]:.1f}°, beta {BLOB_BETA[i]:+.2f}")
    print(bar + "\n")

    # Spectral index on Haslam's own grid, so the power law is applied before any averaging.
    beta = beta_map(haslam["nside"], haslam["frame"], BLOB_RA_DEG, dec_c, BLOB_FWHM_DEG, BLOB_BETA, BETA0)
    t408 = haslam["sky"]

    # Ring positions in Haslam's frame, for the no-beam (pencil) diagnostic.
    ring_vec = np.array(hp.ang2vec(lst, np.full(lst.size, stripe_dec), lonlat=True)).T
    if haslam["frame"] != "C":
        ring_vec = hp.Rotator(coord=["C", haslam["frame"]])(ring_vec)
    ring_th, ring_ph = hp.vec2ang(ring_vec.T)

    tod, pencil = {}, {}
    for case, b in (("blob", beta), ("control", np.full_like(beta, BETA0))):
        for nu in freqs:
            t_nu = t408 * (nu / NU_REF_MHZ) ** b
            m = dict(sky=t_nu, nside=haslam["nside"], frame=haslam["frame"])
            tod[case, nu], coverage = Inputs._native_tod(obs, m)
            assert np.all(coverage > 0)
            pencil[case, nu] = hp.get_interp_val(t_nu, ring_th, ring_ph)
        print(f"observed: {case}")
    print()

    norm = (lambda t: t / np.mean(t)) if NORMALISATION == "mean" else None
    S = {k: norm(v) for k, v in tod.items()}
    P = {k: norm(v) for k, v in pencil.items()}

    # 4. Control: the normalised curve must be identical at every frequency.
    ctrl_dev = {nu: float(np.max(np.abs(S["control", nu] / S["control", NU_REF_MHZ] - 1))) for nu in freqs}
    worst = max(ctrl_dev.values())
    print(f"CONTROL (uniform beta = {BETA0}): max |S_nu / S_408 − 1| over LST and frequency = {worst:.2e}")
    print("  -> the normalised shape IS identical at all frequencies (to rounding)." if worst < 1e-10 else
          "  -> the normalised shape is NOT identical across frequencies.")
    print()

    # 5. Blob case against the control.
    rows = []
    print(f"{'nu [MHz]':>9}{'max shape Δ':>13}{'at LST':>9}{'norm Δ':>10}{'ratio blob min–max':>22}"
          f"{'ratio ctrl':>12}{'no-beam shape Δ':>17}")
    for nu in freqs:
        shape = 100 * (S["blob", nu] / S["control", nu] - 1)
        k = int(np.argmax(np.abs(shape)))
        pshape = 100 * (P["blob", nu] / P["control", nu] - 1)
        r_blob = tod["blob", nu] / tod["blob", NU_REF_MHZ]
        r_ctrl = tod["control", nu] / tod["control", NU_REF_MHZ]
        row = dict(freq_mhz=nu,
                   max_abs_shape_change_pct=float(np.max(np.abs(shape))), at_lst_deg=float(lst[k]),
                   rms_shape_change_pct=float(np.sqrt(np.mean(shape ** 2))),
                   normalisation_change_pct=float(100 * (np.mean(tod["blob", nu]) / np.mean(tod["control", nu]) - 1)),
                   ratio_blob_min=float(r_blob.min()), ratio_blob_max=float(r_blob.max()),
                   ratio_control=float(np.mean(r_ctrl)),
                   control_shape_vs_408_max_abs=ctrl_dev[nu],
                   no_beam_max_abs_shape_change_pct=float(np.max(np.abs(pshape))))
        rows.append(row)
        print(f"{nu:9.1f}{row['max_abs_shape_change_pct']:12.3f}%{row['at_lst_deg']:8.2f}°"
              f"{row['normalisation_change_pct']:+9.3f}%{row['ratio_blob_min']:11.4f}–{row['ratio_blob_max']:.4f}"
              f"{row['ratio_control']:12.4f}{row['no_beam_max_abs_shape_change_pct']:16.3f}%")
    print("\nshape Δ = 100 (S_blob / S_control − 1), S = T / mean_LST(T); norm Δ = mean T_blob / mean T_control − 1;"
          "\nratio = T_nu / T_408 along the ring; no-beam = the same shape Δ sampled at the ring points without the beam.")

    # CSVs
    with open(OUT_DIR / "summary.csv", "w") as f:
        f.write(f"# map = {haslam['file']}; shape change = 100*(S_blob/S_control - 1), S = T/mean_LST(T);"
                " beam = limTOD TRIS cut beam\n")
        f.write(",".join(rows[0]) + "\n")
        for r in rows:
            f.write(",".join(f"{v:.8g}" for v in r.values()) + "\n")
    with open(OUT_DIR / "t_vs_lst.csv", "w") as f:
        cols = ["lst_deg"]
        for nu in freqs:
            cols += [f"T_blob_{nu:g}_K", f"T_control_{nu:g}_K", f"S_blob_{nu:g}", f"S_control_{nu:g}",
                     f"shape_change_pct_{nu:g}", f"T_blob_minus_control_{nu:g}_K"]
        f.write(",".join(cols) + "\n")
        for j in range(lst.size):
            vals = [lst[j]]
            for nu in freqs:
                vals += [tod["blob", nu][j], tod["control", nu][j], S["blob", nu][j], S["control", nu][j],
                         100 * (S["blob", nu][j] / S["control", nu][j] - 1),
                         tod["blob", nu][j] - tod["control", nu][j]]
            f.write(",".join(f"{v:.8g}" for v in vals) + "\n")

    # Plot: one column per frequency; rows = normalised T, T_nu / T_408, shape change.
    blob_col, ctrl_col = "#b2432a", "#256abf"
    # Rows 0 and 2 share y across frequencies; the ratio row does not (its level differs by
    # frequency and the blob variation is ~0.2% of it).
    fig, axes = plt.subplots(3, len(freqs), figsize=(4.6 * len(freqs), 9.5), sharex=True)
    for c in range(1, len(freqs)):
        axes[0, c].sharey(axes[0, 0]); axes[2, c].sharey(axes[2, 0])
    for c, nu in enumerate(freqs):
        ax = axes[0, c]
        ax.plot(lst, S["control", nu], color=ctrl_col, lw=2, label=f"control (β = {BETA0})")
        ax.plot(lst, S["blob", nu], color=blob_col, lw=1.4, ls="--", label="blobs")
        ax.set_title(f"{nu:g} MHz", loc="left", fontsize=10)
        ax = axes[1, c]
        ax.plot(lst, tod["control", nu] / tod["control", NU_REF_MHZ], color=ctrl_col, lw=2)
        ax.plot(lst, tod["blob", nu] / tod["blob", NU_REF_MHZ], color=blob_col, lw=1.4, ls="--")
        ax = axes[2, c]
        ax.plot(lst, 100 * (S["blob", nu] / S["control", nu] - 1), color=blob_col, lw=1.8)
        ax.axhline(0, color="0.4", lw=0.8)
        for a in axes[:, c]:
            for i in range(N_BLOBS):
                a.axvline(BLOB_RA_DEG[i], color="0.82", lw=0.9, zorder=0)
            a.grid(axis="y", alpha=0.2)
        for i in range(N_BLOBS):
            axes[2, c].annotate(f"β {BLOB_BETA[i]:+.1f}", (BLOB_RA_DEG[i], 1), xycoords=("data", "axes fraction"),
                                xytext=(2, -10), textcoords="offset points", fontsize=7, color="0.35")
        axes[2, c].set_xlabel("LST [deg]")
        axes[2, c].set_xlim(0, 360); axes[2, c].set_xticks(np.arange(0, 361, 60))
    axes[0, 0].set_ylabel("T / mean(T)")
    axes[1, 0].set_ylabel("T_ν / T_408")
    axes[2, 0].set_ylabel("shape change: S_blob / S_control − 1 [%]")
    axes[0, 0].legend(fontsize=8, frameon=False)
    fig.suptitle(f"{haslam['name']} ({haslam['file']}) extrapolated with blob β map, observed with the TRIS beam "
                 "along the ring (grey lines: blob LSTs)", x=0.01, ha="left", fontsize=10)
    plt.tight_layout()
    fig.savefig(OUT_DIR / "t_vs_lst_comparison.png", dpi=150, facecolor="white")

    # The beta map itself, in equatorial, with the stripe (for checking where the blobs are).
    beta_eq = beta_map(64, "C", BLOB_RA_DEG, dec_c, BLOB_FWHM_DEG, BLOB_BETA, BETA0)
    plt.figure(figsize=(9, 5))
    hp.mollview(beta_eq, fig=plt.gcf().number, cmap="RdBu_r", min=-2.9, max=-2.5, unit="β",
                title=f"β map (equatorial); dots: TRIS ring at Dec {stripe_dec:.2f}°", notext=True)
    hp.projscatter(lst, np.full(lst.size, stripe_dec), lonlat=True, s=2, color="k")
    hp.graticule(dpar=30, dmer=60, color="#9aa0a6", lw=0.3)
    plt.savefig(OUT_DIR / "beta_map.png", dpi=150, facecolor="white")

    print(f"\nwrote {OUT_DIR / 'summary.csv'}\nwrote {OUT_DIR / 't_vs_lst.csv'}"
          f"\nwrote {OUT_DIR / 't_vs_lst_comparison.png'}\nwrote {OUT_DIR / 'beta_map.png'}")


if __name__ == "__main__":
    main()
