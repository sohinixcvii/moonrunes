"""Bayesian comparison of three spectral models with UltraNest (skeleton).

Models: power law, curved power law, moment expansion (models.py).  Data: the SED
points of the existing stage 3 prototype (``moonrunes.sed_prototype.Inputs``), per TRIS
ring sample (TOD route) or per TRIS-stripe pixel (map route).  Error bars come from
``error_bars.csv``, which marks placeholders.  All choices are in ``config.yaml``.

For each unit (ring sample or pixel) and model, UltraNest returns ln Z.  The units are
treated as independent, so a model's total ln Z is the sum over units.

    python model_comparison/run_comparison.py [--config model_comparison/config.yaml]

Needs ``ultranest`` importable (installed in the limtod environment since 2026-10-21).
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("TQDM_DISABLE", "1")

import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(HERE))

from moonrunes import sky_maps  # noqa: E402
from moonrunes.config import load_config  # noqa: E402
from moonrunes.sed_prototype import Inputs, ring_sigma  # noqa: E402
from models import build_models, prior_transform  # noqa: E402

BANNER = "!" * 78


# ---------------------------------------------------------------------------
# error bars
# ---------------------------------------------------------------------------
def read_error_bars(path):
    rows = {}
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            r["is_placeholder"] = r["is_placeholder"].strip().lower() in ("true", "1", "yes")
            rows[r["key"]] = r
    return rows


def sigma_from_row(row, t_excess):
    """Error bar for a non-TRIS point from its CSV row."""
    if row["kind"] == "fraction":
        return float(row["value"]) * abs(float(t_excess))
    if row["kind"] == "absolute_k":
        return float(row["value"])
    if row["kind"] == "published":   # scale x |T|, zero level and noise, in quadrature
        return float(np.sqrt((float(row["scale_frac"]) * float(t_excess)) ** 2
                             + float(row["zero_level_k"]) ** 2 + float(row["noise_k"]) ** 2))
    raise ValueError("error_bars.csv: kind {!r} not handled for {}".format(row["kind"], row["key"]))


# ---------------------------------------------------------------------------
# SED points per unit, from the existing loader
# ---------------------------------------------------------------------------
def tris_key(f):
    return "TRIS_{}".format(int(f))


def tod_units(inp, errors):
    """One SED per ring sample: as sed_prototype.tod_route builds it (sed_prototype.py:253-264),
    with the non-TRIS error bars taken from the CSV."""
    for k in range(inp.ra.size):
        pts = []
        for f, ring in inp.rings.items():
            nu = float(ring.effective_frequency_mhz)
            pts.append((tris_key(f), nu, ring.temperature_k[k] - sky_maps.cmb_rj_k(nu), ring_sigma(ring)[k]))
        for m in inp.fit_maps:
            d = inp.tods[m["name"]]
            if d["good"][k]:
                t_ex = float(sky_maps.rj_excess(d["tod"][k], m))
                pts.append((m["file"], float(m["freq_mhz"]), t_ex, sigma_from_row(errors[m["file"]], t_ex)))
        yield int(k), float(inp.ra[k]), pts


def map_units(inp, errors, add_zhat):
    """One SED per stripe pixel: as sed_prototype.map_route (sed_prototype.py:140-147, 202-213),
    with z-hat added to TRIS only if ``add_zhat``."""
    st = inp.stage1
    tris = {}
    for f in (600, 820):
        i = st.index(f)
        nu = float(st["effective_freq_mhz"][i])
        t, s = st.band_map("sky_k", f), st.band_map("sky_sigma_k", f)
        if add_zhat:
            t = t + float(st["zero_level_k"][i])
            s = np.hypot(s, float(st["zero_level_sigma_k"][i]))
        tris[f] = (nu, t, s)
    for p in np.flatnonzero(inp.band):
        pts = [(tris_key(f), nu, t[p] - sky_maps.cmb_rj_k(nu), s[p]) for f, (nu, t, s) in tris.items()]
        for m in inp.fit_maps:
            v = inp.grid[m["name"]][p]
            if np.isfinite(v):
                t_ex = float(sky_maps.rj_excess(v, m))
                pts.append((m["file"], float(m["freq_mhz"]), t_ex, sigma_from_row(errors[m["file"]], t_ex)))
        yield int(p), np.nan, pts


def usable(pts):
    """The same rule as sed_prototype.fit_power_law (sed_prototype.py:54): finite, s > 0, t > 0."""
    return [q for q in pts if np.isfinite(q[2]) and np.isfinite(q[3]) and q[3] > 0 and q[2] > 0]


# ---------------------------------------------------------------------------
# evidence
# ---------------------------------------------------------------------------
def log_evidence(names, model_fn, transform, nu, t, s, scfg):
    import ultranest

    norm = -np.sum(np.log(s * np.sqrt(2 * np.pi)))

    def loglike(params):
        r = (t - model_fn(nu, params)) / s
        return -0.5 * np.sum(r * r, axis=-1) + norm

    sampler = ultranest.ReactiveNestedSampler(names, loglike, transform, vectorized=True,
                                              log_dir=scfg.get("log_dir"))
    logging.getLogger("ultranest").setLevel(logging.WARNING)   # the sampler resets it to INFO
    res = sampler.run(min_num_live_points=int(scfg["min_num_live_points"]), dlogz=float(scfg["dlogz"]),
                      show_status=False, viz_callback=False)
    return float(res["logz"]), float(res["logzerr"])


def placeholder_banner(used, errors, n_points):
    n_ph = sum(used.values())
    print("\n" + BANNER)
    print("!!  WARNING: PLACEHOLDER ERROR BARS WERE USED  ".ljust(76) + "!!")
    print("!!  {} of {} data points ({:.0f}%) carry a placeholder error bar:".format(
        n_ph, n_points, 100 * n_ph / max(n_points, 1)).ljust(76) + "!!")
    for key, n in sorted(used.items(), key=lambda kv: -kv[1]):
        row = errors[key]
        what = "{:g}% of T".format(100 * float(row["value"])) if row["kind"] == "fraction" else row["kind"]
        print("!!    {:<34} {:>4} points   ({})".format(row["name"], n, what).ljust(76) + "!!")
    print("!!  The log-evidences below depend on these assumed errors.".ljust(76) + "!!")
    print(BANNER + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", default=str(HERE / "config.yaml"))
    args = ap.parse_args()
    cfg = yaml.safe_load(open(args.config))
    logging.getLogger("ultranest").setLevel(logging.WARNING)

    route = cfg["data"]["route"]
    add_zhat = bool(cfg["zero_level"]["add_zhat_to_tris"])
    errors = read_error_bars(HERE / cfg["data"]["error_bars_csv"])
    models = build_models(cfg)

    print("=" * 78)
    print("Spectral model comparison (UltraNest)")
    print("  route:        {}".format(route))
    print("  zero level:   add z-hat to TRIS = {}{}".format(
        add_zhat, "  (no effect on the TOD route: ring data contain no z-hat)" if route == "tod" else ""))
    print("  models:       {}".format(", ".join(models)))
    bfix = cfg["models"].get("moment_beta_fixed")
    print("  moment exp.:  {}".format("beta free, omega1 = 0" if bfix is None
                                      else "beta fixed at {:g}, omega1 free".format(float(bfix))))
    print("  error bars:   {}".format(cfg["data"]["error_bars_csv"]))
    print("=" * 78)

    config = load_config(REPO / "configs" / "run_config.yaml")
    inp = Inputs(config, REPO / "res", target_file=cfg["data"]["target_file"], verbose=True)
    missing = [m["file"] for m in inp.fit_maps if m["file"] not in errors]
    if missing:
        raise KeyError("no error-bar row in the CSV for: {}".format(missing))

    units = tod_units(inp, errors) if route == "tod" else map_units(inp, errors, add_zhat)
    wanted = cfg["data"]["units"]
    min_points = int(cfg["data"]["min_points"])

    rows, used_ph, n_points, skipped = [], {}, 0, 0
    t0 = time.time()
    for uid, ra, pts in units:
        if wanted is not None and uid not in wanted:
            continue
        pts = usable(pts)
        if len(pts) < min_points:
            skipped += 1
            continue
        keys = [q[0] for q in pts]
        nu, t, s = (np.array([q[i] for q in pts], dtype=float) for i in (1, 2, 3))
        n_points += len(pts)
        for key in keys:
            if errors[key]["is_placeholder"]:
                used_ph[key] = used_ph.get(key, 0) + 1
        row = dict(unit=uid, ra_deg=ra, n_points=len(pts), inputs=";".join(keys))
        for name, (pnames, fn) in models.items():
            lz, lze = log_evidence(pnames, fn, prior_transform(pnames, cfg["priors"]), nu, t, s, cfg["sampler"])
            row["logZ " + name], row["logZerr " + name] = lz, lze
        rows.append(row)
        print("  unit {:>5}: {} points, {} -- {:.0f} s".format(
            uid, len(pts), "  ".join("{:+.2f}".format(row["logZ " + n]) for n in models), time.time() - t0))

    if not rows:
        raise RuntimeError("no unit had >= {} usable points".format(min_points))

    out_dir = HERE / "results"
    out_dir.mkdir(exist_ok=True)
    with open(out_dir / "evidence_per_unit_{}.csv".format(route), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    total = {n: sum(r["logZ " + n] for r in rows) for n in models}
    total_err = {n: float(np.sqrt(sum(r["logZerr " + n] ** 2 for r in rows))) for n in models}
    wins = {n: sum(1 for r in rows if max(models, key=lambda m: r["logZ " + m]) == n) for n in models}
    ref = "power law"

    if used_ph:
        placeholder_banner(used_ph, errors, n_points)
    print("{} units fitted ({} skipped with < {} points), {} data points, {:.0f} s".format(
        len(rows), skipped, min_points, n_points, time.time() - t0))
    print("\n{:<30}{:>14}{:>10}{:>18}{:>14}".format("model", "sum ln Z", "± err", "ln Z − power law", "units won"))
    for n in models:
        print("{:<30}{:>14.2f}{:>10.2f}{:>18.2f}{:>14d}".format(n, total[n], total_err[n], total[n] - total[ref], wins[n]))
    with open(out_dir / "evidence_summary_{}.csv".format(route), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["model", "sum_logZ", "sum_logZ_err", "logZ_minus_power_law", "units_won", "n_units",
                    "route", "add_zhat_to_tris", "placeholder_points", "n_points", "moment_beta_fixed"])
        for n in models:
            w.writerow([n, total[n], total_err[n], total[n] - total[ref], wins[n], len(rows), route, add_zhat,
                        sum(used_ph.values()), n_points, "free" if bfix is None else float(bfix)])
    if used_ph:
        placeholder_banner(used_ph, errors, n_points)
    print("wrote {}".format(out_dir))


if __name__ == "__main__":
    main()
