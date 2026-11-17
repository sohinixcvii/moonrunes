"""Bayesian comparison of three spectral models with UltraNest (skeleton).

Models: power law, curved power law, and the moment expansion at each order in
``models.moment_orders`` (models.py).  Data: the SED
points of the existing stage 3 prototype (``moonrunes.sed_prototype.Inputs``), per TRIS
ring sample (TOD route) or per TRIS-stripe pixel (map route).  Error bars come from
``configs/error_bars.csv``, which marks placeholders.  The model choices are in
``configs/model_comp_config.yaml``; the data (maps, target, beam, ring) in ``configs/run_pipeline.yaml``.

The data sets are the pipeline's: every map ``sed.exclude_maps`` leaves in, plus TRIS 600 and 820 unless
``comparison.exclude_tris`` (true: the comparison is made without TRIS, as pipeline step 3).  Pipeline step 3
(``python run_pipeline.py --step compare``) is the paper's comparison; this script is the standalone version, which
also runs the map route and writes to ``paper_outputs/model_comparison/results/``.

``comparison.statistic`` in the config picks what the models are compared on:

* ``max_likelihood``: for each unit and model, the best fit inside the prior box (minimum
  chi^2, so ln L_max = -chi2/2 + const).  No Occam penalty.  Writes ``maxlike_*.csv``.
* ``evidence``: UltraNest's ln Z, which includes the Occam penalty.  The units are
  treated as independent, so a model's total ln Z is the sum over units.  Writes
  ``evidence_*.csv``.

    python scripts/model_comparison/run_comparison.py [--config configs/model_comp_config.yaml]

Needs ``ultranest`` importable (installed in the limtod environment).
"""

from __future__ import annotations

import argparse
import copy
import csv
import logging
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("TQDM_DISABLE", "1")

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]                 # scripts/model_comparison -> the repo
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(HERE))

from moonrunes import sky_maps  # noqa: E402
from moonrunes.config import (DEFAULT_MODEL_COMP_CONFIG_PATH, Config, error_bars_path, load_config,  # noqa: E402
                              load_model_comp_config)
from moonrunes.sed_prototype import Inputs, ring_sigma, tod_beam as sed_tod_beam  # noqa: E402
from models import build_models, moment_orders, prior_transform  # noqa: E402

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
TRIS_BANDS = (600, 820)


def select_datasets(datasets, config, res_dir, exclude_tris=False):
    """Which data sets go into every SED.

    ``datasets`` is None -- every fit map that run_pipeline.yaml's ``sed.exclude_maps`` leaves in, plus TRIS 600
    and 820 unless ``exclude_tris`` (model_comp_config.yaml ``comparison.exclude_tris``): what the pipeline,
    run_comparison.py and test_results Part C use -- or a list of keys as in error_bars.csv: ``TRIS_600``,
    ``TRIS_820`` and map file names.  A list is used exactly (``exclude_tris`` is ignored), overriding
    ``sed.exclude_maps`` for that caller only (test_results Part J; run_pipeline.yaml is not edited).
    Returns (run config to build ``Inputs`` with, TRIS bands to include, map files included)."""
    fit_files = [m["file"] for m in sky_maps.discover(res_dir) if m["usable"] and m["role"] == "fit"]
    if datasets is None:
        excluded = set(config.get_path("sed.exclude_maps", []) or [])
        return config, (() if exclude_tris else TRIS_BANDS), [f for f in fit_files if f not in excluded]
    datasets = [str(d) for d in datasets]
    valid = [tris_key(f) for f in TRIS_BANDS] + fit_files
    bad = [d for d in datasets if d not in valid]
    if bad:
        raise ValueError("datasets: unknown or unusable key(s) {}; choose from {}".format(bad, valid))
    if len(set(datasets)) != len(datasets):
        raise ValueError("datasets lists a key twice: {}".format(datasets))
    tris = tuple(f for f in TRIS_BANDS if tris_key(f) in datasets)
    maps = [f for f in fit_files if f in datasets]
    data = copy.deepcopy(dict(config))
    data.setdefault("sed", {})["exclude_maps"] = [f for f in fit_files if f not in maps]
    return Config(data, config.path), tris, maps


def tris_key(f):
    return "TRIS_{}".format(int(f))


def tod_units(inp, errors, tris=TRIS_BANDS):
    """One SED per ring sample: as sed_prototype.tod_route builds it (sed_prototype.py:253-264),
    with the non-TRIS error bars taken from the CSV.  ``tris``: the TRIS bands to include
    (600, 820); the maps are ``inp.fit_maps`` (see :func:`select_datasets`)."""
    for k in range(inp.ra.size):
        pts = []
        for f, ring in inp.rings.items():
            if int(f) not in tris:
                continue
            nu = float(ring.effective_frequency_mhz)
            pts.append((tris_key(f), nu, ring.temperature_k[k] - sky_maps.cmb_rj_k(nu), ring_sigma(ring)[k]))
        for m in inp.fit_maps:
            d = inp.tods[m["name"]]
            if d["good"][k]:
                t_ex = float(sky_maps.rj_excess(d["tod"][k], m))
                pts.append((m["file"], float(m["freq_mhz"]), t_ex, sigma_from_row(errors[m["file"]], t_ex)))
        yield int(k), float(inp.ra[k]), pts


def map_units(inp, errors, add_zhat, tris_bands=TRIS_BANDS):
    """One SED per stripe pixel: as sed_prototype.map_route (sed_prototype.py:140-147, 202-213),
    with z-hat added to TRIS only if ``add_zhat``.  ``tris_bands``: the TRIS bands to include."""
    st = inp.stage1
    tris = {}
    for f in tris_bands:
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


# ---------------------------------------------------------------------------
# maximum likelihood (no Occam penalty)
# ---------------------------------------------------------------------------
def max_likelihood(pnames, model_fn, transform, nu, t, s, n_starts, rng):
    """Minimum chi^2 inside the prior box, searched in the unit cube of ``transform``
    (L-BFGS-B from the box centre and ``n_starts - 1`` random points).
    Returns (ln L_max, chi2_min, best-fit parameters)."""
    from scipy.optimize import minimize

    norm = -np.sum(np.log(s * np.sqrt(2 * np.pi)))

    def chi2(u):
        return float(np.sum(((t - model_fn(nu, transform(np.atleast_2d(u)))[0]) / s) ** 2))

    d = len(pnames)
    starts = np.vstack([np.full(d, 0.5), rng.uniform(0.02, 0.98, (max(n_starts, 1) - 1, d))])
    best = None
    for u0 in starts:
        r = minimize(chi2, u0, method="L-BFGS-B", bounds=[(0.0, 1.0)] * d)
        if best is None or r.fun < best.fun:
            best = r
    c2 = float(best.fun)
    return float(-0.5 * c2 + norm), c2, transform(np.atleast_2d(best.x))[0]


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


def write_maxlike(rows, models, out_dir, route, meta, tie=1.0):
    """Per-unit and summary CSVs for the maximum-likelihood comparison.  No sums over units:
    neighbouring units share data, so only per-unit numbers and counts are reported."""
    with open(out_dir / "maxlike_per_unit_{}.csv".format(route), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    ref = "power law"
    lnl = {n: np.array([r["lnLmax " + n] for r in rows]) for n in models}
    c2 = {n: np.array([r["chi2 " + n] for r in rows]) for n in models}
    dof = {n: np.array([r["dof " + n] for r in rows], dtype=float) for n in models}
    best = np.argmax(np.stack([lnl[n] for n in models]), axis=0)
    hdr = ["model", "parameters", "median_chi2", "mean_chi2", "median_chi2_per_dof", "median_dlnLmax_vs_power_law",
           "units_dlnLmax_gt_1", "units_tie", "units_dlnLmax_lt_-1", "units_best", "n_units", "route",
           "add_zhat_to_tris", "placeholder_points", "n_points", "moment_beta_fixed", "tod_beam", "datasets"]
    print("\nmaximum likelihood (no Occam penalty): dlnLmax = lnLmax(model) - lnLmax(power law) = -dchi2/2")
    print("{:<28}{:>8}{:>12}{:>12}{:>12}{:>12}{:>9}{:>7}{:>9}{:>7}".format(
        "model", "params", "med chi2", "mean chi2", "med chi2/nu", "med dlnL", "dlnL>1", "tie", "dlnL<-1", "best"))
    with open(out_dir / "maxlike_summary_{}.csv".format(route), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(hdr)
        for i, n in enumerate(models):
            d = lnl[n] - lnl[ref]
            k = int(round(np.median([r["n_points"] for r in rows]) - np.median(dof[n])))
            with np.errstate(divide="ignore", invalid="ignore"):
                red = np.median(c2[n] / dof[n])
            vals = [n, k, float(np.median(c2[n])), float(np.mean(c2[n])), float(red), float(np.median(d)), int((d > tie).sum()),
                    int((np.abs(d) <= tie).sum()), int((d < -tie).sum()), int((best == i).sum()), len(rows)]
            print("{:<28}{:>8d}{:>12.2f}{:>12.2f}{:>12.2f}{:>12.2f}{:>9d}{:>7d}{:>9d}{:>7d}".format(*vals[:-1]))
            w.writerow(vals + [route, meta["add_zhat_to_tris"], meta["placeholder_points"], meta["n_points"],
                               meta["moment_beta_fixed"], meta["tod_beam"], meta["datasets"]])


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", default=str(DEFAULT_MODEL_COMP_CONFIG_PATH),
                    help="the model-comparison config (default configs/model_comp_config.yaml)")
    ap.add_argument("--run-config", default=None,
                    help="the run config (default configs/run_pipeline.yaml): data sets, target, beam, ring")
    args = ap.parse_args()
    cfg = load_model_comp_config(args.config)
    logging.getLogger("ultranest").setLevel(logging.WARNING)

    route = cfg["data"]["route"]
    add_zhat = bool(cfg["zero_level"]["add_zhat_to_tris"])
    errors = read_error_bars(error_bars_path(cfg, args.config))
    models = build_models(cfg)
    ccfg = cfg.get("comparison") or {}
    stat = ccfg.get("statistic", "evidence")
    if stat not in ("max_likelihood", "evidence"):
        raise ValueError("comparison.statistic must be max_likelihood or evidence, got {!r}".format(stat))
    rng = np.random.default_rng(0)

    print("=" * 78)
    print("Spectral model comparison ({})".format(
        "maximum likelihood, no Occam penalty" if stat == "max_likelihood" else "UltraNest evidence"))
    print("  route:        {}".format(route))
    print("  zero level:   add z-hat to TRIS = {}{}".format(
        add_zhat, "  (no effect on the TOD route: ring data contain no z-hat)" if route == "tod" else ""))
    print("  models:       {}".format(", ".join(models)))
    bfix = cfg["models"].get("moment_beta_fixed")
    print("  moment exp.:  orders {}, {}".format(
        ", ".join(str(o) for o in moment_orders(cfg)),
        "beta free, omega1 = 0" if bfix is None else "beta fixed at {:g}, omega1 free".format(float(bfix))))
    print("  error bars:   {}".format(cfg["data"]["error_bars_csv"]))
    config = load_config(args.run_config)
    exclude_tris = bool(ccfg.get("exclude_tris", False))
    config, tris_bands, map_files = select_datasets(None, config, REPO / "res", exclude_tris=exclude_tris)
    used_keys = [tris_key(f) for f in tris_bands] + map_files
    print("  data sets:    {}  (maps left in by run_pipeline.yaml sed.exclude_maps; TRIS {})".format(
          ", ".join(used_keys) or "none", "excluded (comparison.exclude_tris)" if exclude_tris else "included"))
    tod_beam = sed_tod_beam(config)
    print("  ring TODs:    {} (run_pipeline.yaml sed.tod_beam){}".format(
        tod_beam, "" if route == "tod" else "  (not used by the map route)"))
    print("=" * 78)

    inp = Inputs(config, REPO / "res", target_file=str(config.require("pipeline.target_file")), verbose=True)
    missing = [m["file"] for m in inp.fit_maps if m["file"] not in errors]
    if missing:
        raise KeyError("no error-bar row in the CSV for: {}".format(missing))

    units = (tod_units(inp, errors, tris_bands) if route == "tod"
             else map_units(inp, errors, add_zhat, tris_bands))
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
            tr = prior_transform(pnames, cfg["priors"])
            if stat == "max_likelihood":
                lnl, c2, p = max_likelihood(pnames, fn, tr, nu, t, s, int(ccfg.get("ml_starts", 24)), rng)
                row["chi2 " + name], row["lnLmax " + name] = c2, lnl
                row["dof " + name] = len(pts) - len(pnames)
                for pn, pv in zip(pnames, p):
                    row["{} {}".format(pn, name)] = float(pv)
            else:
                lz, lze = log_evidence(pnames, fn, tr, nu, t, s, cfg["sampler"])
                row["logZ " + name], row["logZerr " + name] = lz, lze
        rows.append(row)
        col = "chi2 " if stat == "max_likelihood" else "logZ "
        print("  unit {:>5}: {} points, {} {} -- {:.0f} s".format(
            uid, len(pts), col.strip(), "  ".join("{:+.2f}".format(row[col + n]) for n in models),
            time.time() - t0))

    if not rows:
        raise RuntimeError("no unit had >= {} usable points".format(min_points))

    out_dir = REPO / "paper_outputs" / "model_comparison" / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    if stat == "max_likelihood":
        write_maxlike(rows, models, out_dir, route, dict(
            n_points=n_points, add_zhat_to_tris=add_zhat, placeholder_points=sum(used_ph.values()),
            moment_beta_fixed="free" if bfix is None else float(bfix), tod_beam=tod_beam if route == "tod" else "",
            datasets=";".join(used_keys)))
        if used_ph:
            placeholder_banner(used_ph, errors, n_points)
        print("{} units fitted ({} skipped with < {} points), {} data points, {:.0f} s".format(
            len(rows), skipped, min_points, n_points, time.time() - t0))
        print("wrote {}".format(out_dir))
        return
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
                    "route", "add_zhat_to_tris", "placeholder_points", "n_points", "moment_beta_fixed", "tod_beam",
                    "datasets"])
        for n in models:
            w.writerow([n, total[n], total_err[n], total[n] - total[ref], wins[n], len(rows), route, add_zhat,
                        sum(used_ph.values()), n_points, "free" if bfix is None else float(bfix),
                        tod_beam if route == "tod" else "", ";".join(used_keys)])
    if used_ph:
        placeholder_banner(used_ph, errors, n_points)
    print("wrote {}".format(out_dir))


if __name__ == "__main__":
    main()
