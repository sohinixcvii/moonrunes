"""The spectral models (power law, curved power law, moment expansion) and their priors.

Every model is a CMB-removed Rayleigh-Jeans temperature T(nu), with x = ln(nu / nu0):

* power law:          T = A exp(beta x)
* curved power law:   T = A exp(beta x + C x^2)
* moment expansion:   T = A exp(beta x) [1 + omega1 x + omega2 x^2 / 2 (+ omega3 x^3 / 6)]
  (Chluba et al. 2017), one model per order in ``models.moment_orders`` (2: up to omega2;
  3: up to omega3).  With ``models.moment_beta_fixed`` set, beta is fixed at that value
  and omega1 is free; with it null, beta is free and omega1 = 0.  Both orders use the same
  beta setting.

Each model is (parameter names, prior transform, model function).  The prior transform
maps UltraNest's unit cube onto the parameters using the config's ``priors`` block.
"""

from __future__ import annotations

from typing import Callable, Dict, List, Tuple

import numpy as np


def _transform_1d(spec):
    lo, hi = float(spec["low"]), float(spec["high"])
    if spec["type"] == "uniform":
        return lambda u: lo + u * (hi - lo)
    if spec["type"] == "loguniform":
        llo, lhi = np.log(lo), np.log(hi)
        return lambda u: np.exp(llo + u * (lhi - llo))
    raise ValueError("unknown prior type {!r}".format(spec["type"]))


def prior_transform(names: List[str], priors: Dict[str, dict]) -> Callable:
    parts = [_transform_1d(priors[n]) for n in names]

    def transform(cube):
        cube = np.asarray(cube, dtype=float)
        out = np.empty_like(cube)
        for i, f in enumerate(parts):
            out[..., i] = f(cube[..., i])
        return out

    return transform


def moment_orders(cfg) -> List[int]:
    """The moment-expansion orders to compare: ``models.moment_orders`` (a list), or the
    older single ``models.moment_order``."""
    m = cfg["models"]
    orders = m["moment_orders"] if "moment_orders" in m else m["moment_order"]
    orders = [int(orders)] if np.isscalar(orders) else [int(o) for o in orders]
    if not orders or len(set(orders)) != len(orders) or any(o not in (2, 3) for o in orders):
        raise ValueError("moment_orders must be a list of distinct values from (2, 3), got {!r}".format(orders))
    return orders


def build_models(cfg) -> Dict[str, Tuple[List[str], Callable]]:
    """{model name: (parameter names, model(nu, params) -> T)}; params is (..., npar).
    One moment-expansion model per order in :func:`moment_orders`, named
    ``"moment expansion (order N)"``."""
    nu0 = float(cfg["models"]["nu0_mhz"])
    beta_fixed = cfg["models"].get("moment_beta_fixed")

    def x_of(nu):
        return np.log(np.asarray(nu, dtype=float) / nu0)

    def power_law(nu, p):
        x = x_of(nu)
        return p[..., 0:1] * np.exp(p[..., 1:2] * x)

    def curved(nu, p):
        x = x_of(nu)
        return p[..., 0:1] * np.exp(p[..., 1:2] * x + p[..., 2:3] * x ** 2)

    def make_moment(order):
        def moment(nu, p):
            # p = (A, beta, omega2[, omega3]) if beta is free, (A, omega1, omega2[, omega3]) if fixed
            x = x_of(nu)
            poly = 1.0 + 0.5 * p[..., 2:3] * x ** 2
            if order == 3:
                poly = poly + p[..., 3:4] * x ** 3 / 6.0
            if beta_fixed is None:
                return p[..., 0:1] * np.exp(p[..., 1:2] * x) * poly
            return p[..., 0:1] * np.exp(float(beta_fixed) * x) * (poly + p[..., 1:2] * x)
        return moment

    models = {
        "power law": (["A_K", "beta"], power_law),
        "curved power law": (["A_K", "beta", "C"], curved),
    }
    for order in moment_orders(cfg):
        names = (["A_K", "beta" if beta_fixed is None else "omega1", "omega2"]
                 + (["omega3"] if order == 3 else []))
        models["moment expansion (order {})".format(order)] = (names, make_moment(order))
    return models
