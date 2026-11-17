"""Loaders for the two configs in ``configs/``:

* ``configs/run_pipeline.yaml`` (:func:`load_config`): the data, beam, ring, corrections, outputs and the
  pipeline's steps -- the run config;
* ``configs/model_comp_config.yaml`` (:func:`load_model_comp_config`): the spectral models, priors, sampler,
  error bars and what the model comparison is decided on.  The run config points to it
  (``pipeline.model_comparison_config``); no setting is in both.

Two rules the stages rely on:

* a value the config leaves ``null`` is an *open decision*, and the stage that
  needs it raises rather than defaulting -- :func:`require` is that raise;
* every stage records the config block it consumed in its manifest, so a
  product can always be traced back to the decisions that produced it.
"""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any, Mapping, Optional

import yaml

# src/moonrunes/config.py -> src/moonrunes -> src -> repo root
REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG_PATH = REPO_ROOT / "configs" / "run_pipeline.yaml"
DEFAULT_MODEL_COMP_CONFIG_PATH = REPO_ROOT / "configs" / "model_comp_config.yaml"


class ConfigError(ValueError):
    """A config value is missing, still open, or of the wrong shape."""


class Config(Mapping):
    """Read-only view of the run config with dotted lookup."""

    def __init__(self, data: Mapping[str, Any], path: Optional[Path] = None) -> None:
        self._data = copy.deepcopy(dict(data))
        self.path = Path(path) if path is not None else None

    # -- Mapping protocol ---------------------------------------------------
    def __getitem__(self, key: str) -> Any:
        return self._data[key]

    def __iter__(self):
        return iter(self._data)

    def __len__(self) -> int:
        return len(self._data)

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return "Config(path={!r}, sections={})".format(
            str(self.path), sorted(self._data)
        )

    # -- lookup -------------------------------------------------------------
    def get_path(self, dotted: str, default: Any = None) -> Any:
        """``config.get_path("tris.prior.relative_sigma")``, ``default`` if absent."""
        node: Any = self._data
        for part in dotted.split("."):
            if not isinstance(node, Mapping) or part not in node:
                return default
            node = node[part]
        return copy.deepcopy(node)

    def require(self, dotted: str) -> Any:
        """Like :meth:`get_path`, but raise if the key is absent or ``None``.

        ``None`` in this config means "still an open decision", so a stage that
        needs the value must stop rather than invent one.
        """
        sentinel = object()
        value = self.get_path(dotted, sentinel)
        if value is sentinel:
            raise ConfigError(
                "{} is not in {}".format(dotted, self.path or "the run config")
            )
        if value is None:
            raise ConfigError(
                "{} is still an open decision (null) in {} -- set it before "
                "running this stage".format(dotted, self.path or "the run config")
            )
        return value

    def section(self, name: str) -> dict:
        """Return a deep copy of one top-level block, for manifest recording."""
        value = self.get_path(name)
        if not isinstance(value, dict):
            raise ConfigError("{!r} is not a config section".format(name))
        return value

    def resolve_path(self, dotted: str) -> Path:
        """Resolve a config path value against the repo root."""
        raw = Path(str(self.require(dotted)))
        return raw if raw.is_absolute() else (REPO_ROOT / raw).resolve()


def load_config(path: Optional[Path] = None) -> Config:
    """Load the run config (``configs/run_pipeline.yaml`` by default)."""
    config_path = Path(path) if path is not None else DEFAULT_CONFIG_PATH
    if not config_path.exists():
        raise ConfigError("no run config at {}".format(config_path))
    with open(config_path, "r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict):
        raise ConfigError("{} did not parse to a mapping".format(config_path))
    return Config(data, path=config_path)


def model_comp_config_path(config: Optional[Config] = None) -> Path:
    """The model-comparison config: ``pipeline.model_comparison_config`` of ``config`` if given, else
    ``configs/model_comp_config.yaml``."""
    if config is not None and config.get_path("pipeline.model_comparison_config") is not None:
        return config.resolve_path("pipeline.model_comparison_config")
    return DEFAULT_MODEL_COMP_CONFIG_PATH


def load_model_comp_config(path: Optional[Path] = None, config: Optional[Config] = None) -> dict:
    """The model-comparison config as a plain dict (``path``, else :func:`model_comp_config_path`)."""
    p = Path(path) if path is not None else model_comp_config_path(config)
    if not p.exists():
        raise ConfigError("no model-comparison config at {}".format(p))
    with open(p, "r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict):
        raise ConfigError("{} did not parse to a mapping".format(p))
    return data


def error_bars_path(mc: Mapping[str, Any], mc_path: Optional[Path] = None) -> Path:
    """``data.error_bars_csv`` of a model-comparison config, relative to that config's folder (``configs/``)."""
    folder = Path(mc_path).parent if mc_path is not None else DEFAULT_MODEL_COMP_CONFIG_PATH.parent
    raw = Path(str(mc["data"]["error_bars_csv"]))
    return raw if raw.is_absolute() else (folder / raw).resolve()
