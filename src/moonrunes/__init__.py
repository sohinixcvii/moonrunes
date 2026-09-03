"""TRIS map-making (limTOD) feeding a Gibbs recalibration of Haslam.

The stages are deliberately importable one at a time: each one reads
``configs/run_config.yaml``, writes its products under ``outputs/<stage>/`` and
drops a JSON manifest recording the exact config block it used.  Nothing
downstream may re-derive a decision the config already records.
"""

__all__ = ["config"]
