#!/usr/bin/env python
"""Single entry point for the TRIS -> Haslam calibration pipeline.

    python run_pipeline.py --stage 1

Every stage reads ``configs/run_config.yaml``, writes under ``outputs/<stage>/``
and records the config block it used in a JSON manifest next to its products.
Stages 2-5 are not implemented yet and say so rather than half-running.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional, Sequence

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from moonrunes import stage1_tris_maps  # noqa: E402

_UNIMPLEMENTED = {
    2: "stage 2 -- Haslam prep, gain operator and the Berkhuijsen prior",
    3: "stage 3 -- assemble the bayesian_skymap .npz",
    4: "stage 4 -- run the Gibbs sampler",
    5: "stage 5 -- validate",
}


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(prog="run_pipeline")
    parser.add_argument(
        "--stage", type=int, required=True, choices=sorted({1, *_UNIMPLEMENTED})
    )
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--archive-dir", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)

    if args.stage != 1:
        raise SystemExit("{} is not implemented yet".format(_UNIMPLEMENTED[args.stage]))

    stage1_tris_maps.run_stage1(
        config_path=args.config,
        archive_dir=args.archive_dir,
        output_dir=args.output_dir,
        overwrite=True if args.overwrite else None,
        verbose=not args.quiet,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
