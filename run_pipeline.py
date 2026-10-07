#!/usr/bin/env python
"""Single entry point for the TRIS -> Haslam calibration pipeline.

    python run_pipeline.py --stage 1

Every stage reads ``configs/run_config.yaml``, writes under ``outputs/<stage>/``
and records the config block it used in a JSON manifest next to its products.
Stages 3-5 are not implemented yet and say so rather than half-running.

``--nside`` and ``--run-name`` override the config for one run without editing
it, which is what comparing grids or keeping several stage 1 products side by
side needs.  Whatever they are set to is what the manifest records, so a
product still says exactly what produced it.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional, Sequence

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from moonrunes.config import load_config  # noqa: E402
from moonrunes import stage1_tris_maps, stage2_beam_matching  # noqa: E402

_UNIMPLEMENTED = {
    3: "stage 3 -- assemble per-pixel SEDs and fit",
    4: "stage 4 -- derive and apply the calibration correction",
    5: "stage 5 -- validate",
}


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(prog="run_pipeline")
    parser.add_argument(
        "--stage", type=int, required=True, choices=sorted({1, 2, *_UNIMPLEMENTED})
    )
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--archive-dir", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument(
        "--nside",
        type=int,
        default=None,
        help="override tris.nside for this run (stage 1). Must be a power of 2.",
    )
    parser.add_argument(
        "--run-name",
        type=str,
        default=None,
        help="override run.name, which names the product file "
        "(tris_maps_<run name>.npz).",
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)

    if args.stage in _UNIMPLEMENTED:
        raise SystemExit("{} is not implemented yet".format(_UNIMPLEMENTED[args.stage]))

    config = load_config(args.config)

    if args.nside is not None:
        if args.nside < 1 or args.nside & (args.nside - 1):
            raise SystemExit(
                "--nside must be a power of 2 (HEALPix), got {}".format(args.nside)
            )
        config["tris"]["nside"] = args.nside
        # nside_hires is the grid the beam is rotated on before being sampled
        # and must not fall below the working grid; lift it if the override
        # would invert them.
        if int(config["tris"]["nside_hires"]) < args.nside:
            config["tris"]["nside_hires"] = args.nside

    if args.run_name is not None:
        config["run"]["name"] = args.run_name

    if args.stage == 1:
        stage1_tris_maps.run_stage1(
            config=config,
            archive_dir=args.archive_dir,
            output_dir=args.output_dir,
            overwrite=True if args.overwrite else None,
            verbose=not args.quiet,
        )
    else:
        if args.nside is not None or args.run_name is not None:
            raise SystemExit(
                "--nside and --run-name apply to stage 1; stage 2 takes its grid "
                "from beam_matching in the config"
            )
        stage2_beam_matching.run_stage2(
            config=config,
            output_dir=args.output_dir,
            overwrite=True if args.overwrite else None,
            verbose=not args.quiet,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
