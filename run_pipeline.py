#!/usr/bin/env python
"""Single entry point for the TRIS -> Haslam pipeline.

    python run_pipeline.py --step all                 # steps 0-7, then the report
    python run_pipeline.py --step fit                 # one step (it reads the products of the steps before it)
    python run_pipeline.py --from-step compare        # this step and every step after it
    python run_pipeline.py --step report              # PIPELINE_REPORT.md from the steps already run

Steps (``moonrunes.pipeline.STEPS``; ``PIPELINE.md`` describes each one):

    0 tris_maps    TRIS maps from the archive rings (stage 1)          -> outputs/stage1/
    1 beam_match   every allowed map beam-matched to the TRIS beam      -> paper_outputs/pipeline/step1_*/
    2 ring_tods    ring TODs of every data set and comparison model     -> .../step2_ring_tods/
    3 compare      model comparison, TRIS excluded (UltraNest)          -> .../step3_model_comparison/
    4 select       the best model                                       -> .../step4_model_selection/
    5 fit          the best model with TRIS: our spectral model         -> .../step5_fit/
    6 predict_408  our model at 408 MHz against Haslam                  -> .../step6_predict_408/
    7 rhino        our model at RHINO frequencies, through both beams   -> .../step7_rhino/

Everything is set in ``configs/run_pipeline.yaml`` (the default ``--config``; data, beam, ring, ``pipeline:``) and
``configs/model_comp_config.yaml`` (models, priors, sampler, error bars, the comparison statistic).

``--stage 1`` and ``--stage 2`` are the older interface: stage 1 is step 0 (``--nside`` and
``--run-name`` apply to it), stage 2 is the older beam matching of Haslam and ARCADE 2 alone
(``moonrunes.stage2_beam_matching``), which nothing downstream reads.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from typing import Optional, Sequence

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from moonrunes.config import load_config  # noqa: E402
from moonrunes import pipeline, stage1_tris_maps, stage2_beam_matching  # noqa: E402


def main(argv: Optional[Sequence[str]] = None) -> int:
    t0 = time.time()
    try:
        return _main(argv)
    finally:
        dt = time.time() - t0
        print("total run time: {:.0f} s ({}h {:02d}m {:02d}s)".format(
            dt, int(dt // 3600), int(dt % 3600 // 60), int(dt % 60)))


def _main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(prog="run_pipeline", description=__doc__.splitlines()[0])
    what = parser.add_mutually_exclusive_group(required=True)
    what.add_argument("--step", choices=["all", "report", *pipeline.STEPS],
                      help="run one step, all of them, or only rewrite the report")
    what.add_argument("--from-step", choices=pipeline.STEPS, help="run this step and every step after it")
    what.add_argument("--stage", type=int, choices=[1, 2, 3, 4, 5],
                      help="the older interface: 1 = step tris_maps, 2 = the older Haslam/ARCADE 2 beam matching")
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--archive-dir", type=Path, default=None, help="stage 1 / step tris_maps only")
    parser.add_argument("--output-dir", type=Path, default=None, help="--stage 1 and 2 only")
    parser.add_argument("--nside", type=int, default=None,
                        help="override tris.nside for this run (stage 1 / step tris_maps). Must be a power of 2.")
    parser.add_argument("--run-name", type=str, default=None,
                        help="override run.name, which names the stage 1 product (tris_maps_<run name>.npz)")
    parser.add_argument("--overwrite", action="store_true",
                        help="redo products that already exist (the fit cache of step 5 included)")
    parser.add_argument("--no-report", action="store_true", help="do not rewrite PIPELINE_REPORT.md at the end")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)

    config = load_config(args.config)

    if args.nside is not None:
        if args.nside < 1 or args.nside & (args.nside - 1):
            raise SystemExit("--nside must be a power of 2 (HEALPix), got {}".format(args.nside))
        config["tris"]["nside"] = args.nside
        # nside_hires is the grid the beam is rotated on before being sampled and must not fall below the
        # working grid; lift it if the override would invert them.
        if int(config["tris"]["nside_hires"]) < args.nside:
            config["tris"]["nside_hires"] = args.nside
    if args.run_name is not None:
        config["run"]["name"] = args.run_name
    stage1_only = args.stage == 1 or args.step == "tris_maps"
    if (args.nside is not None or args.run_name is not None) and not stage1_only:
        raise SystemExit("--nside and --run-name apply to stage 1 (step tris_maps) only: the later steps read "
                         "stage 1's product as configured, so a one-off grid there would not match them")

    if args.stage is not None:
        if args.stage >= 3:
            raise SystemExit("stages 3-5 are pipeline steps: use --step (see --help)")
        if args.stage == 1:
            stage1_tris_maps.run_stage1(config=config, archive_dir=args.archive_dir, output_dir=args.output_dir,
                                        overwrite=True if args.overwrite else None, verbose=not args.quiet)
        else:
            stage2_beam_matching.run_stage2(config=config, output_dir=args.output_dir,
                                            overwrite=True if args.overwrite else None, verbose=not args.quiet)
        return 0

    if args.archive_dir is not None and not stage1_only:
        raise SystemExit("--archive-dir applies to stage 1 (step tris_maps) only")
    if args.output_dir is not None:
        raise SystemExit("--output-dir applies to --stage 1 and 2; pipeline steps write to pipeline.output_dir")

    if args.step == "report":
        print(pipeline.write_report(config))
        return 0
    if args.step == "tris_maps" and args.archive_dir is not None:
        stage1_tris_maps.run_stage1(config=config, archive_dir=args.archive_dir,
                                    overwrite=True if args.overwrite else None, verbose=not args.quiet)
        return 0
    if args.step == "all":
        steps = list(pipeline.STEPS)
    elif args.from_step is not None:
        steps = list(pipeline.STEPS[pipeline.STEPS.index(args.from_step):])
    else:
        steps = [args.step]
    pipeline.run(steps, config=config, verbose=not args.quiet, overwrite=args.overwrite, report=not args.no_report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
