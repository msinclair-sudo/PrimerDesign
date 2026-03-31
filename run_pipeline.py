#!/usr/bin/env python3
"""
Primer Design Pipeline — Top-level entry point.

Usage:
    python run_pipeline.py -i alignment.fasta
    python run_pipeline.py -i alignment.fasta -d path/to/database
    python run_pipeline.py -p my_primers.csv -i alignment.fasta

Debug mode (run individual steps):
    python run_pipeline.py debug analyse -i alignment.fasta -d path/to/database
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

APP_DIR = Path(__file__).parent / "app"
CONFIG = Path(__file__).parent / "data" / "config.yaml"
DEFAULT_OUTPUT = Path("data/output")

STEPS = ["design", "validate", "expand", "analyse", "rank", "report"]


def run(cmd, desc, env=None):
    print(f"\n{'='*60}")
    print(f"  {desc}")
    print(f"{'='*60}\n")
    result = subprocess.run(cmd, text=True, env=env)
    if result.returncode != 0:
        print(f"\nERROR: {desc} failed (exit {result.returncode})")
        sys.exit(result.returncode)


def main():
    parser = argparse.ArgumentParser(description="eDNA Primer Design Pipeline")
    parser.add_argument("mode", nargs="?", default="run",
                        help="'debug' for individual step control (default: run all steps)")
    parser.add_argument("step", nargs="?", default=None,
                        help="Step to run (debug mode only): design, validate, expand, analyse, rank, report")
    parser.add_argument("-i", "--input", default="data/example/input/alignment.fasta",
                        help="Input MSA FASTA")
    parser.add_argument("-o", "--output", default=str(DEFAULT_OUTPUT),
                        help="Output directory (default: data/output)")
    parser.add_argument("-d", "--database", default=None,
                        help="FASTA directory or file for off-target screening")
    parser.add_argument("-g", "--genbank", default=None,
                        help="GenBank file (.gb) for gene annotation mapping")
    parser.add_argument("-p", "--primers", default=None,
                        help="CSV of primers to evaluate (skips design step)")
    parser.add_argument("-c", "--config", default=str(CONFIG),
                        help="Config YAML (default: primer_config.yaml)")
    parser.add_argument("--workers", type=int, default=None,
                        help="Number of parallel workers for Steps 2-4 (default: auto-detect)")
    args = parser.parse_args()

    # Determine run mode
    if args.mode == "debug":
        if args.step not in STEPS:
            parser.error(f"Debug mode requires a step: {', '.join(STEPS)}")
        run_steps = [args.step]
    else:
        # Normal mode — run all steps
        # If mode isn't "debug" or "run", it might be a flag misparse
        if args.mode not in ("run", "debug"):
            parser.error(f"Unknown mode '{args.mode}'. Use 'debug <step>' for individual steps.")
        run_steps = STEPS

    output_dir = Path(args.output)

    # Build environment with optional worker count
    run_env = os.environ.copy()
    if args.workers is not None:
        run_env["PRIMER_WORKERS"] = str(args.workers)

    # Auto-detect genbank if not specified: look for .gb files next to the input
    if args.genbank is None:
        input_dir = Path(args.input).parent
        gb_files = list(input_dir.glob("*.gb"))
        if gb_files:
            args.genbank = str(gb_files[0])
            print(f"Auto-detected GenBank file: {args.genbank}")

    # Ensure output dirs exist
    for sub in ["primers", "analysis"]:
        (output_dir / sub).mkdir(parents=True, exist_ok=True)

    primers_csv = output_dir / "primers" / "primers.csv"
    primers_filtered = output_dir / "primers" / "primers_filtered.csv"
    primers_expanded = output_dir / "primers" / "primers_expanded.csv"
    analysis_json = output_dir / "analysis" / "results.json"
    ranked_csv = output_dir / "analysis" / "primers_ranked.csv"
    report_html = output_dir / "report.html"

    def run_design():
        if args.primers:
            import shutil, csv
            shutil.copy2(args.primers, str(primers_csv))
            with open(primers_csv, newline="") as f:
                count = sum(1 for _ in csv.DictReader(f))
            print(f"\n{'='*60}")
            print(f"  Step 1: Using {count} supplied primer(s) from {args.primers}")
            print(f"{'='*60}\n")
        else:
            run(
                [sys.executable, str(APP_DIR / "design_primers.py"),
                 args.input, "-o", str(primers_csv), "--config", args.config],
                "Step 1: Design primers from MSA",
                env=run_env)

    steps = {
        "design": run_design,

        "validate": lambda: run(
            [sys.executable, str(APP_DIR / "validate_thermodynamics.py"),
             str(primers_csv), "-o", str(primers_filtered), "--config", args.config],
            "Step 2: Thermodynamic validation",
            env=run_env),

        "expand": lambda: run(
            [sys.executable, str(APP_DIR / "expand_combinations.py"),
             str(primers_filtered), "-i", args.input,
             "-o", str(primers_expanded), "--config", args.config],
            "Step 2b: Expand viable FWD+REV combinations",
            env=run_env),

        "analyse": lambda: run(
            [sys.executable, str(APP_DIR / "analyse_primers.py"),
             args.input, "--primers", str(primers_expanded),
             "-o", str(analysis_json), "--config", args.config]
            + (["--database", args.database] if args.database else [])
            + (["--genbank", args.genbank] if args.genbank else []),
            "Steps 3+4: Binding analysis + ecoPCR + annotations",
            env=run_env),

        "rank": lambda: run(
            [sys.executable, str(APP_DIR / "rank_primers.py"),
             str(analysis_json), "-o", str(ranked_csv), "--config", args.config],
            "Step 5: Rank primers",
            env=run_env),

        "report": lambda: run(
            [sys.executable, str(APP_DIR / "generate_report.py"),
             str(analysis_json), "-o", str(report_html), "--config", args.config],
            "Step 6: Generate HTML report",
            env=run_env),
    }

    for step_name in run_steps:
        steps[step_name]()

    if len(run_steps) == len(STEPS):
        print(f"\n{'='*60}")
        print(f"  Pipeline complete!")
        print(f"  Ranked primers: {ranked_csv}")
        print(f"  Report: {report_html}")
        print(f"{'='*60}")


if __name__ == "__main__":
    main()
