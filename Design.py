#!/usr/bin/env python3
"""
eDNA Primer Design Pipeline — Top-level entry point.

Usage:
    # Full primer design pipeline
    python Design.py -i alignment.fasta -g reference.gb

    # With off-target screening
    python Design.py -i alignment.fasta -g reference.gb -d path/to/database

    # Evaluate existing primers (skips design step)
    python Design.py -p my_primers.csv -i alignment.fasta -g reference.gb

    # Prepare alignment first (circular genome rotate + trim)
    python Design.py align -i raw_sequences.fasta -g reference.gb --from cytb --to 16S

    # Debug mode (run individual pipeline step)
    # Steps: design, validate, expand, analyse, report
    python Design.py debug <step> -i alignment.fasta
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

APP_DIR = Path(__file__).parent / "app"
CONFIG = Path(__file__).parent / "data" / "config.yaml"
DEFAULT_OUTPUT = Path("data/output")

STEPS = ["design", "validate", "expand", "analyse", "report"]


def run(cmd, desc, env=None):
    print(f"\n{'='*60}")
    print(f"  {desc}")
    print(f"{'='*60}\n")
    result = subprocess.run(cmd, text=True, env=env)
    if result.returncode != 0:
        print(f"\nERROR: {desc} failed (exit {result.returncode})")
        sys.exit(result.returncode)


# ---------------------------------------------------------------------------
# Subcommand: align
# ---------------------------------------------------------------------------
def cmd_align(args):
    """Run alignment preparation (rotate + realign + trim for circular genomes)."""
    cmd = [
        sys.executable, str(APP_DIR / "prepare_alignment.py"),
        args.input,
        "-g", args.genbank,
        "--from", args.from_gene,
        "--to", args.to_gene,
        "-o", args.output,
    ]
    run(cmd, "Alignment preparation (rotate + realign + trim)")


# ---------------------------------------------------------------------------
# Subcommand: run (default) / debug
# ---------------------------------------------------------------------------
def cmd_pipeline(args):
    """Run the primer design pipeline (all steps or a single debug step)."""
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
    primers_validated = output_dir / "primers" / "primers_validated.csv"
    primers_expanded = output_dir / "primers" / "primers_expanded.csv"
    analysis_json = output_dir / "analysis" / "results.json"
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
             str(primers_csv), "-o", str(primers_validated), "--config", args.config],
            "Step 2: Individual primer validation",
            env=run_env),

        "expand": lambda: run(
            [sys.executable, str(APP_DIR / "expand_combinations.py"),
             str(primers_validated), "-i", args.input,
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

        "report": lambda: run(
            [sys.executable, str(APP_DIR / "generate_report.py"),
             str(analysis_json), "-o", str(report_html), "--config", args.config],
            "Step 5: Generate HTML report",
            env=run_env),
    }

    for step_name in args.run_steps:
        steps[step_name]()

    if len(args.run_steps) == len(STEPS):
        print(f"\n{'='*60}")
        print(f"  Pipeline complete!")
        print(f"  Report: {report_html}")
        print(f"{'='*60}")


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="eDNA Primer Design Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command")

    # --- align subcommand ---------------------------------------------------
    align_parser = subparsers.add_parser(
        "align",
        help="Prepare alignment: rotate circular genomes, realign, and trim to target region",
    )
    align_parser.add_argument("-i", "--input", required=True,
                              help="Input FASTA (unaligned or aligned sequences)")
    align_parser.add_argument("-g", "--genbank", required=True,
                              help="GenBank file (.gb) for the reference sequence")
    align_parser.add_argument("--from", dest="from_gene", required=True,
                              help="Start gene name (e.g. cytb, 12S, Dloop)")
    align_parser.add_argument("--to", dest="to_gene", required=True,
                              help="End gene name (e.g. 16S, cytb)")
    align_parser.add_argument("-o", "--output", default="aligned_trimmed.fasta",
                              help="Output FASTA (default: aligned_trimmed.fasta)")
    align_parser.set_defaults(func=cmd_align)

    # --- debug subcommand ---------------------------------------------------
    debug_parser = subparsers.add_parser(
        "debug",
        help="Run a single pipeline step",
    )
    debug_parser.add_argument("step", choices=STEPS,
                              help="Step to run")
    debug_parser.add_argument("-i", "--input", default="data/example/input/alignment.fasta",
                              help="Input MSA FASTA")
    debug_parser.add_argument("-o", "--output", default=str(DEFAULT_OUTPUT),
                              help="Output directory (default: data/output)")
    debug_parser.add_argument("-d", "--database", default=None,
                              help="FASTA directory or file for off-target screening")
    debug_parser.add_argument("-g", "--genbank", default=None,
                              help="GenBank file (.gb) for gene annotation mapping")
    debug_parser.add_argument("-p", "--primers", default=None,
                              help="CSV of primers to evaluate (skips design step)")
    debug_parser.add_argument("-c", "--config", default=str(CONFIG),
                              help="Config YAML")
    debug_parser.add_argument("--workers", type=int, default=None,
                              help="Number of parallel workers")
    debug_parser.set_defaults(func=cmd_pipeline)

    # --- default: run full pipeline (no subcommand) -------------------------
    # argparse doesn't natively support optional subcommands, so we handle
    # the "no subcommand" case by checking dest="command"
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
                        help="Config YAML")
    parser.add_argument("--workers", type=int, default=None,
                        help="Number of parallel workers")

    args = parser.parse_args()

    if args.command is None:
        # No subcommand — run full pipeline
        args.run_steps = STEPS
        args.func = cmd_pipeline
    elif args.command == "debug":
        args.run_steps = [args.step]

    args.func(args)


if __name__ == "__main__":
    main()
