#!/usr/bin/env python3
"""
Step 2b: Expand all viable FWD+REV combinations from validated primers.

Reads primers_validated.csv (individual primers from validate_thermodynamics.py),
generates all thermodynamically viable FWD+REV combinations filtered by delta Tm
and heterodimer stability. Actual binding and amplicon viability are determined
by obipcr in the analysis step.

Usage:
    python scripts/expand_combinations.py primers_validated.csv -o primers_expanded.csv
    python scripts/expand_combinations.py primers_validated.csv -o primers_expanded.csv --config primer_config.yaml
"""

import argparse
import csv
import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    yaml = None

# Import heterodimer calculation from sibling script
sys.path.insert(0, str(Path(__file__).parent))
from validate_thermodynamics import calc_heterodimer, _init_backend


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
def _get_value(entry):
    """Extract value from a config entry (supports both {value:, action:} and bare value)."""
    if isinstance(entry, dict):
        return entry["value"]
    return entry


def load_config(config_path):
    """Load config for thresholds. Exits if config is missing."""
    if config_path is None:
        print("ERROR: --config is required. No config file specified.", file=sys.stderr)
        sys.exit(1)
    p = Path(config_path)
    if not p.exists():
        print(f"ERROR: Config not found at {p}", file=sys.stderr)
        sys.exit(1)
    if yaml is None:
        print("ERROR: PyYAML is required but not installed.", file=sys.stderr)
        sys.exit(1)
    with open(p) as f:
        raw = yaml.safe_load(f)
    thermo = raw.get("thermodynamics")
    if not thermo:
        print(f"ERROR: Config {p} has no 'thermodynamics' section.", file=sys.stderr)
        sys.exit(1)
    required = ["max_delta_tm", "max_heterodimer_dg"]
    missing = [k for k in required if k not in thermo]
    if missing:
        print(f"ERROR: Config missing required keys: {', '.join(missing)}", file=sys.stderr)
        sys.exit(1)
    return {k: _get_value(thermo[k]) for k in required}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Expand all viable FWD+REV combinations from validated primers."
    )
    parser.add_argument("input", help="Validated primers CSV (from validate_thermodynamics.py)")
    parser.add_argument("-i", "--msa", default=None,
                        help="(unused, kept for pipeline compatibility)")
    parser.add_argument("-o", "--output", default="primers_expanded.csv", help="Output CSV")
    parser.add_argument("--config", required=True, help="Path to config YAML (required)")
    args = parser.parse_args()

    cfg = load_config(args.config)
    max_delta_tm = cfg["max_delta_tm"]
    max_heterodimer_dg = cfg["max_heterodimer_dg"]

    # Initialise thermodynamic backend for heterodimer calculations
    _init_backend()

    # Read validated individual primers
    fwd_list = []  # [{seq, tm, gc, len}, ...]
    rev_list = []

    with open(args.input, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("status") == "REJECT":
                continue

            seq = row["sequence"].upper()
            direction = row["direction"]
            entry = {
                "seq": seq,
                "tm": float(row.get("tm", 0)),
                "gc": float(row.get("gc_pct", 0)),
                "len": len(seq),
            }

            if direction == "FWD":
                fwd_list.append(entry)
            elif direction == "REV":
                rev_list.append(entry)

    print(f"Validated FWD primers: {len(fwd_list)}")
    print(f"Validated REV primers: {len(rev_list)}")
    print(f"Total possible combinations: {len(fwd_list) * len(rev_list)}")

    # Generate all viable combinations (delta Tm + heterodimer only)
    # Actual amplicon viability is determined by obipcr in the analysis step
    viable = []
    n_tm_rejected = 0
    n_heterodimer_rejected = 0

    for fi, fwd in enumerate(fwd_list):
        for ri, rev in enumerate(rev_list):
            # Delta Tm filter
            delta_tm = abs(fwd["tm"] - rev["tm"])
            if delta_tm > max_delta_tm:
                n_tm_rejected += 1
                continue

            # Heterodimer filter
            hetero_dg = calc_heterodimer(fwd["seq"], rev["seq"])
            if hetero_dg < max_heterodimer_dg:
                n_heterodimer_rejected += 1
                continue

            name = f"PAIR_F{fi+1:03d}_R{ri+1:03d}"
            viable.append({
                "name": name,
                "forward": fwd["seq"],
                "reverse": rev["seq"],
            })

    print(f"Rejected: {n_tm_rejected} (delta Tm), {n_heterodimer_rejected} (heterodimer)")
    print(f"Viable combinations: {len(viable)}")

    # Write output
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["name", "forward", "reverse"])
        writer.writeheader()
        writer.writerows(viable)

    print(f"Wrote {len(viable)} primer pairs to {args.output}")


if __name__ == "__main__":
    main()
