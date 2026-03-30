#!/usr/bin/env python3
"""
Step 2b: Expand all viable FWD+REV combinations from filtered primers.

Reads primers_filtered.csv (output of validate_thermodynamics.py), extracts
unique forward and reverse sequences, and generates all thermodynamically
viable combinations. This ensures analyse_primers.py has data for every
FWD+REV pair the user might select in the interactive report.

Usage:
    python scripts/expand_combinations.py primers_filtered.csv -o primers_expanded.csv
    python scripts/expand_combinations.py primers_filtered.csv -o primers_expanded.csv --config primer_config.yaml
"""

import argparse
import csv
import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    yaml = None


def load_config(config_path):
    """Load config for thresholds."""
    defaults = {
        "max_delta_tm": 5.0,
        "min_amplicon_length": 90,
        "max_amplicon_length": 180,
    }
    if config_path is None or yaml is None:
        return defaults
    p = Path(config_path)
    if not p.exists():
        return defaults
    with open(p) as f:
        raw = yaml.safe_load(f)
    thermo = raw.get("thermodynamics", {})
    design = raw.get("design", {})
    return {
        "max_delta_tm": thermo.get("max_delta_tm", defaults["max_delta_tm"]),
        "min_amplicon_length": design.get("min_amplicon_length", defaults["min_amplicon_length"]),
        "max_amplicon_length": design.get("max_amplicon_length", defaults["max_amplicon_length"]),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Expand all viable FWD+REV combinations from filtered primers."
    )
    parser.add_argument("input", help="Filtered primers CSV (from validate_thermodynamics.py)")
    parser.add_argument("-o", "--output", default="primers_expanded.csv", help="Output CSV")
    parser.add_argument("--config", default=None, help="Path to primer_config.yaml")
    args = parser.parse_args()

    cfg = load_config(args.config)
    max_delta_tm = cfg["max_delta_tm"]

    # Read filtered primers, collect unique FWDs and REVs with Tm
    fwd_map = {}  # seq -> {tm, gc, len, ...}
    rev_map = {}
    original_pairs = set()

    with open(args.input, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Skip rejected primers
            if row.get("status") == "REJECT":
                continue

            fwd_seq = row["forward"].upper()
            rev_seq = row["reverse"].upper()
            original_pairs.add((fwd_seq, rev_seq))

            if fwd_seq not in fwd_map:
                fwd_map[fwd_seq] = {
                    "seq": fwd_seq,
                    "tm": float(row.get("tm_fwd", 0)),
                    "gc": float(row.get("gc_fwd", 0)),
                    "len": len(fwd_seq),
                }
            if rev_seq not in rev_map:
                rev_map[rev_seq] = {
                    "seq": rev_seq,
                    "tm": float(row.get("tm_rev", 0)),
                    "gc": float(row.get("gc_rev", 0)),
                    "len": len(rev_seq),
                }

    fwds = list(fwd_map.values())
    revs = list(rev_map.values())

    print(f"Unique FWD primers: {len(fwds)}")
    print(f"Unique REV primers: {len(revs)}")
    print(f"Original pairs: {len(original_pairs)}")
    print(f"Total possible combinations: {len(fwds) * len(revs)}")

    # Generate all viable combinations
    viable = []
    for fi, fwd in enumerate(fwds):
        for ri, rev in enumerate(revs):
            delta_tm = abs(fwd["tm"] - rev["tm"])
            if delta_tm > max_delta_tm:
                continue
            name_suffix = f"F{fi+1:02d}_R{ri+1:02d}"
            is_original = (fwd["seq"], rev["seq"]) in original_pairs
            viable.append({
                "name": f"COMBO_{name_suffix}" if not is_original else f"ORIG_{name_suffix}",
                "forward": fwd["seq"],
                "reverse": rev["seq"],
            })

    print(f"Viable combinations (delta Tm <= {max_delta_tm}°C): {len(viable)}")
    new_count = len(viable) - sum(1 for v in viable if v["name"].startswith("ORIG_"))
    print(f"  Original: {len(viable) - new_count}, New: {new_count}")

    # Write output
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["name", "forward", "reverse"])
        writer.writeheader()
        writer.writerows(viable)

    print(f"Wrote {len(viable)} primer pairs to {args.output}")


if __name__ == "__main__":
    main()
