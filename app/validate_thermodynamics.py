#!/usr/bin/env python3
"""Step 2: Thermodynamic validation of individual primers.

Extracts unique forward and reverse primer sequences from a pairs CSV,
calculates Tm, hairpin/homodimer ΔG, GC%, GC clamp, and homopolymer runs
for each unique primer, then applies hard-reject and soft-flag logic.

Pair-level checks (delta Tm, heterodimer) are deferred to the expansion
step where all viable FWD+REV combinations are generated.

Requires primer3-py (preferred) or ntthal CLI as fallback.
"""

import argparse
import csv
import multiprocessing
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

try:
    import yaml
except ImportError:
    yaml = None

# ---------------------------------------------------------------------------
# Required config keys (thermodynamics section)
REQUIRED_KEYS = ["tm_range", "gc_range", "max_hairpin_dg", "max_homodimer_dg",
                 "max_homopolymer", "gc_clamp_3prime"]

# Map from config key to the flag name used in validation results
CONFIG_KEY_TO_FLAG = {
    "tm_range": "tm_out_of_range",
    "gc_range": "gc_out_of_range",
    "max_hairpin_dg": "hairpin",
    "max_homodimer_dg": "homodimer",
    "max_homopolymer": "homopolymer",
    "gc_clamp_3prime": "gc_clamp",
}


# ---------------------------------------------------------------------------
# Config loading
# ---------------------------------------------------------------------------
def _parse_thermo_entry(key, entry):
    """Parse a thermodynamics config entry. Supports both formats:
    - New: {value: ..., action: reject|flag}
    - Legacy: bare value (defaults to reject)
    """
    if isinstance(entry, dict):
        if "value" not in entry:
            print(f"ERROR: thermodynamics.{key} is missing 'value'", file=sys.stderr)
            sys.exit(1)
        action = entry.get("action", "reject").lower()
        if action not in ("reject", "flag"):
            print(f"ERROR: thermodynamics.{key}.action must be 'reject' or 'flag', got '{action}'",
                  file=sys.stderr)
            sys.exit(1)
        return entry["value"], action
    # Legacy bare value — default to reject
    return entry, "reject"


def load_config(path):
    """Load thermodynamic thresholds from YAML config. Exits if config is missing.
    Returns (values_dict, hard_reject_set) where hard_reject_set contains flag names
    for checks configured as 'reject'."""
    if path is None:
        print("ERROR: --config is required. No config file specified.", file=sys.stderr)
        sys.exit(1)
    p = Path(path)
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
    missing = [k for k in REQUIRED_KEYS if k not in thermo]
    if missing:
        print(f"ERROR: Config missing required keys: {', '.join(missing)}", file=sys.stderr)
        sys.exit(1)

    values = {}
    hard_reject = set()
    for key in REQUIRED_KEYS:
        val, action = _parse_thermo_entry(key, thermo[key])
        values[key] = val
        if action == "reject" and key in CONFIG_KEY_TO_FLAG:
            hard_reject.add(CONFIG_KEY_TO_FLAG[key])

    return values, hard_reject


# ---------------------------------------------------------------------------
# Thermodynamic backend: primer3-py or ntthal CLI
# ---------------------------------------------------------------------------
_BACKEND = None


def _init_backend():
    global _BACKEND
    if _BACKEND is not None:
        return
    try:
        import primer3 as _p3
        _BACKEND = "primer3"
        return
    except ImportError:
        pass
    # Try ntthal CLI
    try:
        subprocess.run(["ntthal", "-h"], capture_output=True, check=False)
        _BACKEND = "ntthal"
        return
    except FileNotFoundError:
        pass
    print(
        "ERROR: Neither primer3-py nor ntthal CLI found.\n"
        "Install primer3-py:  pip install primer3-py\n"
        "Or install Primer3 and ensure ntthal is on PATH.",
        file=sys.stderr,
    )
    sys.exit(1)


# --- primer3-py wrappers ---------------------------------------------------
def _p3_calc_tm(seq):
    import primer3
    return primer3.calc_tm(seq)


def _p3_calc_hairpin(seq):
    import primer3
    return primer3.calc_hairpin(seq).dg / 1000.0  # cal -> kcal


def _p3_calc_homodimer(seq):
    import primer3
    return primer3.calc_homodimer(seq).dg / 1000.0


def _p3_calc_heterodimer(seq1, seq2):
    import primer3
    return primer3.calc_heterodimer(seq1, seq2).dg / 1000.0


# --- ntthal CLI wrappers ---------------------------------------------------
def _ntthal_run(mode, s1, s2=None):
    """Run ntthal and parse ΔG from output."""
    cmd = ["ntthal", "-a", mode, "-s1", s1]
    if s2:
        cmd += ["-s2", s2]
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    # ntthal prints lines; last numerical token is ΔG in cal/mol
    for line in result.stdout.strip().splitlines():
        parts = line.split()
        for part in reversed(parts):
            try:
                return float(part) / 1000.0  # cal -> kcal
            except ValueError:
                continue
    return 0.0


def _ntthal_tm(seq):
    """Approximate Tm via ntthal HAIRPIN mode (complement binding)."""
    cmd = ["ntthal", "-a", "ANY", "-s1", seq, "-s2", seq]
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    for line in result.stdout.strip().splitlines():
        parts = line.split()
        for part in reversed(parts):
            try:
                return float(part)
            except ValueError:
                continue
    # Fallback: basic nearest-neighbor approximation
    return _basic_tm(seq)


def _basic_tm(seq):
    """Wallace rule fallback (only used if ntthal also fails)."""
    seq = seq.upper()
    gc = sum(1 for b in seq if b in "GC")
    at = sum(1 for b in seq if b in "AT")
    if len(seq) < 14:
        return 2 * at + 4 * gc
    return 64.9 + 41 * (gc - 16.4) / len(seq)


# --- IUPAC degenerate base handling ----------------------------------------
IUPAC_TO_FIRST = {
    "R": "A", "Y": "C", "S": "G", "W": "A", "K": "G", "M": "A",
    "B": "C", "D": "A", "H": "A", "V": "A", "N": "A",
}
IUPAC_EXPAND = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}

def resolve_iupac(seq: str) -> str:
    """Replace IUPAC ambiguity codes with a concrete base for primer3."""
    return "".join(IUPAC_TO_FIRST.get(b, b) for b in seq.upper())

def count_degeneracy(seq: str) -> int:
    """Count total degeneracy (product of ambiguity at each position)."""
    d = 1
    for b in seq.upper():
        d *= len(IUPAC_EXPAND.get(b, b))
    return d


# --- Unified interface -----------------------------------------------------
def calc_tm(seq):
    seq = resolve_iupac(seq)
    if _BACKEND == "primer3":
        return _p3_calc_tm(seq)
    return _ntthal_tm(seq)


def calc_hairpin(seq):
    seq = resolve_iupac(seq)
    if _BACKEND == "primer3":
        return _p3_calc_hairpin(seq)
    return _ntthal_run("HAIRPIN", seq)


def calc_homodimer(seq):
    seq = resolve_iupac(seq)
    if _BACKEND == "primer3":
        return _p3_calc_homodimer(seq)
    return _ntthal_run("ANY", seq, seq)


def calc_heterodimer(seq1, seq2):
    seq1, seq2 = resolve_iupac(seq1), resolve_iupac(seq2)
    if _BACKEND == "primer3":
        return _p3_calc_heterodimer(seq1, seq2)
    return _ntthal_run("ANY", seq1, seq2)


# ---------------------------------------------------------------------------
# Sequence-level calculations (no external tool needed)
# ---------------------------------------------------------------------------
def gc_content(seq):
    seq = seq.upper()
    gc = sum(1 for b in seq if b in "GC")
    return 100.0 * gc / len(seq) if seq else 0.0


def gc_clamp_count(seq, window=5):
    """Count G/C bases in the last `window` bases of seq."""
    tail = seq.upper()[-window:]
    return sum(1 for b in tail if b in "GC")


def max_homopolymer(seq):
    """Longest run of identical bases."""
    if not seq:
        return 0
    seq = seq.upper()
    max_run = 1
    run = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            run += 1
            max_run = max(max_run, run)
        else:
            run = 1
    return max_run


# ---------------------------------------------------------------------------
# Individual primer validation
# ---------------------------------------------------------------------------
def validate_primer(seq, cfg, hard_reject):
    """Validate a single primer sequence. Returns dict of metrics + flags.
    hard_reject: set of flag names that cause rejection (from config actions)."""
    seq = seq.strip().upper()

    tm = calc_tm(seq)
    hairpin = calc_hairpin(seq)
    homodimer = calc_homodimer(seq)
    gc = gc_content(seq)
    clamp = gc_clamp_count(seq)
    homopoly = max_homopolymer(seq)

    row = {
        "sequence": seq,
        "length": len(seq),
        "tm": round(tm, 1),
        "hairpin_dg": round(hairpin, 2),
        "homodimer_dg": round(homodimer, 2),
        "gc_pct": round(gc, 1),
        "gc_clamp": clamp,
        "homopolymer": homopoly,
    }

    # --- Flag logic ---
    flags = []
    tm_lo, tm_hi = cfg["tm_range"]
    gc_lo, gc_hi = cfg["gc_range"]
    clamp_lo, clamp_hi = cfg["gc_clamp_3prime"]

    if not (tm_lo <= tm <= tm_hi):
        flags.append("tm_out_of_range")
    if hairpin < cfg["max_hairpin_dg"]:
        flags.append("hairpin")
    if homodimer < cfg["max_homodimer_dg"]:
        flags.append("homodimer")
    if not (gc_lo <= gc <= gc_hi):
        flags.append("gc_out_of_range")
    if not (clamp_lo <= clamp <= clamp_hi):
        flags.append("gc_clamp")
    if homopoly > cfg["max_homopolymer"]:
        flags.append("homopolymer")

    # Hard reject vs soft flag — determined by config actions
    hard = hard_reject & set(flags)

    row["flags"] = ";".join(flags) if flags else ""
    row["status"] = "REJECT" if hard else ("FLAG" if flags else "PASS")

    return row


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------
OUTPUT_COLUMNS = [
    "sequence", "direction", "length",
    "tm", "hairpin_dg", "homodimer_dg",
    "gc_pct", "gc_clamp", "homopolymer",
    "flags", "status",
]


def run(input_csv, output_csv, config_path, keep_rejected):
    cfg, hard_reject = load_config(config_path)
    _init_backend()
    print(f"[info] Backend: {_BACKEND}", file=sys.stderr)
    print(f"[info] Hard-reject checks: {', '.join(sorted(hard_reject)) or 'none'}", file=sys.stderr)

    # Extract unique FWD and REV sequences from pairs CSV
    unique_fwd = {}  # seq -> True
    unique_rev = {}
    with open(input_csv, newline="") as f:
        reader = csv.DictReader(f)
        if not {"name", "forward", "reverse"}.issubset(reader.fieldnames):
            print("ERROR: Input CSV must have columns: name, forward, reverse",
                  file=sys.stderr)
            sys.exit(1)
        for row in reader:
            fwd = row["forward"].strip().upper()
            rev = row["reverse"].strip().upper()
            unique_fwd[fwd] = True
            unique_rev[rev] = True

    all_primers = (
        [(seq, "FWD") for seq in unique_fwd]
        + [(seq, "REV") for seq in unique_rev]
    )
    print(
        f"[info] {len(unique_fwd)} unique FWD + {len(unique_rev)} unique REV "
        f"= {len(all_primers)} individual primers to validate",
        file=sys.stderr,
    )

    workers = int(os.environ.get("PRIMER_WORKERS", max(1, multiprocessing.cpu_count() - 1)))
    print(f"[info] Using {workers} workers", file=sys.stderr)

    def _validate_one(item):
        seq, direction = item
        result = validate_primer(seq, cfg, hard_reject)
        result["direction"] = direction
        return result

    if workers > 1:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            results_all = list(pool.map(_validate_one, all_primers))
    else:
        results_all = [_validate_one(item) for item in all_primers]

    results = []
    n_pass, n_flag, n_reject = 0, 0, 0
    for row in results_all:
        if row["status"] == "PASS":
            n_pass += 1
        elif row["status"] == "FLAG":
            n_flag += 1
        else:
            n_reject += 1
        if keep_rejected or row["status"] != "REJECT":
            results.append(row)

    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        writer.writerows(results)

    print(
        f"[info] {len(all_primers)} primers evaluated: "
        f"{n_pass} PASS, {n_flag} FLAG, {n_reject} REJECT",
        file=sys.stderr,
    )
    print(f"[info] Wrote {len(results)} primers to {output_csv}", file=sys.stderr)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Thermodynamic validation of individual primers (Step 2)."
    )
    parser.add_argument(
        "input", help="Input CSV with columns: name, forward, reverse"
    )
    parser.add_argument(
        "-o", "--output", default=None,
        help="Output CSV (default: primers_validated.csv in same dir as input)"
    )
    parser.add_argument(
        "--config", required=True,
        help="Path to config YAML (required)"
    )
    parser.add_argument(
        "--keep-rejected", action="store_true",
        help="Include rejected primers in output (marked as REJECT)"
    )
    args = parser.parse_args()

    if args.output is None:
        inp = Path(args.input)
        args.output = str(inp.parent / "primers_validated.csv")

    run(args.input, args.output, args.config, args.keep_rejected)


if __name__ == "__main__":
    main()
