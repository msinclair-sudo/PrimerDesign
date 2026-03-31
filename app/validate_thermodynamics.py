#!/usr/bin/env python3
"""Step 2: Thermodynamic validation of primer pairs.

Calculates Tm, hairpin/homodimer/heterodimer ΔG, GC%, GC clamp,
and homopolymer runs. Applies hard-reject and soft-flag logic.

Requires primer3-py (preferred) or ntthal CLI as fallback.
"""

import argparse
import csv
import json
import multiprocessing
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

try:
    import yaml
except ImportError:
    yaml = None

# ---------------------------------------------------------------------------
# Defaults (match primer_config.yaml / spec)
# ---------------------------------------------------------------------------
DEFAULTS = {
    "tm_range": [58, 64],
    "max_delta_tm": 5,
    "gc_range": [40, 60],
    "max_hairpin_dg": -3.0,
    "max_homodimer_dg": -9.0,
    "max_heterodimer_dg": -9.0,
    "max_homopolymer": 4,
    "gc_clamp_3prime": [1, 3],
}

# Hard-reject thresholds: if ANY of these fail, the pair is rejected.
# Everything else is a soft flag.
HARD_REJECT = {"tm_range", "max_delta_tm", "max_hairpin_dg"}


# ---------------------------------------------------------------------------
# Config loading
# ---------------------------------------------------------------------------
def load_config(path):
    """Load thermodynamic thresholds from YAML config, falling back to DEFAULTS."""
    cfg = dict(DEFAULTS)
    if path is None:
        return cfg
    p = Path(path)
    if not p.exists():
        print(f"[warn] Config not found at {p}, using defaults", file=sys.stderr)
        return cfg
    if yaml is None:
        print("[warn] PyYAML not installed, using defaults", file=sys.stderr)
        return cfg
    with open(p) as f:
        raw = yaml.safe_load(f)
    thermo = raw.get("thermodynamics", {})
    for key in DEFAULTS:
        if key in thermo:
            cfg[key] = thermo[key]
    return cfg


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
# Validation
# ---------------------------------------------------------------------------
def validate_pair(name, fwd, rev, cfg):
    """Validate a single primer pair. Returns dict of metrics + flags."""
    fwd, rev = fwd.strip().upper(), rev.strip().upper()

    tm_fwd = calc_tm(fwd)
    tm_rev = calc_tm(rev)
    delta_tm = abs(tm_fwd - tm_rev)
    hairpin_fwd = calc_hairpin(fwd)
    hairpin_rev = calc_hairpin(rev)
    homodimer_fwd = calc_homodimer(fwd)
    homodimer_rev = calc_homodimer(rev)
    heterodimer = calc_heterodimer(fwd, rev)
    gc_fwd = gc_content(fwd)
    gc_rev = gc_content(rev)
    gc_clamp_fwd = gc_clamp_count(fwd)
    gc_clamp_rev = gc_clamp_count(rev)
    homopoly_fwd = max_homopolymer(fwd)
    homopoly_rev = max_homopolymer(rev)

    row = {
        "name": name,
        "forward": fwd,
        "reverse": rev,
        "tm_fwd": round(tm_fwd, 1),
        "tm_rev": round(tm_rev, 1),
        "delta_tm": round(delta_tm, 1),
        "hairpin_dg_fwd": round(hairpin_fwd, 2),
        "hairpin_dg_rev": round(hairpin_rev, 2),
        "homodimer_dg_fwd": round(homodimer_fwd, 2),
        "homodimer_dg_rev": round(homodimer_rev, 2),
        "heterodimer_dg": round(heterodimer, 2),
        "gc_fwd": round(gc_fwd, 1),
        "gc_rev": round(gc_rev, 1),
        "gc_clamp_fwd": gc_clamp_fwd,
        "gc_clamp_rev": gc_clamp_rev,
        "homopolymer_fwd": homopoly_fwd,
        "homopolymer_rev": homopoly_rev,
    }

    # --- Flag logic ---
    flags = []
    tm_lo, tm_hi = cfg["tm_range"]
    gc_lo, gc_hi = cfg["gc_range"]
    clamp_lo, clamp_hi = cfg["gc_clamp_3prime"]

    if not (tm_lo <= tm_fwd <= tm_hi):
        flags.append("tm_fwd_out_of_range")
    if not (tm_lo <= tm_rev <= tm_hi):
        flags.append("tm_rev_out_of_range")
    if delta_tm > cfg["max_delta_tm"]:
        flags.append("delta_tm_too_high")
    if hairpin_fwd < cfg["max_hairpin_dg"]:
        flags.append("hairpin_fwd")
    if hairpin_rev < cfg["max_hairpin_dg"]:
        flags.append("hairpin_rev")
    if homodimer_fwd < cfg["max_homodimer_dg"]:
        flags.append("homodimer_fwd")
    if homodimer_rev < cfg["max_homodimer_dg"]:
        flags.append("homodimer_rev")
    if heterodimer < cfg["max_heterodimer_dg"]:
        flags.append("heterodimer")
    if not (gc_lo <= gc_fwd <= gc_hi):
        flags.append("gc_fwd_out_of_range")
    if not (gc_lo <= gc_rev <= gc_hi):
        flags.append("gc_rev_out_of_range")
    if not (clamp_lo <= gc_clamp_fwd <= clamp_hi):
        flags.append("gc_clamp_fwd")
    if not (clamp_lo <= gc_clamp_rev <= clamp_hi):
        flags.append("gc_clamp_rev")
    if homopoly_fwd > cfg["max_homopolymer"]:
        flags.append("homopolymer_fwd")
    if homopoly_rev > cfg["max_homopolymer"]:
        flags.append("homopolymer_rev")

    # Determine hard reject vs soft flag
    hard_flags = {f for f in flags if any(f.startswith(h.replace("max_", "").replace("_range", ""))
                                        for h in HARD_REJECT)}
    # More precise: check each flag against hard-reject categories
    hard = set()
    for f in flags:
        if f.startswith("tm_") or f == "delta_tm_too_high":
            hard.add(f)
        if f.startswith("hairpin_"):
            hard.add(f)

    row["flags"] = ";".join(flags) if flags else ""
    row["status"] = "REJECT" if hard else ("FLAG" if flags else "PASS")

    return row


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------
OUTPUT_COLUMNS = [
    "name", "forward", "reverse",
    "tm_fwd", "tm_rev", "delta_tm",
    "hairpin_dg_fwd", "hairpin_dg_rev",
    "homodimer_dg_fwd", "homodimer_dg_rev",
    "heterodimer_dg",
    "gc_fwd", "gc_rev",
    "gc_clamp_fwd", "gc_clamp_rev",
    "homopolymer_fwd", "homopolymer_rev",
    "flags", "status",
]


def run(input_csv, output_csv, config_path, keep_rejected):
    cfg = load_config(config_path)
    _init_backend()
    print(f"[info] Backend: {_BACKEND}", file=sys.stderr)

    with open(input_csv, newline="") as f:
        reader = csv.DictReader(f)
        if not {"name", "forward", "reverse"}.issubset(reader.fieldnames):
            print("ERROR: Input CSV must have columns: name, forward, reverse",
                  file=sys.stderr)
            sys.exit(1)
        pairs = list(reader)

    workers = int(os.environ.get("PRIMER_WORKERS", max(1, multiprocessing.cpu_count() - 1)))
    print(f"[info] Using {workers} workers", file=sys.stderr)

    def _validate_one(pair):
        return validate_pair(pair["name"], pair["forward"], pair["reverse"], cfg)

    if workers > 1:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            results_all = list(pool.map(_validate_one, pairs))
    else:
        results_all = [validate_pair(p["name"], p["forward"], p["reverse"], cfg) for p in pairs]

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
        f"[info] {len(pairs)} pairs evaluated: "
        f"{n_pass} PASS, {n_flag} FLAG, {n_reject} REJECT",
        file=sys.stderr,
    )
    print(f"[info] Wrote {len(results)} pairs to {output_csv}", file=sys.stderr)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Thermodynamic validation of primer pairs (Step 2)."
    )
    parser.add_argument(
        "input", help="Input CSV with columns: name, forward, reverse"
    )
    parser.add_argument(
        "-o", "--output", default=None,
        help="Output CSV (default: primers_filtered.csv in same dir as input)"
    )
    parser.add_argument(
        "--config", default=None,
        help="Path to primer_config.yaml (default: primer_config.yaml)"
    )
    parser.add_argument(
        "--keep-rejected", action="store_true",
        help="Include rejected primers in output (marked as REJECT)"
    )
    args = parser.parse_args()

    if args.output is None:
        inp = Path(args.input)
        args.output = str(inp.parent / "primers_filtered.csv")

    # Auto-discover config if not specified
    if args.config is None:
        candidates = [
            Path("primer_config.yaml"),
        ]
        for c in candidates:
            if c.exists():
                args.config = str(c)
                break

    run(args.input, args.output, args.config, args.keep_rejected)


if __name__ == "__main__":
    main()
