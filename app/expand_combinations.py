#!/usr/bin/env python3
"""
Step 2b: Expand all viable FWD+REV combinations from filtered primers.

Reads primers_filtered.csv (output of validate_thermodynamics.py), extracts
unique forward and reverse sequences, and generates all thermodynamically
viable combinations. Combinations are filtered by both delta Tm and estimated
amplicon distance on the reference sequence, so only spatially compatible
pairs reach the analysis step.

Usage:
    python scripts/expand_combinations.py primers_filtered.csv -i alignment.fasta -o primers_expanded.csv
    python scripts/expand_combinations.py primers_filtered.csv -i alignment.fasta -o primers_expanded.csv --config primer_config.yaml
"""

import argparse
import csv
import re
import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    yaml = None


# ---------------------------------------------------------------------------
# IUPAC fuzzy matching
# ---------------------------------------------------------------------------
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}


def iupac_to_regex(seq: str) -> re.Pattern:
    """Convert an IUPAC primer sequence to a regex pattern."""
    return re.compile("".join(IUPAC.get(b, b) for b in seq.upper()))


def revcomp(seq: str) -> str:
    """Reverse complement supporting IUPAC codes."""
    comp = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn",
                         "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")
    return seq.translate(comp)[::-1]


# ---------------------------------------------------------------------------
# MSA reading and primer location
# ---------------------------------------------------------------------------
def read_msa_reference(fasta_path: str) -> tuple[str, str]:
    """Read first non-annotation sequence from MSA. Returns (name, gapped_seq)."""
    name = None
    seq_parts = []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    # Return the first real sequence (skip annotation track)
                    seq = "".join(seq_parts).upper()
                    if not all(c in "-." for c in seq):
                        return name, seq
                header = line[1:]
                # Skip annotation tracks (contain gene names like "tGlu", "12S" etc.)
                if any(tag in header for tag in ["anno", "tGlu", "cytb", "12S", "16S"]):
                    name = None
                    seq_parts = []
                else:
                    name = header.split()[0]
                    seq_parts = []
            elif name is not None and line:
                seq_parts.append(line.upper())
    if name is not None:
        seq = "".join(seq_parts).upper()
        if not all(c in "-." for c in seq):
            return name, seq
    raise ValueError(f"No valid reference sequence found in {fasta_path}")


def locate_primer_on_ref(primer_seq: str, ref_ungapped: str, is_reverse: bool = False) -> int | None:
    """Find the start position of a primer on the ungapped reference sequence.

    For reverse primers, searches the reverse complement against the forward strand.
    Returns the ungapped start position, or None if not found.
    """
    search_seq = revcomp(primer_seq) if is_reverse else primer_seq
    pattern = iupac_to_regex(search_seq)
    match = pattern.search(ref_ungapped)
    if match:
        return match.start()
    return None


def build_ungapped_map(gapped_seq: str) -> dict[int, int]:
    """Build mapping from ungapped position -> gapped position."""
    mapping = {}
    ug_pos = 0
    for g_pos, base in enumerate(gapped_seq):
        if base not in "-. ":
            mapping[ug_pos] = g_pos
            ug_pos += 1
    return mapping


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Expand all viable FWD+REV combinations from filtered primers."
    )
    parser.add_argument("input", help="Filtered primers CSV (from validate_thermodynamics.py)")
    parser.add_argument("-i", "--msa", default=None,
                        help="Input MSA FASTA — used to filter by amplicon distance on the reference")
    parser.add_argument("-o", "--output", default="primers_expanded.csv", help="Output CSV")
    parser.add_argument("--config", default=None, help="Path to primer_config.yaml")
    args = parser.parse_args()

    cfg = load_config(args.config)
    max_delta_tm = cfg["max_delta_tm"]
    min_amp_len = cfg["min_amplicon_length"]
    max_amp_len = cfg["max_amplicon_length"]

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

    # Locate primers on the reference sequence for spatial filtering
    fwd_positions = {}  # seq -> ungapped start position
    rev_positions = {}  # seq -> ungapped start position (of revcomp match)
    spatial_filter = False

    if args.msa:
        try:
            ref_name, ref_gapped = read_msa_reference(args.msa)
            ref_ungapped = ref_gapped.replace("-", "").replace(".", "")
            print(f"Reference: {ref_name} ({len(ref_ungapped)} bp ungapped)")

            for fwd in fwds:
                pos = locate_primer_on_ref(fwd["seq"], ref_ungapped, is_reverse=False)
                if pos is not None:
                    fwd_positions[fwd["seq"]] = pos

            for rev in revs:
                pos = locate_primer_on_ref(rev["seq"], ref_ungapped, is_reverse=True)
                if pos is not None:
                    rev_positions[rev["seq"]] = pos

            n_fwd_located = len(fwd_positions)
            n_rev_located = len(rev_positions)
            print(f"Located on reference: {n_fwd_located}/{len(fwds)} FWD, {n_rev_located}/{len(revs)} REV")

            if n_fwd_located > 0 and n_rev_located > 0:
                spatial_filter = True
                print(f"Spatial filter: amplicon {min_amp_len}–{max_amp_len} bp")
            else:
                print("WARNING: Could not locate enough primers on reference — skipping spatial filter")
        except Exception as e:
            print(f"WARNING: Could not read MSA for spatial filtering: {e}")

    # Generate all viable combinations
    viable = []
    n_tm_rejected = 0
    n_spatial_rejected = 0

    for fi, fwd in enumerate(fwds):
        for ri, rev in enumerate(revs):
            # Delta Tm filter
            delta_tm = abs(fwd["tm"] - rev["tm"])
            if delta_tm > max_delta_tm:
                n_tm_rejected += 1
                continue

            # Spatial filter: check amplicon distance on reference
            is_original = (fwd["seq"], rev["seq"]) in original_pairs
            if spatial_filter:
                fwd_pos = fwd_positions.get(fwd["seq"])
                rev_pos = rev_positions.get(rev["seq"])
                if fwd_pos is not None and rev_pos is not None:
                    # rev_pos is where the revcomp starts on the fwd strand,
                    # so the amplicon runs from fwd_pos to rev_pos + len(rev)
                    amp_len = (rev_pos + len(rev["seq"])) - fwd_pos
                    if amp_len < min_amp_len or amp_len > max_amp_len:
                        n_spatial_rejected += 1
                        continue
                elif not is_original:
                    # Novel combo where we can't verify distance — skip it
                    n_spatial_rejected += 1
                    continue

            name_suffix = f"F{fi+1:02d}_R{ri+1:02d}"
            viable.append({
                "name": f"COMBO_{name_suffix}" if not is_original else f"ORIG_{name_suffix}",
                "forward": fwd["seq"],
                "reverse": rev["seq"],
            })

    orig_count = sum(1 for v in viable if v["name"].startswith("ORIG_"))
    print(f"Rejected: {n_tm_rejected} (delta Tm), {n_spatial_rejected} (amplicon distance)")
    print(f"Viable combinations: {len(viable)}")
    print(f"  Original: {orig_count}, New: {len(viable) - orig_count}")

    # Write output
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["name", "forward", "reverse"])
        writer.writeheader()
        writer.writerows(viable)

    print(f"Wrote {len(viable)} primer pairs to {args.output}")


if __name__ == "__main__":
    main()
