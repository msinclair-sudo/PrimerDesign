#!/usr/bin/env python3
"""
Alignment Preparation for Circular Genomes
===========================================
Prepares a trimmed, correctly-oriented MSA from raw mitogenome sequences.

Handles circular genome artifacts: sequences may start at different positions
in the genome, causing MAFFT to misalign at the linear ends. This script
aligns, rotates to center the target region, realigns, then trims.

Usage:
    python app/prepare_alignment.py input.fasta -g reference.gb --from cytb --to 16S -o output.fasta

Steps:
    1. Initial MAFFT alignment
    2. Identify gene regions from GenBank annotation
    3. Determine rotation point (midpoint of target gene range)
    4. Rotate all ungapped sequences to center the target region
    5. Realign with MAFFT
    6. Trim to target gene range
"""

import argparse
import hashlib
import math
import re
import subprocess
import sys
import tempfile
from pathlib import Path

from Bio import SeqIO

# Import reusable functions from sibling scripts
sys.path.insert(0, str(Path(__file__).parent))
from analyse_primers import (
    FEATURE_SHORT_NAMES,
    _blast_find_offset,
    _extract_and_map_features,
)


# ---------------------------------------------------------------------------
# FASTA I/O
# ---------------------------------------------------------------------------
def read_fasta(path: str) -> tuple[list[tuple[str, str]], dict[str, str]]:
    """Read FASTA file. Returns (records, header_map).
    records: list of (accession, sequence) tuples.
    header_map: accession → full header string (for restoring after MAFFT)."""
    records = []
    header_map = {}
    name = None
    full_header = None
    seq_parts = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq_parts).upper()))
                full_header = line[1:]
                name = full_header.split()[0]
                header_map[name] = full_header
                seq_parts = []
            elif name is not None:
                seq_parts.append(line)
    if name is not None:
        records.append((name, "".join(seq_parts).upper()))
    return records, header_map


def write_fasta(path: str, records: list[tuple[str, str]],
                header_map: dict[str, str] | None = None, line_width: int = 70):
    """Write FASTA file. If header_map provided, restores full headers from accessions."""
    with open(path, "w") as f:
        for name, seq in records:
            header = header_map.get(name, name) if header_map else name
            f.write(f">{header}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i + line_width] + "\n")


# ---------------------------------------------------------------------------
# N-masking: replace runs of N with gaps
# ---------------------------------------------------------------------------
def mask_n_runs(records: list[tuple[str, str]], min_run: int = 3) -> list[tuple[str, str]]:
    """Replace consecutive N bases (runs >= min_run) with gap characters.

    Short isolated Ns (< min_run) are left as-is since they may represent
    genuine single-base ambiguities rather than missing data. Long N-runs
    typically indicate low-coverage assembly regions that would distort
    the alignment if treated as real bases.
    """
    pattern = re.compile(f"[Nn]{{{min_run},}}")
    masked = []
    total_replaced = 0
    seqs_affected = 0
    for name, seq in records:
        new_seq, count = pattern.subn(lambda m: "-" * len(m.group()), seq)
        if count > 0:
            seqs_affected += 1
            n_bases = sum(len(m.group()) for m in pattern.finditer(seq))
            total_replaced += n_bases
        masked.append((name, new_seq))
    if total_replaced > 0:
        print(f"  Masked {total_replaced} N bases → gaps "
              f"across {seqs_affected} sequence(s) (runs >= {min_run})")
    else:
        print(f"  No N-runs >= {min_run} found")
    return masked


# ---------------------------------------------------------------------------
# N-restoration: put Ns back after alignment
# ---------------------------------------------------------------------------
_COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def _masked_positions(seq: str, min_run: int) -> set[int]:
    """Return set of positions that are part of N-runs >= min_run."""
    positions = set()
    for m in re.finditer(f"[Nn]{{{min_run},}}", seq):
        positions.update(range(m.start(), m.end()))
    return positions


def _rotation_in_original(seq: str, non_masked_rot: int, min_run: int) -> int:
    """Convert a rotation count (in non-masked characters) to full sequence position.

    The pipeline computes rotation in terms of characters that MAFFT sees
    (i.e. everything except masked N-runs). This converts that count back
    to a position in the original sequence that includes masked Ns.
    """
    masked = _masked_positions(seq, min_run)
    count = 0
    for i in range(len(seq)):
        if i not in masked:
            if count == non_masked_rot:
                return i
            count += 1
    return len(seq)


def restore_ns(realigned: list[tuple[str, str]], original_seqs: dict[str, str],
               rotations: dict[str, int], flipped: set[str], min_run: int = 3
               ) -> list[tuple[str, str]]:
    """Restore masked N-runs in realigned sequences.

    After N-masking, MAFFT alignment, rotation, and realignment, the masked
    Ns appear as alignment gaps. This function converts those gaps back to N
    by tracking where the Ns were in the original sequences, accounting for
    any reverse-complementation and rotation applied by the pipeline.
    """
    restored = []
    n_total = 0
    n_seqs = 0

    for name, aligned_seq in realigned:
        orig = original_seqs.get(name)
        if not orig:
            restored.append((name, aligned_seq))
            continue

        masked = _masked_positions(orig, min_run)
        if not masked:
            restored.append((name, aligned_seq))
            continue

        # Reverse complement if MAFFT flipped this sequence
        if name in flipped:
            orig = _revcomp(orig)

        # Rotate original to match the pipeline rotation
        rot = rotations.get(name, 0)
        full_rot = _rotation_in_original(orig, rot, min_run)
        rotated_orig = rotate_sequence(orig, full_rot)

        # Restore Ns in the aligned sequence
        new_seq = _restore_ns_in_sequence(aligned_seq, rotated_orig, min_run)

        changes = sum(1 for a, b in zip(aligned_seq, new_seq) if a != b)
        if changes:
            n_total += changes
            n_seqs += 1

        restored.append((name, new_seq))

    if n_total:
        print(f"  Restored {n_total} N bases across {n_seqs} sequence(s)")
    else:
        print(f"  No N bases to restore")

    return restored


def _restore_ns_in_sequence(aligned_seq: str, rotated_original: str,
                            min_run: int) -> str:
    """Convert alignment gaps back to Ns where the original had masked N-runs.

    The aligned sequence's non-gap characters correspond (in order) to the
    non-masked characters in rotated_original. Gaps in the aligned sequence
    that fall between two non-masked characters which had masked Ns between
    them in the original are converted back to N.
    """
    masked = _masked_positions(rotated_original, min_run)
    if not masked:
        return aligned_seq

    # Indices of non-masked characters in the rotated original
    non_masked_idx = [i for i in range(len(rotated_original)) if i not in masked]

    if not non_masked_idx:
        # Entire sequence was masked — convert all gaps to N
        return "".join("N" if c in "-. " else c for c in aligned_seq)

    # Count masked Ns before first non-masked char, between consecutive
    # non-masked chars, and after last non-masked char
    n_before = sum(1 for p in range(non_masked_idx[0]) if p in masked)
    n_between = []
    for i in range(len(non_masked_idx) - 1):
        count = sum(1 for p in range(non_masked_idx[i] + 1, non_masked_idx[i + 1])
                    if p in masked)
        n_between.append(count)
    n_after = sum(1 for p in range(non_masked_idx[-1] + 1, len(rotated_original))
                  if p in masked)

    result = list(aligned_seq)
    non_gap_pos = [i for i, c in enumerate(result) if c not in "-. "]

    if not non_gap_pos:
        # All gaps — convert up to total masked count
        for i in range(min(len(masked), len(result))):
            result[i] = "N"
        return "".join(result)

    # Leading Ns: convert gaps immediately before the first base, right-to-left
    placed = 0
    for i in range(non_gap_pos[0] - 1, -1, -1):
        if placed >= n_before:
            break
        if result[i] in "-. ":
            result[i] = "N"
            placed += 1

    # Ns between consecutive bases: convert gaps left-to-right
    for idx in range(len(non_gap_pos) - 1):
        ns_needed = n_between[idx] if idx < len(n_between) else 0
        if ns_needed == 0:
            continue
        placed = 0
        for i in range(non_gap_pos[idx] + 1, non_gap_pos[idx + 1]):
            if placed >= ns_needed:
                break
            if result[i] in "-. ":
                result[i] = "N"
                placed += 1

    # Trailing Ns: convert gaps after the last base, left-to-right
    placed = 0
    for i in range(non_gap_pos[-1] + 1, len(result)):
        if placed >= n_after:
            break
        if result[i] in "-. ":
            result[i] = "N"
            placed += 1

    return "".join(result)


# ---------------------------------------------------------------------------
# MAFFT
# ---------------------------------------------------------------------------
def run_mafft(records: list[tuple[str, str]], adjust_direction: bool = False
              ) -> tuple[list[tuple[str, str]], set[str]]:
    """Run MAFFT alignment. Returns (aligned_records, flipped_names).
    If adjust_direction=True, uses --adjustdirection to auto-revcomp misoriented sequences."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tmp:
        for name, seq in records:
            tmp.write(f">{name}\n{seq}\n")
        tmp_path = tmp.name

    cmd = ["mafft", "--auto", "--thread", "-1"]
    if adjust_direction:
        cmd.append("--adjustdirection")
    cmd.append(tmp_path)

    result = subprocess.run(cmd, capture_output=True, text=True)
    Path(tmp_path).unlink(missing_ok=True)

    if result.returncode != 0:
        print(f"MAFFT stderr: {result.stderr}", file=sys.stderr)
        sys.exit("ERROR: MAFFT alignment failed")

    # Parse output — MAFFT prepends _R_ to names of sequences it reverse complemented
    aligned = {}
    flipped = []
    current_name = None
    current_seq = []
    for line in result.stdout.splitlines():
        if line.startswith(">"):
            if current_name:
                aligned[current_name] = "".join(current_seq)
            raw_name = line[1:].split()[0]
            if raw_name.startswith("_R_"):
                current_name = raw_name[3:]  # strip _R_ prefix
                flipped.append(current_name)
            else:
                current_name = raw_name
            current_seq = []
        else:
            current_seq.append(line.strip())
    if current_name:
        aligned[current_name] = "".join(current_seq)

    if flipped:
        print(f"  MAFFT reverse-complemented {len(flipped)} sequence(s): {', '.join(flipped)}")

    # Return in input order, plus set of flipped names
    ordered = [(name, aligned[name]) for name in [r[0] for r in records] if name in aligned]
    return ordered, set(flipped)


# ---------------------------------------------------------------------------
# GenBank feature extraction
# ---------------------------------------------------------------------------
def get_genome_features(gb_path: str) -> tuple[str, int, str, list[dict]]:
    """Parse GenBank file. Returns (accession, genome_length, sequence, features).
    Features are sorted by genome position with standardized names."""
    gb = SeqIO.read(gb_path, "genbank")
    gb_seq = str(gb.seq).upper()
    genome_len = len(gb_seq)
    acc = gb.id.split(".")[0]

    features = []
    for feat in gb.features:
        if feat.type in ("source", "gene"):
            continue
        name = feat.qualifiers.get("product", feat.qualifiers.get("gene", ["?"]))[0]
        if name == "?" and feat.type == "D-loop":
            name = "Dloop"
        name = FEATURE_SHORT_NAMES.get(name, name)
        if name.startswith("tRNA-"):
            name = "t" + name[5:]

        features.append({
            "name": name,
            "type": feat.type,
            "start": int(feat.location.start),  # 0-based
            "end": int(feat.location.end),       # 0-based exclusive
        })

    features.sort(key=lambda f: f["start"])
    return acc, genome_len, gb_seq, features


def find_ref_in_alignment(gb_acc: str, gb_seq: str, aligned_records: list[tuple[str, str]]) -> str | None:
    """Find which alignment sequence matches the GenBank record."""
    doubled = gb_seq + gb_seq
    for name, seq in aligned_records:
        if gb_acc in name:
            return name
        ungapped = seq.replace("-", "").replace(".", "")
        if ungapped in doubled:
            return name
    return None


# ---------------------------------------------------------------------------
# Gene range resolution (circular-aware)
# ---------------------------------------------------------------------------
def resolve_gene_name(query: str, features: list[dict]) -> str | None:
    """Match a user-supplied gene name to a feature name, case-insensitive."""
    q = query.lower().strip()
    for feat in features:
        if feat["name"].lower() == q:
            return feat["name"]
    # Try partial match
    for feat in features:
        if q in feat["name"].lower() or feat["name"].lower() in q:
            return feat["name"]
    return None


def get_gene_range_circular(from_name: str, to_name: str,
                            features: list[dict]) -> list[dict]:
    """Get all features between from_name and to_name in circular genome order.
    Wraps around if to_name comes before from_name in linear order."""
    # Get unique major features in genome order (skip tRNAs for counting,
    # but include them in the output)
    names = [f["name"] for f in features]

    from_idx = None
    to_idx = None
    for i, f in enumerate(features):
        if f["name"] == from_name and from_idx is None:
            from_idx = i
        if f["name"] == to_name:
            to_idx = i

    if from_idx is None or to_idx is None:
        return []

    # Collect features in circular order from from_idx to to_idx (inclusive)
    n = len(features)
    result = []
    i = from_idx
    while True:
        result.append(features[i])
        if i == to_idx:
            break
        i = (i + 1) % n
        if len(result) > n:  # safety
            break

    return result


# ---------------------------------------------------------------------------
# Rotation
# ---------------------------------------------------------------------------
def find_rotation_point(features_in_range: list[dict], genome_len: int) -> int:
    """Determine the genome position to start the rotated alignment.

    We want the entire from→to range to be contiguous after rotation.
    So we rotate to start just before the first gene in the target range.
    This places the target region at the beginning of the linear sequence,
    and the non-target region (which we'll trim away) at the end.
    """
    # Start a small distance before the first gene in the range
    first_start = features_in_range[0]["start"]
    # Back up by 100 bp to give some flanking context
    return (first_start - 100) % genome_len


def rotate_sequence(seq: str, amount: int) -> str:
    """Rotate a sequence so position `amount` becomes the start."""
    n = len(seq)
    amount = amount % n
    return seq[amount:] + seq[:amount]


def compute_rotation_per_sequence(
    aligned_records: list[tuple[str, str]],
    ref_name: str,
    ref_genome_rotation: int,
    gb_seq: str,
    genome_len: int,
) -> dict[str, int]:
    """For each sequence, compute how many ungapped bases to rotate.

    The reference sequence's rotation is known from the GenBank coordinates.
    For other sequences, we use the alignment to find the corresponding position.
    """
    # Find the reference's ungapped position corresponding to the genome rotation point
    ref_aligned = None
    for name, seq in aligned_records:
        if name == ref_name:
            ref_aligned = seq
            break

    if ref_aligned is None:
        sys.exit(f"ERROR: Reference {ref_name} not found in alignment")

    ref_ungapped = ref_aligned.replace("-", "").replace(".", "")

    # Find offset: where does the ref ungapped sequence start in the genome?
    offset = _blast_find_offset(ref_ungapped, gb_seq, genome_len)
    if offset is None:
        sys.exit("ERROR: Could not determine reference offset via BLAST")

    # The rotation point in ungapped ref coordinates
    ref_rotation_ungapped = (ref_genome_rotation - offset) % genome_len

    # Find the alignment column that corresponds to this ungapped position in the ref
    ungap_count = 0
    rotation_column = None
    for col, base in enumerate(ref_aligned):
        if base not in "-. ":
            if ungap_count == ref_rotation_ungapped:
                rotation_column = col
                break
            ungap_count += 1

    if rotation_column is None:
        sys.exit("ERROR: Could not map rotation point to alignment column")

    # For each sequence, find the ungapped position at that alignment column
    rotations = {}
    for name, seq in aligned_records:
        ungap_pos = 0
        found = False
        for col, base in enumerate(seq):
            if col == rotation_column:
                rotations[name] = ungap_pos
                found = True
                break
            if base not in "-. ":
                ungap_pos += 1
        if not found:
            # Column is past the end of this sequence; rotate to end
            rotations[name] = ungap_pos

    return rotations


# ---------------------------------------------------------------------------
# Trimming
# ---------------------------------------------------------------------------
def map_features_to_alignment(
    gb_record, offset: int, genome_len: int,
    ref_gapped: str, msa_ungapped_len: int,
) -> list[dict]:
    """Map GenBank features to gapped alignment coordinates."""
    return _extract_and_map_features(
        gb_record, offset, genome_len, msa_ungapped_len, ref_gapped
    )


def trim_alignment(
    aligned_records: list[tuple[str, str]],
    features: list[dict],
    from_name: str,
    to_name: str,
) -> list[tuple[str, str]]:
    """Trim alignment to the region between from_name and to_name features."""
    from_start = None
    to_end = None

    for f in features:
        if f["name"] == from_name and from_start is None:
            from_start = f["start"]
        if f["name"] == to_name:
            to_end = f["end"]

    if from_start is None or to_end is None:
        print(f"WARNING: Could not find {from_name} and/or {to_name} in mapped features",
              file=sys.stderr)
        print(f"Available features: {[f['name'] for f in features]}", file=sys.stderr)
        return aligned_records

    if to_end <= from_start:
        print(f"WARNING: Target region wraps around — this shouldn't happen after rotation",
              file=sys.stderr)
        return aligned_records

    print(f"  Trimming to gapped positions {from_start}–{to_end} ({to_end - from_start} bp)")

    trimmed = []
    for name, seq in aligned_records:
        region = seq[from_start:to_end]
        trimmed.append((name, region))

    # Remove columns that are 100% gaps
    if not trimmed:
        return trimmed

    aln_len = len(trimmed[0][1])
    keep_cols = []
    for col in range(aln_len):
        if any(rec[1][col] not in "-. " for rec in trimmed):
            keep_cols.append(col)

    if len(keep_cols) < aln_len:
        print(f"  Removed {aln_len - len(keep_cols)} all-gap columns")

    result = []
    for name, seq in trimmed:
        result.append((name, "".join(seq[c] for c in keep_cols)))

    return result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Prepare a trimmed, correctly-oriented MSA from raw mitogenome sequences."
    )
    parser.add_argument("input", help="Input FASTA (unaligned or aligned sequences)")
    parser.add_argument("-g", "--genbank", required=True,
                        help="GenBank file for the reference sequence")
    parser.add_argument("--from", dest="from_gene", required=True,
                        help="Start gene name (e.g. cytb, 12S, Dloop)")
    parser.add_argument("--to", dest="to_gene", required=True,
                        help="End gene name (e.g. 16S, cytb)")
    parser.add_argument("-o", "--output", default="aligned_trimmed.fasta",
                        help="Output FASTA (default: aligned_trimmed.fasta)")
    parser.add_argument("--n-mask", type=int, default=3, metavar="N",
                        help="Replace N-runs >= this length with gaps (default: 3, 0 to disable)")
    args = parser.parse_args()

    out_dir = Path(args.output).parent
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── Step 1: Read input and initial alignment ─────────────────────────────
    print("=" * 60)
    print("  Step 1: Initial MAFFT alignment")
    print("=" * 60)
    raw_records, header_map = read_fasta(args.input)
    print(f"  {len(raw_records)} sequences read")

    # Save originals before masking for N-restoration after alignment
    original_seqs = {name: seq for name, seq in raw_records}
    n_mask_min = args.n_mask

    if n_mask_min > 0:
        raw_records = mask_n_runs(raw_records, min_run=n_mask_min)

    aligned, flipped = run_mafft(raw_records, adjust_direction=True)
    aln_len = len(aligned[0][1])
    print(f"  Aligned: {aln_len} bp")

    # ── Step 2: Identify gene regions ────────────────────────────────────────
    print("\n" + "=" * 60)
    print("  Step 2: Identify gene regions from GenBank")
    print("=" * 60)
    gb_acc, genome_len, gb_seq, genome_features = get_genome_features(args.genbank)
    print(f"  GenBank: {gb_acc} ({genome_len} bp, {len(genome_features)} features)")

    from_name = resolve_gene_name(args.from_gene, genome_features)
    to_name = resolve_gene_name(args.to_gene, genome_features)
    if from_name is None:
        sys.exit(f"ERROR: Gene '{args.from_gene}' not found. Available: "
                 + ", ".join(f["name"] for f in genome_features if f["type"] != "tRNA"))
    if to_name is None:
        sys.exit(f"ERROR: Gene '{args.to_gene}' not found. Available: "
                 + ", ".join(f["name"] for f in genome_features if f["type"] != "tRNA"))
    print(f"  Target range: {from_name} → {to_name}")

    genes_in_range = get_gene_range_circular(from_name, to_name, genome_features)
    print(f"  Genes in range ({len(genes_in_range)}): "
          + " → ".join(g["name"] for g in genes_in_range))

    # ── Step 3: Determine rotation point ─────────────────────────────────────
    print("\n" + "=" * 60)
    print("  Step 3: Determine rotation point")
    print("=" * 60)
    rotation_genome_pos = find_rotation_point(genes_in_range, genome_len)
    print(f"  Rotating to place {from_name} at the start of the alignment")
    print(f"  Rotation genome position: {rotation_genome_pos}")

    # ── Step 4: Rotate sequences ─────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("  Step 4: Rotate sequences")
    print("=" * 60)

    ref_name = find_ref_in_alignment(gb_acc, gb_seq, aligned)
    if ref_name is None:
        sys.exit("ERROR: Could not find GenBank reference in alignment")
    print(f"  Reference: {ref_name}")

    rotations = compute_rotation_per_sequence(
        aligned, ref_name, rotation_genome_pos, gb_seq, genome_len
    )

    rotated_records = []
    for name, seq in aligned:
        ungapped = seq.replace("-", "").replace(".", "")
        rot_amount = rotations.get(name, 0)
        rotated = rotate_sequence(ungapped, rot_amount)
        rotated_records.append((name, rotated))
        if name == ref_name:
            print(f"  {name}: rotated by {rot_amount} bases")

    # ── Step 5: Realign ──────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("  Step 5: Realign with MAFFT")
    print("=" * 60)
    realigned, _ = run_mafft(rotated_records)
    aln_len2 = len(realigned[0][1])
    print(f"  Realigned: {aln_len2} bp")

    # Restore masked Ns
    if n_mask_min > 0:
        realigned = restore_ns(realigned, original_seqs, rotations, flipped, n_mask_min)

    # ── Step 6: Trim to target region ────────────────────────────────────────
    print("\n" + "=" * 60)
    print("  Step 6: Trim to target region")
    print("=" * 60)

    ref_realigned = None
    for name, seq in realigned:
        if name == ref_name:
            ref_realigned = seq
            break

    ref_ungapped = ref_realigned.replace("-", "").replace(".", "")
    offset = _blast_find_offset(ref_ungapped, gb_seq, genome_len)
    if offset is None:
        sys.exit("ERROR: Could not determine offset after realignment")
    print(f"  Genome offset after rotation: {offset}")

    gb_record = SeqIO.read(args.genbank, "genbank")
    mapped_features = map_features_to_alignment(
        gb_record, offset, genome_len, ref_realigned, len(ref_ungapped)
    )
    print(f"  Mapped {len(mapped_features)} features to realigned MSA")

    for f in mapped_features:
        f["name"] = FEATURE_SHORT_NAMES.get(f["name"], f["name"])
        if f["name"].startswith("tRNA-"):
            f["name"] = "t" + f["name"][5:]
        if f.get("type") == "D-loop" or (f["name"] == "?" and "loop" in str(f.get("type", "")).lower()):
            f["name"] = "Dloop"

    trimmed = trim_alignment(realigned, mapped_features, from_name, to_name)

    if not trimmed:
        sys.exit("ERROR: Trimming produced no output")

    final_len = len(trimmed[0][1])
    print(f"  Final alignment: {len(trimmed)} sequences × {final_len} bp")

    # Write output with full original headers restored
    write_fasta(args.output, trimmed, header_map=header_map)
    print(f"\n  Output: {args.output}")
    print("  Done.")


if __name__ == "__main__":
    main()
