#!/usr/bin/env python3
"""
Build a FASTA MSA from UGENE BLAST annotation exports and downloaded GenBank files.

Extracts hit subsequences using BLAST coordinates, then aligns with MAFFT.

Usage:
    python blast_to_msa.py [directory] [--output OUTPUT.fasta]
"""

import argparse
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO


def find_annotations_file(directory):
    """Find the UGENE BLAST annotations .gb file (no ORIGIN block)."""
    for path in directory.glob("*.gb"):
        text = path.read_text()
        if "ORIGIN" not in text and "misc_feature" in text:
            return path
    return None


def parse_annotations(path):
    """Parse UGENE BLAST annotation file. Returns (query_length, list of HSP dicts)."""
    text = path.read_text()

    # Query length from LOCUS line: "region [START END]"
    locus_m = re.search(r"region\s*\[(\d+)\s+(\d+)\]", text)
    if not locus_m:
        sys.exit(f"ERROR: Could not parse region from LOCUS line in {path}")
    query_length = int(locus_m.group(2)) - int(locus_m.group(1)) + 1

    # Split into misc_feature blocks
    blocks = re.split(r"(?=     misc_feature\s)", text)
    blocks = [b for b in blocks if b.strip().startswith("misc_feature")]

    hsps = []
    for block in blocks:
        loc_m = re.search(r"misc_feature\s+(\d+)\.\.(\d+)", block)
        acc_m = re.search(r'/accession="(\S+)"', block)
        hf_m = re.search(r"/hit-from=(\d+)", block)
        ht_m = re.search(r"/hit-to=(\d+)", block)
        ident_m = re.search(r'/identities="([^"]+)"', block)
        frame_m = re.search(r'/source_frame="(\w+)"', block)
        hitlen_m = re.search(r'/hit_len=(\d+)', block)
        # /def may span multiple lines
        def_m = re.search(r'/def="(.*?)"', block, re.DOTALL)

        if not (loc_m and acc_m and hf_m and ht_m):
            continue

        # Species: first two words of /def, cleaned up
        species = ""
        if def_m:
            def_text = re.sub(r"\s+", " ", def_m.group(1)).strip()
            words = def_text.split()
            if len(words) >= 2:
                species = f"{words[0]}_{words[1]}"

        hit_from = int(hf_m.group(1))
        hit_to = int(ht_m.group(1))
        is_reverse = (frame_m and frame_m.group(1) != "direct") or (hit_from > hit_to)

        hsps.append({
            "accession": acc_m.group(1),
            "query_start": int(loc_m.group(1)),
            "query_end": int(loc_m.group(2)),
            "hit_from": hit_from,
            "hit_to": hit_to,
            "identity_str": ident_m.group(1) if ident_m else "?",
            "species": species,
            "is_reverse": is_reverse,
            "hit_len": int(hitlen_m.group(1)) if hitlen_m else None,
        })

    return query_length, hsps


def group_by_accession(hsps):
    """Group HSPs by accession. Returns dict: accession -> {species, hsps}."""
    grouped = defaultdict(lambda: {"species": "", "hsps": []})
    for h in hsps:
        acc = h["accession"]
        grouped[acc]["species"] = h["species"]
        grouped[acc]["hsps"].append(h)
    return grouped


def extract_circular(seq, start_0based, length):
    """Extract `length` bases starting at `start_0based`, wrapping if circular."""
    n = len(seq)
    return "".join(str(seq[(start_0based + i) % n]) for i in range(length))


def rev_comp(seq_str):
    """Return the reverse complement of a DNA string."""
    comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq_str.translate(comp)[::-1]


def build_reference(ref_accession, ref_group, query_length, gb_records):
    """Build the ungapped reference sequence using self-hit coordinates + circular extraction.
    Also returns the genome offset for annotation mapping."""
    seq = gb_records[ref_accession].seq
    first_hsp = sorted(ref_group["hsps"], key=lambda h: h["query_start"])[0]
    # Both hit_from and query_start are 1-based, so the difference is the 0-based offset
    offset = first_hsp["hit_from"] - first_hsp["query_start"]
    return extract_circular(seq, offset, query_length), offset


def build_annotation_track(ref_accession, gb_records, query_length, genome_offset, aligned_ref):
    """Build an annotation row showing gene/feature boundaries mapped through the alignment.

    Returns a string the same length as the aligned reference, with feature names
    bracketed like >12S-------12S< and dashes elsewhere.
    """
    rec = gb_records[ref_accession]
    genome_len = len(rec.seq)

    # Collect features in query space (0-based query coords)
    feature_spans = []
    for feat in rec.features:
        if feat.type in ("source", "gene"):
            continue
        fs = int(feat.location.start)   # 0-based genome start
        fe = int(feat.location.end)     # 0-based exclusive genome end

        # Convert to 0-based query coordinates (circular)
        qs = (fs - genome_offset) % genome_len
        qe = (fe - 1 - genome_offset) % genome_len  # inclusive end

        # Handle features that wrap around the circular origin
        # (qs > qe means the feature crosses the origin boundary)
        if qs > qe:
            # Only include if at least part falls within query region
            if qs < query_length:
                # Portion from qs to end of query
                qe_clamped = query_length - 1
                feature_end = qe_clamped
            elif qe < query_length:
                # Portion from start to qe
                qs = 0
                feature_end = qe
            else:
                continue
        else:
            feature_end = qe

        # Skip features entirely outside the query region
        if qs >= query_length and feature_end >= query_length:
            continue
        # Clamp to query bounds
        if feature_end >= query_length:
            feature_end = query_length - 1

        name = feat.qualifiers.get("product", feat.qualifiers.get("gene", ["?"]))[0]
        # Use feature type as fallback for unlabeled features
        if name == "?" and feat.type == "D-loop":
            name = "Dloop"
        # Shorten common names
        short = {
            "s-rRNA": "12S", "l-rRNA": "16S",
            "cytochrome b": "cytb", "CYTB": "cytb",
        }
        name = short.get(name, name)
        # tRNAs: use 3-letter amino acid
        if name.startswith("tRNA-"):
            name = name.replace("tRNA-", "t")

        feature_spans.append((qs, feature_end, name))

    # Build annotation in query space (pre-alignment, length = query_length)
    anno_query = ["-"] * query_length
    for qs, qe, name in feature_spans:
        span_len = qe - qs + 1
        if span_len < 1:
            continue
        # Place >NAME at start, NAME< at end, fill middle with -
        tag_open = f">{name}"
        tag_close = f"{name}<"

        if span_len >= len(tag_open) + len(tag_close):
            # Full tags fit
            for i, c in enumerate(tag_open):
                anno_query[qs + i] = c
            for i, c in enumerate(tag_close):
                anno_query[qe - len(tag_close) + 1 + i] = c
        elif span_len >= len(name):
            # Just the name
            for i, c in enumerate(name):
                anno_query[qs + i] = c
        else:
            # Too short, just mark with name chars that fit
            for i in range(span_len):
                if i < len(name):
                    anno_query[qs + i] = name[i]

    # Map query positions through the alignment gaps
    # aligned_ref has gaps inserted by MAFFT; we map each query base to its
    # alignment column, then place annotation characters accordingly
    aln_len = len(aligned_ref)
    anno_aln = ["-"] * aln_len

    qi = 0  # index into query (ungapped)
    for ai in range(aln_len):
        if aligned_ref[ai] != "-":
            if qi < query_length:
                anno_aln[ai] = anno_query[qi]
            qi += 1

    return "".join(anno_aln)


def extract_hit_sequence(acc_group, gb_records):
    """Extract and concatenate hit subsequences for one accession.

    For multi-HSP hits, the HSPs cover different parts of the query region.
    We concatenate them in query-coordinate order so MAFFT can align the
    combined sequence against the reference (MAFFT naturally inserts gaps
    where data is missing between HSPs).

    Handles circular genomes and reverse-complement hits.
    Returns (sequence_str, gap_bp) where gap_bp is the total unsequenced
    gap between HSPs in query coordinates.
    """
    acc = acc_group["hsps"][0]["accession"]
    if acc not in gb_records:
        return None, 0
    seq = gb_records[acc].seq

    sorted_hsps = sorted(acc_group["hsps"], key=lambda h: h["query_start"])

    fragments = []
    gap_bp = 0
    prev_query_end = None
    for hsp in sorted_hsps:
        h_from = hsp["hit_from"]
        h_to = hsp["hit_to"]
        is_reverse = hsp["is_reverse"]

        # Track gaps between HSPs in query space
        if prev_query_end is not None:
            gap = hsp["query_start"] - prev_query_end - 1
            if gap > 0:
                gap_bp += gap
        prev_query_end = hsp["query_end"]

        if is_reverse:
            # Reverse-complement hit: hit_from > hit_to
            start_0 = h_to - 1  # 0-based, smaller coordinate
            length = h_from - h_to + 1
            fragment = extract_circular(seq, start_0, length)
            fragment = rev_comp(fragment)
        else:
            # Forward hit: use circular extraction for wrap-around safety
            start_0 = h_from - 1  # 0-based
            length = h_to - h_from + 1
            fragment = extract_circular(seq, start_0, length)

        fragments.append(fragment)

    return "".join(fragments), gap_bp


def run_mafft(sequences, names):
    """Run MAFFT on a list of (name, sequence) pairs. Returns aligned sequences."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tmp_in:
        for name, seq in zip(names, sequences):
            tmp_in.write(f">{name}\n{seq}\n")
        tmp_in_path = tmp_in.name

    result = subprocess.run(
        ["mafft", "--auto", "--thread", "-1", tmp_in_path],
        capture_output=True, text=True
    )

    if result.returncode != 0:
        print("MAFFT stderr:", result.stderr, file=sys.stderr)
        sys.exit("ERROR: MAFFT failed")

    # Parse MAFFT output
    aligned = {}
    current_name = None
    current_seq = []
    for line in result.stdout.splitlines():
        if line.startswith(">"):
            if current_name:
                aligned[current_name] = "".join(current_seq)
            current_name = line[1:].strip()
            current_seq = []
        else:
            current_seq.append(line.strip())
    if current_name:
        aligned[current_name] = "".join(current_seq)

    # Clean up
    Path(tmp_in_path).unlink(missing_ok=True)

    return aligned


def compute_identity(ref_seq, hit_seq):
    """Compute identity between ref and hit in aligned (non-gap) positions."""
    matches = 0
    aligned = 0
    for r, h in zip(ref_seq, hit_seq):
        if r != "-" and h != "-":
            aligned += 1
            if r == h:
                matches += 1
    return matches, aligned


def identify_reference(grouped):
    """Find the reference accession (100% identity self-hit).

    Parses the identity fraction (e.g. '1633/1633 (100%)') to check
    for exact matches rather than relying on substring search.
    """
    for acc, grp in grouped.items():
        for hsp in grp["hsps"]:
            ident = hsp["identity_str"]
            # Try to parse "N/N (100%)" format
            m = re.match(r"(\d+)/(\d+)", ident)
            if m and m.group(1) == m.group(2):
                return acc
            # Fallback: check for "100%" anywhere
            if "100%" in ident:
                return acc
    return None


def main():
    parser = argparse.ArgumentParser(description="UGENE BLAST annotations to FASTA MSA")
    parser.add_argument("directory", nargs="?", default=".", help="Directory with .gb files")
    parser.add_argument("--output", "-o", default="blast_msa_output.fasta", help="Output FASTA file")
    args = parser.parse_args()

    directory = Path(args.directory)
    output_path = Path(args.output)

    # Find annotations file
    ann_path = find_annotations_file(directory)
    if not ann_path:
        sys.exit("ERROR: No BLAST annotations file found (expected .gb without ORIGIN block)")
    print(f"Annotations file: {ann_path.name}")

    # Parse annotations
    query_length, hsps = parse_annotations(ann_path)
    print(f"Query length: {query_length} bp")
    print(f"Total HSPs parsed: {len(hsps)}")

    # Report any reverse-complement hits
    rev_count = sum(1 for h in hsps if h["is_reverse"])
    if rev_count:
        print(f"  Note: {rev_count} reverse-complement HSP(s) detected")

    # Group by accession
    grouped = group_by_accession(hsps)
    print(f"Unique accessions: {len(grouped)}")

    # Load all GenBank sequence files
    gb_records = {}
    for path in directory.glob("*.gb"):
        if path == ann_path:
            continue
        try:
            rec = SeqIO.read(path, "genbank")
            acc = rec.id.split(".")[0]
            gb_records[acc] = rec
        except ValueError as e:
            print(f"  WARNING: Skipping {path.name} — {e}")
        except Exception as e:
            print(f"  WARNING: Could not parse {path.name} — {type(e).__name__}: {e}")
    print(f"GenBank records loaded: {len(gb_records)}")

    # Identify reference (100% identity self-hit)
    ref_accession = identify_reference(grouped)

    if not ref_accession:
        # Show available identities to help debug
        print("\n  Available identity values:")
        for acc, grp in grouped.items():
            for hsp in grp["hsps"]:
                print(f"    {acc}: {hsp['identity_str']}")
        sys.exit("ERROR: No 100% identity self-hit found — cannot identify reference")
    if ref_accession not in gb_records:
        sys.exit(f"ERROR: Reference accession {ref_accession} has no GenBank file")

    ref_species = grouped[ref_accession]["species"]
    print(f"Reference: {ref_accession} ({ref_species})")

    # Build reference sequence (circular extraction, ungapped)
    ref_seq, genome_offset = build_reference(ref_accession, grouped[ref_accession], query_length, gb_records)
    assert len(ref_seq) == query_length, f"Reference length {len(ref_seq)} != {query_length}"

    # Extract unaligned hit sequences
    names = [f"{ref_accession}_{ref_species}_REF"]
    sequences = [ref_seq]

    for acc, grp in sorted(grouped.items()):
        if acc == ref_accession:
            continue
        if acc not in gb_records:
            print(f"  WARNING: No GenBank file for {acc} — skipping")
            continue
        hit_seq, gap_bp = extract_hit_sequence(grp, gb_records)
        if hit_seq is None:
            continue
        label = f"{acc}_{grp['species']}"
        names.append(label)
        sequences.append(hit_seq)
        total_hsps = len(grp["hsps"])
        gap_info = f" ({gap_bp} bp gap between HSPs)" if gap_bp > 0 else ""
        print(f"  {label}: {len(hit_seq)} bp extracted ({total_hsps} HSP{'s' if total_hsps > 1 else ''}){gap_info}")

    # Align with MAFFT
    print(f"\nRunning MAFFT on {len(sequences)} sequences...")
    aligned = run_mafft(sequences, names)

    # Build annotation track mapped through the alignment
    ref_name = names[0]
    ref_aligned = aligned[ref_name]
    anno_track = build_annotation_track(
        ref_accession, gb_records, query_length, genome_offset, ref_aligned
    )

    # Write output FASTA (annotation first, then reference, then hits)
    with open(output_path, "w") as f:
        f.write(f">anno\n{anno_track}\n")
        f.write(f">{ref_name}\n{aligned[ref_name]}\n")
        for name in names[1:]:
            if name in aligned:
                f.write(f">{name}\n{aligned[name]}\n")

    # Summary
    aln_length = len(ref_aligned)

    print(f"\n{'='*65}")
    print(f"Output: {output_path}")
    print(f"Sequences: {len(aligned)}")
    print(f"Alignment length: {aln_length} bp")
    print(f"{'='*65}")
    print(f"{'Sequence':<45} {'Bases':>7} {'Gaps':>7} {'Identity':>10}")
    print(f"{'-'*45} {'-'*7} {'-'*7} {'-'*10}")

    for name in [ref_name] + names[1:]:
        seq = aligned.get(name, "")
        bases = sum(1 for c in seq if c != "-")
        gaps = aln_length - bases
        matches, compared = compute_identity(ref_aligned, seq)
        pct = f"{matches/compared*100:.1f}%" if compared > 0 else "N/A"
        print(f"{name:<45} {bases:>7} {gaps:>7} {pct:>10}")


if __name__ == "__main__":
    main()
