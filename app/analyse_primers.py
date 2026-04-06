#!/usr/bin/env python3
"""
analyse_primers.py
==================
Analyses primer pairs against an MSA and optional off-target database.

All primer binding analysis uses obipcr (OBITools4) — no custom fuzzy search.
All thermodynamic properties use primer3-py — no Wallace rule approximations.

Computes:
- Per-sequence primer binding via obipcr (in silico PCR)
- Pairwise identity matrix and sliding-window identity (alignment-level)
- Per-position amplicon variability
- Off-target screening via obipcr against broader database

Usage:
    python scripts/analyse_primers.py alignment.fasta --primers primers.csv
    python scripts/analyse_primers.py alignment.fasta --primers primers.csv --database seqs/mam_groups/
"""

import argparse
import csv
import hashlib
import json
import multiprocessing
import os
import re
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import primer3
from Bio import SeqIO

# ── CONFIG ────────────────────────────────────────────────────────────────────
INPUT_FASTA = "blast_msa_output.fasta"
OUTPUT_JSON = "blast_vis_data.json"

WINDOW = 400
STEP = 100
MIN_SITES = 30

# Feature type → colour mapping for gene track
FEATURE_COLORS = {
    "tRNA": "#e9c46a",
    "rRNA_12S": "#57cc99",
    "rRNA_16S": "#a8dadc",
    "rRNA": "#57cc99",
    "CDS": "#4cc9f0",
    "D-loop": "#f4a261",
    "misc_feature": "#b0b0b0",
}

# Standard feature name shortening
FEATURE_SHORT_NAMES = {
    "s-rRNA": "12S", "l-rRNA": "16S",
    "cytochrome b": "cytb", "CYTB": "cytb",
    "cytochrome c oxidase subunit I": "COX1",
    "cytochrome c oxidase subunit II": "COX2",
    "cytochrome c oxidase subunit III": "COX3",
    "NADH dehydrogenase subunit 1": "ND1",
    "NADH dehydrogenase subunit 2": "ND2",
    "NADH dehydrogenase subunit 3": "ND3",
    "NADH dehydrogenase subunit 4": "ND4",
    "NADH dehydrogenase subunit 4L": "ND4L",
    "NADH dehydrogenase subunit 5": "ND5",
    "NADH dehydrogenase subunit 6": "ND6",
    "ATP synthase F0 subunit 6": "ATP6",
    "ATP synthase F0 subunit 8": "ATP8",
}


# ── GENBANK ANNOTATION MAPPING ──────────────────────────────────────────────
def map_annotations(genbank_path: str, seqs_ungapped: dict, seqs_gapped: dict,
                    seq_names: list, cache_path: str | None = None) -> list[dict]:
    """Map GenBank features onto MSA coordinates using BLAST with doubled reference.

    Returns a list of feature dicts: [{name, start, end, color}, ...]
    where start/end are gapped MSA coordinates.
    """
    # Check cache
    if cache_path:
        cache = Path(cache_path)
        if cache.exists():
            cache_data = json.loads(cache.read_text())
            # Validate cache by checking checksums
            gb_hash = _file_hash(genbank_path)
            if cache_data.get("genbank_hash") == gb_hash:
                print("      Using cached annotations")
                return cache_data["features"]

    # Parse GenBank
    gb_record = SeqIO.read(genbank_path, "genbank")
    gb_seq = str(gb_record.seq).upper()
    genome_len = len(gb_seq)
    gb_acc = gb_record.id.split(".")[0]
    print(f"      GenBank: {gb_acc} ({genome_len} bp)")

    # Find matching MSA sequence
    ref_name = _find_matching_sequence(gb_acc, gb_seq, seqs_ungapped, seq_names)
    if not ref_name:
        print("      WARNING: No matching sequence found in MSA for GenBank record")
        return []
    print(f"      Matched to MSA sequence: {ref_name}")

    ref_ungapped = seqs_ungapped[ref_name]
    ref_gapped = seqs_gapped[ref_name]

    # BLAST ungapped MSA sequence against doubled GenBank reference
    offset = _blast_find_offset(ref_ungapped, gb_seq, genome_len)
    if offset is None:
        print("      WARNING: BLAST could not map MSA to genome")
        return []
    print(f"      Genome offset: {offset} (MSA starts at genome pos {offset})")

    # Extract features and map to MSA coordinates
    features = _extract_and_map_features(
        gb_record, offset, genome_len, len(ref_ungapped), ref_gapped
    )
    print(f"      Mapped {len(features)} features to MSA coordinates")

    # Cache
    if cache_path:
        cache_data = {
            "genbank_hash": _file_hash(genbank_path),
            "genbank_accession": gb_acc,
            "genome_length": genome_len,
            "msa_offset": offset,
            "matched_sequence": ref_name,
            "features": features,
        }
        Path(cache_path).parent.mkdir(parents=True, exist_ok=True)
        Path(cache_path).write_text(json.dumps(cache_data, indent=2))

    return features


def _file_hash(path: str) -> str:
    """MD5 hash of file contents for cache validation."""
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def _find_matching_sequence(gb_acc: str, gb_seq: str, seqs_ungapped: dict,
                            seq_names: list) -> str | None:
    """Find which MSA sequence matches the GenBank record."""
    # Try matching by accession in the sequence name
    for name in seq_names:
        if gb_acc in name:
            return name

    # Fallback: find the MSA sequence that is a subsequence of the GenBank sequence
    for name in seq_names:
        ungapped = seqs_ungapped[name].upper()
        if ungapped in gb_seq or ungapped in (gb_seq + gb_seq):
            return name

    return None


def _blast_find_offset(msa_seq: str, gb_seq: str, genome_len: int) -> int | None:
    """Use BLAST to find where the MSA sequence sits in the genome.
    Uses doubled reference to handle circular genomes."""
    doubled = gb_seq + gb_seq

    with tempfile.TemporaryDirectory(prefix="blast_map_") as tmpdir:
        tmpdir = Path(tmpdir)

        # Write doubled reference
        ref_path = tmpdir / "ref.fasta"
        ref_path.write_text(f">doubled_ref\n{doubled}\n")

        # Write query (ungapped MSA sequence)
        query_path = tmpdir / "query.fasta"
        query_path.write_text(f">msa_query\n{msa_seq}\n")

        # Make BLAST database
        subprocess.run(
            ["makeblastdb", "-in", str(ref_path), "-dbtype", "nucl",
             "-out", str(tmpdir / "ref_db")],
            capture_output=True, check=True,
        )

        # Run BLAST
        result = subprocess.run(
            ["blastn", "-query", str(query_path), "-db", str(tmpdir / "ref_db"),
             "-outfmt", "6 sstart send qstart qend length pident",
             "-max_target_seqs", "1", "-evalue", "1e-10"],
            capture_output=True, text=True, check=True,
        )

    if not result.stdout.strip():
        return None

    # Parse best hit — take the one with longest alignment
    best = None
    for line in result.stdout.strip().splitlines():
        parts = line.split("\t")
        sstart, send = int(parts[0]), int(parts[1])
        length = int(parts[4])
        if best is None or length > best[2]:
            best = (sstart, send, length)

    if best is None:
        return None

    # Convert doubled-reference position back to real genome position
    offset = (best[0] - 1) % genome_len  # 0-based

    # Verify: extract region from genome and compare to MSA sequence
    extracted = "".join(
        gb_seq[(offset + i) % genome_len] for i in range(len(msa_seq))
    )
    mismatches = sum(a != b for a, b in zip(msa_seq.upper(), extracted.upper()))
    if mismatches > 0:
        pct = round(mismatches / len(msa_seq) * 100, 2)
        if pct > 1.0:
            print(f"      WARNING: {mismatches} mismatches ({pct}%) between MSA sequence and GenBank at offset {offset}")
            print(f"      Annotation positions may be inaccurate")
        else:
            print(f"      Verified: {mismatches} mismatches ({pct}%) — minor, likely sequencing variants")
    else:
        print(f"      Verified: exact match at offset {offset}")

    return offset


def _extract_and_map_features(gb_record, offset: int, genome_len: int,
                               msa_ungapped_len: int, ref_gapped: str) -> list[dict]:
    """Extract GenBank features and map to gapped MSA coordinates."""
    # Build ungapped→gapped coordinate map
    ungap_to_gap = {}
    ungap_pos = 0
    for gap_pos, ch in enumerate(ref_gapped):
        if ch != "-":
            ungap_to_gap[ungap_pos] = gap_pos
            ungap_pos += 1

    features = []
    for feat in gb_record.features:
        if feat.type in ("source", "gene"):
            continue

        feat_start = int(feat.location.start)  # 0-based
        feat_end = int(feat.location.end)       # 0-based exclusive

        # Get feature name
        name = feat.qualifiers.get("product", feat.qualifiers.get("gene", ["?"]))[0]
        if name == "?" and feat.type == "D-loop":
            name = "Dloop"
        name = FEATURE_SHORT_NAMES.get(name, name)
        if name.startswith("tRNA-"):
            name = "t" + name[5:]

        # Determine colour
        if "12S" in name or "s-rRNA" in name:
            color = FEATURE_COLORS["rRNA_12S"]
        elif "16S" in name or "l-rRNA" in name:
            color = FEATURE_COLORS["rRNA_16S"]
        elif feat.type == "tRNA":
            color = FEATURE_COLORS["tRNA"]
        elif feat.type == "CDS":
            color = FEATURE_COLORS["CDS"]
        elif feat.type == "D-loop" or "loop" in name.lower():
            color = FEATURE_COLORS["D-loop"]
        elif feat.type == "rRNA":
            color = FEATURE_COLORS["rRNA"]
        else:
            color = FEATURE_COLORS.get(feat.type, "#b0b0b0")

        # Convert genome coordinates to MSA-region coordinates (circular-aware)
        msa_start = (feat_start - offset) % genome_len
        msa_end = (feat_end - offset) % genome_len

        # For trimmed alignments (MSA covers only part of the genome), reject
        # features that map into range via modular wrap but aren't actually present.
        # A feature is truly in-range if its linear distance from the offset
        # (unwrapped) falls within the MSA ungapped length.
        if msa_ungapped_len < genome_len * 0.95:
            linear_dist = (feat_start - offset) % genome_len
            if linear_dist >= msa_ungapped_len and (feat_end - offset) % genome_len >= msa_ungapped_len:
                continue

        # Handle feature wrapping around in MSA space
        if msa_start > msa_end:
            # Feature crosses the boundary — check if either part is in our region
            if msa_start < msa_ungapped_len:
                msa_end_clamped = min(msa_ungapped_len, genome_len)
                # Map to gapped coordinates
                gap_s = ungap_to_gap.get(msa_start, None)
                gap_e = ungap_to_gap.get(min(msa_end_clamped - 1, msa_ungapped_len - 1), None)
                if gap_s is not None and gap_e is not None:
                    features.append({"name": name, "start": gap_s, "end": gap_e + 1, "color": color})
            if msa_end < msa_ungapped_len and msa_end > 0:
                gap_s = ungap_to_gap.get(0, 0)
                gap_e = ungap_to_gap.get(min(msa_end - 1, msa_ungapped_len - 1), None)
                if gap_e is not None:
                    features.append({"name": name, "start": gap_s, "end": gap_e + 1, "color": color})
            continue

        # Normal case — feature doesn't wrap
        if msa_start >= msa_ungapped_len or msa_end <= 0:
            continue  # Feature outside MSA region

        # Clamp to MSA bounds
        msa_start = max(0, msa_start)
        msa_end = min(msa_end, msa_ungapped_len)

        # Map to gapped coordinates
        gap_s = ungap_to_gap.get(msa_start, None)
        gap_e = ungap_to_gap.get(msa_end - 1, None)
        if gap_s is not None and gap_e is not None:
            features.append({"name": name, "start": gap_s, "end": gap_e + 1, "color": color})

    # Sort by start position
    features.sort(key=lambda f: f["start"])
    return features


# ── IUPAC resolution for primer3-py (which requires ACGT only) ───────────────
IUPAC_TO_FIRST = {
    "R": "A", "Y": "C", "S": "G", "W": "A", "K": "G", "M": "A",
    "B": "C", "D": "A", "H": "A", "V": "A", "N": "A",
}


def resolve_iupac(seq: str) -> str:
    """Replace IUPAC ambiguity codes with a concrete base for primer3."""
    return "".join(IUPAC_TO_FIRST.get(b, b) for b in seq.upper())


# ── PRIMER PROPERTIES VIA PRIMER3-PY ────────────────────────────────────────
def compute_primer_properties(fwd: str, rev: str) -> dict:
    """Compute thermodynamic properties using primer3-py (nearest-neighbor)."""
    fwd_resolved = resolve_iupac(fwd)
    rev_resolved = resolve_iupac(rev)

    fwd_tm = primer3.calc_tm(fwd_resolved)
    rev_tm = primer3.calc_tm(rev_resolved)
    fwd_hairpin = primer3.calc_hairpin(fwd_resolved).dg / 1000.0
    rev_hairpin = primer3.calc_hairpin(rev_resolved).dg / 1000.0
    fwd_homodimer = primer3.calc_homodimer(fwd_resolved).dg / 1000.0
    rev_homodimer = primer3.calc_homodimer(rev_resolved).dg / 1000.0
    heterodimer = primer3.calc_heterodimer(fwd_resolved, rev_resolved).dg / 1000.0

    def gc_pct(s):
        s = s.upper()
        return round((s.count("G") + s.count("C")) / len(s) * 100, 1) if s else 0.0

    return {
        "fwd": {
            "seq": fwd,
            "length": len(fwd),
            "gc_pct": gc_pct(fwd),
            "tm": round(fwd_tm, 1),
            "hairpin_dg": round(fwd_hairpin, 2),
            "homodimer_dg": round(fwd_homodimer, 2),
        },
        "rev": {
            "seq": rev,
            "length": len(rev),
            "gc_pct": gc_pct(rev),
            "tm": round(rev_tm, 1),
            "hairpin_dg": round(rev_hairpin, 2),
            "homodimer_dg": round(rev_homodimer, 2),
        },
        "delta_tm": round(abs(fwd_tm - rev_tm), 1),
        "heterodimer_dg": round(heterodimer, 2),
    }


def compute_flags(props: dict, thermo_cfg: dict) -> list[str]:
    """Compute thermodynamic flags for a primer pair, matching validate_thermodynamics logic."""
    flags = []
    fwd = props["fwd"]
    rev = props["rev"]
    tm_lo, tm_hi = thermo_cfg.get("tm_range", [58, 64])
    gc_lo, gc_hi = thermo_cfg.get("gc_range", [40, 60])
    clamp_lo, clamp_hi = thermo_cfg.get("gc_clamp_3prime", [1, 3])
    max_delta_tm = thermo_cfg.get("max_delta_tm", 5)
    max_hairpin = thermo_cfg.get("max_hairpin_dg", -3.0)
    max_homodimer = thermo_cfg.get("max_homodimer_dg", -9.0)
    max_heterodimer = thermo_cfg.get("max_heterodimer_dg", -9.0)
    max_homopoly = thermo_cfg.get("max_homopolymer", 4)

    def gc_clamp_count(seq, window=5):
        return sum(1 for b in seq.upper()[-window:] if b in "GC")

    def max_homopolymer(seq):
        if not seq:
            return 0
        max_run = run = 1
        for i in range(1, len(seq)):
            if seq[i].upper() == seq[i-1].upper():
                run += 1
                if run > max_run:
                    max_run = run
            else:
                run = 1
        return max_run

    if not (tm_lo <= fwd["tm"] <= tm_hi):
        flags.append("tm_fwd_out_of_range")
    if not (tm_lo <= rev["tm"] <= tm_hi):
        flags.append("tm_rev_out_of_range")
    if props["delta_tm"] > max_delta_tm:
        flags.append("delta_tm_too_high")
    if fwd["hairpin_dg"] < max_hairpin:
        flags.append("hairpin_fwd")
    if rev["hairpin_dg"] < max_hairpin:
        flags.append("hairpin_rev")
    if fwd["homodimer_dg"] < max_homodimer:
        flags.append("homodimer_fwd")
    if rev["homodimer_dg"] < max_homodimer:
        flags.append("homodimer_rev")
    if props["heterodimer_dg"] < max_heterodimer:
        flags.append("heterodimer")
    if not (gc_lo <= fwd["gc_pct"] <= gc_hi):
        flags.append("gc_fwd_out_of_range")
    if not (gc_lo <= rev["gc_pct"] <= gc_hi):
        flags.append("gc_rev_out_of_range")

    fwd_clamp = gc_clamp_count(fwd["seq"])
    rev_clamp = gc_clamp_count(rev["seq"])
    if not (clamp_lo <= fwd_clamp <= clamp_hi):
        flags.append("gc_clamp_fwd")
    if not (clamp_lo <= rev_clamp <= clamp_hi):
        flags.append("gc_clamp_rev")

    fwd_hp = max_homopolymer(fwd["seq"])
    rev_hp = max_homopolymer(rev["seq"])
    if fwd_hp > max_homopoly:
        flags.append("homopolymer_fwd")
    if rev_hp > max_homopoly:
        flags.append("homopolymer_rev")

    return flags


# ── ALIGNMENT-LEVEL FUNCTIONS (no primers involved) ─────────────────────────
def pairwise_id(s1: str, s2: str, min_sites: int = 1) -> float | None:
    m = t = 0
    for a, b in zip(s1, s2):
        if a != "-" and b != "-":
            t += 1
            if a.upper() == b.upper():
                m += 1
    if t < min_sites:
        return None
    return round(m / t * 100, 2)


def sliding_window_id(ref: str, query: str, window: int, step: int, min_sites: int):
    results = []
    x_pos = []
    for i in range(0, len(ref) - window + 1, step):
        pid = pairwise_id(ref[i:i + window], query[i:i + window], min_sites)
        results.append(pid)
        x_pos.append(i + window // 2)
    return x_pos, results


def per_position_variability(ref_amp: str, other_amps: list[str]) -> list[float]:
    scores = []
    for i in range(len(ref_amp)):
        diffs = sum(1 for s in other_amps if i < len(s) and s[i].upper() != ref_amp[i].upper())
        scores.append(round(diffs / len(other_amps) * 100, 1))
    return scores


# ── PARSE FASTA ──────────────────────────────────────────────────────────────
def parse_fasta(path: str):
    seqs_gapped = {}
    seqs_ungapped = {}
    current = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">anno") or line.startswith(">tGlu"):
                current = None
                continue
            if line.startswith(">"):
                current = line[1:]
                seqs_gapped[current] = ""
                seqs_ungapped[current] = ""
            elif current:
                seqs_gapped[current] += line
                seqs_ungapped[current] += line.replace("-", "")
    return seqs_gapped, seqs_ungapped


def build_coord_maps(gapped_seq: str):
    ungap_to_gap = {}
    gap_to_ungap = {}
    ungap_pos = 0
    for g_pos, ch in enumerate(gapped_seq):
        if ch != "-":
            ungap_to_gap[ungap_pos] = g_pos
            gap_to_ungap[g_pos] = ungap_pos
            ungap_pos += 1
    return ungap_to_gap, gap_to_ungap


# ── CSV PRIMER LOADING ──────────────────────────────────────────────────────
def load_primers_csv(path: str) -> list[dict]:
    primers = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            primers.append({
                "name": row["name"].strip(),
                "forward": row["forward"].strip().upper(),
                "reverse": row["reverse"].strip().upper(),
            })
    return primers


# ── OBIPCR WRAPPER ──────────────────────────────────────────────────────────
def run_obipcr(fwd: str, rev: str, fasta_input, max_mm: int = 3,
               min_len: int = 50, max_len: int = 500,
               circular: bool = False) -> dict | None:
    """Run obipcr against FASTA file(s). fasta_input can be a file path,
    directory, or list of file paths. Returns parsed results dict."""
    if not shutil.which("obipcr"):
        print("      WARNING: obipcr not found")
        return None

    # Resolve input files
    if isinstance(fasta_input, (list, tuple)):
        fasta_files = [Path(f) for f in fasta_input]
    else:
        db = Path(fasta_input)
        if db.is_dir():
            fasta_files = sorted(db.glob("*.fasta")) + sorted(db.glob("*.fa")) + sorted(db.glob("*.fas"))
        else:
            fasta_files = [db]

    if not fasta_files:
        return None

    # Count total sequences
    total_seqs = 0
    for f in fasta_files:
        with open(f) as fh:
            total_seqs += sum(1 for line in fh if line.startswith(">"))

    cmd = [
        "obipcr", "--fasta",
        "--forward", fwd,
        "--reverse", rev,
        "-e", str(max_mm),
        "-l", str(min_len),
        "-L", str(max_len),
        "--fasta-output",
        "--no-progressbar",
    ]
    if circular:
        cmd.append("--circular")

    try:
        cat_proc = subprocess.Popen(
            ["cat"] + [str(f) for f in fasta_files],
            stdout=subprocess.PIPE,
        )
        result = subprocess.run(
            cmd, stdin=cat_proc.stdout,
            capture_output=True, text=True, timeout=300,
        )
        cat_proc.wait()
    except (subprocess.TimeoutExpired, Exception) as e:
        print(f"      WARNING: obipcr failed — {e}")
        return None

    # Parse output — both headers and sequences
    hits = []
    current_header = None
    current_seq_parts = []

    for line in result.stdout.splitlines():
        if line.startswith(">"):
            # Save previous entry
            if current_header is not None:
                _save_hit(current_header, "".join(current_seq_parts), hits)
            current_header = line
            current_seq_parts = []
        elif current_header is not None:
            current_seq_parts.append(line.strip())

    # Save last entry
    if current_header is not None:
        _save_hit(current_header, "".join(current_seq_parts), hits)

    # Build summary
    mm_counts = {}
    amp_lengths = []
    for h in hits:
        key = f"{h['fwd_mm']}+{h['rev_mm']}"
        mm_counts[key] = mm_counts.get(key, 0) + 1
        amp_lengths.append(h["amp_len"])

    return {
        "db_sequences": total_seqs,
        "db_files": [f.name for f in fasta_files],
        "total_hits": len(hits),
        "mismatch_counts": mm_counts,
        "amplicon_length_min": min(amp_lengths) if amp_lengths else 0,
        "amplicon_length_max": max(amp_lengths) if amp_lengths else 0,
        "amplicon_length_median": sorted(amp_lengths)[len(amp_lengths) // 2] if amp_lengths else 0,
        "hits": hits,
    }


def _save_hit(header: str, sequence: str, hits: list):
    """Parse an obipcr FASTA entry and append to hits list."""
    json_m = re.search(r"\{.*\}", header)
    id_m = re.match(r">(\S+?)_sub\[(\d+)\.\.(\d+)\]", header)
    if not json_m:
        return

    meta = json.loads(json_m.group())
    seq_id = id_m.group(1) if id_m else header.split()[0][1:]
    start = int(id_m.group(2)) if id_m else 0
    end = int(id_m.group(3)) if id_m else 0

    hits.append({
        "id": seq_id,
        "definition": meta.get("definition", ""),
        "fwd_mm": meta.get("forward_error", 0),
        "rev_mm": meta.get("reverse_error", 0),
        "amp_len": end - start,
        "start": start,
        "end": end,
        "direction": meta.get("direction", "forward"),
        "sequence": sequence,
        "fwd_match": meta.get("forward_match", ""),
        "rev_match": meta.get("reverse_match", ""),
    })


# ── ANALYSE ONE PRIMER SET VIA OBIPCR ───────────────────────────────────────
def analyse_primer_set(
    primer_name: str,
    fwd_primer: str,
    rev_primer: str,
    msa_fasta_path: str,
    seq_names: list,
    ref_name: str,
    ungap_to_gap: dict,
    all_ungap_to_gap: dict | None = None,
    on_target_mm: int = 2,
    amp_range: tuple = (10, 500),
    thermo_cfg: dict | None = None,
) -> dict:
    """Run primer analysis using obipcr for binding + primer3-py for properties."""

    # Primer properties via primer3-py (nearest-neighbor thermodynamics)
    primer_props = compute_primer_properties(fwd_primer, rev_primer)

    # Compute flags if thermo config provided
    if thermo_cfg:
        flags = compute_flags(primer_props, thermo_cfg)
        primer_props["flags"] = flags
        primer_props["status"] = "FLAG" if flags else "PASS"
    else:
        primer_props["flags"] = []
        primer_props["status"] = "PASS"

    # In silico PCR via obipcr against MSA sequences
    obipcr_result = run_obipcr(
        fwd_primer, rev_primer, msa_fasta_path,
        max_mm=on_target_mm, min_len=amp_range[0], max_len=amp_range[1], circular=False,
    )

    # Build per-sequence hit data from obipcr results
    hits = {}
    amplicons_dict = {}
    for name in seq_names:
        hits[name] = {"fwd_hits": [], "rev_hits": [], "amplicons": []}
        amplicons_dict[name] = None

    if obipcr_result:
        for h in obipcr_result["hits"]:
            # Match obipcr hit ID back to MSA sequence name
            matched_name = None
            for name in seq_names:
                if h["id"] in name or name.startswith(h["id"]):
                    matched_name = name
                    break
            if not matched_name:
                continue

            amp = {
                "start": h["start"],
                "end": h["end"],
                "length": h["amp_len"],
                "sequence": h["sequence"],
                "inner": h["sequence"][len(fwd_primer):-len(rev_primer)] if len(h["sequence"]) > len(fwd_primer) + len(rev_primer) else "",
                "fwd_mm": h["fwd_mm"],
                "rev_mm": h["rev_mm"],
            }
            hits[matched_name]["amplicons"].append(amp)
            if amplicons_dict[matched_name] is None:
                amplicons_dict[matched_name] = amp

    # Map primer coords to gapped alignment
    # Prefer reference; fall back to any amplified sequence if reference didn't amplify
    coord_amp = amplicons_dict.get(ref_name)
    coord_gap_map = ungap_to_gap
    coord_source = ref_name

    if coord_amp is None and all_ungap_to_gap:
        # Find the sequence with the fewest total mismatches
        best_name, best_mm = None, 999
        for name in seq_names:
            amp = amplicons_dict.get(name)
            if amp is not None and name in all_ungap_to_gap:
                mm = amp["fwd_mm"] + amp["rev_mm"]
                if mm < best_mm:
                    best_mm = mm
                    best_name = name
        if best_name:
            coord_amp = amplicons_dict[best_name]
            coord_gap_map = all_ungap_to_gap[best_name]
            coord_source = best_name

    if coord_amp:
        fwd_ug_s = coord_amp["start"]
        rev_ug_s = coord_amp["end"] - len(rev_primer)

        fwd_gap_s = coord_gap_map.get(fwd_ug_s, fwd_ug_s)
        fwd_gap_e = coord_gap_map.get(fwd_ug_s + len(fwd_primer) - 1, fwd_gap_s + len(fwd_primer)) + 1
        rev_gap_s = coord_gap_map.get(rev_ug_s, rev_ug_s)
        rev_gap_e = coord_gap_map.get(rev_ug_s + len(rev_primer) - 1, rev_gap_s + len(rev_primer)) + 1

        primer_props["gapped_coords"] = {
            "fwd_start": fwd_gap_s,
            "fwd_end": fwd_gap_e,
            "rev_start": rev_gap_s,
            "rev_end": rev_gap_e,
            "amp_start": fwd_gap_s,
            "amp_end": rev_gap_e,
        }
        if coord_source != ref_name:
            print(f"      Coords via: {coord_source} (ref did not amplify)")
        print(f"      FWD gapped: {fwd_gap_s}–{fwd_gap_e}")
        print(f"      REV gapped: {rev_gap_s}–{rev_gap_e}")
        print(f"      Amplicon:   {rev_gap_e - fwd_gap_s} gapped bp")
    else:
        print("      WARNING: No amplicon found in any sequence — check primers")
        primer_props["gapped_coords"] = {}

    # Per-position variability across amplicons
    ref_amp_seq = (amplicons_dict.get(ref_name) or {}).get("sequence", "")
    other_amp_seqs = [
        (amplicons_dict.get(n) or {}).get("sequence", "")
        for n in seq_names[1:]
        if amplicons_dict.get(n)
    ]
    variability = (
        per_position_variability(ref_amp_seq.upper(), [s.upper() for s in other_amp_seqs])
        if ref_amp_seq and other_amp_seqs
        else []
    )

    return {
        "name": primer_name,
        "fwd_primer": fwd_primer,
        "rev_primer": rev_primer,
        "properties": primer_props,
        "hits": hits,
        "amplicons": amplicons_dict,
        "variability": variability,
    }


# ── PARALLEL WORKER FUNCTIONS (module-level for pickling) ─────────────────────
def _analyse_worker(args):
    """Worker for parallel MSA binding analysis."""
    name, fwd, rev, msa_path, seq_names, ref_name, ungap_to_gap, all_ungap_to_gap, on_target_mm, amp_range, thermo_cfg = args
    return analyse_primer_set(name, fwd, rev, msa_path, seq_names, ref_name, ungap_to_gap, all_ungap_to_gap, on_target_mm, amp_range, thermo_cfg)


def _ecopcr_worker(args):
    """Worker for parallel ecoPCR off-target screening."""
    fwd_primer, rev_primer, database, off_target_mm, eco_min, eco_max, eco_circular = args
    return run_obipcr(
        fwd_primer, rev_primer, database,
        max_mm=off_target_mm, min_len=eco_min, max_len=eco_max, circular=eco_circular,
    )


# ── MAIN ──────────────────────────────────────────────────────────────────────
def load_config(config_path: str | None) -> dict:
    """Load analysis config from YAML. Exits if config is missing."""
    import yaml
    if config_path is None:
        print("ERROR: --config is required. No config file specified.", file=sys.stderr)
        sys.exit(1)
    p = Path(config_path)
    if not p.exists():
        print(f"ERROR: Config not found at {p}", file=sys.stderr)
        sys.exit(1)
    with open(p) as f:
        raw = yaml.safe_load(f)
    mm = raw.get("mismatches", {})
    eco = raw.get("ecopcr", {})
    design = raw.get("design", {})
    thermo = raw.get("thermodynamics", {})
    return {
        "on_target_mm": mm.get("on_target", 2),
        "off_target_mm": mm.get("off_target", 3),
        "amp_min": design.get("min_amplicon_length", 90),
        "amp_max": design.get("max_amplicon_length", 300),
        "eco_min": eco.get("min_amplicon_length", 10),
        "eco_max": eco.get("max_amplicon_length", 500),
        "eco_circular": eco.get("circular", True),
        "thermo": thermo,
    }


def main(fasta_path: str, primers_csv: str | None = None,
         output_path: str = OUTPUT_JSON, database: str | None = None,
         genbank_path: str | None = None, config_path: str | None = None):
    workers = int(os.environ.get("PRIMER_WORKERS", max(1, multiprocessing.cpu_count() - 1)))
    print(f"Using {workers} workers for binding analysis")

    cfg = load_config(config_path)
    on_target_mm = cfg["on_target_mm"]
    off_target_mm = cfg["off_target_mm"]
    amp_range = (cfg["amp_min"], cfg["amp_max"])
    eco_min = cfg["eco_min"]
    eco_max = cfg["eco_max"]
    eco_circular = cfg["eco_circular"]
    thermo_cfg = cfg.get("thermo", {})

    has_db = database is not None
    has_gb = genbank_path is not None
    n_steps = 5 + (1 if has_db else 0) + (1 if has_gb else 0)
    step = 0

    step += 1
    print(f"[{step}/{n_steps}] Parsing FASTA: {fasta_path}")
    seqs_gapped, seqs_ungapped = parse_fasta(fasta_path)
    seq_names = list(seqs_gapped.keys())
    align_len = len(next(iter(seqs_gapped.values())))
    print(f"      {len(seq_names)} sequences · {align_len} bp alignment")

    ref_name = seq_names[0]
    ref_gapped = seqs_gapped[ref_name]
    ungap_to_gap, _ = build_coord_maps(ref_gapped)

    # Build gap maps for all sequences (used as fallback when ref doesn't amplify)
    all_ungap_to_gap = {}
    for name, gapped in seqs_gapped.items():
        all_ungap_to_gap[name], _ = build_coord_maps(gapped)

    # ── Gene annotations from GenBank
    genes = []
    if has_gb:
        step += 1
        print(f"[{step}/{n_steps}] Mapping GenBank annotations ...")
        cache_path = str(Path(output_path).parent / "annotations.json")
        genes = map_annotations(
            genbank_path, seqs_ungapped, seqs_gapped, seq_names, cache_path
        )
    if not genes:
        print("      No gene annotations (provide --genbank for annotation track)")

    # Write ungapped MSA sequences to temp file for obipcr
    msa_tmp = tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=False, prefix="msa_ungapped_"
    )
    for name, seq in seqs_ungapped.items():
        msa_tmp.write(f">{name}\n{seq}\n")
    msa_tmp.close()
    msa_fasta_path = msa_tmp.name
    print(f"      Wrote {len(seqs_ungapped)} ungapped sequences to temp file")

    # ── Pairwise identity matrix
    step += 1
    print(f"[{step}/{n_steps}] Computing pairwise identity matrix ...")
    matrix = []
    for n1 in seq_names:
        row = []
        for n2 in seq_names:
            row.append(pairwise_id(seqs_gapped[n1], seqs_gapped[n2]) or 0.0)
        matrix.append(row)

    id_vs_ref = {n: pairwise_id(ref_gapped, seqs_gapped[n]) for n in seq_names}

    # ── Sliding window
    step += 1
    print(f"[{step}/{n_steps}] Sliding window identity vs reference ...")
    sw_x, _ = sliding_window_id(ref_gapped, ref_gapped, WINDOW, STEP, MIN_SITES)
    sw_results = {}
    for n in seq_names[1:]:
        _, vals = sliding_window_id(ref_gapped, seqs_gapped[n], WINDOW, STEP, MIN_SITES)
        sw_results[n] = vals
    print(f"      {len(sw_x)} windows")

    # ── Load primer pairs
    if primers_csv:
        primer_list = load_primers_csv(primers_csv)
        step += 1
        print(f"[{step}/{n_steps}] Loaded {len(primer_list)} primer set(s) from {primers_csv}")
    else:
        sys.exit("ERROR: --primers CSV is required. See primers.csv for format.")

    # ── Analyse each primer set via obipcr
    if workers > 1 and len(primer_list) > 1:
        print(f"\n  Running {len(primer_list)} primer sets in parallel ({workers} workers) ...")
        worker_args = [
            (p["name"], p["forward"], p["reverse"], msa_fasta_path, seq_names, ref_name, ungap_to_gap, all_ungap_to_gap, on_target_mm, amp_range, thermo_cfg)
            for p in primer_list
        ]
        with ProcessPoolExecutor(max_workers=workers) as pool:
            primer_sets = list(pool.map(_analyse_worker, worker_args))
        for i, (p, result) in enumerate(zip(primer_list, primer_sets), 1):
            n_amp = sum(1 for a in result["amplicons"].values() if a is not None)
            print(f"  ── Primer set {i}/{len(primer_list)}: {p['name']} — {n_amp}/{len(seq_names)} amplify")
    else:
        primer_sets = []
        for i, p in enumerate(primer_list, 1):
            print(f"\n  ── Primer set {i}/{len(primer_list)}: {p['name']} ──")
            print(f"      FWD: {p['forward']}")
            print(f"      REV: {p['reverse']}")
            result = analyse_primer_set(
                p["name"], p["forward"], p["reverse"],
                msa_fasta_path, seq_names, ref_name, ungap_to_gap,
                all_ungap_to_gap, on_target_mm, amp_range, thermo_cfg,
            )
            primer_sets.append(result)

    # Clean up temp file
    Path(msa_fasta_path).unlink(missing_ok=True)

    # ── ecoPCR against broader database
    if has_db:
        viable = [(i, ps) for i, ps in enumerate(primer_sets)
                  if any(a is not None for a in ps["amplicons"].values())]
        skipped = len(primer_sets) - len(viable)
        step += 1
        print(f"\n[{step}/{n_steps}] Running ecoPCR (obipcr) against {database} ...")
        print(f"      {len(viable)} pairs with MSA amplicons, skipping {skipped} non-amplifying pairs")

        if workers > 1 and len(viable) > 1:
            print(f"      Running {len(viable)} ecoPCR jobs in parallel ({workers} workers) ...")
            ecopcr_args = [
                (ps["fwd_primer"], ps["rev_primer"], database, off_target_mm, eco_min, eco_max, eco_circular)
                for _, ps in viable
            ]
            with ProcessPoolExecutor(max_workers=workers) as pool:
                ecopcr_results = list(pool.map(_ecopcr_worker, ecopcr_args))
            for vi, ((i, ps), ecopcr_result) in enumerate(zip(viable, ecopcr_results)):
                ps["ecopcr"] = ecopcr_result
                if ecopcr_result:
                    print(f"  ── ecoPCR {vi+1}/{len(viable)}: {ps['name']} — "
                          f"{ecopcr_result['total_hits']}/{ecopcr_result['db_sequences']} sequences amplify")
                else:
                    print(f"  ── ecoPCR {vi+1}/{len(viable)}: {ps['name']} — not available or failed")
        else:
            for vi, (i, ps) in enumerate(viable):
                print(f"\n  ── ecoPCR {vi+1}/{len(viable)}: {ps['name']} ──")
                ecopcr_result = run_obipcr(
                    ps["fwd_primer"], ps["rev_primer"], database,
                    max_mm=off_target_mm, min_len=eco_min, max_len=eco_max, circular=eco_circular,
                )
                ps["ecopcr"] = ecopcr_result
                if ecopcr_result:
                    print(f"      {ecopcr_result['total_hits']}/{ecopcr_result['db_sequences']} sequences amplify")
                    if ecopcr_result["total_hits"] > 0:
                        print(f"      Amplicon range: {ecopcr_result['amplicon_length_min']}–{ecopcr_result['amplicon_length_max']} bp")
                        for mm_key, count in sorted(ecopcr_result["mismatch_counts"].items(),
                                                     key=lambda x: -x[1])[:5]:
                            print(f"      {mm_key} mm: {count} hits")
                else:
                    print("      ecoPCR not available or failed")

    # ── Assemble output
    step += 1
    print(f"\n[{step}/{n_steps}] Writing output JSON ...")
    output = {
        "meta": {
            "fasta_file": Path(fasta_path).name,
            "n_sequences": len(seq_names),
            "alignment_length": align_len,
            "reference": ref_name,
            "n_primer_sets": len(primer_sets),
            "window": WINDOW,
            "step": STEP,
            "ecopcr_database": database if has_db else None,
        },
        "seq_names": seq_names,
        "gapped_sequences": {name: seqs_gapped[name] for name in seq_names},
        "genes": genes,
        "matrix": matrix,
        "id_vs_ref": id_vs_ref,
        "sliding_window": {
            "x_positions": sw_x,
            "data": sw_results,
        },
        "primer_sets": primer_sets,
    }

    with open(output_path, "w") as fh:
        json.dump(output, fh, indent=2)

    print(f"\nDone. Output written to {output_path}")

    # Quick terminal summary
    for ps in primer_sets:
        n_amp = sum(1 for a in ps["amplicons"].values() if a is not None)
        print(f"\n── {ps['name']} ({n_amp}/{len(seq_names)} amplify) ──")
        print(f"   FWD: {ps['fwd_primer']}  REV: {ps['rev_primer']}")
        print(f"   Tm: {ps['properties']['fwd']['tm']}°C / {ps['properties']['rev']['tm']}°C  ΔTm: {ps['properties']['delta_tm']}°C")

    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="MSA primer binding analysis (uses obipcr + primer3-py)")
    parser.add_argument("fasta", nargs="?", default=INPUT_FASTA, help="Input MSA FASTA")
    parser.add_argument("--primers", "-p", required=True, help="CSV file with primer pairs")
    parser.add_argument("--output", "-o", default=OUTPUT_JSON, help="Output JSON file")
    parser.add_argument("--database", "-d", default=None, help="FASTA dir/file for ecoPCR off-target screening")
    parser.add_argument("--genbank", "-g", default=None,
                        help="GenBank file (.gb) for gene annotation mapping. Must match at least one MSA sequence.")
    parser.add_argument("--config", "-c", default=None,
                        help="Config YAML for mismatch tolerances and ecoPCR parameters")
    args = parser.parse_args()
    main(args.fasta, args.primers, args.output, args.database, args.genbank, args.config)
