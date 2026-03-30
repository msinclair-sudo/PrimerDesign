#!/usr/bin/env python3
"""Step 5: Rank primer pairs by composite score from analysis results.

Reads the JSON output of analyse_primers.py, computes normalised scores
for binding universality, thermodynamic quality, species resolution,
off-target specificity, and amplicon suitability, then writes a ranked CSV.

Usage:
    python scripts/rank_primers.py results.json -o primers_ranked.csv
"""

import argparse
import csv
import json
import sys
from pathlib import Path

import yaml


def load_config(config_path):
    with open(config_path) as f:
        return yaml.safe_load(f)


def load_analysis(json_path):
    with open(json_path) as f:
        return json.load(f)


# ── Scoring functions (each returns 0-1, higher = better) ────────────────────

def score_binding_universality(ps):
    """Fraction of MSA sequences that produce an amplicon."""
    amplicons = ps.get("amplicons", {})
    if not amplicons:
        return 0.0
    total = len(amplicons)
    amplified = sum(1 for v in amplicons.values() if v is not None)
    return amplified / total


def score_thermodynamic_quality(ps, config):
    """Score from primer properties (Tm, GC). Proper NN thermo is in Step 2;
    this uses the Wallace Tm stored in the JSON as a proxy."""
    thermo = config["thermodynamics"]
    tm_min, tm_max = thermo["tm_range"]
    tm_optimal = (tm_min + tm_max) / 2

    props = ps.get("properties", {})
    scores = []

    for direction in ("fwd", "rev"):
        info = props.get(direction, {})
        # Use Wallace Tm as proxy (proper NN Tm is in filtered CSV)
        tm = info.get("tm_wallace")
        if tm is not None:
            if tm_min <= tm <= tm_max:
                scores.append(1.0 - abs(tm - tm_optimal) / ((tm_max - tm_min) / 2))
            else:
                # Penalise proportionally to how far out of range
                dist = min(abs(tm - tm_min), abs(tm - tm_max))
                scores.append(max(0.0, 1.0 - dist / 10))

        gc = info.get("gc_pct")
        if gc is not None:
            gc_min, gc_max = thermo["gc_range"]
            if gc_min <= gc <= gc_max:
                scores.append(1.0)
            else:
                scores.append(0.5)

    # Delta Tm
    fwd_tm = props.get("fwd", {}).get("tm_wallace")
    rev_tm = props.get("rev", {}).get("tm_wallace")
    if fwd_tm is not None and rev_tm is not None:
        delta = abs(fwd_tm - rev_tm)
        scores.append(max(0.0, 1.0 - delta / thermo["max_delta_tm"]))

    # Secondary structure penalties
    ss = props.get("secondary_structure", {})
    for key in ("fwd_hairpin", "rev_hairpin", "fwd_self_dimer", "rev_self_dimer", "hetero_dimer"):
        val = ss.get(key)
        if val is None:
            scores.append(1.0)
        else:
            # Longer complementary runs = worse
            scores.append(max(0.0, 1.0 - len(val) / 8))

    return sum(scores) / len(scores) if scores else 0.5


def score_species_resolution(ps):
    """Score based on amplicon sequence variability across species."""
    variability = ps.get("variability", [])
    if not variability:
        return 0.5

    # Fraction of positions that are variable (>0% difference)
    variable_positions = sum(1 for v in variability if v > 0)
    total = len(variability)
    if total == 0:
        return 0.5

    # More variable positions = better species resolution
    # Cap at 50% variable for full score (typical for good metabarcoding amplicons)
    ratio = variable_positions / total
    return min(1.0, ratio / 0.5)


def score_off_target(ps):
    """Score based on ecoPCR results. For metabarcoding we WANT broad amplification
    of target taxa, so more hits in the target database = better.
    But too many hits with high mismatches = poor specificity."""
    ecopcr = ps.get("ecopcr")
    if not ecopcr:
        return 0.5

    total_seqs = ecopcr.get("db_sequences", 0)
    total_hits = ecopcr.get("total_hits", 0)

    if total_seqs == 0:
        return 0.5

    # Hit ratio — for metabarcoding, higher coverage is better
    hit_ratio = total_hits / total_seqs

    # But penalise if most hits require many mismatches
    mm_counts = ecopcr.get("mismatch_counts", {})
    if mm_counts and total_hits > 0:
        perfect = mm_counts.get("0+0", 0)
        low_mm = sum(v for k, v in mm_counts.items()
                     if sum(int(x) for x in k.split("+")) <= 1)
        quality_ratio = low_mm / total_hits
    else:
        quality_ratio = 0.5

    # Composite: coverage * quality
    return min(1.0, hit_ratio * 0.6 + quality_ratio * 0.4)


def score_amplicon_suitability(ps, config):
    """Score by closeness of amplicon size to optimal."""
    optimal = config["design"]["optimal_amplicon_length"]
    amp_min = config["design"]["min_amplicon_length"]
    amp_max = config["design"]["max_amplicon_length"]
    max_deviation = max(optimal - amp_min, amp_max - optimal)

    amplicons = ps.get("amplicons", {})
    lengths = [a["length"] for a in amplicons.values() if a is not None]

    if not lengths:
        return 0.0

    mean_len = sum(lengths) / len(lengths)

    # Out of range entirely
    if mean_len < amp_min or mean_len > amp_max:
        return 0.0

    deviation = abs(mean_len - optimal)
    return max(0.0, 1.0 - deviation / max_deviation)


# ── Main ranking logic ───────────────────────────────────────────────────────

def rank_primers(analysis_data, config):
    weights = config["ranking"]
    primer_sets = analysis_data.get("primer_sets", [])

    results = []
    for ps in primer_sets:
        s_bind = score_binding_universality(ps)
        s_thermo = score_thermodynamic_quality(ps, config)
        s_resolution = score_species_resolution(ps)
        s_offtarget = score_off_target(ps)
        s_amplicon = score_amplicon_suitability(ps, config)

        composite = (
            weights["binding_universality"] * s_bind
            + weights["thermodynamic_quality"] * s_thermo
            + weights["species_resolution"] * s_resolution
            + weights["off_target_specificity"] * s_offtarget
            + weights["amplicon_suitability"] * s_amplicon
        )

        results.append({
            "name": ps.get("name", "unknown"),
            "forward": ps.get("fwd_primer", ""),
            "reverse": ps.get("rev_primer", ""),
            "composite_score": round(composite, 4),
            "binding_universality": round(s_bind, 4),
            "thermodynamic_quality": round(s_thermo, 4),
            "species_resolution": round(s_resolution, 4),
            "off_target_specificity": round(s_offtarget, 4),
            "amplicon_suitability": round(s_amplicon, 4),
        })

    results.sort(key=lambda r: r["composite_score"], reverse=True)

    for i, r in enumerate(results, 1):
        r["rank"] = i

    return results


def write_ranked_csv(results, output_path):
    if not results:
        print("No primer sets to rank.", file=sys.stderr)
        return

    fieldnames = [
        "rank", "name", "forward", "reverse", "composite_score",
        "binding_universality", "thermodynamic_quality",
        "species_resolution", "off_target_specificity", "amplicon_suitability",
    ]

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    print(f"Wrote {len(results)} ranked primers to {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Rank primer pairs by composite score.")
    parser.add_argument("analysis_json", help="Path to analysis JSON")
    parser.add_argument("-o", "--output", default="primers_ranked.csv", help="Output CSV")
    parser.add_argument("--config", default="primer_config.yaml", type=Path)
    args = parser.parse_args()

    config = load_config(args.config)
    analysis = load_analysis(args.analysis_json)
    results = rank_primers(analysis, config)
    write_ranked_csv(results, args.output)

    # Show top 10
    print(f"\nTop 10:")
    for r in results[:10]:
        print(f"  #{r['rank']} {r['name']}: {r['composite_score']:.3f}  bind={r['binding_universality']:.2f}  thermo={r['thermodynamic_quality']:.2f}  res={r['species_resolution']:.2f}  spec={r['off_target_specificity']:.2f}  amp={r['amplicon_suitability']:.2f}")


if __name__ == "__main__":
    main()
