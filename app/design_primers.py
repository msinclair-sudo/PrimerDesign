#!/usr/bin/env python3
"""
Step 1: Design primers from MSA using DECIPHER (R) or Primer3 (Python fallback).

Reads an MSA FASTA, designs candidate primer pairs, and writes primers.csv.
DECIPHER DesignSignatures is preferred; falls back to primer3-py if R/DECIPHER
is not available.

Usage:
    python scripts/design_primers.py alignment.fasta -o primers.csv
    python scripts/design_primers.py alignment.fasta -o primers.csv --append
    python scripts/design_primers.py alignment.fasta -o primers.csv --method primer3
"""

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import yaml


def load_config(config_path: str) -> dict:
    """Load primer design config from YAML."""
    with open(config_path) as f:
        return yaml.safe_load(f)


def read_msa_fasta(fasta_path: str) -> list[tuple[str, str]]:
    """Read aligned sequences from FASTA. Returns list of (name, sequence)."""
    sequences = []
    name = None
    seq_parts = []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    sequences.append((name, "".join(seq_parts)))
                name = line[1:].split()[0]
                seq_parts = []
            elif line:
                seq_parts.append(line.upper())
    if name is not None:
        sequences.append((name, "".join(seq_parts)))
    return sequences


def consensus_from_msa(sequences: list[tuple[str, str]]) -> str:
    """Build a simple majority-rule consensus from aligned sequences, removing gap columns."""
    if not sequences:
        raise ValueError("No sequences in MSA")
    aln_len = len(sequences[0][1])
    consensus = []
    for i in range(aln_len):
        counts: dict[str, int] = {}
        for _, seq in sequences:
            base = seq[i] if i < len(seq) else "-"
            if base not in ("-", "."):
                counts[base] = counts.get(base, 0) + 1
        if counts:
            consensus.append(max(counts, key=counts.get))
    return "".join(consensus)


# ---------------------------------------------------------------------------
# DECIPHER (R subprocess)
# ---------------------------------------------------------------------------

def check_decipher_available() -> bool:
    """Check if R and DECIPHER are available."""
    if not shutil.which("Rscript"):
        return False
    result = subprocess.run(
        ["Rscript", "-e", "library(DECIPHER); cat('ok')"],
        capture_output=True, text=True, timeout=30,
    )
    return result.returncode == 0 and "ok" in result.stdout


R_DESIGN_SCRIPT = r"""
suppressPackageStartupMessages({
  library(DECIPHER)
  library(jsonlite)
  library(RSQLite)
})

args <- commandArgs(trailingOnly = TRUE)
fasta_file   <- args[1]
config_json  <- args[2]
output_csv   <- args[3]

config <- fromJSON(config_json)

# Load MSA into in-memory DECIPHER database
dbConn <- dbConnect(SQLite(), ":memory:")
Seqs2DB(fasta_file, "FASTA", dbConn, "target_seqs")
n_seqs <- dbGetQuery(dbConn, "SELECT COUNT(*) FROM Seqs")[[1]]
cat("Loaded", n_seqs, "sequences\n")

# DesignSignatures for maximal amplicon diversity
n_cores <- max(1, min(parallel::detectCores() - 1, 4))
cat("Using", n_cores, "cores\n")

sigs <- DesignSignatures(
  dbConn,
  type        = "sequence",
  minProductSize = config$min_amplicon_length,
  maxProductSize = config$max_amplicon_length,
  resolution  = config$resolution,
  levels      = 5,
  processors  = n_cores
)

dbDisconnect(dbConn)

if (nrow(sigs) == 0) {
  cat("DesignSignatures returned 0 results\n")
  write.csv(data.frame(name=character(), forward=character(), reverse=character()),
            output_csv, row.names=FALSE)
  quit(status = 0)
}

cat("DesignSignatures returned", nrow(sigs), "primer pairs\n")

# Build output dataframe
out <- data.frame(
  name    = paste0("DECIPHER_", sprintf("%03d", seq_len(nrow(sigs)))),
  forward = sigs$forward_primer,
  reverse = sigs$reverse_primer,
  stringsAsFactors = FALSE
)

write.csv(out, output_csv, row.names = FALSE)
cat("Wrote", nrow(out), "primers to", output_csv, "\n")
"""


def design_with_decipher(fasta_path: str, config: dict) -> list[dict]:
    """Run DECIPHER DesignSignatures via R subprocess."""
    design_cfg = config.get("design", {})
    decipher_cfg = config.get("decipher", {})

    r_config = {
        "min_amplicon_length": design_cfg.get("min_amplicon_length", 90),
        "max_amplicon_length": design_cfg.get("max_amplicon_length", 180),
        "resolution": decipher_cfg.get("resolution", 5),
    }

    with tempfile.TemporaryDirectory() as tmpdir:
        r_script = os.path.join(tmpdir, "design.R")
        config_json = os.path.join(tmpdir, "config.json")
        output_csv = os.path.join(tmpdir, "primers.csv")

        with open(r_script, "w") as f:
            f.write(R_DESIGN_SCRIPT)
        with open(config_json, "w") as f:
            json.dump(r_config, f)

        result = subprocess.run(
            ["Rscript", r_script, os.path.abspath(fasta_path), config_json, output_csv],
            capture_output=True, text=True, timeout=600,
        )

        if result.returncode != 0:
            print(f"DECIPHER stderr:\n{result.stderr}", file=sys.stderr)
            raise RuntimeError(f"DECIPHER R script failed (exit {result.returncode})")

        if result.stdout:
            print(result.stdout, end="")

        primers = []
        if os.path.exists(output_csv):
            with open(output_csv) as f:
                reader = csv.DictReader(f)
                for row in reader:
                    primers.append({
                        "name": row["name"],
                        "forward": row["forward"],
                        "reverse": row["reverse"],
                    })
    return primers


# ---------------------------------------------------------------------------
# Primer3 fallback
# ---------------------------------------------------------------------------

def design_with_primer3(fasta_path: str, config: dict) -> list[dict]:
    """Design primers using primer3-py on the consensus sequence."""
    import primer3 as p3

    sequences = read_msa_fasta(fasta_path)
    consensus = consensus_from_msa(sequences)

    if len(consensus) < 50:
        raise ValueError(f"Consensus too short ({len(consensus)} bp) for primer design")

    design_cfg = config.get("design", {})
    thermo_cfg = config.get("thermodynamics", {})
    tm_range = thermo_cfg.get("tm_range", [58, 64])
    gc_range = thermo_cfg.get("gc_range", [40, 60])

    size_range = [
        design_cfg.get("min_amplicon_length", 90),
        design_cfg.get("max_amplicon_length", 180),
    ]

    result = p3.bindings.design_primers(
        seq_args={"SEQUENCE_TEMPLATE": consensus},
        global_args={
            "PRIMER_NUM_RETURN": 20,
            "PRIMER_PRODUCT_SIZE_RANGE": size_range,
            "PRIMER_MIN_SIZE": design_cfg.get("min_primer_length", 18),
            "PRIMER_MAX_SIZE": design_cfg.get("max_primer_length", 25),
            "PRIMER_MIN_TM": tm_range[0],
            "PRIMER_OPT_TM": (tm_range[0] + tm_range[1]) / 2,
            "PRIMER_MAX_TM": tm_range[1],
            "PRIMER_MIN_GC": gc_range[0],
            "PRIMER_MAX_GC": gc_range[1],
        },
    )

    n_returned = result.get("PRIMER_PAIR_NUM_RETURNED", 0)
    primers = []
    for i in range(n_returned):
        fwd = result.get(f"PRIMER_LEFT_{i}_SEQUENCE", "")
        rev = result.get(f"PRIMER_RIGHT_{i}_SEQUENCE", "")
        if fwd and rev:
            primers.append({
                "name": f"Primer3_{i + 1:03d}",
                "forward": fwd,
                "reverse": rev,
            })

    return primers


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def write_primers_csv(primers: list[dict], output_path: str, append: bool = False):
    """Write primers to CSV. If append=True, add to existing file."""
    existing = []
    if append and os.path.exists(output_path):
        with open(output_path) as f:
            reader = csv.DictReader(f)
            existing = list(reader)

    all_primers = existing + primers

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["name", "forward", "reverse"])
        writer.writeheader()
        writer.writerows(all_primers)

    return len(all_primers)


def main():
    parser = argparse.ArgumentParser(
        description="Design primers from MSA (DECIPHER or Primer3)"
    )
    parser.add_argument("msa", help="MSA FASTA file")
    parser.add_argument("-o", "--output", default="primers.csv",
                        help="Output CSV path (default: primers.csv)")
    parser.add_argument("--config", default="primer_config.yaml",
                        help="Config YAML path")
    parser.add_argument("--method", choices=["decipher", "primer3", "auto"],
                        default="auto",
                        help="Design method (default: auto — tries DECIPHER first)")
    parser.add_argument("--append", action="store_true",
                        help="Append to existing primers.csv instead of overwriting")
    args = parser.parse_args()

    if not os.path.exists(args.msa):
        sys.exit(f"Error: MSA file not found: {args.msa}")

    config = load_config(args.config)

    method = args.method
    if method == "auto":
        if check_decipher_available():
            method = "decipher"
            print("Using DECIPHER (R) for primer design")
        else:
            method = "primer3"
            print("DECIPHER not available, falling back to Primer3")

    if method == "decipher":
        primers = design_with_decipher(args.msa, config)
    else:
        primers = design_with_primer3(args.msa, config)

    if not primers:
        print("Warning: no primers were designed", file=sys.stderr)

    total = write_primers_csv(primers, args.output, append=args.append)
    print(f"Wrote {len(primers)} new primers ({total} total) to {args.output}")


if __name__ == "__main__":
    main()
