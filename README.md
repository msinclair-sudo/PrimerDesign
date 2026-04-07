# eDNA Primer Design Pipeline

A pipeline for designing and validating metabarcoding primers from a multiple sequence alignment. It identifies candidate primer binding sites, validates their thermodynamic properties, tests binding across all sequences in the alignment, screens for off-target amplification against a broader sequence database, and ranks the results — producing an interactive report for primer selection.

---

## Usage

Design new primers from an alignment:

```bash
python Design.py -i alignment.fasta -g reference.gb
```

With off-target screening against a sequence database:

```bash
python Design.py -i alignment.fasta -g reference.gb -d path/to/fasta_files
```

Evaluate existing primers instead of designing new ones:

```bash
python Design.py -i alignment.fasta -g reference.gb -p my_primers.csv
```

When `-p` is supplied, the design step is skipped and the supplied primers are used as the starting point. They pass through the same validation, analysis, and ranking as designed candidates. This makes the tool dual-purpose — it can design new primers or evaluate existing ones through the same analytical pipeline. The primers CSV must have columns `name`, `forward`, `reverse`.

Custom output directory:

```bash
python Design.py -i alignment.fasta -g reference.gb -o results/my_run
```

Run with the included example data (11-species rodent mitogenome panel):

```bash
python Design.py -i data/example/input/rodent_aligned.fasta -g data/example/input/NC_006914.1.gb
```

### Arguments

| Flag              | Description                                                         |
| ----------------- | ------------------------------------------------------------------- |
| `-i / --input`    | Input MSA FASTA file                                                |
| `-o / --output`   | Output directory (default: `data/output`)                           |
| `-d / --database` | FASTA directory or file for off-target screening (optional)         |
| `-g / --genbank`  | GenBank file for gene annotation mapping (auto-detected if omitted) |
| `-p / --primers`  | CSV of primers to evaluate — skips design step (optional)           |
| `-c / --config`   | Config YAML file (default: `data/config.yaml`)                      |
| `--workers`       | Number of parallel workers (default: auto, env: PRIMER_WORKERS)     |

The database flag accepts a directory containing multiple FASTA files, or a single FASTA file. If `--genbank` is not specified, the pipeline looks for a `.gb` file in the same directory as the input alignment and uses it automatically.

### Alignment Preparation

Prepare a trimmed, correctly-oriented MSA from raw mitogenome sequences. Handles circular genome artifacts by rotating all sequences to center the target region, realigning with MAFFT, then trimming to the specified gene range. Long runs of N characters (from low-coverage assembly regions) are masked to gap characters before alignment and restored afterwards — see N-masking below.

```bash
python Design.py align -i raw_sequences.fasta -g reference.gb --from cytb --to 16S -o aligned_trimmed.fasta
```

| Flag              | Description                                              |
| ----------------- | -------------------------------------------------------- |
| `-i / --input`    | Input FASTA (unaligned or aligned sequences)             |
| `-g / --genbank`  | GenBank file for the reference sequence (required)       |
| `--from`          | Start gene name (e.g. cytb, 12S, Dloop)                 |
| `--to`            | End gene name (e.g. 16S, cytb)                           |
| `-o / --output`   | Output FASTA (default: `aligned_trimmed.fasta`)          |
| `--n-mask N`      | Replace N-runs >= N with gaps before alignment (default: 3, 0 to disable) |

### Debug Mode

Individual pipeline steps can be invoked separately using debug mode:

```bash
python Design.py debug analyse -i alignment.fasta -d path/to/fasta_files
python Design.py debug rank
python Design.py debug report
```

Available steps: `design`, `validate`, `expand`, `analyse`, `rank`, `report`.

---

## Inputs

The pipeline requires one file, and optionally accepts up to three more:

- **MSA FASTA** (required) — A multiple sequence alignment of the target region. This can come from any alignment tool (MAFFT, MUSCLE, ClustalO, or manual curation). The pipeline does not construct alignments — it accepts them as-is. The first sequence in the file is treated as the reference.

- **GenBank file (.gb)** (recommended) — A GenBank-format file containing the full genome (or full mitogenome) for at least one of the sequences in the MSA. This is used to annotate the gene track in the report — mapping features like coding regions, rRNAs, tRNAs, and the D-loop onto the alignment. The GenBank record must match at least one MSA sequence by accession or sequence identity. If placed in the same directory as the input alignment, it is detected automatically.

- **Off-target database** (optional) — A directory containing FASTA files of sequences from non-target taxa, or a single FASTA file. This is used for specificity screening to identify whether the designed primers would amplify unintended organisms.

- **Pre-made primers** (optional) — A CSV file with columns `name`, `forward`, `reverse` containing existing primer pairs to evaluate alongside the computationally designed candidates.

## Outputs

All outputs are written to the output directory (default `data/output/`):

- **primers/** — CSV files at each stage: initial candidates, individually validated primers, expanded combinations.
- **analysis/** — JSON file containing the complete analysis results (binding data, identity matrices, ecoPCR hits, annotations). Also contains the cached annotation mapping and the final ranked primers CSV.
- **report.html** — An interactive HTML dashboard for visually comparing and selecting primer pairs.

---

## Configuration

All pipeline parameters are controlled through `data/config.yaml`. A custom config can be specified with `-c/--config`. The config file is required — the pipeline will not run without it.

**header** — Report title fields: main text, coloured species name, subtitle, and optional suffix. These populate the header bar in the HTML report.

**design** — Primer length range, amplicon size range, and optimal amplicon size. The amplicon size range controls which primer combinations are considered viable in the expansion and analysis steps.

**decipher** — Parameters for the DECIPHER DesignSignatures algorithm: minimum coverage (fraction of MSA sequences a primer must bind), resolution (k-mer size for amplicon differentiation), maximum degeneracy (permutations per primer; 1 = no degenerate bases), `num_primer_sets` (number of final forward-reverse pairs returned per window — the top-scoring combinations), `search_primers` (number of individual candidate primer sequences DECIPHER evaluates per window — controls how broadly it explores binding sites before forming pairs), and sliding window parameters (`window_size`, `window_overlap`) for scanning long alignments.

**thermodynamics** — Tm range, maximum delta Tm, GC content range, free energy thresholds for hairpin, homodimer, and heterodimer formation, maximum homopolymer run length, and 3-prime GC clamp requirements. Tm and hairpin thresholds are hard-reject (primers removed); all others are soft-flag (kept with warning).

**mismatches** — Separate mismatch tolerances for on-target binding (how forgiving the MSA analysis is) and off-target screening (how broadly the ecoPCR specificity search looks).

**ecopcr** — Parameters for off-target screening: amplicon size range and whether to treat database sequences as circular.

**ranking** — Weights for each of the five scoring components. These can be adjusted to prioritise different aspects of primer quality depending on the application. Weights must sum to 1.0.

---

## Dependencies

The pipeline requires the following tools, all installable via conda:

- **DECIPHER** (R/Bioconductor) — MSA-aware primer design using DesignSignatures
- **Primer3 / primer3-py** — Nearest-neighbour thermodynamic calculations
- **OBITools4 (obipcr)** — In silico PCR for binding analysis and off-target screening
- **BLAST+** — Sequence alignment for GenBank annotation mapping
- **Biopython** — GenBank file parsing
- **MAFFT** — Multiple sequence alignment (for use outside this pipeline)
- **PyYAML** — Configuration file parsing

Install all dependencies:

```bash
conda env create -f environment.yml
conda activate primer
```

---

## Pipeline Methods

### Alignment Preparation — N-Masking and Restoration

Genome assemblies from low-coverage sequencing often contain long stretches of N characters representing regions where bases could not be called. If left in the sequences, these N-runs cause MAFFT to treat them as real characters and attempt to align them, distorting the surrounding alignment and introducing artificial gaps.

The alignment preparation step handles this with a mask-align-restore strategy:

1. **Mask:** Before alignment, runs of N characters at or above a configurable threshold (default: 3 consecutive Ns) are replaced with gap characters. Short isolated Ns (below the threshold) are left as-is since they may represent genuine single-base ambiguities rather than missing data.

2. **Align:** MAFFT strips gap characters from input sequences before aligning, so the masked N-runs are effectively removed. MAFFT aligns only the real bases (plus any short Ns below the threshold), producing a cleaner alignment without N-driven distortion.

3. **Restore:** After the final alignment, masked Ns are restored in their correct positions. The pipeline tracks which ungapped positions in each original sequence contained masked N-runs, accounting for any reverse-complementation (from MAFFT's `--adjustdirection`) and circular rotation applied during the pipeline. For each sequence, alignment gaps that correspond to originally-masked N positions are converted back to N characters. This preserves the alignment structure while correctly representing ambiguous regions as N rather than as absence (gaps).

The threshold is configurable via `--n-mask` (default 3, set to 0 to disable). The restoration step reports how many N bases were restored and in how many sequences.

### Step 1 — Primer Design

The pipeline identifies candidate primer binding sites within conserved regions of the input alignment using DECIPHER's DesignSignatures algorithm (Wright, 2016). DECIPHER is an R/Bioconductor package specifically designed for designing PCR primers from multiple sequence alignments. It uses a signature-based approach that maximises the diversity of amplicon sequences (to allow species discrimination) while ensuring the primer binding sites are conserved across the aligned sequences (to allow universal amplification).

The alignment is imported into a DECIPHER sequence database, and DesignSignatures is called with configurable constraints on primer length, amplicon size, GC content, and maximum degeneracy. The number of primer sets returned and the breadth of the initial candidate search are also configurable — increasing these values causes DECIPHER to explore primer binding sites across all gene regions in the alignment rather than concentrating only on the single highest-scoring region.

For alignments longer than a single window, the pipeline uses a **sliding window approach** to ensure primer candidates are discovered across the full length of the alignment. The MSA is divided into overlapping windows (default: 1000 bp windows with 200 bp overlap), and DECIPHER's DesignSignatures is run independently on each window. Primers designed in each window are collected and deduplicated — if the same forward-reverse pair is found in adjacent overlapping windows, only one copy is retained. This strategy prevents DECIPHER from concentrating all candidates in a single high-scoring region and ensures that viable primer sites in lower-scoring but still conserved regions are represented. Window size and overlap are configurable via `decipher.window_size` and `decipher.window_overlap` in the config file. If the alignment is short enough to fit in a single window, the windowing step is skipped and DECIPHER runs on the full alignment directly.

If DECIPHER is not available in the environment, the pipeline falls back to Primer3 via its Python bindings. In this mode, a majority-rule consensus sequence is computed from the alignment, and Primer3 designs primers against that consensus using nearest-neighbour thermodynamic models.

If pre-made primers are supplied via the `--primers` flag, the design step is skipped entirely and the supplied primers are used as the starting point.

The output is a CSV file of candidate primer pairs (name, forward sequence, reverse sequence).

### Step 2 — Individual Primer Validation

Each unique primer sequence (forward and reverse extracted separately from the design output) is validated individually using the Primer3 thermodynamic engine through its Python bindings.

For each primer, the following properties are calculated:

- **Melting temperature (Tm)** — Primers with Tm outside the configured range are rejected.
- **Hairpin stability** — Primers with hairpin delta-G more negative than the configured threshold are rejected.
- **Homodimer stability** — Primers with homodimer delta-G more negative than the configured threshold are flagged.
- **GC content** — Primers outside the configured range are flagged.
- **3-prime GC clamp** — The number of G or C bases in the last five positions at the 3-prime end. Outside the configured range is flagged.
- **Homopolymer runs** — Primers with runs exceeding the configured maximum are flagged.

Primers are classified as PASS (all checks pass), FLAG (soft threshold violations annotated but retained), or REJECT (hard threshold failures removed). Hard-reject criteria are Tm range and hairpin stability. All other checks are soft flags.

Degenerate primers containing IUPAC ambiguity codes are resolved to their most common concrete base before thermodynamic calculation, as the Primer3 engine requires unambiguous sequences.

This individual validation approach ensures that a good primer is never lost because it was originally paired with a bad partner. All pair-level checks (delta Tm, heterodimer) are deferred to the expansion step.

### Step 2b — Combination Expansion

All unique forward and reverse primers that passed individual validation are combined into every possible forward-reverse pairing. Each combination is filtered by:

- **Delta Tm** — The absolute Tm difference between forward and reverse must be within the configured maximum.
- **Heterodimer stability** — The free energy of the most stable heterodimer structure must be above the configured threshold.

Actual amplicon viability (whether the pair produces a product of the right size) is determined by obipcr in the analysis step, not by spatial filtering at this stage.

### Steps 3 and 4 — Binding Analysis and Off-Target Screening

Both binding analysis and off-target screening use obipcr from OBITools4, which performs in silico PCR — simulating the physical process of PCR amplification on a sequence database.

**In silico PCR on the MSA (Step 3):** For each primer pair, obipcr searches for both primer binding sites allowing a configurable number of mismatches per primer. When both primers bind with an intervening distance within the configured amplicon size range, obipcr reports a successful amplification.

**Pairwise identity matrix:** A full pairwise identity matrix is computed between all sequences in the alignment, counting matching non-gap positions as a fraction of total comparable sites.

**Sliding window identity:** Per-sequence identity relative to the reference is computed using a sliding window across the full alignment. A window of configurable size (default 400 bp) advances in configurable steps (default 100 bp) along the gapped alignment. At each position, pairwise identity between the reference and each other sequence is calculated over the window, requiring a minimum number of non-gap sites for a valid measurement. This produces a per-position identity profile for each sequence, which is displayed as the variability chart in the report.

**GenBank annotation mapping:** When a GenBank file is provided, the pipeline maps genomic features onto the alignment coordinates using BLAST against a doubled copy of the GenBank sequence. Doubling the reference sequence linearises the circular genome, so that regions spanning the origin of replication appear as a single contiguous stretch.

**Off-target screening (Step 4):** For every primer pair that successfully amplifies at least one sequence in the MSA, obipcr is run against the off-target sequence database with a separate mismatch tolerance.

### Step 5 — Ranking and Selection

Each primer pair is scored on five components, each normalised to a 0-1 scale:

- **Binding universality** (configurable weight) — Fraction of MSA sequences amplified.
- **Thermodynamic quality** (configurable weight) — Tm optimality, delta Tm, and secondary structure stability.
- **Species resolution** (configurable weight) — Proportion of variable positions in the amplicon.
- **Off-target specificity** (configurable weight) — Quality of ecoPCR hits.
- **Amplicon size suitability** (configurable weight) — Closeness to optimal amplicon length.

### Step 6 — Interactive Report

The analysis results are rendered as a self-contained HTML dashboard. The report includes:

- **Per-position variability chart** — Variability across the full alignment, rendered on the same canvas as the gene annotation track for pixel-perfect alignment.
- **Gene track** — Coloured gene segments with labels for major regions (coding, rRNA, D-loop). tRNAs are visible on hover.
- **Primer tracks** — Forward and reverse primers displayed at their binding positions using percentage-based coordinates that scale with window resize. No primers are pre-selected on load. Clicking a primer greys out incompatible partners — those where the resulting amplicon would fall outside the configured size range (`design.min_amplicon_length` / `design.max_amplicon_length`) or where the delta Tm exceeds 5 °C.
- **Selected primer detail** — Sequences, Tm, hairpin, homodimer, homopolymer, GC clamp, delta Tm, amplicon length, and amplification count.
- **Amplification prediction table** — Per-species mismatch data for the selected primer combination. Shows confirmed results (YES/NO) for tested pairs and mismatch estimates for untested combinations.
- **Amplicon alignment viewer** — Gapped MSA region between selected primers with FWD/REV regions coloured. Scrolls horizontally and vertically independently. Per-position variability bar chart aligned to sequence columns below.
- **Amplicon length distribution** — Histogram of ungapped amplicon lengths across species.
- **ecoPCR specificity table** — Off-target hits with mismatch counts and amplicon lengths.
