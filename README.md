# eDNA Primer Design Pipeline

A pipeline for designing and validating metabarcoding primers from a multiple sequence alignment. It identifies candidate primer binding sites, validates their thermodynamic properties, tests binding across all sequences in the alignment, screens for off-target amplification against a broader sequence database, and ranks the results — producing an interactive report for primer selection.

---

## Usage

Design new primers from an alignment:

```bash
python run_pipeline.py -i alignment.fasta -g reference.gb
```

With off-target screening against a sequence database:

```bash
python run_pipeline.py -i alignment.fasta -g reference.gb -d path/to/fasta_files
```

Evaluate existing primers instead of designing new ones:

```bash
python run_pipeline.py -i alignment.fasta -g reference.gb -p my_primers.csv
```

When `-p` is supplied, the design step is skipped and the supplied primers are used as the starting point. They pass through the same validation, analysis, and ranking as designed candidates. The primers CSV must have columns `name`, `forward`, `reverse`.

Custom output directory:

```bash
python run_pipeline.py -i alignment.fasta -g reference.gb -o results/my_run
```

Run with the included example data:

```bash
python run_pipeline.py -i data/example/input/alignment.fasta -g data/example/input/NC_008135.gb
```

### Arguments

| Flag             | Description                                                         |
| ---------------- | ------------------------------------------------------------------- |
| `-i / --input`   | Input MSA FASTA file                                                |
| `-o / --output`  | Output directory (default: `data/output`)                           |
| `-d / --database`| FASTA directory or file for off-target screening (optional)         |
| `-g / --genbank` | GenBank file for gene annotation mapping (auto-detected if omitted) |
| `-p / --primers` | CSV of primers to evaluate — skips design step (optional)           |
| `-c / --config`  | Config YAML file (default: `data/config.yaml`)                      |

The database flag accepts a directory containing multiple FASTA files, or a single FASTA file. If `--genbank` is not specified, the pipeline looks for a `.gb` file in the same directory as the input alignment and uses it automatically.

### Debug Mode

Individual pipeline steps can be invoked separately using debug mode:

```bash
python run_pipeline.py debug analyse -i alignment.fasta -d path/to/fasta_files
python run_pipeline.py debug rank
python run_pipeline.py debug report
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

- **primers/** — CSV files at each stage: initial candidates, thermodynamically filtered, expanded combinations.
- **analysis/** — JSON file containing the complete analysis results (binding data, identity matrices, ecoPCR hits, annotations). Also contains the cached annotation mapping and the final ranked primers CSV.
- **report.html** — An interactive HTML dashboard for visually comparing and selecting primer pairs.

---

## Configuration

All pipeline parameters are controlled through `data/config.yaml`. A custom config can be specified with `-c/--config`. Key sections:

**design** — Primer length range, amplicon size range, and optimal amplicon size. Adjusting the amplicon size range changes which primer combinations are considered viable. For eDNA work with degraded samples, short amplicons (90–180 bp) are standard. For long-read sequencing (e.g. Nanopore), this can be extended.

**thermodynamics** — Tm range, maximum delta Tm, GC content range, free energy thresholds for hairpin, homodimer, and heterodimer formation, maximum homopolymer run length, and 3' GC clamp requirements.

**binding** — Maximum number of mismatches allowed when scanning for primer binding sites, and maximum distance between forward and reverse primer hits to be considered an amplicon.

**ecopcr** — Parameters passed to obipcr for off-target screening: mismatch tolerance, amplicon size range, and whether to treat sequences as circular.

**ranking** — Weights for each of the five scoring components. These can be adjusted to prioritise different aspects of primer quality depending on the application.

**decipher** — Parameters for the DECIPHER DesignSignatures algorithm: minimum coverage, resolution, and maximum degeneracy.

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

### Step 1 — Primer Design

The pipeline identifies candidate primer binding sites within conserved regions of the input alignment using DECIPHER's DesignSignatures algorithm (Wright, 2016). DECIPHER is an R/Bioconductor package specifically designed for designing PCR primers from multiple sequence alignments. It uses a signature-based approach that maximises the diversity of amplicon sequences (to allow species discrimination) while ensuring the primer binding sites are conserved across the aligned sequences (to allow universal amplification).

The alignment is imported into a DECIPHER sequence database, and DesignSignatures is called with constraints on primer length (18-25 bp by default), amplicon size (90-180 bp), GC content (40-60%), and maximum degeneracy (128-fold). DECIPHER evaluates millions of potential primer combinations and returns the top-scoring pairs ranked by their predicted ability to produce diverse, species-discriminating amplicons.

If DECIPHER is not available in the environment, the pipeline falls back to Primer3 via its Python bindings. In this mode, a majority-rule consensus sequence is computed from the alignment, and Primer3 designs primers against that consensus using nearest-neighbour thermodynamic models.

If pre-made primers are supplied via the `--primers` flag, the design step is skipped entirely and the supplied primers are used as the starting point. This allows the pipeline to be used purely for evaluation of existing primers — they pass through the same thermodynamic validation, binding analysis, off-target screening, and ranking as computationally designed candidates.

The output is a CSV file of candidate primer pairs (name, forward sequence, reverse sequence).

### Step 2 — Thermodynamic Validation

Each candidate primer pair is evaluated for thermodynamic viability using the Primer3 thermodynamic engine (Untergasser et al., 2012) through its Python bindings (primer3-py). Primer3 uses the SantaLucia unified nearest-neighbour model for DNA thermodynamics, which accounts for base stacking interactions, salt concentration effects, and strand concentration.

For each primer, the following properties are calculated:

**Melting temperature (Tm)** — The temperature at which 50% of primer-template duplexes are dissociated. Calculated using nearest-neighbour thermodynamics, which considers not just base composition but the identity of each adjacent base pair. Primers with Tm outside the configured range (default 58-64 degrees C) are rejected.

**Delta Tm** — The absolute difference in melting temperature between the forward and reverse primers. A large delta Tm means one primer will bind much more strongly than the other at a given annealing temperature, leading to asymmetric amplification. Pairs exceeding the configured threshold (default 5 degrees C) are rejected.

**Hairpin stability** — The free energy of the most stable intramolecular hairpin structure the primer can form by folding back on itself. A primer with a stable hairpin will self-anneal instead of binding the template, reducing PCR efficiency. Primers with hairpin delta-G more negative than the threshold (default -3.0 kcal/mol) are rejected.

**Self-dimer stability** — The free energy of the most stable homodimer (two copies of the same primer annealing to each other). Self-dimers consume primer molecules and produce artefact bands. Primers with homodimer delta-G more negative than the threshold (default -9.0 kcal/mol) are flagged.

**Hetero-dimer stability** — The free energy of the most stable structure formed between the forward and reverse primers annealing to each other. Hetero-dimers are particularly problematic because both primers are present in every PCR reaction. Pairs with hetero-dimer delta-G more negative than the threshold (default -9.0 kcal/mol) are flagged.

**GC content** — The percentage of guanine and cytosine bases. GC content affects Tm and binding strength. Primers outside the configured range (default 40-60%) are flagged.

**3-prime GC clamp** — The number of G or C bases in the last five positions at the 3-prime end. A moderate GC clamp (1-3 G/C in the last 5 bases) stabilises the 3-prime end where polymerase extension begins, improving priming efficiency. Too few or too many are flagged.

**Homopolymer runs** — The maximum number of consecutive identical bases. Long runs (e.g. AAAA) can cause polymerase slippage and mispriming. Primers with runs exceeding the configured maximum (default 4) are flagged.

Primers are classified as PASS (all checks pass), FLAG (soft threshold violations annotated but retained), or REJECT (hard threshold failures removed). Degenerate primers containing IUPAC ambiguity codes are resolved to their most common concrete base before thermodynamic calculation, as the Primer3 engine requires unambiguous sequences.

### Step 2b — Combination Expansion

DECIPHER designs specific forward-reverse pairs, but a forward primer designed with one reverse may also work well with a different reverse primer at a nearby position. This step extracts all unique forward and all unique reverse sequences from the filtered candidates, then generates every possible forward-reverse combination that passes a delta Tm compatibility check (within the configured maximum, default 5 degrees C).

This expansion ensures the interactive report can display all viable pairings, not just the ones the design algorithm happened to output. For example, 4 unique forwards and 75 unique reverses might produce 300 possible combinations, of which a subset will be thermodynamically compatible.

### Steps 3 and 4 — Binding Analysis and Off-Target Screening

Both binding analysis and off-target screening use obipcr from OBITools4 (Boyer et al., 2016), which performs in silico PCR — simulating the physical process of PCR amplification on a sequence database.

**In silico PCR on the MSA (Step 3):**

The ungapped sequences from the MSA are written to a temporary FASTA file. For each primer pair, obipcr searches for the forward primer binding site on one strand and the reverse primer binding site on the other strand, allowing up to 3 mismatches per primer (configurable). When both primers bind with an intervening distance within the configured amplicon size range, obipcr reports a successful amplification — returning the amplicon sequence, the number of mismatches at each primer binding site, and the binding positions.

This replaces traditional approaches that use simple string matching or custom mismatch counting. obipcr correctly handles the directionality of primer binding (forward primer on the sense strand, reverse primer on the antisense strand) and applies the mismatch tolerance uniformly using the Agrep pattern-matching algorithm.

The results tell us: which sequences in the alignment each primer pair would amplify, with how many mismatches, and what the resulting amplicon sequences look like. Amplicon sequences are then compared position-by-position to calculate per-site variability — positions where species differ are the ones that provide taxonomic resolution.

**Primer properties** are calculated for each pair using primer3-py nearest-neighbour thermodynamics, replacing the less accurate Wallace rule (Tm = 2(A+T) + 4(G+C)) that is sometimes used as an approximation.

**GenBank annotation mapping:**

When a GenBank file is provided (or auto-detected), the pipeline maps genomic features (genes, tRNAs, rRNAs, D-loop) onto the alignment coordinates. This is done by BLASTing the ungapped reference MSA sequence against a doubled copy of the GenBank sequence. Doubling the reference sequence (concatenating it with itself) linearises the circular genome, so that regions spanning the origin of replication appear as a single contiguous stretch rather than being split across the two ends. BLAST returns the position in the doubled reference where the MSA sequence aligns, which is then converted back to real genome coordinates using modular arithmetic.

After finding the offset, the pipeline verifies the mapping by extracting the corresponding region from the genome and comparing it base-by-base to the ungapped MSA sequence. An exact match confirms the mapping is correct. Any mismatches are reported with a percentage — minor variants (under 1%) are noted but accepted, while larger discrepancies trigger a warning that annotation positions may be inaccurate.

Each GenBank feature is then converted from genome coordinates to MSA-region coordinates (using the offset and modular arithmetic for circular genome handling), and then from ungapped MSA coordinates to gapped alignment coordinates (using a lookup table built from the reference sequence's gap pattern in the alignment). Features that fall outside the MSA region are discarded. Features that partially overlap are clipped to the MSA boundaries.

The annotation mapping is cached to an annotations.json file with an MD5 checksum of the GenBank file. On subsequent runs, if the GenBank file has not changed, the cached annotations are reused without repeating the BLAST step.

**Off-target screening (Step 4):**

For every primer pair that successfully amplifies at least one sequence in the MSA, obipcr is run again against the off-target sequence database. This is identical to the MSA analysis — same tool, same parameters — but applied to a broader set of sequences to assess specificity.

The results show how many non-target sequences each primer pair would amplify, the distribution of mismatches (pairs with many zero-mismatch off-target hits are less specific), and the range of amplicon lengths produced. Primer pairs that do not produce any amplicon in the MSA are skipped for off-target screening, since there is no point checking the specificity of primers that do not work on the target.

### Step 5 — Ranking and Selection

Each primer pair is scored on five components, each normalised to a 0-1 scale:

**Binding universality (weight: 0.25)** — The fraction of MSA sequences that the primer pair successfully amplifies. A score of 1.0 means the primers work on every sequence in the alignment. This measures how universal the primers are across the target taxa.

**Thermodynamic quality (weight: 0.20)** — A composite of how close the Tm values are to the optimal midpoint of the configured range, how small the delta Tm is, and how stable the hairpin and dimer structures are (less stable is better, meaning the primer prefers binding to the template over forming secondary structures).

**Species resolution (weight: 0.20)** — The fraction of positions within the amplicon that are variable across species. A higher proportion of variable positions means the amplicon sequences differ more between species, which is essential for metabarcoding — the amplicon must contain enough variation to distinguish species after sequencing.

**Off-target specificity (weight: 0.20)** — Derived from the ecoPCR results. For metabarcoding primers, broader taxonomic coverage within the target group is desirable (more hits is better), but the quality of those hits matters — pairs where most hits are perfect matches score higher than pairs where hits require multiple mismatches.

**Amplicon size suitability (weight: 0.15)** — How close the mean amplicon length is to the configured optimal size (default 130 bp). For eDNA applications, shorter amplicons are preferred because environmental DNA is typically fragmented. This can be adjusted in the configuration for other sequencing platforms.

The five component scores are combined using a weighted sum to produce a composite score. Primer pairs are ranked by composite score in descending order.

### Step 6 — Interactive Report

The analysis results are rendered as a self-contained HTML dashboard that runs entirely in the browser with no server required. The report includes:

**Gene track** — A visual representation of the aligned region with gene annotations (from the GenBank file) shown as coloured blocks. The track provides spatial context for where the primers bind relative to known genomic features.

**Primer location tracks** — Separate tracks for forward and reverse primers, displayed below the gene track at their binding positions. Forward primers are shown in amber, reverse primers in pink. Clicking a forward primer greys out incompatible reverse primers (those with too large a Tm difference or that would produce an amplicon outside the configured size range), and vice versa. Clicking an already-selected primer deselects it. Incompatible primers retain a coloured border so they remain identifiable.

**Selected primer detail** — Shows the sequences, lengths, Tm values, delta Tm, amplicon length, and amplification count for the currently selected forward-reverse combination.

**Amplification prediction table** — For the selected primer pair, shows each sequence in the alignment with the number of forward and reverse primer mismatches and whether amplification is predicted.

**ecoPCR specificity panel** — Mismatch distribution chart showing how many off-target sequences amplify at each mismatch level, amplicon length distribution histogram, and a scrollable table of all ecoPCR hits.

**Amplicon alignment viewer** — Base-by-base display of the amplicon sequences across all species, with primer regions highlighted, mismatches marked in red, and the inner variable region shown in cyan.

**Variability charts** — Per-position variability across the amplicon (showing which positions differ between species) and inner-region identity versus the reference.

**Alignment overview** — Sliding-window identity plot across the full alignment, pairwise identity matrix heatmap, and per-sequence identity bars.
