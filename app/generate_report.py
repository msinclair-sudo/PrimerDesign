#!/usr/bin/env python3
"""
generate_viz.py
===============
Reads blast_vis_data.json (output of blast_msa_analysis.py) and produces a
self-contained HTML dashboard with interactive primer-set switching.

Usage:
    python scripts/generate_viz.py blast_vis_data.json [-o output.html]
"""

import argparse
import json
import sys
from pathlib import Path

# ── Colour palette for sequences ─────────────────────────────────────────────
PALETTE = [
    "#0da888", "#d97706", "#0284c7", "#7c3aed", "#ea580c",
    "#dc2626", "#16a34a", "#0a8870", "#9333ea", "#c2410c",
    "#2563eb", "#b91c1c", "#059669", "#6d28d9", "#ca8a04",
    "#0369a1", "#a21caf", "#15803d", "#be123c", "#4f46e5",
]


import re as _re

# Regex for NCBI accessions: e.g. NC_014696.1, GU570660.1, OR840808.1
_ACCESSION_RE = _re.compile(r'^[A-Z]{1,3}_?\d{5,}(?:\.\d+)?$')

_METADATA_WORDS = {"mitochondrion,", "mitochondrion", "isolate", "strain",
                   "voucher", "breed", "complete", "partial", "genome",
                   "genome,", "sequence"}


def _is_accession(token: str) -> bool:
    """Check if a token looks like an NCBI accession number."""
    return bool(_ACCESSION_RE.match(token))


def short_name(full: str) -> str:
    """Shorten a FASTA header to a readable species label.

    Handles formats:
      NC_014696.1 Leggadina lakedownensis mitochondrion, complete genome
      GU570660.1 Rattus leucopus isolate RleuPN66 mitochondrion, complete genome
      NC_008135_Petaurus_breviceps_REF
      Conilurus_penicillatus
      OR840808.1
    """
    # Try space-separated NCBI-style: "ACCESSION Species name ... mitochondrion"
    if " " in full:
        parts = full.split()
        acc = parts[0]
        words = parts[1:]
        species_words = []
        for w in words:
            if w.lower() in _METADATA_WORDS:
                break
            species_words.append(w)
        if species_words:
            return f"{' '.join(species_words)} ({acc})"
        return acc

    # Underscore-separated: try to identify accession prefix
    parts = full.split("_")
    is_ref = parts[-1] == "REF"
    if is_ref:
        parts = parts[:-1]

    # Try to reconstruct an accession from leading parts
    # e.g. ["NC", "014696.1", "Species", "name"] -> acc = "NC_014696.1"
    acc = None
    sp_parts = parts
    for i in range(1, min(len(parts), 3)):
        candidate = "_".join(parts[:i + 1])
        if _is_accession(candidate):
            acc = candidate
            sp_parts = parts[i + 1:]
            break
    # Single-part accession: e.g. ["OR840808.1", ...]
    if acc is None and _is_accession(parts[0]):
        acc = parts[0]
        sp_parts = parts[1:]

    # No accession found — treat entire string as species name
    if acc is None:
        sp = " ".join(parts)
        return sp

    tag = f"{acc} · REF" if is_ref else acc
    if sp_parts:
        return f"{' '.join(sp_parts)} ({tag})"
    return tag


def build_js_data(data: dict, short_names: dict, colors: dict,
                   max_amp_len: int = 300, min_amp_len: int = 90) -> str:
    """Build the JavaScript data constants block."""
    seq_names = data["seq_names"]
    meta = data["meta"]
    genes = data["genes"]

    # Build primer_sets JS array
    ps_js = []
    for ps in data["primer_sets"]:
        gc = ps["properties"].get("gapped_coords", {})
        amplicons = {}
        for name in seq_names:
            amps = ps["hits"][name]["amplicons"]
            if amps:
                amplicons[name] = amps[0]["sequence"]

        # Per-sequence hit summary for the table
        hit_summary = {}
        for name in seq_names:
            h = ps["hits"][name]
            if h["amplicons"]:
                a = h["amplicons"][0]
                hit_summary[name] = {
                    "has_amp": True,
                    "fwd_mm": a["fwd_mm"],
                    "rev_mm": a["rev_mm"],
                    "amp_len": a["length"],
                }
            else:
                fwd_best = min(h["fwd_hits"] or [{"mm": 99}], key=lambda x: x["mm"])
                rev_best = min(h["rev_hits"] or [{"mm": 99}], key=lambda x: x["mm"])
                hit_summary[name] = {
                    "has_amp": False,
                    "fwd_mm": fwd_best["mm"],
                    "rev_mm": rev_best["mm"],
                }

        props = ps["properties"]
        ss = props.get("secondary_structure", {})

        # ecoPCR data (if available)
        ecopcr = ps.get("ecopcr")
        ecopcr_js = None
        if ecopcr:
            ecopcr_js = {
                "db_sequences": ecopcr["db_sequences"],
                "total_hits": ecopcr["total_hits"],
                "mismatch_counts": ecopcr["mismatch_counts"],
                "amp_len_min": ecopcr["amplicon_length_min"],
                "amp_len_max": ecopcr["amplicon_length_max"],
                "amp_len_median": ecopcr["amplicon_length_median"],
                "amp_lengths": [h["amp_len"] for h in ecopcr["hits"]],
                "hits": [
                    {"id": h["id"], "def": h["definition"][:60],
                     "fwd_mm": h["fwd_mm"], "rev_mm": h["rev_mm"],
                     "amp_len": h["amp_len"]}
                    for h in ecopcr["hits"]
                ],
            }

        ps_obj = {
            "name": ps["name"],
            "fwd": ps["fwd_primer"],
            "rev": ps["rev_primer"],
            "fwd_len": len(ps["fwd_primer"]),
            "rev_len": len(ps["rev_primer"]),
            "fwd_gc": props["fwd"]["gc_pct"],
            "rev_gc": props["rev"]["gc_pct"],
            "fwd_tm": props["fwd"].get("tm", props["fwd"].get("tm_wallace", 0)),
            "rev_tm": props["rev"].get("tm", props["rev"].get("tm_wallace", 0)),
            "fwd_hairpin": props["fwd"].get("hairpin_dg", 0),
            "rev_hairpin": props["rev"].get("hairpin_dg", 0),
            "fwd_homodimer": props["fwd"].get("homodimer_dg", 0),
            "rev_homodimer": props["rev"].get("homodimer_dg", 0),
            "heterodimer": props.get("heterodimer_dg", 0),
            "delta_tm": props["delta_tm"],
            "ss": ss,
            "gapped": gc,
            "amplicons": amplicons,
            "hit_summary": hit_summary,
            "ecopcr": ecopcr_js,
            "properties": {
                "flags": props.get("flags", []),
                "status": props.get("status", "PASS"),
            },
        }
        ps_js.append(ps_obj)

    lines = [
        f"const SEQ_NAMES = {json.dumps(seq_names)};",
        f"const SHORT = {json.dumps(short_names)};",
        f"const COLORS = {json.dumps(colors)};",
        f"const GENES = {json.dumps(genes)};",
        f"const ALIGN_LEN = {meta['alignment_length']};",
        f"const GAPPED = {json.dumps(data.get('gapped_sequences', {}))};",
        f"const MATRIX = {json.dumps(data['matrix'])};",
        f"const ID_VS_REF = {json.dumps(data['id_vs_ref'])};",
        f"const X_POS = {json.dumps(data['sliding_window']['x_positions'])};",
        f"const SLIDING = {json.dumps(data['sliding_window']['data'])};",
        f"const PRIMER_SETS = {json.dumps(ps_js)};",
        f"const MAX_AMP_LEN = {max_amp_len};",
        f"const MIN_AMP_LEN = {min_amp_len};",
        "let activePSIndex = -1;",
        "function PS() { return activePSIndex >= 0 ? PRIMER_SETS[activePSIndex] : null; }",
    ]
    return "\n".join(lines)


def generate_html(data: dict, header_cfg: dict | None = None,
                   full_cfg: dict | None = None) -> str:
    """Generate the full HTML dashboard."""
    seq_names = data["seq_names"]
    meta = data["meta"]
    ref_name = meta["reference"]

    # Extract amplicon length limits from config
    design_cfg = (full_cfg or {}).get("design", {})
    max_amp_len = design_cfg.get("max_amplicon_length", 300)
    min_amp_len = design_cfg.get("min_amplicon_length", 90)

    # Embed Chart.js from local file (fully self-contained HTML)
    chartjs_path = Path(__file__).parent / "chartjs_4.4.1.min.js"
    chartjs_src = chartjs_path.read_text(encoding="utf-8")

    sn = {name: short_name(name) for name in seq_names}
    colors = {name: PALETTE[i % len(PALETTE)] for i, name in enumerate(seq_names)}

    js_data = build_js_data(data, sn, colors, max_amp_len, min_amp_len)

    # Reference species for title fallback
    parts = ref_name.split("_")
    if len(parts) >= 3 and parts[1].isdigit():
        ref_species = " ".join(parts[2:]).replace("REF", "").strip()
        ref_acc = f"{parts[0]}_{parts[1]}"
    else:
        ref_species = " ".join(parts[1:]).replace("REF", "").strip()
        ref_acc = parts[0]

    # Header from config or fallback to auto-generated
    if header_cfg:
        h_main = header_cfg.get("main", "Primers for")
        h_coloured = header_cfg.get("coloured", ref_species)
        h_last = header_cfg.get("last", "")
        h_subtitle = header_cfg.get("subtitle", "Primer binding assessment")
    else:
        h_main = "Primers for"
        h_coloured = ref_species
        h_last = ""
        h_subtitle = "Primer binding assessment"

    align_len = meta["alignment_length"]
    n_primer_sets = meta.get("n_primer_sets", len(data.get("primer_sets", [])))

    genes = data["genes"]

    axis_ticks = [0]
    pos = 500
    while pos < align_len:
        axis_ticks.append(pos)
        pos += 500
    axis_ticks.append(align_len)
    axis_spans = "".join(f"<span>{t}</span>" for t in axis_ticks)

    html = f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{h_main} {h_coloured}{h_last}</title>
<script>{chartjs_src}</script>
<link href="https://fonts.googleapis.com/css2?family=Playfair+Display:ital,wght@0,600;0,700;1,400&family=IBM+Plex+Sans:wght@300;400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap" rel="stylesheet">
<style>
  :root {{
    /* ── Surfaces ── */
    --bg:      #e8ede6;
    --bg2:     #f2f5ef;
    --bg3:     #dae1d4;
    --bg4:     #cdd6c6;
    --border:  #b4c0ac;
    /* ── Text ── */
    --text:    #0e1a0d;
    --dim:     #374832;
    --muted:   #5f7358;
    /* ── Semantic colours ── */
    --accent:  #1d6457;
    --teal:    #1d6457;
    --gold:    #c98a18;
    --cyan:    #1a7a9e;
    --green:   #16a34a;
    --coral:   #b83228;
    --purple:  #6e3ab8;
    --orange:  #b85c14;
    /* ── Primer pair colours ── */
    --primer-fwd: #c98a18;
    --primer-rev: #1d6457;
    /* ── Transparent variants ── */
    --green-bg:   rgba(22,163,74,0.14);
    --gold-bg:    rgba(201,138,24,0.14);
    --coral-bg:   rgba(184,50,40,0.12);
    --cyan-bg:    rgba(26,122,158,0.10);
    --purple-bg:  rgba(110,58,184,0.10);
    --fwd-bg:     rgba(201,138,24,0.08);
    --rev-bg:     rgba(29,100,87,0.08);
    --amplicon:   rgba(29,100,87,0.10);
    --flag-bg:    rgba(201,138,24,0.14);
    --diff-bg:    rgba(184,50,40,0.12);
    --sel-sp:     rgba(22,163,74,0.10);
    --sel-sp-hov: rgba(22,163,74,0.18);
    /* ── Typography ── */
    --font-mono: 'IBM Plex Mono', monospace;
    --font-head: 'Playfair Display', serif;
    --font-body: 'IBM Plex Sans', sans-serif;
  }}
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{
    background: var(--bg);
    color: var(--text);
    font-family: var(--font-body);
    font-weight: 300;
    min-height: 100vh;
    padding: 28px 28px 40px;
  }}
  .header {{
    display: flex; align-items: flex-start; justify-content: space-between;
    margin-bottom: 28px; padding-bottom: 20px;
    border-bottom: 1px solid var(--border); gap: 16px; flex-wrap: wrap;
  }}
  .header h1 {{
    font-family: var(--font-head); font-size: 20px; font-weight: 800;
    color: var(--text); letter-spacing: -0.3px;
  }}
  .header .sub {{ font-size: 11px; color: var(--muted); margin-top: 5px; }}
  .chip-row {{ display: flex; gap: 8px; flex-wrap: wrap; align-items: flex-start; }}
  .chip {{
    background: var(--bg3); border: 1px solid var(--border); border-radius: 5px;
    padding: 5px 10px; font-size: 10px; color: var(--muted); white-space: nowrap;
  }}
  .chip b {{ color: var(--accent); }}
  .section-label {{
    font-family: var(--font-head); font-size: 11px; font-weight: 700;
    letter-spacing: 2px; text-transform: uppercase; color: var(--muted);
    margin: 28px 0 12px; display: flex; align-items: center; gap: 10px;
  }}
  .section-label::after {{ content: ''; flex: 1; height: 1px; background: var(--border); }}
  .grid {{ display: grid; gap: 14px; }}
  .grid-2 {{ grid-template-columns: 1fr 1fr; }}
  .card {{
    background: var(--bg2); border: 1px solid var(--border);
    border-radius: 10px; padding: 18px 20px;
  }}
  .card.accent-border {{ border-color: var(--gold); }}
  .card-title {{
    font-size: 10px; font-weight: 600; text-transform: uppercase;
    letter-spacing: 1.2px; color: var(--muted); margin-bottom: 14px;
    display: flex; align-items: center; gap: 7px;
  }}
  .dot {{ width: 6px; height: 6px; border-radius: 50%; flex-shrink: 0; }}

  /* ── PRIMER SELECTOR ── */
  .ps-list {{ display: flex; flex-direction: column; gap: 10px; }}
  .ps-card {{
    background: var(--bg3); border: 2px solid var(--border); border-radius: 8px;
    padding: 12px 14px; cursor: pointer; transition: border-color 0.2s, background 0.2s;
  }}
  .ps-card:hover {{ border-color: var(--accent); background: var(--bg4); }}
  .ps-card.active {{ border-color: var(--primer-fwd); background: var(--fwd-bg); }}
  .ps-card .ps-name {{
    font-size: 11px; font-weight: 700; color: var(--text); margin-bottom: 8px;
    letter-spacing: 0.5px;
  }}
  .ps-card .ps-seq {{
    font-size: 9px; letter-spacing: 1px; margin-bottom: 3px; word-break: break-all;
  }}
  .ps-card .ps-stats {{
    display: flex; gap: 10px; flex-wrap: wrap; margin-top: 8px;
  }}

  /* ── GENE TRACK ── */
  .gene-track {{
    display: flex; height: 26px; border-radius: 4px; overflow: visible; position: relative;
  }}
  .gene-seg {{
    display: flex; align-items: center; justify-content: center;
    font-size: 8.5px; font-weight: 700; letter-spacing: 0.3px;
    color: rgba(0,0,0,0.8); overflow: hidden; position: relative;
    cursor: default; transition: filter 0.15s;
  }}
  .gene-seg:hover {{ filter: brightness(1.15); z-index: 5; }}
  .gene-axis {{
    display: flex; justify-content: space-between; margin-top: 5px;
    font-size: 8.5px; color: var(--muted);
  }}
  .gene-legend {{ display: flex; flex-wrap: wrap; gap: 10px; margin-top: 12px; }}
  .leg-item {{ display: flex; align-items: center; gap: 5px; font-size: 10px; color: var(--muted); }}
  .leg-swatch {{ width: 10px; height: 10px; border-radius: 2px; flex-shrink: 0; }}
  .track-container {{ position: relative; }}

  /* ── PRIMER SPECS (inline in card) ── */
  .primer-specs {{ display: flex; flex-direction: column; gap: 16px; }}
  .primer-block {{ border-left: 3px solid var(--primer-fwd); padding-left: 12px; }}
  .primer-block.rev {{ border-color: var(--primer-rev); }}
  .primer-tag {{ font-size: 9px; font-weight: 700; text-transform: uppercase; letter-spacing: 1px; margin-bottom: 6px; }}
  .primer-seq {{ font-size: 13px; font-weight: 600; letter-spacing: 1.5px; margin-bottom: 8px; word-break: break-all; }}
  .primer-stats {{ display: flex; gap: 14px; flex-wrap: wrap; }}
  .pstat {{ font-size: 10px; color: var(--muted); }}
  .pstat b {{ color: var(--text); }}

  /* ── ALIGNMENT ── */
  .aln-viewer {{
    overflow-x: auto; overflow-y: hidden;
    background: var(--bg); border-radius: 6px;
    border: 1px solid var(--border); padding: 12px;
  }}
  .aln-row {{ display: flex; gap: 0; line-height: 1.6; align-items: center; }}
  .aln-label {{
    width: 200px; flex-shrink: 0; font-size: 10px; color: var(--muted);
    padding-right: 10px; text-align: right; white-space: nowrap;
    overflow: hidden; text-overflow: ellipsis;
  }}
  .aln-label.ref {{ color: var(--gold); font-weight: 600; }}
  .aln-seq {{ font-size: 10.5px; letter-spacing: 1.2px; white-space: nowrap; line-height: 1.7; }}
  .aln-primer-fwd {{ color: var(--primer-fwd); font-weight: 700; }}
  .aln-primer-rev {{ color: var(--primer-rev); font-weight: 700; }}
  .aln-match  {{ color: var(--muted); }}
  .aln-diff   {{ color: var(--coral); font-weight: 700; background: var(--diff-bg); border-radius: 2px; }}
  .aln-inner  {{ color: var(--cyan); }}
  .aln-ruler-row {{ margin-bottom: 4px; }}

  /* ── HIT TABLE ── */
  .hit-table {{ width: 100%; border-collapse: collapse; font-size: 11px; }}
  .hit-table th {{
    text-align: left; padding: 5px 10px; font-size: 9px;
    text-transform: uppercase; letter-spacing: 1px; color: var(--muted);
    border-bottom: 1px solid var(--border); font-weight: 600;
  }}
  .hit-table td {{ padding: 7px 10px; border-bottom: 1px solid var(--bg3); vertical-align: middle; }}
  .hit-table tr:last-child td {{ border-bottom: none; }}
  .hit-table tr:hover td {{ background: var(--bg3); }}
  .hit-table tr.species-selected td {{ background: var(--sel-sp); }}
  .hit-table tr.species-selected:hover td {{ background: var(--sel-sp-hov); }}
  .hit-table tr {{ cursor: pointer; user-select: none; }}
  .mm-badge {{
    display: inline-flex; align-items: center; justify-content: center;
    width: 22px; height: 22px; border-radius: 50%; font-size: 11px; font-weight: 700;
  }}
  .mm-0 {{ background: var(--green-bg); color: var(--green); }}
  .mm-1 {{ background: var(--gold-bg); color: var(--gold); }}
  .mm-bad {{ background: var(--coral-bg); color: var(--coral); }}
  .sp-name {{ font-style: italic; }}
  .accession {{ font-size: 9px; color: var(--muted); display: block; margin-top: 1px; }}

  /* ── IDENTITY BARS ── */
  .id-bars {{ display: flex; flex-direction: column; gap: 5px; }}
  .id-bar-row {{ display: flex; align-items: center; gap: 8px; }}
  .id-bar-label {{ width: 185px; flex-shrink: 0; text-align: right; color: var(--muted); font-size: 10px; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }}
  .id-bar-track {{ flex: 1; height: 12px; background: var(--bg3); border-radius: 2px; overflow: hidden; position: relative; }}
  .id-bar-fill {{ height: 100%; border-radius: 2px; transition: width 1s ease; }}
  .id-bar-val {{ width: 46px; flex-shrink: 0; font-size: 10px; font-weight: 700; }}

  /* ── HEATMAP ── */
  .heatmap-table {{ border-collapse: collapse; font-size: 9px; width: 100%; }}
  .heatmap-table th {{ padding: 3px; color: var(--muted); font-weight: 400; text-align: center; }}
  .heatmap-table th.row-head {{ text-align: right; padding-right: 6px; max-width: 100px; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; font-size: 9px; }}
  .heatmap-table td {{ width: 38px; height: 28px; text-align: center; font-size: 8.5px; font-weight: 700; border: 1px solid var(--bg); border-radius: 2px; cursor: default; }}

  /* ── CHARTS ── */
  .chart-h240 {{ position: relative; height: 240px; }}
  .chart-h200 {{ position: relative; height: 200px; }}
  .sp-legend {{ display: flex; flex-wrap: wrap; gap: 8px; margin-top: 10px; }}
  .sp-leg-item {{ display: flex; align-items: center; gap: 5px; font-size: 9.5px; color: var(--muted); cursor: pointer; }}
  .sp-leg-item:hover {{ color: var(--text); }}
  .sp-leg-item.inactive {{ opacity: 0.3; }}
  .sp-dot {{ width: 7px; height: 7px; border-radius: 50%; flex-shrink: 0; }}

  .footnote {{
    margin-top: 28px; font-size: 9.5px; color: var(--muted); text-align: center;
    border-top: 1px solid var(--border); padding-top: 14px;
  }}

  /* ── PRIMER TRACK ── */
  .primer-track-row {{
    position: relative;
    height: 14px;
    margin-bottom: 2px;
  }}
  .primer-bar {{
    position: absolute;
    height: 12px;
    border-radius: 3px;
    cursor: pointer;
    transition: filter 0.15s, transform 0.1s;
    display: flex;
    align-items: center;
    justify-content: center;
    font-size: 7px;
    font-weight: 700;
    overflow: hidden;
    white-space: nowrap;
  }}
  .primer-bar:hover {{ filter: brightness(1.3); transform: scaleY(1.3); z-index: 10; }}
  .primer-bar.active {{ outline: 2px solid var(--text); outline-offset: 1px; z-index: 11; }}

  /* ── SELECTED PRIMER DETAIL ── */
  .sel-primer-seq {{
    font-size: 12px; font-weight: 600; letter-spacing: 1.2px;
    margin-bottom: 4px; word-break: break-all;
  }}
  .sel-primer-stats {{
    display: flex; gap: 12px; flex-wrap: wrap; margin-top: 8px;
  }}

  /* ── NO-AMPLICON message ── */
  .no-amp-msg {{
    padding: 20px; text-align: center; color: var(--muted); font-size: 12px;
  }}
  .no-amp-msg b {{ color: var(--coral); }}
</style>
</head>
<body>

<!-- HEADER -->
<div class="header">
  <div>
    <h1>{h_main} <em style="font-style:italic;color:var(--accent)">{h_coloured}</em>{h_last}</h1>
    <div class="sub">{h_subtitle}</div>
  </div>
  <div class="chip-row">
    <div class="chip">Sequences: <b>{len(seq_names)}</b></div>
    <div class="chip">Alignment: <b>{align_len:,} bp</b></div>
    <label class="chip" style="cursor:pointer;user-select:none;" title="Hide primers that amplify any off-target sequence">
      <input type="checkbox" id="hideOffTargetToggle" style="margin:0 4px 0 0;vertical-align:middle;accent-color:var(--teal);"> Hide off-target hits
    </label>
  </div>
</div>

<!-- MAIN DASHBOARD GRID — single grid so all 4fr/1fr columns align -->
<div class="grid" style="grid-template-columns:7fr 3fr;margin-top:14px;align-items:start;">

  <!-- Row 1 left: Variability + Gene Track + Primers -->
  <div class="card">
    <div class="card-title"><div class="dot" style="background:var(--coral)"></div>Per-Position Variability &amp; Gene Track ({align_len:,} bp)</div>
    <div style="height:180px;position:relative;"><canvas id="varFullChart"></canvas></div>
    <div class="gene-axis">{axis_spans}</div>
    <div style="margin-top:14px;">
      <div style="font-size:9px;color:var(--primer-fwd);margin-bottom:4px;text-transform:uppercase;letter-spacing:1px;font-weight:600;">▶ Forward primers</div>
      <div id="fwdTrack" style="position:relative;"></div>
    </div>
    <div style="margin-top:8px;">
      <div style="font-size:9px;color:var(--primer-rev);margin-bottom:4px;text-transform:uppercase;letter-spacing:1px;font-weight:600;">◀ Reverse primers</div>
      <div id="revTrack" style="position:relative;min-height:16px;"></div>
    </div>
  </div>

  <!-- Row 1 right: Selected Primer + Amplification Prediction -->
  <div style="display:flex;flex-direction:column;gap:14px;">
    <div class="card accent-border" id="selectedPrimerCard">
      <div class="card-title"><div class="dot" style="background:var(--primer-fwd)"></div><span id="selectedPrimerTitle">Selected Primer</span></div>
      <div id="selectedPrimerDetail"></div>
    </div>
    <div class="card">
      <div class="card-title"><div class="dot" style="background:var(--green)"></div><span id="ampPredTitle">Amplification Prediction · click species to filter primers</span></div>
      <div style="max-height:400px;overflow-y:auto;">
        <table class="hit-table">
          <thead>
            <tr>
              <th>Species / Accession</th>
              <th style="text-align:center">FWD mm</th>
              <th style="text-align:center">REV mm</th>
              <th style="text-align:center">Amplicon</th>
              <th style="text-align:center">Call</th>
            </tr>
          </thead>
          <tbody id="hitTableBody"></tbody>
        </table>
      </div>
    </div>
  </div>

  <!-- Row 2 left: Amplicon Alignment -->
  <div class="card" style="overflow:hidden;">
    <div class="card-title"><div class="dot" style="background:var(--cyan)"></div>Amplicon Alignment</div>
    <div style="display:flex;gap:12px;margin-bottom:8px;font-size:9px;flex-wrap:wrap;">
      <span><span style="color:var(--primer-fwd);font-weight:700">■</span> FWD primer</span>
      <span><span style="color:var(--primer-rev);font-weight:700">■</span> REV primer</span>
    </div>
    <div id="alnHScroll" style="overflow-x:auto;overflow-y:hidden;">
      <div style="display:flex;">
        <div id="varYAxis" style="height:60px;display:flex;flex-direction:column;justify-content:space-between;padding-right:4px;text-align:right;font-family:var(--font-mono);font-size:8px;color:var(--muted);width:140px;flex-shrink:0;position:sticky;left:0;z-index:2;background:var(--bg2);">
          <span>100</span><span>0</span>
        </div>
        <canvas id="varChart" style="display:block;height:60px;"></canvas>
      </div>
      <div id="alnRows"></div>
    </div>
  </div>

  <!-- Row 2 right: Amplicon Length Distribution -->
  <div class="card">
    <div class="card-title"><div class="dot" style="background:var(--purple)"></div><span id="ampLenTitle">Amplicon Length Distribution</span></div>
    <div class="chart-h200"><canvas id="ampLenChart"></canvas></div>
  </div>

</div>

<!-- ROW 4: ecoPCR Hits (full width) -->
<div style="margin-top:14px;">
  <div class="card" id="ecopcrSection">
    <div class="card-title"><div class="dot" style="background:var(--teal)"></div><span id="ecopcrTableTitle">ecoPCR Hits</span></div>
    <div style="max-height:300px;overflow-y:auto;">
      <table class="hit-table">
        <thead>
          <tr>
            <th>Accession</th>
            <th>Description</th>
            <th style="text-align:center">FWD mm</th>
            <th style="text-align:center">REV mm</th>
            <th style="text-align:center">Amplicon</th>
          </tr>
        </thead>
        <tbody id="ecopcrTableBody"></tbody>
      </table>
    </div>
  </div>
</div>

<div class="footnote">
  <em>{ref_species}</em> · {len(seq_names)} sequences
</div>

<!-- JAVASCRIPT -->
<script>
// ── DATA ──────────────────────────────────────────────────────────────────────
{js_data}

// ── THEME (single source of truth for JS colours) ───────────────────────────
const S = getComputedStyle(document.documentElement);
const THEME = {{
  bg: S.getPropertyValue('--bg').trim(),
  bg2: S.getPropertyValue('--bg2').trim(),
  bg3: S.getPropertyValue('--bg3').trim(),
  bg4: S.getPropertyValue('--bg4').trim(),
  border: S.getPropertyValue('--border').trim(),
  text: S.getPropertyValue('--text').trim(),
  muted: S.getPropertyValue('--muted').trim(),
  accent: S.getPropertyValue('--accent').trim(),
  gold: S.getPropertyValue('--gold').trim(),
  cyan: S.getPropertyValue('--cyan').trim(),
  green: S.getPropertyValue('--green').trim(),
  coral: S.getPropertyValue('--coral').trim(),
  purple: S.getPropertyValue('--purple').trim(),
  teal: S.getPropertyValue('--teal').trim(),
  fwd: S.getPropertyValue('--primer-fwd').trim(),
  rev: S.getPropertyValue('--primer-rev').trim(),
  // Chart defaults
  gridColor: S.getPropertyValue('--bg3').trim(),
  tickColor: S.getPropertyValue('--muted').trim(),
  tooltipBg: S.getPropertyValue('--bg2').trim(),
  tooltipBorder: S.getPropertyValue('--border').trim(),
  tooltipTitle: S.getPropertyValue('--text').trim(),
  tooltipBody: S.getPropertyValue('--dim').trim(),
  font: 'IBM Plex Mono',
}};
// Chart.js defaults
const CHART_TICK = {{ color: THEME.tickColor, font: {{ family: THEME.font, size: 9 }} }};
const CHART_TICK_SM = {{ color: THEME.tickColor, font: {{ family: THEME.font, size: 8 }} }};
const CHART_GRID = {{ color: THEME.gridColor }};
const CHART_TOOLTIP = {{ backgroundColor: THEME.tooltipBg, borderColor: THEME.tooltipBorder, borderWidth: 1,
  titleColor: THEME.tooltipTitle, bodyColor: THEME.tooltipBody,
  titleFont: {{ family: THEME.font, size: 10 }}, bodyFont: {{ family: THEME.font, size: 10 }} }};

// ── HELPERS ───────────────────────────────────────────────────────────────────
function speciesLabel(name) {{
  // Space-separated NCBI header: "NC_014696.1 Leggadina lakedownensis ..."
  const accRe = /^[A-Z]{{1,3}}_?\d{{5,}}(\.\d+)?$/;
  if (name.includes(' ')) {{
    const parts = name.split(' ');
    const acc = parts[0];
    const isRef = parts[parts.length-1] === 'REF';
    const meta = ['mitochondrion,','mitochondrion','isolate','strain','voucher',
                  'breed','complete','partial','genome','genome,','sequence'];
    const words = [];
    for (let i = 1; i < parts.length; i++) {{
      if (parts[i] === 'REF') continue;
      if (meta.includes(parts[i].toLowerCase())) break;
      words.push(parts[i]);
    }}
    const sp = words.length ? words.join(' ') : acc;
    return {{sp, acc: acc + (isRef ? ' · REF' : ''), isRef}};
  }}
  // Underscore-separated
  const p = name.split('_');
  const isRef = p[p.length-1] === 'REF';
  const parts = isRef ? p.slice(0, -1) : p;
  // Try to find accession from leading parts (e.g. NC_014696.1)
  let acc = null, spParts = parts;
  for (let i = 1; i < Math.min(parts.length, 3); i++) {{
    const cand = parts.slice(0, i+1).join('_');
    if (accRe.test(cand)) {{ acc = cand; spParts = parts.slice(i+1); break; }}
  }}
  if (!acc && accRe.test(parts[0])) {{ acc = parts[0]; spParts = parts.slice(1); }}
  // No accession — treat as species name
  if (!acc) {{
    const sp = spParts.length ? spParts.join(' ') : parts.join(' ');
    return {{sp, acc: '', isRef}};
  }}
  const tag = acc + (isRef ? ' · REF' : '');
  const sp = spParts.length ? spParts.join(' ') : acc;
  return {{sp, acc: tag, isRef}};
}}

function ssDisplay(val) {{
  if (val) return `<b style="color:var(--gold)">⚠ ${{val.length}}-nt (${{val}})</b>`;
  return '<b style="color:var(--green)">✓ None</b>';
}}

function gcClamp(seq) {{
  const tail = seq.slice(-5).toUpperCase();
  const gc = [...tail].filter(c => c==='G'||c==='C').length;
  return `${{gc}}/5 (${{tail}})`;
}}

// ── OFF-TARGET FILTER ───────────────────────────────────────────────────────
let hideOffTarget = false;

function primerHasOffTarget(primerSeq, isFwd) {{
  return PRIMER_SETS.some(ps => {{
    if (isFwd ? ps.fwd !== primerSeq : ps.rev !== primerSeq) return false;
    return ps.ecopcr && ps.ecopcr.total_hits > 0;
  }});
}}

document.getElementById('hideOffTargetToggle').addEventListener('change', function() {{
  hideOffTarget = this.checked;
  buildPrimerIndex();
  buildPrimerTrack();
  updateAll();
}});

// ── SPECIES FILTER ──────────────────────────────────────────────────────────
let selectedSpecies = new Set();

function primerAmplifiesSelected(primerSeq, isFwd) {{
  if (selectedSpecies.size === 0) return true;
  return PRIMER_SETS.some(ps => {{
    if (isFwd ? ps.fwd !== primerSeq : ps.rev !== primerSeq) return false;
    for (const sp of selectedSpecies) {{
      const h = ps.hit_summary[sp];
      if (!h || !h.has_amp) return false;
    }}
    return true;
  }});
}}

// ── FWD/REV PRIMER TRACKS ────────────────────────────────────────────────────
// Build unique FWD and REV lists with their gapped positions
let uniqueFwds = [];  // {{ seq, indices, gapped_start, gapped_end }}
let uniqueRevs = [];
let selectedFwdIdx = null;  // index into uniqueFwds
let selectedRevIdx = null;  // index into uniqueRevs

function buildPrimerIndex() {{
  // First pass: collect all known positions for each unique FWD and REV sequence
  const fwdPositions = {{}};  // seq -> {{ start, end }}
  const revPositions = {{}};
  PRIMER_SETS.forEach(ps => {{
    const gc = ps.gapped;
    if (!gc || gc.fwd_start === undefined || gc.fwd_start === null) return;
    if (!fwdPositions[ps.fwd]) fwdPositions[ps.fwd] = {{ start: gc.fwd_start, end: gc.fwd_end }};
    if (!revPositions[ps.rev]) revPositions[ps.rev] = {{ start: gc.rev_start, end: gc.rev_end }};
  }});

  // Second pass: build unique lists using positions from any primer set
  const fwdMap = {{}};
  const revMap = {{}};
  PRIMER_SETS.forEach((ps, i) => {{
    const fk = ps.fwd;
    const fp = fwdPositions[fk];
    if (fp && !fwdMap[fk]) {{
      fwdMap[fk] = {{ seq: fk, indices: [], start: fp.start, end: fp.end,
                      tm: ps.fwd_tm, gc_pct: ps.fwd_gc, len: ps.fwd_len }};
    }}
    if (fwdMap[fk]) fwdMap[fk].indices.push(i);

    const rk = ps.rev;
    const rp = revPositions[rk];
    if (rp && !revMap[rk]) {{
      revMap[rk] = {{ seq: rk, indices: [], start: rp.start, end: rp.end,
                      tm: ps.rev_tm, gc_pct: ps.rev_gc, len: ps.rev_len }};
    }}
    if (revMap[rk]) revMap[rk].indices.push(i);
  }});

  uniqueFwds = Object.values(fwdMap).sort((a,b) => a.start - b.start);
  uniqueRevs = Object.values(revMap).sort((a,b) => a.start - b.start);

  syncActivePrimerSet();
}}

function syncActivePrimerSet() {{
  if (selectedFwdIdx === null || selectedRevIdx === null) return;
  const fwd = uniqueFwds[selectedFwdIdx];
  const rev = uniqueRevs[selectedRevIdx];
  // Find primer set that matches both
  const match = PRIMER_SETS.findIndex(ps => ps.fwd === fwd.seq && ps.rev === rev.seq);
  if (match >= 0) activePSIndex = match;
}}

function currentPairMatches() {{
  if (selectedFwdIdx === null || selectedRevIdx === null) return false;
  const ps = PS();
  if (!ps) return false;
  return ps.fwd === uniqueFwds[selectedFwdIdx].seq && ps.rev === uniqueRevs[selectedRevIdx].seq;
}}

// Find the best hit_summary data for a given fwd or rev primer across all primer sets
function bestHitsForPrimer(seq, isFwd) {{
  // Returns a map of {{ seqName: {{ fwd_mm, rev_mm }} }} from any PS using this primer
  const result = {{}};
  PRIMER_SETS.forEach(ps => {{
    if (isFwd ? ps.fwd !== seq : ps.rev !== seq) return;
    SEQ_NAMES.forEach(name => {{
      const h = ps.hit_summary[name];
      if (!h) return;
      const prev = result[name];
      if (isFwd) {{
        if (!prev || h.fwd_mm < prev.fwd_mm) result[name] = {{ ...result[name], fwd_mm: h.fwd_mm, fwd_has_amp: h.has_amp }};
      }} else {{
        if (!prev || h.rev_mm < prev.rev_mm) result[name] = {{ ...result[name], rev_mm: h.rev_mm, rev_has_amp: h.has_amp }};
      }}
    }});
  }});
  return result;
}}

function isCompatible(fwdItem, revItem) {{
  // REV must start after FWD ends (no negative amplicons)
  if (revItem.start <= fwdItem.end) return false;
  const ampLen = revItem.end - fwdItem.start;
  if (ampLen > MAX_AMP_LEN || ampLen < MIN_AMP_LEN) return false;
  const deltaTm = Math.abs(fwdItem.tm - revItem.tm);
  return deltaTm <= 5;
}}

function hasPrimerSet(fwdItem, revItem) {{
  return PRIMER_SETS.some(ps => ps.fwd === fwdItem.seq && ps.rev === revItem.seq);
}}

function buildPrimerTrack() {{
  const fwdContainer = document.getElementById('fwdTrack');
  const revContainer = document.getElementById('revTrack');
  fwdContainer.innerHTML = '';
  revContainer.innerHTML = '';
  const trackEl = document.getElementById('geneTrack');

  const selFwd = selectedFwdIdx !== null ? uniqueFwds[selectedFwdIdx] : null;
  const selRev = selectedRevIdx !== null ? uniqueRevs[selectedRevIdx] : null;

  // FWD track — row pack using percentages so positions scale with the gene track
  const fwdSorted = uniqueFwds.map((f, fi) => ({{ ...f, fi }})).sort((a,b) => a.start - b.start);
  const fwdRows = [];
  fwdSorted.forEach(item => {{
    let placed = false;
    for (let row = 0; row < fwdRows.length; row++) {{
      if (fwdRows[row] <= item.start) {{ fwdRows[row] = item.end; item.row = row; placed = true; break; }}
    }}
    if (!placed) {{ item.row = fwdRows.length; fwdRows.push(item.end); }}
  }});
  fwdContainer.style.height = (Math.max(fwdRows.length, 1) * 16 + 2) + 'px';
  fwdContainer.style.position = 'relative';

  fwdSorted.forEach(item => {{
    const fi = item.fi;
    const f = uniqueFwds[fi];
    if (hideOffTarget && primerHasOffTarget(f.seq, true)) return;
    const leftPct = (f.start / ALIGN_LEN) * 100;
    const widthPct = Math.max(0.5, (f.len / ALIGN_LEN) * 100);
    const isSelected = fi === selectedFwdIdx;
    const sameTrackGrey = selFwd && !isSelected;
    const compatible = sameTrackGrey ? false : (selRev ? isCompatible(f, selRev) : true);
    const ampOk = primerAmplifiesSelected(f.seq, true);

    const bar = document.createElement('div');
    bar.className = 'primer-bar' + (isSelected ? ' active' : '');
    bar.style.cssText = `left:${{leftPct}}%;width:${{widthPct}}%;min-width:6px;top:${{item.row * 16}}px;height:14px;`
      + (compatible && ampOk ? `background:${{isSelected ? 'var(--primer-fwd)' : THEME.gold+'99'}};color:#fff;`
                             : `background:var(--bg4);color:var(--muted);border:1px solid ${{THEME.gold}}66;`);
    bar.title = `FWD: ${{f.seq}}\nTm: ${{f.tm}}°C · GC: ${{f.gc_pct}}% · ${{f.len}} bp\nUsed by ${{f.indices.length}} pair(s)`;
    bar.addEventListener('click', () => {{
      selectedFwdIdx = selectedFwdIdx === fi ? null : fi;
      syncActivePrimerSet();
      updateAll();
    }});
    fwdContainer.appendChild(bar);
  }});

  // REV track — row pack using percentages
  const revSorted = uniqueRevs.map((r, ri) => ({{ ...r, ri }})).sort((a,b) => a.start - b.start);
  const revRows = [];
  revSorted.forEach(item => {{
    let placed = false;
    for (let row = 0; row < revRows.length; row++) {{
      if (revRows[row] <= item.start) {{
        revRows[row] = item.end;
        item.row = row;
        placed = true;
        break;
      }}
    }}
    if (!placed) {{
      item.row = revRows.length;
      revRows.push(item.end);
    }}
  }});
  revContainer.style.height = (Math.max(revRows.length, 1) * 16 + 2) + 'px';
  revContainer.style.position = 'relative';

  revSorted.forEach(item => {{
    const ri = item.ri;
    const r = uniqueRevs[ri];
    if (hideOffTarget && primerHasOffTarget(r.seq, false)) return;
    const leftPct = (r.start / ALIGN_LEN) * 100;
    const widthPct = Math.max(0.5, (r.len / ALIGN_LEN) * 100);
    const isSelected = ri === selectedRevIdx;
    const sameTrackGrey = selRev && !isSelected;
    const compatible = sameTrackGrey ? false : (selFwd ? isCompatible(selFwd, r) : true);
    const ampOk = primerAmplifiesSelected(r.seq, false);

    const bar = document.createElement('div');
    bar.className = 'primer-bar' + (isSelected ? ' active' : '');
    bar.style.cssText = `left:${{leftPct}}%;width:${{widthPct}}%;min-width:6px;top:${{item.row * 16}}px;height:14px;`
      + (compatible && ampOk ? `background:${{isSelected ? 'var(--primer-rev)' : THEME.teal+'99'}};color:#fff;`
                             : `background:var(--bg4);color:var(--muted);border:1px solid ${{THEME.teal}}66;`);
    bar.title = `REV: ${{r.seq}}\nTm: ${{r.tm}}°C · GC: ${{r.gc_pct}}% · ${{r.len}} bp\nUsed by ${{r.indices.length}} pair(s)`;
    bar.addEventListener('click', () => {{
      selectedRevIdx = selectedRevIdx === ri ? null : ri;
      syncActivePrimerSet();
      updateAll();
    }});
    revContainer.appendChild(bar);
  }});
}}

function gcClamp(seq) {{
  const tail = seq.toUpperCase().slice(-5);
  return [...tail].filter(b => b === 'G' || b === 'C').length;
}}
function maxHomopoly(seq) {{
  let max = 1, run = 1;
  for (let i = 1; i < seq.length; i++) {{
    if (seq[i].toUpperCase() === seq[i-1].toUpperCase()) {{ run++; if (run > max) max = run; }}
    else run = 1;
  }}
  return max;
}}
function flagColor(val, lo, hi, invert) {{
  // invert=true means lower is worse (like dG: more negative = worse)
  if (invert) return val < lo ? 'var(--coral)' : val < hi ? 'var(--gold)' : 'var(--green)';
  return val >= lo && val <= hi ? 'var(--green)' : 'var(--gold)';
}}

function buildSelectedDetail() {{
  const detail = document.getElementById('selectedPrimerDetail');
  const title = document.getElementById('selectedPrimerTitle');

  const fwd = selectedFwdIdx !== null ? uniqueFwds[selectedFwdIdx] : null;
  const rev = selectedRevIdx !== null ? uniqueRevs[selectedRevIdx] : null;

  if (!fwd || !rev) {{
    title.textContent = 'Select a forward and reverse primer';
    detail.innerHTML = '<div class="no-amp-msg">Click primers on the tracks above</div>';
    return;
  }}

  const paired = currentPairMatches() && hasPrimerSet(fwd, rev);
  const compatible = isCompatible(fwd, rev);
  const ampLen = rev.end - fwd.start;
  const deltaTm = Math.abs(fwd.tm - rev.tm).toFixed(1);
  const tmIncompatible = parseFloat(deltaTm) > 5;
  const ps = paired ? PS() : null;
  const tested = ps && Object.values(ps.hit_summary).some(h => h.has_amp || h.fwd_mm < 99);
  const ampCount = tested ? Object.values(ps.hit_summary).filter(h => h.has_amp).length : null;
  const pairName = ps && tested ? ps.name : (tmIncompatible ? 'Incompatible pair' : 'Untested pair');

  // Get thermo values from the PS if available, else compute what we can
  const fHairpin = ps ? ps.fwd_hairpin : '—';
  const rHairpin = ps ? ps.rev_hairpin : '—';
  const fHomodimer = ps ? ps.fwd_homodimer : '—';
  const rHomodimer = ps ? ps.rev_homodimer : '—';
  const hetero = ps ? ps.heterodimer : '—';
  const fClamp = gcClamp(fwd.seq);
  const rClamp = gcClamp(rev.seq);
  const fHomopoly = maxHomopoly(fwd.seq);
  const rHomopoly = maxHomopoly(rev.seq);

  const flags = ps ? (ps.properties.flags || []) : [];
  const hasFlag = f => flags.includes(f);

  function sv(val, unit) {{ return typeof val === 'number' ? val.toFixed(1) + unit : val; }}
  function flagBg(flagName) {{ return hasFlag(flagName) ? 'background:var(--flag-bg);border-radius:4px;padding:1px 4px;' : ''; }}

  title.textContent = pairName.replace(/_/g, ' ');

  detail.innerHTML = `
    <div style="margin-bottom:8px;">
      <div class="sel-primer-seq" style="color:var(--primer-fwd)">▶ ${{fwd.seq}}</div>
      <div class="sel-primer-stats" style="margin:4px 0;">
        <span class="pstat" style="${{flagBg('')}}"><b>${{fwd.len}}</b> bp</span>
        <span class="pstat" style="${{flagBg('gc_fwd_out_of_range')}}">GC <b>${{fwd.gc_pct}}%</b></span>
        <span class="pstat" style="${{flagBg('tm_fwd_out_of_range')}}">Tm <b>${{fwd.tm}}°C</b></span>
        <span class="pstat" style="${{flagBg('hairpin_fwd')}}">hairpin <b>${{sv(fHairpin, '')}}</b></span>
        <span class="pstat" style="${{flagBg('homodimer_fwd')}}">homodimer <b>${{sv(fHomodimer, '')}}</b></span>
        <span class="pstat" style="${{flagBg('homopolymer_fwd')}}">homopoly <b>${{fHomopoly}}</b></span>
        <span class="pstat" style="${{flagBg('gc_clamp_fwd')}}">3′ GC <b>${{fClamp}}/5</b></span>
      </div>
    </div>
    <div style="margin-bottom:8px;">
      <div class="sel-primer-seq" style="color:var(--primer-rev)">◀ ${{rev.seq}}</div>
      <div class="sel-primer-stats" style="margin:4px 0;">
        <span class="pstat" style="${{flagBg('')}}"><b>${{rev.len}}</b> bp</span>
        <span class="pstat" style="${{flagBg('gc_rev_out_of_range')}}">GC <b>${{rev.gc_pct}}%</b></span>
        <span class="pstat" style="${{flagBg('tm_rev_out_of_range')}}">Tm <b>${{rev.tm}}°C</b></span>
        <span class="pstat" style="${{flagBg('hairpin_rev')}}">hairpin <b>${{sv(rHairpin, '')}}</b></span>
        <span class="pstat" style="${{flagBg('homodimer_rev')}}">homodimer <b>${{sv(rHomodimer, '')}}</b></span>
        <span class="pstat" style="${{flagBg('homopolymer_rev')}}">homopoly <b>${{rHomopoly}}</b></span>
        <span class="pstat" style="${{flagBg('gc_clamp_rev')}}">3′ GC <b>${{rClamp}}/5</b></span>
      </div>
    </div>
    <div style="background:var(--bg3);border-radius:6px;padding:8px 12px;display:flex;gap:16px;flex-wrap:wrap;">
      <span class="pstat" style="${{flagBg('delta_tm_too_high')}}">ΔTm <b style="color:${{parseFloat(deltaTm) <= 2 ? 'var(--green)' : parseFloat(deltaTm) <= 5 ? 'var(--gold)' : 'var(--coral)'}}">${{deltaTm}}°C</b></span>
      ${{ps ? `<span class="pstat" style="${{flagBg('heterodimer')}}">heterodimer <b>${{sv(hetero, '')}}</b></span>` : ''}}
      <span class="pstat">Amplicon <b style="color:var(--gold)">${{ampLen > 0 ? ampLen + ' bp' : '—'}}</b></span>
      <span class="pstat">Amplifies <b style="color:${{ampCount !== null && ampCount > 0 ? 'var(--green)' : ampCount === null ? 'var(--muted)' : 'var(--coral)'}}">${{ampCount !== null ? ampCount + '/' + SEQ_NAMES.length : 'Not tested'}}</b></span>
      ${{ps && ps.ecopcr ? `<span class="pstat">Off-target <b style="color:var(--teal)">${{ps.ecopcr.total_hits}}/${{ps.ecopcr.db_sequences}}</b></span>` : ''}}
      ${{tmIncompatible ? '<span class="pstat"><b style="color:var(--coral)">ΔTm incompatible</b></span>' : ''}}
    </div>
  `;
}}

function updateAll() {{
  buildPrimerTrack();
  buildSelectedDetail();
  buildHitTable();
  buildEcoPCR();
  buildAlignment();
  buildAmpLenChart();
}}

function selectPrimerSet(idx) {{
  // Find which FWD/REV indices match this primer set
  const ps = PRIMER_SETS[idx];
  const fi = uniqueFwds.findIndex(f => f.seq === ps.fwd);
  const ri = uniqueRevs.findIndex(r => r.seq === ps.rev);
  if (fi >= 0) selectedFwdIdx = fi;
  if (ri >= 0) selectedRevIdx = ri;
  activePSIndex = idx;
  updateAll();
}}

// ── HIT TABLE ────────────────────────────────────────────────────────────────
function buildHitTable() {{
  const tbody = document.getElementById('hitTableBody');
  tbody.innerHTML = '';
  const titleEl = document.getElementById('ampPredTitle');
  titleEl.textContent = selectedSpecies.size > 0
    ? `Amplification Prediction · ${{selectedSpecies.size}} species required`
    : 'Amplification Prediction · click species to filter primers';

  const paired = currentPairMatches();
  const ps = paired ? PS() : null;

  // Always gather per-primer best mismatch data (used for untested pairs AND
  // tested pairs where obipcr found no amplicon but individual primers do bind)
  const fwdHits = selectedFwdIdx !== null ? bestHitsForPrimer(uniqueFwds[selectedFwdIdx].seq, true) : null;
  const revHits = selectedRevIdx !== null ? bestHitsForPrimer(uniqueRevs[selectedRevIdx].seq, false) : null;

  // Estimate amplicon length from gapped coords
  const selF = selectedFwdIdx !== null ? uniqueFwds[selectedFwdIdx] : null;
  const selR = selectedRevIdx !== null ? uniqueRevs[selectedRevIdx] : null;
  const estAmp = (selF && selR && selR.end > selF.start) ? (selR.end - selF.start) : null;

  // Check if obipcr found ANY binding for this pair across all species.
  // If all species have fwd_mm=99 and rev_mm=99, obipcr effectively didn't test it
  // (e.g. amplicon outside the configured size range).
  const pairAllBlank = ps ? SEQ_NAMES.every(n => {{
    const h = ps.hit_summary[n];
    return h && !h.has_amp && h.fwd_mm >= 99 && h.rev_mm >= 99;
  }}) : true;

  SEQ_NAMES.forEach(name => {{
    const {{sp, acc, isRef}} = speciesLabel(name);
    const spStyle = isRef ? ' style="color:var(--gold)"' : '';
    const isSel = selectedSpecies.has(name);
    const tr = document.createElement('tr');
    if (isSel) tr.classList.add('species-selected');
    const selIcon = isSel ? '<span style="color:var(--green);margin-right:4px">●</span>' : '<span style="color:var(--muted);margin-right:4px">○</span>';

    // Per-primer mismatch data — always available regardless of pair
    const fh = fwdHits ? fwdHits[name] : null;
    const rh = revHits ? revHits[name] : null;
    const fmm = fh ? fh.fwd_mm : null;
    const rmm = rh ? rh.rev_mm : null;
    const fDisp = fmm !== null && fmm < 99 ? fmm : '—';
    const rDisp = rmm !== null && rmm < 99 ? rmm : '—';
    const fCls = fmm === null || fmm >= 99 ? 'mm-bad' : (fmm === 0 ? 'mm-0' : fmm <= 1 ? 'mm-1' : 'mm-bad');
    const rCls = rmm === null || rmm >= 99 ? 'mm-bad' : (rmm === 0 ? 'mm-0' : rmm <= 1 ? 'mm-1' : 'mm-bad');

    // Determine call: YES / NO / Not tested
    // A pair is only meaningfully tested if obipcr actually found primer binding
    // (not all mm=99). Pairs outside the amplicon range get all-99 and aren't real tests.
    const pairHit = ps && ps.hit_summary[name];
    const hasAmp = pairHit && pairHit.has_amp;
    const reallyTested = pairHit && !pairAllBlank;

    const ampDisp = hasAmp ? pairHit.amp_len + ' bp'
                   : estAmp ? `~${{estAmp}} bp` : '—';
    const ampStyle = hasAmp ? 'color:var(--gold)' : 'color:var(--muted)';

    let callHtml;
    if (hasAmp) {{
      callHtml = '<span style="color:var(--green);font-weight:700">✓ YES</span>';
    }} else if (reallyTested) {{
      callHtml = '<span style="color:var(--coral);font-weight:700">✗ NO</span>';
    }} else {{
      callHtml = '<span style="color:var(--muted);font-size:9px">Not tested</span>';
    }}

    tr.innerHTML = `
      <td>${{selIcon}}<span class="sp-name"${{spStyle}}>${{sp}}</span><span class="accession">${{acc}}</span></td>
      <td style="text-align:center"><span class="mm-badge ${{fCls}}">${{fDisp}}</span></td>
      <td style="text-align:center"><span class="mm-badge ${{rCls}}">${{rDisp}}</span></td>
      <td style="text-align:center;${{ampStyle}}">${{ampDisp}}</td>
      <td style="text-align:center">${{callHtml}}</td>`;
    tr.addEventListener('click', () => {{
      if (selectedSpecies.has(name)) selectedSpecies.delete(name);
      else selectedSpecies.add(name);
      buildHitTable();
      buildPrimerTrack();
    }});
    tbody.appendChild(tr);
  }});
}}

// ── ecoPCR PANEL (rebuilt on switch) ─────────────────────────────────────────
function buildEcoPCR() {{
  const section = document.getElementById('ecopcrSection');
  const ps = PS();
  const eco = ps ? ps.ecopcr : null;

  const tbody = document.getElementById('ecopcrTableBody');
  tbody.innerHTML = '';

  if (!eco || eco.total_hits === 0) {{
    document.getElementById('ecopcrTableTitle').textContent = 'ecoPCR Hits · No data';
    tbody.innerHTML = '<tr><td colspan="5" style="text-align:center;color:var(--muted);padding:20px;">No off-target database provided or no hits for this primer pair</td></tr>';
    return;
  }}

  document.getElementById('ecopcrTableTitle').textContent =
    `ecoPCR Hits · ${{eco.total_hits}} sequences from ${{eco.db_sequences}} in database`;

  eco.hits.sort((a,b) => (a.fwd_mm + a.rev_mm) - (b.fwd_mm + b.rev_mm)).forEach(h => {{
    const tr = document.createElement('tr');
    const fCls = h.fwd_mm === 0 ? 'mm-0' : h.fwd_mm <= 1 ? 'mm-1' : 'mm-bad';
    const rCls = h.rev_mm === 0 ? 'mm-0' : h.rev_mm <= 1 ? 'mm-1' : 'mm-bad';
    tr.innerHTML = `
      <td style="font-size:10px;white-space:nowrap">${{h.id}}</td>
      <td style="font-size:10px;color:var(--muted);max-width:300px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap">${{h.def}}</td>
      <td style="text-align:center"><span class="mm-badge ${{fCls}}">${{h.fwd_mm}}</span></td>
      <td style="text-align:center"><span class="mm-badge ${{rCls}}">${{h.rev_mm}}</span></td>
      <td style="text-align:center;color:var(--gold)">${{h.amp_len}} bp</td>`;
    tbody.appendChild(tr);
  }});
}}



// ── AMPLICON ALIGNMENT + PER-AMPLICON VARIABILITY ───────────────────────────
// Measure one monospace character width (cached)
let _charW = null;
function charWidth() {{
  if (_charW) return _charW;
  const span = document.createElement('span');
  span.style.cssText = 'font-family:var(--font-mono);font-size:10.5px;letter-spacing:1.2px;position:absolute;visibility:hidden;white-space:pre;';
  span.textContent = 'A';
  document.body.appendChild(span);
  _charW = span.getBoundingClientRect().width;
  document.body.removeChild(span);
  return _charW;
}}

function buildAlignment() {{
  const rows = document.getElementById('alnRows');
  const canvas = document.getElementById('varChart');
  rows.innerHTML = '';

  const fwd = selectedFwdIdx !== null ? uniqueFwds[selectedFwdIdx] : null;
  const rev = selectedRevIdx !== null ? uniqueRevs[selectedRevIdx] : null;

  if (!fwd || !rev || !Object.keys(GAPPED).length) {{
    rows.innerHTML = '<div style="padding:20px;text-align:center;color:var(--muted);font-size:12px;">Select a FWD and REV primer</div>';
    canvas.width = 0;
    return;
  }}

  const start = fwd.start;
  const end = rev.end;
  if (end <= start) {{
    rows.innerHTML = '<div style="padding:20px;text-align:center;color:var(--muted);font-size:12px;">REV is before FWD</div>';
    canvas.width = 0;
    return;
  }}

  const fwdEnd = fwd.end;
  const revStart = rev.start;
  const regionLen = end - start;
  const cw = charWidth();
  const totalSeqWidth = regionLen * cw;

  // Build unified rows — each row contains a sticky label + scrollable sequence
  SEQ_NAMES.forEach(name => {{
    const fullSeq = GAPPED[name];
    if (!fullSeq) return;
    const slice = fullSeq.substring(start, end).toUpperCase();

    const row = document.createElement('div');
    row.style.cssText = 'display:flex;align-items:center;';

    // Label — sticky so it stays visible during horizontal scroll
    const lbl = document.createElement('div');
    lbl.style.cssText = `font-size:10px;color:${{COLORS[name]}};padding-right:8px;text-align:right;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;line-height:1.7;width:140px;flex-shrink:0;position:sticky;left:0;z-index:1;background:var(--bg2);`;
    if (name.includes('REF')) lbl.style.fontWeight = '600';
    lbl.textContent = SHORT[name];
    row.appendChild(lbl);

    // Sequence
    const seqEl = document.createElement('div');
    seqEl.style.cssText = `font-family:var(--font-mono);font-size:10.5px;letter-spacing:1.2px;white-space:pre;line-height:1.7;width:${{totalSeqWidth}}px;flex-shrink:0;`;
    let html = '';
    for (let i = 0; i < slice.length; i++) {{
      const gPos = start + i;
      const c = slice[i];
      if (gPos < fwdEnd) {{
        html += `<span style="color:var(--primer-fwd);font-weight:700">${{c}}</span>`;
      }} else if (gPos >= revStart) {{
        html += `<span style="color:var(--primer-rev);font-weight:700">${{c}}</span>`;
      }} else {{
        html += `<span style="color:var(--text)">${{c}}</span>`;
      }}
    }}
    seqEl.innerHTML = html;
    row.appendChild(seqEl);
    rows.appendChild(row);
  }});

  // Draw variability directly on canvas — one bar per character position
  const allSeqs = SEQ_NAMES.map(n => GAPPED[n]).filter(Boolean);
  const dpr = window.devicePixelRatio || 1;
  canvas.width = Math.round(totalSeqWidth * dpr);
  canvas.height = Math.round(60 * dpr);
  canvas.style.width = totalSeqWidth + 'px';
  canvas.style.height = '60px';
  const ctx = canvas.getContext('2d');
  ctx.scale(dpr, dpr);
  const h = 60;

  for (let i = 0; i < regionLen; i++) {{
    const gPos = start + i;
    const counts = {{}};
    let total = 0;
    allSeqs.forEach(s => {{
      const b = s[gPos];
      if (b && b !== '-' && b !== '.') {{ counts[b] = (counts[b] || 0) + 1; total++; }}
    }});
    let pct = 0;
    if (total > 1) {{
      const maxCount = Math.max(...Object.values(counts));
      pct = ((total - maxCount) / total) * 100;
    }}
    const barH = (pct / 100) * h;
    const x = i * cw;
    // Colour: primer regions use primer colours, inner uses coral
    if (gPos < fwdEnd) ctx.fillStyle = THEME.gold + '99';
    else if (gPos >= revStart) ctx.fillStyle = THEME.teal + '99';
    else ctx.fillStyle = pct > 0 ? THEME.coral + 'cc' : THEME.coral + '18';
    ctx.fillRect(x, h - barH, cw - 0.5, barH);
  }}
}}

// ── FULL-ALIGNMENT VARIABILITY CHART (built once, shares x-axis with gene track) ──
function buildVariabilityChart() {{
  if (!Object.keys(GAPPED).length) return;
  const allSeqs = SEQ_NAMES.map(n => GAPPED[n]).filter(Boolean);
  if (allSeqs.length < 2) return;

  // Compute per-position variability: at each column, what fraction of sequences
  // differ from the most common base (consensus variability, not vs reference)
  const varScore = [];
  for (let i = 0; i < ALIGN_LEN; i++) {{
    const counts = {{}};
    let total = 0;
    allSeqs.forEach(s => {{
      const b = s[i];
      if (b && b !== '-' && b !== '.') {{
        counts[b] = (counts[b] || 0) + 1;
        total++;
      }}
    }});
    if (total <= 1) {{ varScore.push(0); continue; }}
    const maxCount = Math.max(...Object.values(counts));
    varScore.push(((total - maxCount) / total) * 100);
  }}

  // Gene track plugin — draws coloured gene segments in the bottom strip of the chart
  const geneTrackHeight = 28;
  const geneTrackPlugin = {{
    id: 'geneTrack',
    afterDraw(chart) {{
      const {{ ctx, chartArea: {{ left, right, bottom }}, scales: {{ x }} }} = chart;
      const trackTop = bottom + 4;
      const trackBot = trackTop + geneTrackHeight;
      ctx.save();
      GENES.forEach(g => {{
        const x0 = x.getPixelForValue(g.start);
        const x1 = x.getPixelForValue(g.end);
        const w = x1 - x0;
        ctx.fillStyle = g.color;
        ctx.fillRect(x0, trackTop, w, geneTrackHeight);
        // Label major regions (not tRNAs) if wide enough
        const isMajor = !g.name.startsWith('t') || g.name === 'tRNA';
        if (isMajor && w > 30) {{
          ctx.fillStyle = 'rgba(0,0,0,0.75)';
          ctx.font = '700 8.5px IBM Plex Mono, monospace';
          ctx.textAlign = 'center';
          ctx.textBaseline = 'middle';
          ctx.fillText(g.name, x0 + w / 2, trackTop + geneTrackHeight / 2);
        }}
      }});
      ctx.restore();
    }}
  }};

  // Tooltip: include gene name when hovering
  function geneAtPos(pos) {{
    for (const g of GENES) {{ if (pos >= g.start && pos < g.end) return g.name; }}
    return null;
  }}

  const ctx = document.getElementById('varFullChart').getContext('2d');
  new Chart(ctx, {{
    type: 'line',
    data: {{
      datasets: [{{
        data: varScore.map((v, i) => ({{x: i, y: v}})),
        borderColor: THEME.coral,
        backgroundColor: THEME.coral + '18',
        borderWidth: 1,
        pointRadius: 0,
        fill: true,
        tension: 0.1
      }}]
    }},
    plugins: [geneTrackPlugin],
    options: {{
      responsive: true, maintainAspectRatio: false,
      layout: {{ padding: {{ bottom: geneTrackHeight + 8, left: 0, right: 0, top: 0 }} }},
      scales: {{
        x: {{ type: 'linear', min: 0, max: ALIGN_LEN, display: false }},
        y: {{ min: 0, max: 100, display: false }}
      }},
      plugins: {{
        legend: {{ display: false }},
        tooltip: {{ ...CHART_TOOLTIP,
          callbacks: {{
            title: i => {{
              const pos = Math.round(i[0].parsed.x);
              const gene = geneAtPos(pos);
              return gene ? `Position ${{pos}} · ${{gene}}` : `Position ${{pos}}`;
            }},
            label: i => ` ${{i.parsed.y.toFixed(1)}}% variable`
          }}
        }}
      }}
    }}
  }});
}}

// ── AMPLICON LENGTH DISTRIBUTION (rebuilt on primer switch) ──────────────────
let ampLenChartObj = null;
function buildAmpLenChart() {{
  if (ampLenChartObj) {{ ampLenChartObj.destroy(); ampLenChartObj = null; }}
  const titleEl = document.getElementById('ampLenTitle');

  const fwd = selectedFwdIdx !== null ? uniqueFwds[selectedFwdIdx] : null;
  const rev = selectedRevIdx !== null ? uniqueRevs[selectedRevIdx] : null;

  if (!fwd || !rev || !Object.keys(GAPPED).length) {{
    titleEl.textContent = 'Amplicon Length Distribution';
    return;
  }}

  const start = fwd.start;
  const end = rev.end;
  if (end <= start) return;

  // For each sequence, count ungapped bases in the amplicon region
  const lengths = [];
  const labels = [];
  const colors = [];
  SEQ_NAMES.forEach(name => {{
    const fullSeq = GAPPED[name];
    if (!fullSeq) return;
    const slice = fullSeq.substring(start, end);
    const ungapped = slice.replace(/[-.]/g, '').length;
    if (ungapped > 0) {{
      lengths.push(ungapped);
      labels.push(SHORT[name]);
      colors.push(COLORS[name] + 'cc');
    }}
  }});

  if (lengths.length === 0) return;

  const minL = Math.min(...lengths);
  const maxL = Math.max(...lengths);
  titleEl.textContent = `Amplicon Length · ${{minL}}–${{maxL}} bp (${{lengths.length}} seqs)`;

  // Build histogram bins
  const binSize = Math.max(1, Math.ceil((maxL - minL + 1) / 15));
  const bins = {{}};
  lengths.forEach(l => {{
    const bin = Math.floor((l - minL) / binSize) * binSize + minL;
    bins[bin] = (bins[bin] || 0) + 1;
  }});
  const binKeys = Object.keys(bins).map(Number).sort((a,b) => a - b);

  const ctx = document.getElementById('ampLenChart').getContext('2d');
  ampLenChartObj = new Chart(ctx, {{
    type: 'bar',
    data: {{
      labels: binKeys.map(k => binSize === 1 ? `${{k}}` : `${{k}}–${{k+binSize-1}}`),
      datasets: [{{
        data: binKeys.map(k => bins[k]),
        backgroundColor: THEME.accent + 'b3',
        borderWidth: 0,
        borderRadius: 3,
        barPercentage: 1.0,
        categoryPercentage: 0.9
      }}]
    }},
    options: {{
      responsive: true, maintainAspectRatio: false,
      scales: {{
        x: {{ title: {{ display: true, text: 'Amplicon length (bp)', ...CHART_TICK }}, ticks: CHART_TICK_SM, grid: {{ display: false }} }},
        y: {{ title: {{ display: true, text: 'Sequences', ...CHART_TICK }}, ticks: {{ ...CHART_TICK, stepSize: 1 }}, grid: CHART_GRID }}
      }},
      plugins: {{
        legend: {{ display: false }},
        tooltip: {{ ...CHART_TOOLTIP,
          callbacks: {{
            label: i => ` ${{i.parsed.y}} sequence${{i.parsed.y !== 1 ? 's' : ''}}`
          }}
        }}
      }}
    }}
  }});
}}

// ── INIT ──────────────────────────────────────────────────────────────────────
buildVariabilityChart();
buildPrimerIndex();
buildPrimerTrack();
buildSelectedDetail();
buildHitTable();
buildEcoPCR();
buildAlignment();
buildAmpLenChart();
</script>
</body>
</html>'''

    return html


def main():
    parser = argparse.ArgumentParser(
        description="Generate HTML visualisation from blast_vis_data.json"
    )
    parser.add_argument("json_file", nargs="?", default="blast_vis_data.json", help="Input JSON file")
    parser.add_argument("--output", "-o", default="blast_msa_viz.html", help="Output HTML file")
    parser.add_argument("--config", "-c", default=None, help="Config YAML for report header")
    args = parser.parse_args()

    json_path = Path(args.json_file)
    if not json_path.exists():
        sys.exit(f"ERROR: {json_path} not found")

    print(f"Reading {json_path} ...")
    with open(json_path) as f:
        data = json.load(f)

    # Load config if provided
    header_cfg = None
    full_cfg = None
    if args.config:
        try:
            import yaml
            cfg_path = Path(args.config)
            if cfg_path.exists():
                with open(cfg_path) as f:
                    full_cfg = yaml.safe_load(f)
                header_cfg = full_cfg.get("header")
        except ImportError:
            pass

    print(f"  {data['meta']['n_sequences']} sequences · {data['meta']['alignment_length']} bp alignment")
    print(f"  {data['meta'].get('n_primer_sets', '?')} primer set(s)")
    print("Generating HTML ...")

    html = generate_html(data, header_cfg, full_cfg)

    out_path = Path(args.output)
    with open(out_path, "w") as f:
        f.write(html)

    print(f"Done. Output written to {out_path}")


if __name__ == "__main__":
    main()
