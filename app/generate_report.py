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
    "#f0a732", "#ffd166", "#38d9d9", "#ff6b6b", "#2ec4b6",
    "#ff8fa3", "#ff6b9d", "#c77dff", "#e040fb", "#f06292",
    "#ba68c8", "#4fc3f7", "#81c784", "#ffb74d", "#e57373",
    "#aed581", "#4dd0e1", "#9575cd", "#f06292", "#dce775",
]


def short_name(full: str) -> str:
    """NC_008135_Petaurus_breviceps_REF → P. breviceps (REF)"""
    parts = full.split("_")
    if parts[-1] == "REF":
        tag = "REF"
        species_parts = parts[1:-1] if len(parts) > 3 else parts[1:-1]
    else:
        tag = parts[0]
        species_parts = parts[1:]
    if len(parts) >= 3 and parts[1].isdigit():
        tag = "REF" if parts[-1] == "REF" else f"{parts[0]}_{parts[1]}"
        species_parts = parts[2:-1] if parts[-1] == "REF" else parts[2:]
    elif len(parts) >= 2 and not parts[1].isdigit():
        if parts[-1] == "REF":
            tag = "REF"
            species_parts = parts[1:-1]
        else:
            tag = parts[0]
            species_parts = parts[1:]
    if len(species_parts) >= 2:
        return f"{species_parts[0][0]}. {species_parts[1]} ({tag})"
    return full


def build_js_data(data: dict, short_names: dict, colors: dict) -> str:
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
        f"const MATRIX = {json.dumps(data['matrix'])};",
        f"const ID_VS_REF = {json.dumps(data['id_vs_ref'])};",
        f"const X_POS = {json.dumps(data['sliding_window']['x_positions'])};",
        f"const SLIDING = {json.dumps(data['sliding_window']['data'])};",
        f"const PRIMER_SETS = {json.dumps(ps_js)};",
        f"let activePSIndex = {next((i for i, ps in enumerate(ps_js) if any(v is not None for v in ps['amplicons'].values())), 0)};",
        "function PS() { return PRIMER_SETS[activePSIndex]; }",
    ]
    return "\n".join(lines)


def generate_html(data: dict, header_cfg: dict | None = None) -> str:
    """Generate the full HTML dashboard."""
    seq_names = data["seq_names"]
    meta = data["meta"]
    ref_name = meta["reference"]

    sn = {name: short_name(name) for name in seq_names}
    colors = {name: PALETTE[i % len(PALETTE)] for i, name in enumerate(seq_names)}

    js_data = build_js_data(data, sn, colors)

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
<script src="https://cdnjs.cloudflare.com/ajax/libs/Chart.js/4.4.1/chart.umd.min.js"></script>
<link href="https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@300;400;600;700&family=Syne:wght@400;600;700;800&display=swap" rel="stylesheet">
<style>
  :root {{
    --bg:     #0a0e14;
    --bg2:    #111620;
    --bg3:    #1a2030;
    --bg4:    #222b3a;
    --border: #2a3548;
    --text:   #dde4f0;
    --muted:  #6b7d99;
    --accent: #4fa3f7;
    --gold:   #f0b429;
    --cyan:   #38d9d9;
    --green:  #4cd964;
    --coral:  #ff6b6b;
    --purple: #a78bfa;
    --teal:   #2ec4b6;
    --primer-fwd: #f59e0b;
    --primer-rev: #ec4899;
    --amplicon:   rgba(250,204,21,0.12);
    --font-mono: 'JetBrains Mono', monospace;
    --font-head: 'Syne', sans-serif;
  }}
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{
    background: var(--bg);
    color: var(--text);
    font-family: var(--font-mono);
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
  .card.accent-border {{ border-color: rgba(245,158,11,0.4); }}
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
  .ps-card.active {{ border-color: var(--primer-fwd); background: rgba(245,158,11,0.08); }}
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
    overflow-x: auto; background: var(--bg); border-radius: 6px;
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
  .aln-diff   {{ color: var(--coral); font-weight: 700; background: rgba(255,107,107,0.12); border-radius: 2px; }}
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
  .hit-table tr.species-selected td {{ background: rgba(34,197,94,0.12); }}
  .hit-table tr.species-selected:hover td {{ background: rgba(34,197,94,0.2); }}
  .hit-table tr {{ cursor: pointer; user-select: none; }}
  .mm-badge {{
    display: inline-flex; align-items: center; justify-content: center;
    width: 22px; height: 22px; border-radius: 50%; font-size: 11px; font-weight: 700;
  }}
  .mm-0 {{ background: rgba(76,217,100,0.2); color: var(--green); }}
  .mm-1 {{ background: rgba(240,180,41,0.2); color: var(--gold); }}
  .mm-bad {{ background: rgba(255,107,107,0.2); color: var(--coral); }}
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
    <div class="chip">Primer sets: <b>{n_primer_sets}</b></div>
    <div class="chip">Ref: <b>{ref_acc}</b></div>
  </div>
</div>

<!-- SECTION 1: GENE TRACK + PRIMER LOCATIONS -->
<div class="section-label">01 · Mitogenome Annotation &amp; Primer Binding Locations</div>
<div class="card">
  <div class="card-title"><div class="dot" style="background:var(--accent)"></div>Gene Track ({align_len:,} bp alignment) · click a primer pair to select</div>
  <div class="track-container" style="position:relative;">
    <div class="gene-track" id="geneTrack" style="position:relative;overflow:visible;"></div>
  </div>
  <div class="gene-axis">{axis_spans}</div>
  <div style="display:flex;gap:20px;margin-top:8px;flex-wrap:wrap;align-items:center;">
    <div class="gene-legend" id="geneLegend"></div>
  </div>
  <!-- FWD primer track -->
  <div style="margin-top:14px;">
    <div style="font-size:9px;color:var(--primer-fwd);margin-bottom:4px;text-transform:uppercase;letter-spacing:1px;font-weight:600;">▶ Forward primers</div>
    <div id="fwdTrack" style="position:relative;"></div>
  </div>
  <!-- REV primer track -->
  <div style="margin-top:8px;">
    <div style="font-size:9px;color:var(--primer-rev);margin-bottom:4px;text-transform:uppercase;letter-spacing:1px;font-weight:600;">◀ Reverse primers</div>
    <div id="revTrack" style="position:relative;min-height:16px;"></div>
  </div>
</div>

<!-- Selected primer detail + hit table -->
<div class="grid grid-2" style="margin-top:14px;">
  <div class="card accent-border" id="selectedPrimerCard">
    <div class="card-title"><div class="dot" style="background:var(--primer-fwd)"></div><span id="selectedPrimerTitle">Selected Primer</span></div>
    <div id="selectedPrimerDetail"></div>
  </div>
  <div class="card">
    <div class="card-title"><div class="dot" style="background:var(--green)"></div><span id="ampPredTitle">Amplification Prediction · click species to filter primers</span></div>
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

<!-- ecoPCR PANEL (shown if data available) -->
<div id="ecopcrSection" style="display:none;">
  <div class="section-label">01b · ecoPCR In Silico Specificity</div>
  <div class="grid grid-2">
    <div class="card">
      <div class="card-title"><div class="dot" style="background:var(--teal)"></div><span id="ecopcrMmTitle">Mismatch Distribution</span></div>
      <div class="chart-h200"><canvas id="ecopcrMmChart"></canvas></div>
    </div>
    <div class="card">
      <div class="card-title"><div class="dot" style="background:var(--purple)"></div><span id="ecopcrLenTitle">Amplicon Length Distribution</span></div>
      <div class="chart-h200"><canvas id="ecopcrLenChart"></canvas></div>
    </div>
  </div>
  <div class="card" style="margin-top:14px;">
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

<!-- SECTION 2: AMPLICON ALIGNMENT -->
<div class="section-label">02 · Amplicon Sequence Alignment</div>
<div class="card">
  <div class="card-title"><div class="dot" style="background:var(--cyan)"></div>Per-base alignment · primers coloured · mismatches highlighted</div>
  <div class="aln-viewer" id="alnViewer"></div>
  <div style="display:flex;gap:20px;margin-top:10px;font-size:10px;flex-wrap:wrap;">
    <span><span style="color:var(--primer-fwd);font-weight:700">■</span> FWD primer</span>
    <span><span style="color:var(--primer-rev);font-weight:700">■</span> REV primer</span>
    <span><span style="color:var(--cyan);font-weight:700">■</span> Inner variable region</span>
    <span><span style="color:var(--coral);font-weight:700;background:rgba(255,107,107,0.15);padding:0 3px;">N</span> Mismatch vs reference</span>
  </div>
</div>

<div class="grid grid-2" style="margin-top:14px;">
  <div class="card">
    <div class="card-title"><div class="dot" style="background:var(--purple)"></div>Inner Region Identity vs Reference (%)</div>
    <div class="chart-h200"><canvas id="innerIdChart"></canvas></div>
  </div>
  <div class="card">
    <div class="card-title"><div class="dot" style="background:var(--coral)"></div>Per-position variability within amplicon</div>
    <div class="chart-h200"><canvas id="varChart"></canvas></div>
  </div>
</div>

<!-- SECTION 3: FULL MSA OVERVIEW -->
<div class="section-label">03 · Full Alignment Overview</div>
<div class="card">
  <div class="card-title"><div class="dot" style="background:var(--gold)"></div>Sliding Window Identity vs Reference ({meta['window']} bp window · {meta['step']} bp step)</div>
  <div class="chart-h240"><canvas id="slidingChart"></canvas></div>
  <div class="sp-legend" id="speciesLegend"></div>
</div>

<div class="grid grid-2" style="margin-top:14px;">
  <div class="card">
    <div class="card-title"><div class="dot" style="background:var(--green)"></div>Overall Identity vs Reference</div>
    <div class="id-bars" id="idBars"></div>
  </div>
  <div class="card">
    <div class="card-title"><div class="dot" style="background:var(--coral)"></div>Pairwise Identity Matrix (%)</div>
    <div style="overflow-x:auto;"><table class="heatmap-table" id="heatmapTable"></table></div>
  </div>
</div>

<div class="footnote">
  <em>{ref_species}</em> · {ref_acc} · {len(seq_names)} sequences · {n_primer_sets} primer set(s) analysed
</div>

<!-- JAVASCRIPT -->
<script>
// ── DATA ──────────────────────────────────────────────────────────────────────
{js_data}

// ── HELPERS ───────────────────────────────────────────────────────────────────
function speciesLabel(name) {{
  const p = name.split('_');
  const isRef = p[p.length-1] === 'REF';
  let sp, acc;
  if (p[1] && /^\\d+$/.test(p[1])) {{
    acc = p[0]+'_'+p[1] + (isRef ? ' · REF' : '');
    sp = p.slice(2, isRef ? -1 : undefined).join(' ');
  }} else {{
    acc = p[0] + (isRef ? ' · REF' : '');
    sp = p.slice(1, isRef ? -1 : undefined).join(' ');
  }}
  return {{sp, acc, isRef}};
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

  if (selectedFwdIdx === null && uniqueFwds.length > 0) selectedFwdIdx = 0;
  if (selectedRevIdx === null && uniqueRevs.length > 0) selectedRevIdx = 0;
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
  // Compatible if: rev starts after fwd ends, amplicon in range, delta Tm ≤ 5
  if (revItem.start <= fwdItem.end) return false;
  const ampLen = revItem.end - fwdItem.start;
  if (ampLen < 50 || ampLen > 500) return false;
  const deltaTm = Math.abs(fwdItem.tm - revItem.tm);
  if (deltaTm > 5) return false;
  return true;
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

  requestAnimationFrame(() => {{
    const w = trackEl.offsetWidth;

    // FWD track — row pack
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
      const left = (f.start / ALIGN_LEN) * w;
      const width = Math.max(8, ((f.end - f.start) / ALIGN_LEN) * w);
      const isSelected = fi === selectedFwdIdx;
      const compatible = selRev ? isCompatible(f, selRev) : true;
      const ampOk = primerAmplifiesSelected(f.seq, true);

      const bar = document.createElement('div');
      bar.className = 'primer-bar' + (isSelected ? ' active' : '');
      bar.style.cssText = `left:${{left}}px;width:${{width}}px;top:${{item.row * 16}}px;height:14px;`
        + (compatible && ampOk ? `background:${{isSelected ? 'var(--primer-fwd)' : 'rgba(245,158,11,0.6)'}};color:#000;`
                               : `background:var(--bg4);color:var(--muted);border:1px solid rgba(245,158,11,0.5);`);
      bar.textContent = width > 30 ? f.seq.substring(0,6)+'…' : '';
      bar.title = `FWD: ${{f.seq}}\nTm: ${{f.tm}}°C · GC: ${{f.gc_pct}}% · ${{f.len}} bp\nUsed by ${{f.indices.length}} pair(s)`;
      bar.addEventListener('click', () => {{
        selectedFwdIdx = selectedFwdIdx === fi ? null : fi;
        syncActivePrimerSet();
        updateAll();
      }});
      fwdContainer.appendChild(bar);
    }});

    // REV track — row pack since they can overlap
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
      const left = (r.start / ALIGN_LEN) * w;
      const width = Math.max(8, ((r.end - r.start) / ALIGN_LEN) * w);
      const isSelected = ri === selectedRevIdx;
      const compatible = selFwd ? isCompatible(selFwd, r) : true;
      const ampOk = primerAmplifiesSelected(r.seq, false);

      const bar = document.createElement('div');
      bar.className = 'primer-bar' + (isSelected ? ' active' : '');
      bar.style.cssText = `left:${{left}}px;width:${{width}}px;top:${{item.row * 16}}px;height:14px;`
        + (compatible && ampOk ? `background:${{isSelected ? 'var(--primer-rev)' : 'rgba(236,72,153,0.6)'}};color:#fff;`
                               : `background:var(--bg4);color:var(--muted);border:1px solid rgba(236,72,153,0.5);`);
      bar.textContent = width > 30 ? r.seq.substring(0,6)+'…' : '';
      bar.title = `REV: ${{r.seq}}\nTm: ${{r.tm}}°C · GC: ${{r.gc_pct}}% · ${{r.len}} bp\nUsed by ${{r.indices.length}} pair(s)`;
      bar.addEventListener('click', () => {{
        selectedRevIdx = selectedRevIdx === ri ? null : ri;
        syncActivePrimerSet();
        updateAll();
      }});
      revContainer.appendChild(bar);
    }});
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
  const ampCount = ps ? Object.values(ps.hit_summary).filter(h => h.has_amp).length : '?';
  const pairName = ps ? ps.name : (tmIncompatible ? 'Incompatible pair' : 'Untested pair');

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
  function flagBg(flagName) {{ return hasFlag(flagName) ? 'background:rgba(245,158,11,0.15);border-radius:4px;padding:1px 4px;' : ''; }}

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
      <span class="pstat">Amplicon <b style="color:var(--gold)">${{ampLen}} bp</b></span>
      <span class="pstat">Amplifies <b style="color:${{ampCount !== '?' && ampCount > 0 ? 'var(--green)' : ampCount === '?' ? 'var(--muted)' : 'var(--coral)'}}">${{ampCount}}${{ampCount !== '?' ? '/' + SEQ_NAMES.length : ''}}</b></span>
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
  rebuildInnerIdChart();
  rebuildVarChart();
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

  // If no exact pair, gather per-primer mismatch data independently
  const fwdHits = !paired && selectedFwdIdx !== null ? bestHitsForPrimer(uniqueFwds[selectedFwdIdx].seq, true) : null;
  const revHits = !paired && selectedRevIdx !== null ? bestHitsForPrimer(uniqueRevs[selectedRevIdx].seq, false) : null;

  SEQ_NAMES.forEach(name => {{
    const {{sp, acc, isRef}} = speciesLabel(name);
    const spStyle = isRef ? ' style="color:var(--gold)"' : '';
    const isSel = selectedSpecies.has(name);
    const tr = document.createElement('tr');
    if (isSel) tr.classList.add('species-selected');
    const selIcon = isSel ? '<span style="color:var(--green);margin-right:4px">●</span>' : '<span style="color:var(--muted);margin-right:4px">○</span>';

    if (ps) {{
      // Exact pair match — full data
      const h = ps.hit_summary[name];
      if (h.has_amp) {{
        const fCls = h.fwd_mm === 0 ? 'mm-0' : h.fwd_mm <= 1 ? 'mm-1' : 'mm-bad';
        const rCls = h.rev_mm === 0 ? 'mm-0' : h.rev_mm <= 1 ? 'mm-1' : 'mm-bad';
        tr.innerHTML = `
          <td>${{selIcon}}<span class="sp-name"${{spStyle}}>${{sp}}</span><span class="accession">${{acc}}</span></td>
          <td style="text-align:center"><span class="mm-badge ${{fCls}}">${{h.fwd_mm}}</span></td>
          <td style="text-align:center"><span class="mm-badge ${{rCls}}">${{h.rev_mm}}</span></td>
          <td style="text-align:center;color:var(--gold)">${{h.amp_len}} bp</td>
          <td style="text-align:center;color:var(--green);font-weight:700">✓ YES</td>`;
      }} else {{
        const fmm = h.fwd_mm >= 99 ? '—' : h.fwd_mm;
        const rmm = h.rev_mm >= 99 ? '—' : h.rev_mm;
        const fCls = h.fwd_mm >= 99 ? 'mm-bad' : (h.fwd_mm <= 1 ? 'mm-1' : 'mm-bad');
        const rCls = h.rev_mm >= 99 ? 'mm-bad' : (h.rev_mm <= 1 ? 'mm-1' : 'mm-bad');
        tr.innerHTML = `
          <td>${{selIcon}}<span class="sp-name"${{spStyle}}>${{sp}}</span><span class="accession">${{acc}}</span></td>
          <td style="text-align:center"><span class="mm-badge ${{fCls}}">${{fmm}}</span></td>
          <td style="text-align:center"><span class="mm-badge ${{rCls}}">${{rmm}}</span></td>
          <td style="text-align:center;color:var(--muted)">—</td>
          <td style="text-align:center;color:var(--coral);font-weight:700">✗ NO</td>`;
      }}
    }} else {{
      // No exact pair — show individual primer mismatch data
      const fh = fwdHits ? fwdHits[name] : null;
      const rh = revHits ? revHits[name] : null;
      const fmm = fh ? fh.fwd_mm : null;
      const rmm = rh ? rh.rev_mm : null;
      const fDisp = fmm !== null && fmm < 99 ? fmm : '—';
      const rDisp = rmm !== null && rmm < 99 ? rmm : '—';
      const fCls = fmm === null || fmm >= 99 ? 'mm-bad' : (fmm === 0 ? 'mm-0' : fmm <= 1 ? 'mm-1' : 'mm-bad');
      const rCls = rmm === null || rmm >= 99 ? 'mm-bad' : (rmm === 0 ? 'mm-0' : rmm <= 1 ? 'mm-1' : 'mm-bad');
      tr.innerHTML = `
        <td>${{selIcon}}<span class="sp-name"${{spStyle}}>${{sp}}</span><span class="accession">${{acc}}</span></td>
        <td style="text-align:center"><span class="mm-badge ${{fCls}}">${{fDisp}}</span></td>
        <td style="text-align:center"><span class="mm-badge ${{rCls}}">${{rDisp}}</span></td>
        <td style="text-align:center;color:var(--muted)">—</td>
        <td style="text-align:center;color:var(--muted);font-weight:700">?</td>`;
    }}
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
let ecopcrMmChartObj = null;
let ecopcrLenChartObj = null;
function buildEcoPCR() {{
  const section = document.getElementById('ecopcrSection');
  const ps = PS();
  const eco = ps.ecopcr;

  if (!eco || eco.total_hits === 0) {{
    section.style.display = 'none';
    return;
  }}
  section.style.display = '';

  const pctHit = (eco.total_hits / eco.db_sequences * 100).toFixed(1);
  document.getElementById('ecopcrMmTitle').textContent =
    `Mismatch Distribution · ${{eco.total_hits}}/${{eco.db_sequences}} amplify (${{pctHit}}%)`;
  document.getElementById('ecopcrLenTitle').textContent =
    `Amplicon Length Distribution · ${{eco.amp_len_min}}–${{eco.amp_len_max}} bp (median ${{eco.amp_len_median}})`;
  document.getElementById('ecopcrTableTitle').textContent =
    `ecoPCR Hits · ${{eco.total_hits}} sequences from ${{eco.db_sequences}} in database`;

  // ── Mismatch bar chart
  if (ecopcrMmChartObj) {{ ecopcrMmChartObj.destroy(); ecopcrMmChartObj = null; }}
  const mmLabels = Object.keys(eco.mismatch_counts).sort();
  const mmData = mmLabels.map(k => eco.mismatch_counts[k]);
  const mmColors = mmLabels.map(k => {{
    const [f,r] = k.split('+').map(Number);
    if (f+r === 0) return 'rgba(76,217,100,0.8)';
    if (f+r <= 1) return 'rgba(240,180,41,0.8)';
    if (f+r <= 2) return 'rgba(79,163,247,0.8)';
    return 'rgba(255,107,107,0.8)';
  }});
  const mmCtx = document.getElementById('ecopcrMmChart').getContext('2d');
  ecopcrMmChartObj = new Chart(mmCtx, {{
    type: 'bar',
    data: {{
      labels: mmLabels.map(k => k.replace('+', ' + ') + ' mm'),
      datasets: [{{ data: mmData, backgroundColor: mmColors, borderWidth: 0, borderRadius: 3 }}]
    }},
    options: {{
      responsive: true, maintainAspectRatio: false,
      scales: {{
        y: {{ title: {{ display: true, text: 'Sequences', color: '#6b7d99', font: {{ family: 'JetBrains Mono', size: 9 }} }}, ticks: {{ color: '#6b7d99', font: {{ family: 'JetBrains Mono', size: 9 }} }}, grid: {{ color: '#1a2030' }} }},
        x: {{ ticks: {{ color: '#6b7d99', font: {{ family: 'JetBrains Mono', size: 8 }}, maxRotation: 45 }}, grid: {{ display: false }} }}
      }},
      plugins: {{ legend: {{ display: false }},
        tooltip: {{ backgroundColor: '#111620', borderColor: '#2a3548', borderWidth: 1,
          titleFont: {{ family: 'JetBrains Mono', size: 10 }}, bodyFont: {{ family: 'JetBrains Mono', size: 10 }},
          callbacks: {{ label: i => ` ${{i.parsed.y}} sequences (${{(i.parsed.y/eco.total_hits*100).toFixed(1)}}%)` }}
        }}
      }}
    }}
  }});

  // ── Amplicon length histogram
  if (ecopcrLenChartObj) {{ ecopcrLenChartObj.destroy(); ecopcrLenChartObj = null; }}
  const lengths = eco.amp_lengths;
  const minL = Math.min(...lengths), maxL = Math.max(...lengths);
  const binSize = Math.max(1, Math.ceil((maxL - minL + 1) / 20));
  const bins = {{}};
  lengths.forEach(l => {{
    const bin = Math.floor((l - minL) / binSize) * binSize + minL;
    bins[bin] = (bins[bin] || 0) + 1;
  }});
  const binKeys = Object.keys(bins).map(Number).sort((a,b) => a-b);
  const lenCtx = document.getElementById('ecopcrLenChart').getContext('2d');
  ecopcrLenChartObj = new Chart(lenCtx, {{
    type: 'bar',
    data: {{
      labels: binKeys.map(k => binSize === 1 ? `${{k}}` : `${{k}}–${{k+binSize-1}}`),
      datasets: [{{ data: binKeys.map(k => bins[k]), backgroundColor: 'rgba(167,139,250,0.7)', borderWidth: 0, borderRadius: 3, barPercentage: 1.0, categoryPercentage: 0.9 }}]
    }},
    options: {{
      responsive: true, maintainAspectRatio: false,
      scales: {{
        y: {{ title: {{ display: true, text: 'Sequences', color: '#6b7d99', font: {{ family: 'JetBrains Mono', size: 9 }} }}, ticks: {{ color: '#6b7d99', font: {{ family: 'JetBrains Mono', size: 9 }} }}, grid: {{ color: '#1a2030' }} }},
        x: {{ title: {{ display: true, text: 'Amplicon length (bp, primers excluded)', color: '#6b7d99', font: {{ family: 'JetBrains Mono', size: 9 }} }}, ticks: {{ color: '#6b7d99', font: {{ family: 'JetBrains Mono', size: 8 }}, maxRotation: 45 }}, grid: {{ display: false }} }}
      }},
      plugins: {{ legend: {{ display: false }},
        tooltip: {{ backgroundColor: '#111620', borderColor: '#2a3548', borderWidth: 1,
          titleFont: {{ family: 'JetBrains Mono', size: 10 }}, bodyFont: {{ family: 'JetBrains Mono', size: 10 }}
        }}
      }}
    }}
  }});

  // ── Hits table
  const tbody = document.getElementById('ecopcrTableBody');
  tbody.innerHTML = '';
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

// ── GENE TRACK (built once) ──────────────────────────────────────────────────
function buildGeneTrack() {{
  const track = document.getElementById('geneTrack');
  const legend = document.getElementById('geneLegend');
  const seen = new Set();
  GENES.forEach(g => {{
    const pct = (g.end - g.start) / ALIGN_LEN * 100;
    const el = document.createElement('div');
    el.className = 'gene-seg';
    el.style.cssText = `width:${{pct}}%;background:${{g.color}};flex-shrink:0;`;
    el.textContent = pct > 3 ? g.name : '';
    el.title = `${{g.name}}: ${{g.start}}–${{g.end}} (${{g.end-g.start}} bp)`;
    track.appendChild(el);
    if (!seen.has(g.name)) {{
      seen.add(g.name);
      const li = document.createElement('div');
      li.className = 'leg-item';
      li.innerHTML = `<div class="leg-swatch" style="background:${{g.color}}"></div>${{g.name}}`;
      legend.appendChild(li);
    }}
  }});
}}


// ── AMPLICON ALIGNMENT (rebuilt on switch) ───────────────────────────────────
function buildAlignment() {{
  const container = document.getElementById('alnViewer');
  container.innerHTML = '';
  const paired = currentPairMatches();
  const ps = PS();
  const amps = ps.amplicons;
  const refAmp = amps[SEQ_NAMES[0]];
  if (!paired) {{
    container.innerHTML = '<div class="no-amp-msg">No amplicon data for this combination — select a tested FWD+REV pair to view alignment</div>';
    return;
  }}
  if (!refAmp) {{
    container.innerHTML = '<div class="no-amp-msg"><b>No amplicon</b> found in reference for this primer set</div>';
    return;
  }}
  const ref = refAmp.toUpperCase();
  const fwdLen = ps.fwd_len, revLen = ps.rev_len, ampLen = ref.length;

  // Ruler
  const rulerRow = document.createElement('div');
  rulerRow.className = 'aln-row aln-ruler-row';
  const rulerLabel = document.createElement('div');
  rulerLabel.className = 'aln-label';
  rulerLabel.style.cssText = 'font-size:9px;color:var(--muted);text-align:right;padding-right:10px;';
  rulerLabel.textContent = 'pos';
  rulerRow.appendChild(rulerLabel);
  const rulerSeq = document.createElement('div');
  rulerSeq.className = 'aln-seq';
  rulerSeq.style.cssText = 'font-size:9px;color:var(--muted);letter-spacing:1.2px;';
  let ruler = '';
  for (let i = 0; i < ampLen; i++) {{
    if (i % 10 === 0) {{ ruler += `<span style="color:var(--accent)">${{String(i).padStart(2,'0')}}</span>`; i += String(i).length - 1; }}
    else ruler += '·';
  }}
  rulerSeq.innerHTML = ruler;
  rulerRow.appendChild(rulerSeq);
  container.appendChild(rulerRow);

  SEQ_NAMES.forEach(name => {{
    const raw = amps[name];
    if (!raw) return;
    const amp = raw.toUpperCase();
    const row = document.createElement('div');
    row.className = 'aln-row';
    const label = document.createElement('div');
    label.className = 'aln-label' + (name.includes('REF') ? ' ref' : '');
    label.style.color = COLORS[name];
    label.textContent = SHORT[name];
    row.appendChild(label);
    const seqEl = document.createElement('div');
    seqEl.className = 'aln-seq';
    let html = '';
    for (let i = 0; i < amp.length; i++) {{
      const c = amp[i], r = ref[i];
      if (i < fwdLen) {{
        html += `<span class="aln-primer-fwd${{c !== r ? ' aln-diff' : ''}}">${{c}}</span>`;
      }} else if (i >= ampLen - revLen) {{
        html += `<span class="aln-primer-rev${{c !== r ? ' aln-diff' : ''}}">${{c}}</span>`;
      }} else {{
        if (c === r || name.includes('REF')) html += `<span class="aln-inner">${{c}}</span>`;
        else html += `<span class="aln-diff">${{c}}</span>`;
      }}
    }}
    seqEl.innerHTML = html;
    row.appendChild(seqEl);
    container.appendChild(row);
  }});
}}

// ── INNER ID CHART (rebuilt on switch) ───────────────────────────────────────
let innerIdChartObj = null;
function rebuildInnerIdChart() {{
  if (innerIdChartObj) {{ innerIdChartObj.destroy(); innerIdChartObj = null; }}
  if (!currentPairMatches()) {{
    document.getElementById('innerIdChart').getContext('2d').clearRect(0,0,9999,9999);
    return;
  }}
  const ps = PS();
  const refAmp = ps.amplicons[SEQ_NAMES[0]];
  if (!refAmp) {{
    document.getElementById('innerIdChart').getContext('2d').clearRect(0,0,9999,9999);
    return;
  }}
  const ref = refAmp.toUpperCase();
  const fwdLen = ps.fwd_len, revLen = ps.rev_len;
  const innerRef = ref.substring(fwdLen, ref.length - revLen);
  const innerIds = {{}};
  SEQ_NAMES.forEach(name => {{
    const a = ps.amplicons[name];
    if (!a) {{ innerIds[name] = 0; return; }}
    const inner = a.toUpperCase().substring(fwdLen, a.length - revLen);
    let m = 0, t = 0;
    for (let i = 0; i < innerRef.length; i++) {{ t++; if (inner[i] === innerRef[i]) m++; }}
    innerIds[name] = t > 0 ? Math.round(m / t * 1000) / 10 : 0;
  }});
  const ctx = document.getElementById('innerIdChart').getContext('2d');
  innerIdChartObj = new Chart(ctx, {{
    type: 'bar',
    data: {{
      labels: SEQ_NAMES.map(n => SHORT[n]),
      datasets: [{{ data: SEQ_NAMES.map(n => innerIds[n] || 0), backgroundColor: SEQ_NAMES.map(n => COLORS[n] + 'cc'), borderColor: SEQ_NAMES.map(n => COLORS[n]), borderWidth: 1, borderRadius: 3 }}]
    }},
    options: {{
      responsive: true, maintainAspectRatio: false,
      scales: {{
        y: {{ min: 80, max: 101, ticks: {{ color: '#6b7d99', font: {{ family: 'JetBrains Mono', size: 9 }} }}, grid: {{ color: '#1a2030' }} }},
        x: {{ ticks: {{ color: '#6b7d99', font: {{ family: 'JetBrains Mono', size: 8 }}, maxRotation: 40 }}, grid: {{ display: false }} }}
      }},
      plugins: {{ legend: {{ display: false }}, tooltip: {{ backgroundColor: '#111620', borderColor: '#2a3548', borderWidth: 1, titleFont: {{ family: 'JetBrains Mono', size: 10 }}, bodyFont: {{ family: 'JetBrains Mono', size: 10 }}, callbacks: {{ label: i => ` ${{i.parsed.y.toFixed(1)}}% identity (inner region)` }} }} }}
    }}
  }});
}}

// ── PER-POSITION VARIABILITY CHART (rebuilt on switch) ───────────────────────
let varChartObj = null;
function rebuildVarChart() {{
  if (varChartObj) {{ varChartObj.destroy(); varChartObj = null; }}
  if (!currentPairMatches()) {{
    document.getElementById('varChart').getContext('2d').clearRect(0,0,9999,9999);
    return;
  }}
  const ps = PS();
  const refAmp = ps.amplicons[SEQ_NAMES[0]];
  if (!refAmp) {{
    document.getElementById('varChart').getContext('2d').clearRect(0,0,9999,9999);
    return;
  }}
  const ref = refAmp.toUpperCase();
  const fwdLen = ps.fwd_len, revLen = ps.rev_len, ampLen = ref.length;
  const others = SEQ_NAMES.slice(1).map(n => ps.amplicons[n] ? ps.amplicons[n].toUpperCase() : null).filter(Boolean);
  const varScore = [];
  for (let i = 0; i < ampLen; i++) {{
    const diffs = others.filter(s => s[i] && s[i] !== ref[i]).length;
    varScore.push(others.length > 0 ? diffs / others.length * 100 : 0);
  }}
  const positions = Array.from({{length: ampLen}}, (_, i) => i);
  const ctx = document.getElementById('varChart').getContext('2d');
  varChartObj = new Chart(ctx, {{
    type: 'bar',
    data: {{
      labels: positions,
      datasets: [{{
        data: varScore,
        backgroundColor: positions.map((_, i) => {{
          if (i < fwdLen) return 'rgba(245,158,11,0.7)';
          if (i >= ampLen - revLen) return 'rgba(236,72,153,0.7)';
          return varScore[i] > 0 ? 'rgba(56,217,217,0.8)' : 'rgba(56,217,217,0.15)';
        }}),
        borderWidth: 0, barPercentage: 1.0, categoryPercentage: 1.0
      }}]
    }},
    options: {{
      responsive: true, maintainAspectRatio: false,
      scales: {{
        y: {{ min: 0, max: 100, title: {{ display: true, text: '% seqs differ', color: '#6b7d99', font: {{ family: 'JetBrains Mono', size: 9 }} }}, ticks: {{ color: '#6b7d99', font: {{ family: 'JetBrains Mono', size: 9 }} }}, grid: {{ color: '#1a2030' }} }},
        x: {{ ticks: {{ display: false }}, grid: {{ display: false }} }}
      }},
      plugins: {{
        legend: {{ display: false }},
        tooltip: {{ backgroundColor: '#111620', borderColor: '#2a3548', borderWidth: 1, titleFont: {{ family: 'JetBrains Mono', size: 10 }}, bodyFont: {{ family: 'JetBrains Mono', size: 10 }},
          callbacks: {{
            title: i => {{ const pos = i[0].dataIndex; const region = pos < fwdLen ? 'FWD primer' : pos >= ampLen - revLen ? 'REV primer' : 'Variable region'; return `Pos ${{pos}} [${{region}}]`; }},
            label: i => ` ${{i.parsed.y.toFixed(0)}}% of seqs differ vs REF`
          }}
        }}
      }}
    }}
  }});
}}

// ── SLIDING WINDOW CHART (built once) ────────────────────────────────────────
let slidingChart = null;
function buildSlidingChart() {{
  const ctx = document.getElementById('slidingChart').getContext('2d');
  const datasets = Object.entries(SLIDING).map(([name, vals]) => ({{
    label: SHORT[name],
    data: X_POS.map((x,i) => vals[i] !== null ? {{x,y:vals[i]}} : null).filter(Boolean),
    borderColor: COLORS[name], backgroundColor: COLORS[name]+'18',
    borderWidth: 1.5, pointRadius: 0, tension: 0.4, spanGaps: false
  }}));
  slidingChart = new Chart(ctx, {{
    type: 'line', data: {{ datasets }},
    options: {{
      responsive: true, maintainAspectRatio: false,
      interaction: {{ mode: 'index', intersect: false }},
      scales: {{
        x: {{ type: 'linear', min: 0, max: ALIGN_LEN, title: {{ display: true, text: 'Alignment Position (bp)', color: '#6b7d99', font: {{ family: 'JetBrains Mono', size: 9 }} }}, grid: {{ color: '#1a2030' }}, ticks: {{ color: '#6b7d99', font: {{ family: 'JetBrains Mono', size: 9 }}, maxTicksLimit: 12 }} }},
        y: {{ min: 75, max: 100, title: {{ display: true, text: 'Identity (%)', color: '#6b7d99', font: {{ family: 'JetBrains Mono', size: 9 }} }}, grid: {{ color: '#1a2030' }}, ticks: {{ color: '#6b7d99', font: {{ family: 'JetBrains Mono', size: 9 }} }} }}
      }},
      plugins: {{
        legend: {{ display: false }},
        tooltip: {{ backgroundColor: '#111620', borderColor: '#2a3548', borderWidth: 1, titleFont: {{ family: 'JetBrains Mono', size: 10 }}, bodyFont: {{ family: 'JetBrains Mono', size: 9 }},
          callbacks: {{
            title: i => {{ const pos = i[0].parsed.x; const gc = PS().gapped; const inAmp = gc && gc.amp_start && pos >= gc.amp_start && pos <= gc.amp_end; return `Pos: ${{pos}} bp${{inAmp ? ' ← amplicon region' : ''}}`; }},
            label: i => ` ${{i.dataset.label}}: ${{i.parsed.y.toFixed(1)}}%`
          }}
        }}
      }}
    }}
  }});
  const leg = document.getElementById('speciesLegend');
  Object.keys(SLIDING).forEach(name => {{
    const d = document.createElement('div');
    d.className = 'sp-leg-item';
    d.innerHTML = `<div class="sp-dot" style="background:${{COLORS[name]}}"></div>${{SHORT[name]}}`;
    d.addEventListener('click', () => {{
      const ds = slidingChart.data.datasets.find(x => x.label === SHORT[name]);
      if (ds) {{ ds.hidden = !ds.hidden; d.classList.toggle('inactive', ds.hidden); slidingChart.update(); }}
    }});
    leg.appendChild(d);
  }});
}}

// ── IDENTITY BARS (built once) ───────────────────────────────────────────────
function buildIdBars() {{
  const c = document.getElementById('idBars');
  SEQ_NAMES.forEach(n => {{
    const v = ID_VS_REF[n];
    const pct = (v - 83) / 17 * 100;
    const row = document.createElement('div');
    row.className = 'id-bar-row';
    row.innerHTML = `<div class="id-bar-label" style="color:${{COLORS[n]}}">${{SHORT[n]}}</div>
      <div class="id-bar-track"><div class="id-bar-fill" style="width:${{pct}}%;background:${{COLORS[n]}};"></div></div>
      <div class="id-bar-val" style="color:${{COLORS[n]}}">${{v.toFixed(2)}}%</div>`;
    c.appendChild(row);
  }});
}}

// ── HEATMAP (built once) ─────────────────────────────────────────────────────
function idToColor(v) {{
  const stops = [[83,[31,41,55]],[87,[37,99,148]],[91,[56,189,210]],[96,[250,173,60]],[100,[240,167,50]]];
  let i = 0;
  while (i < stops.length-2 && v > stops[i+1][0]) i++;
  const [v0,c0] = stops[i], [v1,c1] = stops[i+1];
  const t = Math.max(0,Math.min(1,(v-v0)/(v1-v0)));
  const r=Math.round(c0[0]+t*(c1[0]-c0[0])), g=Math.round(c0[1]+t*(c1[1]-c0[1])), b=Math.round(c0[2]+t*(c1[2]-c0[2]));
  const lum=(0.299*r+0.587*g+0.114*b)/255;
  return {{bg:`rgb(${{r}},${{g}},${{b}})`,txt:lum>0.5?'#111':'#eee'}};
}}
function buildHeatmap() {{
  const t = document.getElementById('heatmapTable');
  const ax = SEQ_NAMES.map(n => {{
    const s = SHORT[n];
    if (s.includes('REF')) return 'REF';
    return s.split(' ')[0].charAt(0)+'.'+s.split(' ')[1].substring(0,4);
  }});
  let h = '<tr><th></th>' + ax.map(l => `<th title="${{l}}">${{l.substring(0,6)}}</th>`).join('') + '</tr>';
  MATRIX.forEach((row,i) => {{
    h += `<tr><th class="row-head" style="color:${{COLORS[SEQ_NAMES[i]]}}">${{ax[i]}}</th>`;
    row.forEach(v => {{ const c=idToColor(v); h += `<td style="background:${{c.bg}};color:${{c.txt}}" title="${{v}}%">${{v.toFixed(0)}}</td>`; }});
    h += '</tr>';
  }});
  t.innerHTML = h;
}}

// ── INIT ──────────────────────────────────────────────────────────────────────
buildGeneTrack();
buildPrimerIndex();
buildPrimerTrack();
buildSelectedDetail();
buildHitTable();
buildEcoPCR();
buildAlignment();
rebuildInnerIdChart();
rebuildVarChart();
buildSlidingChart();
buildIdBars();
buildHeatmap();
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

    # Load header config if provided
    header_cfg = None
    if args.config:
        try:
            import yaml
            cfg_path = Path(args.config)
            if cfg_path.exists():
                with open(cfg_path) as f:
                    raw = yaml.safe_load(f)
                header_cfg = raw.get("header")
        except ImportError:
            pass

    print(f"  {data['meta']['n_sequences']} sequences · {data['meta']['alignment_length']} bp alignment")
    print(f"  {data['meta'].get('n_primer_sets', '?')} primer set(s)")
    print("Generating HTML ...")

    html = generate_html(data, header_cfg)

    out_path = Path(args.output)
    with open(out_path, "w") as f:
        f.write(html)

    print(f"Done. Output written to {out_path}")


if __name__ == "__main__":
    main()
