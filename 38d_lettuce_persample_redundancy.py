"""
38d_lettuce_persample_redundancy.py
====================================
Per-sample functional redundancy + N-cycle completeness for lettuce.

Mirrors pepper's `14_functional_redundancy.py` + `16_ncycle_pgp_specificity.py`.
Current `38_veggie_lettuce_figures.py` only reports group-level genus counts;
this script computes per-sample distributions enabling MWU / Cliff's δ.

Also computes:
  - N-cycle completeness per sample (fraction of 5 FAPROTAX N-cycle steps
    with ≥1 contributing genus at non-zero RA)
  - Aggregate + mission-stratified
  - Pepper comparable format

Requires lettuce FAPROTAX_report.txt — if missing, the script runs FAPROTAX
on the DADA2-exported genus table.

Output:
  - acta_astronautica/lettuce/lettuce_persample_redundancy.csv
  - acta_astronautica/lettuce/lettuce_redundancy_summary.csv
  - acta_astronautica/lettuce/lettuce_ncycle_completeness_persample.csv
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config

LETTUCE_DIR = config.VEGGIE_LETTUCE_DIR
EXPORTS_DIR = LETTUCE_DIR / "exports"
OUT_DIR = config.AA_DIR / "lettuce"
OUT_DIR.mkdir(parents=True, exist_ok=True)

FAPROTAX_REPORT = OUT_DIR / "lettuce_FAPROTAX_report.txt"
FAPROTAX_TABLE = OUT_DIR / "lettuce_FAPROTAX_functional_table.tsv"

NCYCLE_STEPS = [
    'nitrogen_fixation',
    'nitrification',
    'denitrification',
    'nitrate_reduction',
    'ureolysis',
]

FOCUS_FUNCTIONS = [
    'nitrification',
    'nitrogen_fixation',
    'denitrification',
    'nitrate_reduction',
    'fermentation',
    'aerobic_chemoheterotrophy',
    'ureolysis',
    'cellulolysis',
    'methylotrophy',
    'sulfate_respiration',
]


def parse_genus(taxon_str):
    if pd.isna(taxon_str):
        return "Unassigned"
    for part in reversed(str(taxon_str).split(";")):
        p = part.strip()
        if p.startswith("g__"):
            name = p.split("__", 1)[-1].strip()
            if name and name not in ("uncultured", "metagenome", "unidentified"):
                return name
    return "Unassigned"


def parse_faprotax_report(report_path):
    """Return dict: function_name → set of genus names."""
    func2genera = {}
    current_func = None
    with open(report_path) as fh:
        for line in fh:
            line = line.rstrip()
            m = re.match(r'^# (\w+) \(\d+ records\):', line)
            if m:
                current_func = m.group(1)
                func2genera[current_func] = set()
                continue
            if current_func and line.startswith('    '):
                for part in line.split(';'):
                    part = part.strip()
                    if part.startswith('g__'):
                        name = part.split('__', 1)[-1].strip()
                        if name and name not in ('uncultured', 'metagenome',
                                                 'unidentified', ''):
                            func2genera[current_func].add(name)
    return func2genera


def build_genus_table(ft, tax):
    """Build genus-level relative-abundance (%) table for lettuce."""
    shared = ft.index.intersection(tax.index)
    ft = ft.loc[shared]
    genus = tax.loc[shared, "Taxon"].apply(parse_genus)
    ft = ft.copy()
    ft["Genus"] = genus.values
    g = ft.groupby("Genus").sum()
    g = g.drop(index="Unassigned", errors="ignore")
    col_sums = g.sum(axis=0).replace(0, 1)
    return g.div(col_sums, axis=1) * 100


def build_faprotax_input(genus_ra, tax, out_path):
    """
    Build a FAPROTAX-compatible TSV: first column = genus-level lineage,
    subsequent columns = samples (counts or RA).

    FAPROTAX expects:
      col 1: Taxonomy string (semicolon-delimited lineage with d__/p__/.../g__)
      cols 2..: per-sample abundances
    """
    # Reconstruct lineage string per genus: pick first taxonomy row per genus
    shared = genus_ra.index
    tax_by_genus = {}
    for aid, row in tax.iterrows():
        taxon_str = row.get("Taxon", "")
        if pd.isna(taxon_str):
            continue
        g = parse_genus(taxon_str)
        if g in shared and g not in tax_by_genus:
            tax_by_genus[g] = taxon_str
    df = genus_ra.copy()
    df.index = df.index.map(lambda g: tax_by_genus.get(g, f"d__Bacteria;g__{g}"))
    df.index.name = "OTU_ID"
    df.to_csv(out_path, sep="\t")


def run_faprotax(input_tsv, output_table, report):
    """Invoke FAPROTAX's collapse_table.py to classify genera into functions."""
    script = config.FAPROTAX_COLLAPSE_SCRIPT
    db = config.FAPROTAX_DB
    if not script.exists() or not db.exists():
        raise FileNotFoundError(f"FAPROTAX tool not found at {script} or DB at {db}")
    cmd = [
        sys.executable, str(script),
        "-i", str(input_tsv),
        "-o", str(output_table),
        "-g", str(db),
        "-r", str(report),
        "--row_names_are_in_column", "OTU_ID",
        "-f",
        "-v",
    ]
    print(f"\n  Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def cliffs_delta(a, b):
    a = np.asarray(a); b = np.asarray(b)
    if len(a) == 0 or len(b) == 0:
        return np.nan
    dom = np.sum(a[:, None] > b[None, :]) - np.sum(a[:, None] < b[None, :])
    return dom / (len(a) * len(b))


def hedges_g(a, b):
    a = np.asarray(a); b = np.asarray(b)
    if len(a) < 2 or len(b) < 2:
        return np.nan
    na, nb = len(a), len(b)
    sa2, sb2 = np.var(a, ddof=1), np.var(b, ddof=1)
    sp = np.sqrt(((na - 1) * sa2 + (nb - 1) * sb2) / (na + nb - 2))
    if sp == 0:
        return np.nan
    d = (np.mean(a) - np.mean(b)) / sp
    j = 1 - 3 / (4 * (na + nb) - 9)
    return d * j


def compute_redundancy(genus_ra, func2genera, functions):
    """Returns DataFrame: samples × functions, each cell = # contributing genera."""
    records = {}
    for func in functions:
        assigned = func2genera.get(func, set())
        present = genus_ra.index.intersection(assigned)
        if len(present) == 0:
            records[func] = pd.Series(0, index=genus_ra.columns)
        else:
            records[func] = (genus_ra.loc[present] > 0).sum(axis=0)
    return pd.DataFrame(records)


def compute_completeness(genus_ra, func2genera, steps=NCYCLE_STEPS):
    """Per-sample completeness = fraction of N-cycle steps with ≥1 contributing genus."""
    rows = []
    for s in genus_ra.columns:
        present_steps = []
        for step in steps:
            assigned = func2genera.get(step, set())
            in_table = genus_ra.index.intersection(assigned)
            if len(in_table) == 0:
                continue
            if (genus_ra.loc[in_table, s] > 0).any():
                present_steps.append(step)
        rows.append({
            "sample_id":       s,
            "n_steps_present": len(present_steps),
            "completeness":    round(100 * len(present_steps) / len(steps), 1),
            "present_steps":   ";".join(present_steps),
        })
    return pd.DataFrame(rows)


def compare_groups(redund_df, meta, functions, subset_name):
    sp_idx = meta.index[meta['Study_Group'] == 'Space_Flight']
    tr_idx = meta.index[meta['Study_Group'] == 'Terrestrial_Soil']
    rows = []
    for func in functions:
        if func not in redund_df.columns:
            continue
        sp = redund_df.loc[redund_df.index.intersection(sp_idx), func].dropna().values
        tr = redund_df.loc[redund_df.index.intersection(tr_idx), func].dropna().values
        if len(sp) == 0 or len(tr) == 0:
            continue
        try:
            stat, p = mannwhitneyu(sp, tr, alternative="two-sided")
        except ValueError:
            stat, p = np.nan, np.nan
        rows.append({
            "Subset":        subset_name,
            "Function":      func,
            "Space_n":       len(sp),
            "Terr_n":        len(tr),
            "Space_mean":    round(float(np.mean(sp)), 3),
            "Space_SD":      round(float(np.std(sp, ddof=1)) if len(sp) > 1 else 0.0, 3),
            "Terr_mean":     round(float(np.mean(tr)), 3),
            "Terr_SD":       round(float(np.std(tr, ddof=1)) if len(tr) > 1 else 0.0, 3),
            "MWU_p":         round(p, 4) if not np.isnan(p) else "NA",
            "Cliffs_delta":  round(cliffs_delta(sp, tr), 3),
            "Hedges_g":      round(hedges_g(sp, tr), 3),
        })
    return pd.DataFrame(rows)


def main():
    print("=" * 60)
    print("38d: Per-sample functional redundancy + N-cycle completeness")
    print("=" * 60)

    ft = pd.read_csv(EXPORTS_DIR / "feature-table.tsv", sep="\t", skiprows=1, index_col=0)
    tax = pd.read_csv(EXPORTS_DIR / "taxonomy.tsv", sep="\t", index_col=0)
    meta = pd.read_csv(LETTUCE_DIR / "metadata_harmonized.tsv", sep="\t", index_col=0)

    print(f"\n  Feature table: {ft.shape[0]} ASVs × {ft.shape[1]} samples")
    print(f"  Metadata:      {len(meta)} samples")

    # Restrict to root-zone samples (all 8 are root)
    sample_ids = [s for s in meta.index if s in ft.columns]
    ft = ft[sample_ids]

    # Build genus RA table
    genus_ra = build_genus_table(ft, tax)
    print(f"  Genus-level table: {genus_ra.shape[0]} genera × {genus_ra.shape[1]} samples")

    # Run FAPROTAX if report not present
    if not FAPROTAX_REPORT.exists():
        print("\n[1] Running FAPROTAX on DADA2-based lettuce genus table...")
        fa_input = OUT_DIR / "lettuce_faprotax_input.tsv"
        build_faprotax_input(genus_ra, tax, fa_input)
        run_faprotax(fa_input, FAPROTAX_TABLE, FAPROTAX_REPORT)
    else:
        print("\n[1] FAPROTAX report cached, loading...")

    func2genera = parse_faprotax_report(FAPROTAX_REPORT)
    print(f"  Parsed {len(func2genera)} functions from FAPROTAX report")
    for f in NCYCLE_STEPS:
        print(f"    {f}: {len(func2genera.get(f, set()))} genera")

    # ── Per-sample redundancy ───────────────────────────────────────────────
    print("\n[2] Computing per-sample functional redundancy...")
    redund = compute_redundancy(genus_ra, func2genera, FOCUS_FUNCTIONS)
    redund.insert(0, "Study_Group", meta.loc[redund.index, "Study_Group"].values)
    redund.insert(1, "OSD", meta.loc[redund.index, "OSD"].values)
    redund.insert(2, "Mission", meta.loc[redund.index, "Mission"].values)
    redund.to_csv(OUT_DIR / "lettuce_persample_redundancy.csv")
    print(f"  Saved: {OUT_DIR}/lettuce_persample_redundancy.csv")
    print(redund.to_string())

    # ── N-cycle completeness ────────────────────────────────────────────────
    print("\n[3] Computing per-sample N-cycle completeness...")
    comp = compute_completeness(genus_ra, func2genera)
    comp["Study_Group"] = meta.loc[comp["sample_id"], "Study_Group"].values
    comp["OSD"] = meta.loc[comp["sample_id"], "OSD"].values
    comp["Mission"] = meta.loc[comp["sample_id"], "Mission"].values
    comp.to_csv(OUT_DIR / "lettuce_ncycle_completeness_persample.csv", index=False)
    print(f"  Saved: {OUT_DIR}/lettuce_ncycle_completeness_persample.csv")
    print(comp.to_string(index=False))

    # ── Statistical comparisons: aggregate + mission-stratified ─────────────
    print("\n[4] Statistical comparisons...")
    # Drop metadata cols for redundancy numeric table
    redund_numeric = redund[FOCUS_FUNCTIONS]

    all_summary = []
    summary_agg = compare_groups(redund_numeric, meta, FOCUS_FUNCTIONS, "aggregate_n4x4")
    all_summary.append(summary_agg)
    print(f"\n-- Aggregate --")
    print(summary_agg.to_string(index=False))

    for mission, label in [("OSD-267", "VEG-01A"), ("OSD-269", "VEG-03A")]:
        mis_ids = meta.index[meta["OSD"] == mission]
        summary = compare_groups(redund_numeric.loc[redund_numeric.index.intersection(mis_ids)],
                                  meta.loc[mis_ids], FOCUS_FUNCTIONS,
                                  f"{mission}_{label}")
        all_summary.append(summary)
        print(f"\n-- {mission} ({label}) --")
        print(summary.to_string(index=False))

    summary_df = pd.concat(all_summary, ignore_index=True)
    summary_df.to_csv(OUT_DIR / "lettuce_redundancy_summary.csv", index=False)
    print(f"\nSaved: {OUT_DIR}/lettuce_redundancy_summary.csv")

    # N-cycle completeness summary
    print("\n[5] N-cycle completeness group summary...")
    comp_summary_rows = []
    for subset_name, subset_mask in [
        ("aggregate_n4x4", pd.Series(True, index=meta.index)),
        ("OSD-267_VEG-01A", meta["OSD"] == "OSD-267"),
        ("OSD-269_VEG-03A", meta["OSD"] == "OSD-269"),
    ]:
        ids_in = [s for s in meta.index[subset_mask] if s in comp["sample_id"].values]
        sub_comp = comp[comp["sample_id"].isin(ids_in)]
        sp = sub_comp.loc[sub_comp["Study_Group"] == "Space_Flight", "completeness"].values
        tr = sub_comp.loc[sub_comp["Study_Group"] == "Terrestrial_Soil", "completeness"].values
        if len(sp) > 0 and len(tr) > 0:
            try:
                _, p = mannwhitneyu(sp, tr, alternative="two-sided")
            except ValueError:
                p = np.nan
        else:
            p = np.nan
        comp_summary_rows.append({
            "Subset":        subset_name,
            "Space_n":       len(sp),
            "Terr_n":        len(tr),
            "Space_mean_%":  round(float(np.mean(sp)), 1) if len(sp) else "NA",
            "Terr_mean_%":   round(float(np.mean(tr)), 1) if len(tr) else "NA",
            "MWU_p":         round(p, 4) if not np.isnan(p) else "NA",
            "Cliffs_delta":  round(cliffs_delta(sp, tr), 3) if len(sp) and len(tr) else "NA",
            "Hedges_g":      round(hedges_g(sp, tr), 3) if len(sp) > 1 and len(tr) > 1 else "NA",
        })
    comp_summary = pd.DataFrame(comp_summary_rows)
    comp_summary.to_csv(OUT_DIR / "lettuce_ncycle_completeness_summary.csv", index=False)
    print(comp_summary.to_string(index=False))
    print(f"\nSaved: {OUT_DIR}/lettuce_ncycle_completeness_summary.csv")

    print("\nDone.")


if __name__ == "__main__":
    main()
