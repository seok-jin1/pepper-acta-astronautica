"""
38c_lettuce_nitrifier_direct.py
===============================
Direct ASV-level nitrifier detection for lettuce — mirrors the pepper
signature finding (pepper: 0/20 Space vs 19/20 Terrestrial, Cliff's δ = -0.95).

Uses the DADA2-generated lettuce feature table (from completed 37_veggie_lettuce_reanalysis.sh),
NOT the GeneLab-processed tables, for methodological symmetry with pepper.

For each of the four canonical nitrifier genera (Nitrosomonas, Nitrosospira,
Nitrospira, Nitrobacter), compute:
  - Prevalence (X/4 flight vs Y/4 ground) — composite (any genus)
  - Mean relative abundance
  - Fisher's exact test p-value
  - Cliff's δ effect size

Reports both (a) aggregated n=4+4 and (b) mission-stratified (OSD-267 n=2+2, OSD-269 n=2+2).

Output:
  - acta_astronautica/lettuce/lettuce_nitrifier_direct_aggregate.csv
  - acta_astronautica/lettuce/lettuce_nitrifier_direct_by_mission.csv
  - acta_astronautica/lettuce/lettuce_nitrifier_persample.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config

NITRIFIER_GENERA = ["Nitrosomonas", "Nitrosospira", "Nitrospira", "Nitrobacter"]
LETTUCE_DIR = config.VEGGIE_LETTUCE_DIR
EXPORTS_DIR = LETTUCE_DIR / "exports"
OUT_DIR = config.AA_DIR / "lettuce"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def parse_genus(taxon_str):
    """Extract genus name from a semicolon-delimited SILVA lineage string."""
    if pd.isna(taxon_str):
        return "Unassigned"
    for part in reversed(str(taxon_str).split(";")):
        p = part.strip()
        if p.startswith("g__"):
            name = p.split("__", 1)[-1].strip()
            if name and name not in ("uncultured", "metagenome", "unidentified"):
                return name
    return "Unassigned"


def cliffs_delta(a, b):
    a = np.asarray(a); b = np.asarray(b)
    if len(a) == 0 or len(b) == 0:
        return np.nan
    dom = np.sum(a[:, None] > b[None, :]) - np.sum(a[:, None] < b[None, :])
    return dom / (len(a) * len(b))


def load_lettuce_table():
    """Load DADA2-exported feature table + taxonomy + metadata for lettuce."""
    # Feature table (after biom convert --to-tsv → has extra first line #Constructed from biom)
    ft_path = EXPORTS_DIR / "feature-table.tsv"
    if not ft_path.exists():
        raise FileNotFoundError(
            f"{ft_path} not found. Run 37_veggie_lettuce_reanalysis.sh completely first.")
    ft = pd.read_csv(ft_path, sep="\t", skiprows=1, index_col=0)

    tax_path = EXPORTS_DIR / "taxonomy.tsv"
    tax = pd.read_csv(tax_path, sep="\t", index_col=0)

    meta = pd.read_csv(LETTUCE_DIR / "metadata_harmonized.tsv",
                       sep="\t", index_col=0)
    return ft, tax, meta


def compute_nitrifier_signal(ft, tax, meta, subset_name="aggregate",
                             sample_filter=None):
    """
    Compute per-sample nitrifier RA + per-group prevalence + Fisher's exact
    for a given sample subset.
    """
    # Filter metadata if specified
    if sample_filter is not None:
        meta_sub = meta[sample_filter]
    else:
        meta_sub = meta
    sample_ids = [s for s in meta_sub.index if s in ft.columns]
    ft_sub = ft[sample_ids]

    # Build genus relative-abundance table
    shared = ft_sub.index.intersection(tax.index)
    ft_sub = ft_sub.loc[shared]
    genus = tax.loc[shared, "Taxon"].apply(parse_genus)
    ft_sub = ft_sub.copy()
    ft_sub["Genus"] = genus.values
    genus_ra = ft_sub.groupby("Genus").sum()
    genus_ra = genus_ra.drop(index="Unassigned", errors="ignore")
    # Normalize to relative abundance (%)
    col_sums = genus_ra.sum(axis=0).replace(0, 1)
    genus_ra = genus_ra.div(col_sums, axis=1) * 100

    # Per-sample nitrifier composite (any of 4 genera)
    present_genera = [g for g in NITRIFIER_GENERA if g in genus_ra.index]
    comp = (genus_ra.loc[present_genera].sum(axis=0)
            if present_genera else pd.Series(0.0, index=genus_ra.columns))

    sp_ids = [s for s in meta_sub.index[meta_sub["Study_Group"] == "Space_Flight"]
              if s in comp.index]
    tr_ids = [s for s in meta_sub.index[meta_sub["Study_Group"] == "Terrestrial_Soil"]
              if s in comp.index]

    sp_vals = comp.loc[sp_ids].values
    tr_vals = comp.loc[tr_ids].values

    # Per-genus + composite
    rows = []
    for genus_name in NITRIFIER_GENERA + ["COMPOSITE"]:
        if genus_name == "COMPOSITE":
            sp_v = sp_vals
            tr_v = tr_vals
        elif genus_name in genus_ra.index:
            sp_v = genus_ra.loc[genus_name, sp_ids].values
            tr_v = genus_ra.loc[genus_name, tr_ids].values
        else:
            sp_v = np.zeros(len(sp_ids))
            tr_v = np.zeros(len(tr_ids))

        sp_prev = int((sp_v > 0).sum())
        tr_prev = int((tr_v > 0).sum())
        n_sp = len(sp_ids)
        n_tr = len(tr_ids)

        # Fisher's exact on prevalence
        if n_sp > 0 and n_tr > 0:
            try:
                _, fisher_p = fisher_exact(
                    [[sp_prev, n_sp - sp_prev], [tr_prev, n_tr - tr_prev]])
            except Exception:
                fisher_p = np.nan
        else:
            fisher_p = np.nan

        log2fc = np.log2((sp_v.mean() + 1e-6) / (tr_v.mean() + 1e-6))

        rows.append({
            "Subset":               subset_name,
            "Genus":                genus_name,
            "Space_n":              n_sp,
            "Terr_n":               n_tr,
            "Space_prevalence":     f"{sp_prev}/{n_sp}",
            "Terr_prevalence":      f"{tr_prev}/{n_tr}",
            "Space_prevalence_pct": round(100 * sp_prev / n_sp, 1) if n_sp else np.nan,
            "Terr_prevalence_pct":  round(100 * tr_prev / n_tr, 1) if n_tr else np.nan,
            "Space_mean_RA":        round(float(sp_v.mean()), 4),
            "Terr_mean_RA":         round(float(tr_v.mean()), 4),
            "Log2FC":               round(float(log2fc), 2),
            "Cliffs_delta":         round(cliffs_delta(sp_v, tr_v), 3),
            "Fisher_p":             round(fisher_p, 4) if not np.isnan(fisher_p) else "NA",
        })

    return pd.DataFrame(rows), comp, sp_ids, tr_ids


def main():
    print("=" * 60)
    print("38c: Direct ASV-level nitrifier detection for lettuce")
    print("=" * 60)

    ft, tax, meta = load_lettuce_table()
    print(f"\n  Feature table: {ft.shape[0]} ASVs × {ft.shape[1]} samples")
    print(f"  Taxonomy:      {tax.shape[0]} ASVs")
    print(f"  Metadata:      {len(meta)} samples\n")

    all_results = []
    all_persample = []

    # (a) Aggregate: n=4+4 (both missions combined)
    print("[1] Aggregate (all missions)...")
    agg_df, agg_comp, sp_ids, tr_ids = compute_nitrifier_signal(
        ft, tax, meta, subset_name="aggregate_n4x4")
    print(agg_df.to_string(index=False))
    all_results.append(agg_df)

    # Per-sample table
    persample = pd.DataFrame({
        "sample_id": list(agg_comp.index),
        "composite_nitrifier_RA": agg_comp.values,
    })
    persample["Study_Group"] = meta.loc[persample["sample_id"], "Study_Group"].values
    persample["Mission"] = meta.loc[persample["sample_id"], "Mission"].values
    persample["OSD"] = meta.loc[persample["sample_id"], "OSD"].values
    all_persample.append(persample)

    # (b) Mission-stratified
    for mission, label in [("OSD-267", "VEG-01A"), ("OSD-269", "VEG-03A")]:
        print(f"\n[2] Mission {mission} ({label}) n=2+2...")
        sel = (meta["OSD"] == mission)
        mis_df, _, _, _ = compute_nitrifier_signal(
            ft, tax, meta, subset_name=f"{mission}_{label}",
            sample_filter=sel)
        print(mis_df.to_string(index=False))
        all_results.append(mis_df)

    # Concatenate + save
    summary = pd.concat(all_results, ignore_index=True)
    summary.to_csv(OUT_DIR / "lettuce_nitrifier_direct.csv", index=False)
    print(f"\nSaved: {OUT_DIR}/lettuce_nitrifier_direct.csv")

    persample_df = all_persample[0]
    persample_df.to_csv(OUT_DIR / "lettuce_nitrifier_persample.csv", index=False)
    print(f"Saved: {OUT_DIR}/lettuce_nitrifier_persample.csv")

    # Pepper comparison row (read from existing pepper output if available)
    pepper_file = config.VERSION_DIR / "SuppTable_S1c_rarefaction_sensitivity.csv"
    if pepper_file.exists():
        pep = pd.read_csv(pepper_file)
        print("\n[3] Pepper comparison (from SuppTable_S1c at depth=1000):")
        pep_1k = pep[pep.Rarefaction_depth == 1000].iloc[0]
        print(f"  Pepper: Space prev = {pep_1k.get('Space_prevalence_%', 'NA')}%, "
              f"Terr prev = {pep_1k.get('Terr_prevalence_%', 'NA')}%, "
              f"Log2FC = {pep_1k.get('Nitrifier_composite_Log2FC', 'NA')}")

    print("\nDone.")


if __name__ == "__main__":
    main()
