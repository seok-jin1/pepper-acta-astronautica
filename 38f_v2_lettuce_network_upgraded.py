"""
38f_v2_lettuce_network_upgraded.py
===================================
Upgrade of 38f_lettuce_network_exploratory.py with parameters aligned
to pepper's script 09 + robustness extensions for low-n:

  - TOP_N_GENERA = 15  (same as pepper 09)
  - PREVALENCE_MIN = 0.125  (genus in ≥1/8 samples)
  - R_THRESHOLD = 0.5  (unchanged; tighter helps at n=4)
  - Parallel Spearman, Kendall τ, CLR-Pearson correlations
  - Sign-split degree (positive / negative edges separately)
  - Jaccard binary co-occurrence (already in 38f, re-computed here)
  - Genus × genus correlation heatmaps (Space / Terr / combined)

Genus-level throughout (NOT ASV-level). Lettuce GeneLab-processed data is
NOT used — uses DADA2-reprocessed table at
`version-2_integrated/external/veggie_lettuce/exports/feature-table.tsv`.

Still EXPLORATORY at n=4 per group. Nothing in this script is inferential.

Outputs (acta_astronautica/lettuce/):
  - lettuce_network_top15_spearman.csv
  - lettuce_network_top15_kendall.csv
  - lettuce_network_top15_clr_pearson.csv
  - lettuce_network_sign_split.csv
  - lettuce_network_heatmap_{aggregate,space,terrestrial}.{pdf,png}
"""

from __future__ import annotations

import sys
import warnings
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import kendalltau, pearsonr, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config

LETTUCE_DIR = config.VEGGIE_LETTUCE_DIR
EXPORTS_DIR = LETTUCE_DIR / "exports"
OUT_DIR = config.AA_DIR / "lettuce"
OUT_DIR.mkdir(parents=True, exist_ok=True)

TOP_N_GENERA = 15
R_THRESHOLD = 0.5
PREVALENCE_MIN = 0.125  # ≥1/8 samples

NITRIFIER_GENERA = {"Nitrosomonas", "Nitrosospira", "Nitrospira", "Nitrobacter"}


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


def build_genus_ra(ft, tax, sample_ids):
    ft = ft[[s for s in sample_ids if s in ft.columns]]
    shared = ft.index.intersection(tax.index)
    ft = ft.loc[shared].copy()
    ft["Genus"] = tax.loc[shared, "Taxon"].apply(parse_genus).values
    g = ft.groupby("Genus").sum()
    g = g.drop(index="Unassigned", errors="ignore")
    col_sums = g.sum(axis=0).replace(0, 1)
    return g.div(col_sums, axis=1) * 100


def top_n_genera(genus_ra, n=TOP_N_GENERA, prev_min=PREVALENCE_MIN):
    prev = (genus_ra > 0).sum(axis=1) / genus_ra.shape[1]
    filtered = genus_ra.loc[prev >= prev_min]
    mean_ra = filtered.mean(axis=1).sort_values(ascending=False)
    return filtered.loc[mean_ra.head(n).index]


def clr_transform(X, pseudocount=1e-6):
    """Centered log-ratio transformation on genera (rows) × samples (columns).
    For each sample (column), subtract geometric mean of log-RA.
    """
    X = X + pseudocount
    logX = np.log(X)
    gm = logX.mean(axis=0)  # per-sample geometric mean
    return logX.subtract(gm, axis=1)


def pairwise_corr(mat, method="spearman"):
    """
    Return DataFrame of pairwise correlations.
    method: 'spearman' | 'kendall' | 'pearson'
    """
    genera = mat.index.tolist()
    edges = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for g1, g2 in combinations(genera, 2):
            x = mat.loc[g1].values
            y = mat.loc[g2].values
            if np.std(x) == 0 or np.std(y) == 0:
                continue
            if method == "spearman":
                r, p = spearmanr(x, y)
            elif method == "kendall":
                r, p = kendalltau(x, y)
            else:
                r, p = pearsonr(x, y)
            if np.isnan(r) or abs(r) < R_THRESHOLD:
                continue
            edges.append({
                "genus_1": g1, "genus_2": g2,
                f"{method}_r": round(r, 3),
                f"{method}_p": round(p, 4),
                "sign": "+" if r > 0 else "-",
            })
    return pd.DataFrame(edges)


def sign_split_degree(edges_df, genera, r_col):
    """Return DataFrame with positive/negative degree per genus."""
    pos = {g: 0 for g in genera}
    neg = {g: 0 for g in genera}
    pos_w = {g: 0.0 for g in genera}
    neg_w = {g: 0.0 for g in genera}
    for _, row in edges_df.iterrows():
        g1, g2 = row["genus_1"], row["genus_2"]
        r = row[r_col]
        if g1 in pos and g2 in pos:
            if r > 0:
                pos[g1] += 1; pos[g2] += 1
                pos_w[g1] += r; pos_w[g2] += r
            else:
                neg[g1] += 1; neg[g2] += 1
                neg_w[g1] += abs(r); neg_w[g2] += abs(r)
    return pd.DataFrame({
        "genus":            genera,
        "pos_degree":       [pos[g] for g in genera],
        "neg_degree":       [neg[g] for g in genera],
        "total_degree":     [pos[g] + neg[g] for g in genera],
        "pos_weighted":     [round(pos_w[g], 3) for g in genera],
        "neg_weighted":     [round(neg_w[g], 3) for g in genera],
        "is_nitrifier":     [g in NITRIFIER_GENERA for g in genera],
    })


def heatmap_corr_matrix(top_ra, subset_name, method="spearman"):
    """Compute full correlation matrix (for heatmap) and save figure."""
    genera = top_ra.index.tolist()
    n = len(genera)
    M = np.zeros((n, n))
    for i, g1 in enumerate(genera):
        for j, g2 in enumerate(genera):
            if i == j:
                M[i, j] = 1.0
                continue
            x = top_ra.loc[g1].values
            y = top_ra.loc[g2].values
            if np.std(x) == 0 or np.std(y) == 0:
                M[i, j] = 0.0
                continue
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                r, _ = spearmanr(x, y) if method == "spearman" else kendalltau(x, y)
            M[i, j] = r if not np.isnan(r) else 0.0

    # Highlight nitrifier rows/cols if any
    labels = [f"* {g}" if g in NITRIFIER_GENERA else g for g in genera]
    fig, ax = plt.subplots(figsize=(9, 8))
    sns.heatmap(
        pd.DataFrame(M, index=labels, columns=labels),
        cmap="RdBu_r", vmin=-1, vmax=1, center=0,
        square=True, cbar_kws={"label": f"{method.capitalize()} r"},
        xticklabels=True, yticklabels=True, ax=ax,
    )
    ax.set_title(f"Lettuce top-{TOP_N_GENERA} genus-level {method} correlation "
                 f"— subset: {subset_name}\n(* = canonical nitrifier if present)",
                 fontsize=11)
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)
    plt.tight_layout()
    for ext in ("pdf", "png"):
        out = OUT_DIR / f"lettuce_network_heatmap_{subset_name}.{ext}"
        plt.savefig(out, dpi=300 if ext == "png" else None, bbox_inches="tight")
    plt.close()
    return M


def main():
    print("=" * 60)
    print("38f_v2: Upgraded lettuce network (top-15, Spearman+Kendall+CLR)")
    print("=" * 60)

    # Use RAW (un-rarefied) table to keep all 8 samples — critical for
    # getting a non-empty Terrestrial subset (at depth 1000, OSD-267 ground
    # samples are dropped leaving n=2; raw retains n=4+4).
    # Relative abundance normalization applied per sample, so depth bias is
    # minimized at genus level.
    ft = pd.read_csv(EXPORTS_DIR / "feature-table-raw.tsv",
                     sep="\t", skiprows=1, index_col=0)
    tax = pd.read_csv(EXPORTS_DIR / "taxonomy.tsv", sep="\t", index_col=0)
    meta = pd.read_csv(LETTUCE_DIR / "metadata_harmonized.tsv", sep="\t", index_col=0)
    print(f"  Feature table (raw): {ft.shape[0]} ASVs × {ft.shape[1]} samples")
    print(f"  Sample depths: min={ft.sum().min():,.0f}, "
          f"median={ft.sum().median():,.0f}, max={ft.sum().max():,.0f}\n")

    all_spearman = []
    all_kendall  = []
    all_clr      = []
    all_sign     = []

    for subset_name, subset_ids in [
        ("aggregate",    list(meta.index)),
        ("space",        list(meta.index[meta["Study_Group"] == "Space_Flight"])),
        ("terrestrial",  list(meta.index[meta["Study_Group"] == "Terrestrial_Soil"])),
    ]:
        # Only samples present in rarefied table
        subset_ids = [s for s in subset_ids if s in ft.columns]
        print(f"\n[{subset_name}] n = {len(subset_ids)}")
        if len(subset_ids) < 3:
            print(f"  SKIP (need ≥3 samples; got {len(subset_ids)})")
            continue

        genus_ra = build_genus_ra(ft, tax, subset_ids)
        top_ra = top_n_genera(genus_ra)
        print(f"  Top-{top_ra.shape[0]} genera selected. "
              f"Nitrifier in top: "
              f"{sorted(set(top_ra.index) & NITRIFIER_GENERA) or 'none'}")

        # CLR transform for compositional-aware pearson
        top_clr = clr_transform(top_ra)

        # Spearman
        sp_e = pairwise_corr(top_ra, method="spearman")
        sp_e["Subset"] = subset_name
        all_spearman.append(sp_e)
        print(f"  Spearman |r|>{R_THRESHOLD} edges: {len(sp_e)} "
              f"(+: {(sp_e['sign']=='+').sum() if len(sp_e) else 0}, "
              f"-: {(sp_e['sign']=='-').sum() if len(sp_e) else 0})")

        # Kendall
        kd_e = pairwise_corr(top_ra, method="kendall")
        kd_e["Subset"] = subset_name
        all_kendall.append(kd_e)
        print(f"  Kendall |τ|>{R_THRESHOLD} edges: {len(kd_e)}")

        # CLR-Pearson
        cp_e = pairwise_corr(top_clr, method="pearson")
        cp_e["Subset"] = subset_name
        all_clr.append(cp_e)
        print(f"  CLR-Pearson |r|>{R_THRESHOLD} edges: {len(cp_e)}")

        # Sign-split degree (Spearman-based)
        sign_df = sign_split_degree(sp_e, list(top_ra.index), "spearman_r")
        sign_df["Subset"] = subset_name
        all_sign.append(sign_df)
        print("  Sign-split degree (top 5 by total):")
        print("  " + sign_df.sort_values("total_degree", ascending=False)
              .head(5).to_string(index=False).replace("\n", "\n  "))

        # Heatmap
        heatmap_corr_matrix(top_ra, subset_name, method="spearman")
        print(f"  Saved heatmap: lettuce_network_heatmap_{subset_name}.pdf/.png")

    # Concatenate and save
    if all_spearman:
        pd.concat(all_spearman, ignore_index=True).to_csv(
            OUT_DIR / "lettuce_network_top15_spearman.csv", index=False)
    if all_kendall:
        pd.concat(all_kendall, ignore_index=True).to_csv(
            OUT_DIR / "lettuce_network_top15_kendall.csv", index=False)
    if all_clr:
        pd.concat(all_clr, ignore_index=True).to_csv(
            OUT_DIR / "lettuce_network_top15_clr_pearson.csv", index=False)
    if all_sign:
        pd.concat(all_sign, ignore_index=True).to_csv(
            OUT_DIR / "lettuce_network_sign_split.csv", index=False)

    print("\nSaved:")
    print(f"  {OUT_DIR}/lettuce_network_top15_spearman.csv")
    print(f"  {OUT_DIR}/lettuce_network_top15_kendall.csv")
    print(f"  {OUT_DIR}/lettuce_network_top15_clr_pearson.csv")
    print(f"  {OUT_DIR}/lettuce_network_sign_split.csv")
    print(f"  {OUT_DIR}/lettuce_network_heatmap_*.pdf/.png (3 subsets)")
    print("\nEXPLORATORY — n=4 per group remains underpowered. "
          "Use this for visualization + hypothesis generation only.")


if __name__ == "__main__":
    main()
