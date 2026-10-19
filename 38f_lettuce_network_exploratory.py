"""
38f_lettuce_network_exploratory.py
===================================
EXPLORATORY co-occurrence network for lettuce root-zone — n=4 per condition
is below the threshold for reliable BH-FDR Spearman inference, so this
script reports:

  (a) Raw Spearman |r| > 0.5 edge set with NO FDR correction
  (b) Jaccard presence/absence co-occurrence (more robust at low n)
  (c) Top-8 shared genera metrics (degree, weighted degree)

Flagged clearly as "exploratory, not inferential". Pepper's
`09_network_analysis.py` uses BH-FDR q<0.05 which at n=4 would return
near-zero edges.

Both aggregate and within-group (Space-only / Terrestrial-only) networks
are computed for comparison.

Output:
  - acta_astronautica/lettuce/lettuce_network_edges.csv
  - acta_astronautica/lettuce/lettuce_network_jaccard.csv
  - acta_astronautica/lettuce/lettuce_network_metrics.csv
"""

from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config

LETTUCE_DIR = config.VEGGIE_LETTUCE_DIR
EXPORTS_DIR = LETTUCE_DIR / "exports"
OUT_DIR = config.AA_DIR / "lettuce"
OUT_DIR.mkdir(parents=True, exist_ok=True)

TOP_N_GENERA = 8          # reduced from pepper's 15 due to tiny n
R_THRESHOLD = 0.5         # raw |r| cutoff
PREVALENCE_MIN = 0.25     # a genus must appear in ≥25% of samples to be considered


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
    """Select top-N genera by mean RA, filtered by prevalence."""
    prev = (genus_ra > 0).sum(axis=1) / genus_ra.shape[1]
    filtered = genus_ra.loc[prev >= prev_min]
    mean_ra = filtered.mean(axis=1).sort_values(ascending=False)
    return filtered.loc[mean_ra.head(n).index]


def compute_spearman_edges(genus_mat, r_threshold=R_THRESHOLD):
    """Return edge list DataFrame with |r| above threshold (no FDR)."""
    genera = genus_mat.index.tolist()
    edges = []
    for g1, g2 in combinations(genera, 2):
        x = genus_mat.loc[g1].values
        y = genus_mat.loc[g2].values
        if np.std(x) == 0 or np.std(y) == 0:
            continue
        try:
            r, p = spearmanr(x, y)
        except Exception:
            continue
        if np.isnan(r) or abs(r) < r_threshold:
            continue
        edges.append({
            "genus_1": g1, "genus_2": g2,
            "spearman_r": round(r, 3),
            "p_raw": round(p, 4),
            "sign": "+" if r > 0 else "-",
        })
    return pd.DataFrame(edges)


def compute_jaccard_edges(genus_mat, prev_threshold=0):
    """Binary co-occurrence Jaccard index."""
    binary = (genus_mat > prev_threshold).astype(int)
    genera = binary.index.tolist()
    edges = []
    for g1, g2 in combinations(genera, 2):
        a = binary.loc[g1].values; b = binary.loc[g2].values
        inter = int(np.sum((a == 1) & (b == 1)))
        union = int(np.sum((a == 1) | (b == 1)))
        if union == 0:
            continue
        jacc = inter / union
        edges.append({
            "genus_1": g1, "genus_2": g2,
            "jaccard": round(jacc, 3),
            "cooccurring_samples": inter,
            "union_samples": union,
        })
    return pd.DataFrame(edges)


def network_metrics(edges_df, genera):
    """Degree (unweighted) and weighted degree for each genus."""
    deg = {g: 0 for g in genera}
    wdeg = {g: 0.0 for g in genera}
    for _, row in edges_df.iterrows():
        g1, g2 = row["genus_1"], row["genus_2"]
        w = row.get("spearman_r", row.get("jaccard", 0))
        if g1 in deg:
            deg[g1] += 1; wdeg[g1] += abs(w)
        if g2 in deg:
            deg[g2] += 1; wdeg[g2] += abs(w)
    return pd.DataFrame({
        "genus":           genera,
        "degree":          [deg[g] for g in genera],
        "weighted_degree": [round(wdeg[g], 3) for g in genera],
    }).sort_values("weighted_degree", ascending=False)


def main():
    print("=" * 60)
    print("38f: EXPLORATORY network for lettuce (n=4 per group, no FDR)")
    print("=" * 60)

    ft = pd.read_csv(EXPORTS_DIR / "feature-table.tsv", sep="\t", skiprows=1, index_col=0)
    tax = pd.read_csv(EXPORTS_DIR / "taxonomy.tsv", sep="\t", index_col=0)
    meta = pd.read_csv(LETTUCE_DIR / "metadata_harmonized.tsv", sep="\t", index_col=0)

    all_edges_spearman = []
    all_edges_jaccard = []
    all_metrics = []

    for subset_name, subset_ids in [
        ("all_8_samples",    list(meta.index)),
        ("space_flight",     list(meta.index[meta["Study_Group"] == "Space_Flight"])),
        ("terrestrial_soil", list(meta.index[meta["Study_Group"] == "Terrestrial_Soil"])),
    ]:
        print(f"\n[{subset_name}] n = {len(subset_ids)}")
        genus_ra = build_genus_ra(ft, tax, subset_ids)
        top = top_n_genera(genus_ra)
        print(f"  Top-{top.shape[0]} genera by mean RA: {list(top.index)}")

        # Spearman edges
        sp_edges = compute_spearman_edges(top)
        sp_edges["Subset"] = subset_name
        all_edges_spearman.append(sp_edges)
        print(f"  Spearman |r|>{R_THRESHOLD} edges: {len(sp_edges)}")

        # Jaccard edges
        jac_edges = compute_jaccard_edges(top)
        jac_edges["Subset"] = subset_name
        all_edges_jaccard.append(jac_edges)
        print(f"  Jaccard edges: {len(jac_edges)}")

        # Metrics
        metrics = network_metrics(sp_edges, list(top.index))
        metrics["Subset"] = subset_name
        all_metrics.append(metrics)
        print(metrics.to_string(index=False))

    spearman_df = pd.concat(all_edges_spearman, ignore_index=True)
    jaccard_df = pd.concat(all_edges_jaccard, ignore_index=True)
    metrics_df = pd.concat(all_metrics, ignore_index=True)

    spearman_df.to_csv(OUT_DIR / "lettuce_network_edges_spearman.csv", index=False)
    jaccard_df.to_csv(OUT_DIR / "lettuce_network_edges_jaccard.csv", index=False)
    metrics_df.to_csv(OUT_DIR / "lettuce_network_metrics.csv", index=False)

    print(f"\nSaved:")
    print(f"  {OUT_DIR}/lettuce_network_edges_spearman.csv  ({len(spearman_df)} rows)")
    print(f"  {OUT_DIR}/lettuce_network_edges_jaccard.csv   ({len(jaccard_df)} rows)")
    print(f"  {OUT_DIR}/lettuce_network_metrics.csv         ({len(metrics_df)} rows)")

    print("\nNote: flagged as EXPLORATORY — n=4 per group is below reliable")
    print("BH-FDR threshold. Report as descriptive, not inferential.")
    print("\nDone.")


if __name__ == "__main__":
    main()
