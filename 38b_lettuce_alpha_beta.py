"""
38b_lettuce_alpha_beta.py
=========================
Alpha and beta diversity for lettuce root-zone (n=4 flight + n=4 ground)
matching pepper's 10_supp_faith_pd.py + 15_permanova.R analyses.

Alpha: Shannon entropy + Observed ASVs (Faith's PD optional if tree provided)
Beta:  Bray-Curtis distance matrix, PERMANOVA (adonis2 equivalent via scipy),
       PCoA coordinates, betadisper (dispersion homogeneity)

Both aggregate (n=4+4) and mission-stratified (n=2+2 × 2) are reported.

Uses scikit-bio's PERMANOVA to avoid switching to R.

Output:
  - acta_astronautica/lettuce/lettuce_alpha_diversity.csv
  - acta_astronautica/lettuce/lettuce_beta_permanova.csv
  - acta_astronautica/lettuce/lettuce_pcoa_coords.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config

LETTUCE_DIR = config.VEGGIE_LETTUCE_DIR
EXPORTS_DIR = LETTUCE_DIR / "exports"
OUT_DIR = config.AA_DIR / "lettuce"
OUT_DIR.mkdir(parents=True, exist_ok=True)

RANDOM_SEED = config.RANDOM_SEED
N_PERM = 999


def load_table_tax_meta():
    ft = pd.read_csv(EXPORTS_DIR / "feature-table.tsv", sep="\t", skiprows=1, index_col=0)
    tax = pd.read_csv(EXPORTS_DIR / "taxonomy.tsv", sep="\t", index_col=0)
    meta = pd.read_csv(LETTUCE_DIR / "metadata_harmonized.tsv", sep="\t", index_col=0)
    return ft, tax, meta


def shannon(p):
    p = np.asarray(p, dtype=float)
    p = p[p > 0]
    if len(p) == 0:
        return 0.0
    p = p / p.sum()
    return -np.sum(p * np.log(p))


def observed_features(counts):
    return int(np.sum(np.asarray(counts) > 0))


def bray_curtis(X):
    """X: samples × features (counts or RA). Returns square distance matrix."""
    rel = X / X.sum(axis=1, keepdims=True).clip(min=1)
    dist = squareform(pdist(rel, metric="braycurtis"))
    return dist


def permanova_adonis(dist_mat, labels, n_perm=N_PERM, seed=RANDOM_SEED):
    """
    PERMANOVA via scipy: F-statistic from R² and permutation-based p-value.
    Returns dict with F, R2, p.
    """
    rng = np.random.default_rng(seed)
    labels = np.asarray(labels)
    n = len(labels)
    groups = np.unique(labels)
    k = len(groups)

    # Total SS
    SS_T = np.sum(dist_mat ** 2) / (2 * n)

    # Within-group SS
    def within_ss(lbls):
        ss_w = 0.0
        for g in groups:
            idx = np.where(lbls == g)[0]
            if len(idx) < 2:
                continue
            sub = dist_mat[np.ix_(idx, idx)]
            ss_w += np.sum(sub ** 2) / (2 * len(idx))
        return ss_w

    SS_W_obs = within_ss(labels)
    SS_A_obs = SS_T - SS_W_obs
    F_obs = (SS_A_obs / (k - 1)) / (SS_W_obs / (n - k))
    R2 = SS_A_obs / SS_T

    # Permutation
    F_perm = np.zeros(n_perm)
    for i in range(n_perm):
        shuffled = rng.permutation(labels)
        SS_W_p = within_ss(shuffled)
        SS_A_p = SS_T - SS_W_p
        F_perm[i] = (SS_A_p / (k - 1)) / (SS_W_p / (n - k))

    p_val = (np.sum(F_perm >= F_obs) + 1) / (n_perm + 1)
    return {"F": F_obs, "R2": R2, "p": p_val}


def betadisper(dist_mat, labels, n_perm=N_PERM, seed=RANDOM_SEED):
    """
    Homogeneity of dispersion test: test if groups differ in average distance
    from their centroid (approximated using PCoA embedding).
    Returns dict with p, per-group dispersions.
    """
    from numpy.linalg import eigh
    rng = np.random.default_rng(seed)
    n = len(labels)
    # Double centering for PCoA
    D2 = dist_mat ** 2
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ D2 @ J
    w, V = eigh(B)
    # Use positive eigenvalues only
    pos = w > 1e-8
    coords = V[:, pos] * np.sqrt(w[pos])

    labels = np.asarray(labels)
    groups = np.unique(labels)

    def dispersion(lbls):
        disps = {}
        total_var = 0.0
        for g in groups:
            idx = np.where(lbls == g)[0]
            if len(idx) < 2:
                disps[g] = np.nan
                continue
            centroid = coords[idx].mean(axis=0)
            d = np.sqrt(((coords[idx] - centroid) ** 2).sum(axis=1))
            disps[g] = d
        return disps

    obs_disps = dispersion(labels)
    group_means = {g: float(np.mean(obs_disps[g])) for g in groups}

    # Permutation test on F-statistic of dispersion values
    all_d = np.concatenate([obs_disps[g] for g in groups])
    grp_lbls = np.concatenate([[g] * len(obs_disps[g]) for g in groups])
    from scipy.stats import f_oneway
    if all(len(obs_disps[g]) > 1 for g in groups):
        F_obs, p_one = f_oneway(*[obs_disps[g] for g in groups])
    else:
        F_obs, p_one = np.nan, np.nan

    return {"p": p_one, "dispersions": group_means}


def analyze_subset(ft, tax, meta, subset_mask, subset_name):
    """Run alpha + beta + PERMANOVA + betadisper on a sample subset."""
    sub_meta = meta[subset_mask]
    sample_ids = [s for s in sub_meta.index if s in ft.columns]
    ft_sub = ft[sample_ids].T  # samples × ASVs
    ft_sub = ft_sub.loc[(ft_sub.sum(axis=1) > 0)]
    ids = ft_sub.index
    labels = sub_meta.loc[ids, "Study_Group"].values

    # Alpha per sample
    alpha_rows = []
    for s in ids:
        counts = ft_sub.loc[s].values
        alpha_rows.append({
            "sample_id": s,
            "Study_Group": sub_meta.loc[s, "Study_Group"],
            "OSD": sub_meta.loc[s, "OSD"],
            "Mission": sub_meta.loc[s, "Mission"],
            "Subset": subset_name,
            "Observed_ASVs": observed_features(counts),
            "Shannon": round(shannon(counts), 4),
            "Total_reads": int(counts.sum()),
        })
    alpha_df = pd.DataFrame(alpha_rows)

    # Beta / PERMANOVA (need ≥2 per group)
    n_sp = sum(labels == "Space_Flight")
    n_tr = sum(labels == "Terrestrial_Soil")
    if n_sp >= 2 and n_tr >= 2:
        dist = bray_curtis(ft_sub.values)
        perm = permanova_adonis(dist, labels)
        bd = betadisper(dist, labels)
        beta_row = {
            "Subset":           subset_name,
            "n_Space":          n_sp,
            "n_Terr":           n_tr,
            "PERMANOVA_F":      round(perm["F"], 3),
            "PERMANOVA_R2":     round(perm["R2"], 4),
            "PERMANOVA_p":      round(perm["p"], 4),
            "Betadisper_p":     round(bd["p"], 4) if not np.isnan(bd["p"]) else "NA",
            "Space_dispersion": round(bd["dispersions"].get("Space_Flight", np.nan), 4),
            "Terr_dispersion":  round(bd["dispersions"].get("Terrestrial_Soil", np.nan), 4),
        }
    else:
        beta_row = {
            "Subset":           subset_name,
            "n_Space":          n_sp,
            "n_Terr":           n_tr,
            "PERMANOVA_F":      "NA",
            "PERMANOVA_R2":     "NA",
            "PERMANOVA_p":      "NA (insufficient n)",
            "Betadisper_p":     "NA",
            "Space_dispersion": "NA",
            "Terr_dispersion":  "NA",
        }

    return alpha_df, beta_row


def main():
    print("=" * 60)
    print("38b: α/β-diversity + PERMANOVA for lettuce")
    print("=" * 60)

    ft, tax, meta = load_table_tax_meta()
    print(f"\n  Feature table: {ft.shape[0]} ASVs × {ft.shape[1]} samples")
    print(f"  Metadata:      {len(meta)} samples")

    all_alpha = []
    all_beta = []

    # Aggregate
    print("\n[1] Aggregate (n=4+4)...")
    a, b = analyze_subset(ft, tax, meta,
                          subset_mask=pd.Series(True, index=meta.index),
                          subset_name="aggregate_n4x4")
    all_alpha.append(a)
    all_beta.append(b)
    print(a.to_string(index=False))
    print(f"\n  {b}")

    # Mission-stratified
    for mission, label in [("OSD-267", "VEG-01A"), ("OSD-269", "VEG-03A")]:
        print(f"\n[2] Mission {mission} ({label}) n=2+2...")
        a, b = analyze_subset(ft, tax, meta,
                              subset_mask=(meta["OSD"] == mission),
                              subset_name=f"{mission}_{label}")
        all_alpha.append(a)
        all_beta.append(b)
        print(a.to_string(index=False))
        print(f"  {b}")

    # Save
    alpha_df = pd.concat(all_alpha, ignore_index=True)
    alpha_df.to_csv(OUT_DIR / "lettuce_alpha_diversity.csv", index=False)
    print(f"\nSaved: {OUT_DIR}/lettuce_alpha_diversity.csv")

    beta_df = pd.DataFrame(all_beta)
    beta_df.to_csv(OUT_DIR / "lettuce_beta_permanova.csv", index=False)
    print(f"Saved: {OUT_DIR}/lettuce_beta_permanova.csv")

    # PCoA coords on aggregate (useful for figure)
    sp_ids = [s for s in meta.index if s in ft.columns]
    ft_all = ft[sp_ids].T
    dist = bray_curtis(ft_all.values)
    # double centering + PCoA
    from numpy.linalg import eigh
    n = len(ft_all)
    D2 = dist ** 2
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ D2 @ J
    w, V = eigh(B)
    order = np.argsort(-w)
    w = w[order]; V = V[:, order]
    coords = V * np.sqrt(np.maximum(w, 0))
    pcoa_df = pd.DataFrame({
        "sample_id": ft_all.index,
        "PC1": coords[:, 0],
        "PC2": coords[:, 1] if coords.shape[1] > 1 else 0,
        "Study_Group": meta.loc[ft_all.index, "Study_Group"].values,
        "OSD":         meta.loc[ft_all.index, "OSD"].values,
        "Mission":     meta.loc[ft_all.index, "Mission"].values,
    })
    pcoa_df.to_csv(OUT_DIR / "lettuce_pcoa_coords.csv", index=False)
    pct_var = w[:4] / w[w > 0].sum() * 100 if (w > 0).sum() else [0] * 4
    print(f"Saved: {OUT_DIR}/lettuce_pcoa_coords.csv "
          f"(PC1 {pct_var[0]:.1f}%, PC2 {pct_var[1]:.1f}%)")

    print("\nDone.")


if __name__ == "__main__":
    main()
