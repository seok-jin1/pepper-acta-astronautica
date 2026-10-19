"""
43_nonindependence_check.py
===========================
Non-independence diagnostic for the analytical root-zone PERMANOVA.

Design: Space_Flight n=20 root-zone samples are clustered into 4 quadrants
(Q1-Q4; 5 samples per quadrant: 4 Plant Roots + 1 Rhizosphere).  Samples in
the same quadrant likely share growth-pillow-level environment, violating
the independence assumption of standard PERMANOVA.

We report three diagnostics:

  (A) Within- vs between-group Bray-Curtis distance summary, including
      within-quadrant similarity inside Space_Flight.

  (B) Quadrant-level aggregation PERMANOVA: Space samples collapsed to their
      per-quadrant centroid (n=4), Terrestrial subsampled to n=4 once per
      iteration.  Tests whether the group difference survives strict
      independence enforcement.

  (C) One-sample-per-Quadrant bootstrap PERMANOVA: 1 Space sample per Q
      (n=4 Space) + 4 randomly drawn Terr samples; 1000 iterations;
      reports F/R^2 distribution.

Output: version-2_integrated/SuppTable_S1b_nonindependence.csv
        stdout summary lines ready to be pasted into Methods/Results.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config

RAREFY_SEED = 42
N_ITER = 1000

ROOT_ZONE_TYPES = {"Plant Roots", "Plant_Roots", "Rhizosphere"}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_analytical_root_zone():
    ft = pd.read_csv(config.FEATURE_TABLE_CLEAN, sep="\t", skiprows=1, index_col=0)
    meta = pd.read_csv(config.INTEGRATED_METADATA, sep="\t", index_col=0)
    sp_ids = meta.index[
        (meta["Study_Group"] == "Space_Flight")
        & meta["sample_type"].isin(ROOT_ZONE_TYPES)
    ]
    tr_ids = meta.index[meta["Study_Group"] == "Terrestrial_Soil"]
    keep = [s for s in list(sp_ids) + list(tr_ids) if s in ft.columns]
    rel = ft[keep].div(ft[keep].sum(axis=0), axis=1).fillna(0.0).T
    return rel, meta.loc[keep]


def bc_matrix(rel):
    return pd.DataFrame(
        squareform(pdist(rel.values, metric="braycurtis")),
        index=rel.index, columns=rel.index,
    )


def permanova(rel, groups, n_perm=999, seed=42):
    """Two-group PERMANOVA on Bray-Curtis distances. Returns dict with F, R2, p."""
    d = squareform(pdist(rel.values, metric="braycurtis"))
    labels = np.asarray(groups)
    n = len(labels)
    uniq = np.unique(labels)
    if len(uniq) != 2:
        return dict(F=np.nan, R2=np.nan, p=np.nan)

    def _stat(lab):
        ss_t = (d ** 2).sum() / (2 * n)
        ss_w = 0.0
        for g in uniq:
            idx = np.where(lab == g)[0]
            if len(idx) < 2:
                continue
            sub = d[np.ix_(idx, idx)]
            ss_w += (sub ** 2).sum() / (2 * len(idx))
        ss_a = ss_t - ss_w
        df_a = len(uniq) - 1
        df_w = n - len(uniq)
        if ss_w <= 0 or df_w <= 0:
            return np.nan, np.nan
        F = (ss_a / df_a) / (ss_w / df_w)
        R2 = ss_a / ss_t
        return F, R2

    F_obs, R2_obs = _stat(labels)
    rng = np.random.default_rng(seed)
    F_perm = np.empty(n_perm)
    for i in range(n_perm):
        F_perm[i], _ = _stat(labels[rng.permutation(n)])
    p = (np.sum(F_perm >= F_obs) + 1) / (n_perm + 1)
    return dict(F=F_obs, R2=R2_obs, p=p)


# ---------------------------------------------------------------------------
# (A) Within vs between BC distances
# ---------------------------------------------------------------------------

def within_between_report(rel, meta):
    bc = bc_matrix(rel)
    sp_mask = meta["Study_Group"] == "Space_Flight"
    tr_mask = meta["Study_Group"] == "Terrestrial_Soil"
    sp_ids = meta.index[sp_mask]
    tr_ids = meta.index[tr_mask]

    def mean_off_diag(ids):
        sub = bc.loc[ids, ids].values
        iu = np.triu_indices_from(sub, k=1)
        return float(sub[iu].mean()) if len(iu[0]) else np.nan

    within_sp = mean_off_diag(sp_ids)
    within_tr = mean_off_diag(tr_ids)
    between = float(bc.loc[sp_ids, tr_ids].values.mean())

    # Within-quadrant distances inside Space_Flight
    meta_sp = meta.loc[sp_ids]
    q_rows = []
    for q in sorted(meta_sp["location"].dropna().unique()):
        ids_q = meta_sp.index[meta_sp["location"] == q]
        if len(ids_q) < 2:
            continue
        q_rows.append({"quadrant": q, "n": len(ids_q),
                       "mean_within_BC": mean_off_diag(ids_q)})
    q_df = pd.DataFrame(q_rows)

    print("\n(A) Bray-Curtis distance summary")
    print(f"  Within Space (n=20):       {within_sp:.3f}")
    print(f"  Within Terrestrial (n=20): {within_tr:.3f}")
    print(f"  Between groups:            {between:.3f}")
    print("  Within-quadrant (Space only):")
    print(q_df.to_string(index=False))
    return dict(within_sp=within_sp, within_tr=within_tr, between=between,
                within_quadrants=q_df)


# ---------------------------------------------------------------------------
# (B) Quadrant-centroid PERMANOVA
# ---------------------------------------------------------------------------

def quadrant_centroid_permanova(rel, meta, n_iter=N_ITER, seed=RAREFY_SEED):
    sp_mask = meta["Study_Group"] == "Space_Flight"
    tr_mask = meta["Study_Group"] == "Terrestrial_Soil"
    meta_sp = meta.loc[sp_mask]
    sp_centroids = {}
    for q in sorted(meta_sp["location"].dropna().unique()):
        ids_q = meta_sp.index[meta_sp["location"] == q]
        if len(ids_q) == 0:
            continue
        sp_centroids[q] = rel.loc[ids_q].mean(axis=0)
    sp_mat = pd.DataFrame(sp_centroids).T  # Q × features
    tr_ids_all = meta.index[tr_mask].tolist()
    n_q = sp_mat.shape[0]

    rng = np.random.default_rng(seed)
    records = []
    for i in range(n_iter):
        tr_sub = rng.choice(tr_ids_all, size=n_q, replace=False)
        combined = pd.concat([sp_mat, rel.loc[list(tr_sub)]])
        labels = np.array(["Space"] * n_q + ["Terr"] * n_q)
        res = permanova(combined, labels, n_perm=199, seed=int(rng.integers(1e9)))
        records.append({"iter": i + 1, **res})
    df = pd.DataFrame(records).dropna()
    print("\n(B) Quadrant-centroid PERMANOVA "
          f"(Space n={n_q} quadrant centroids vs Terr n={n_q} subsampled, "
          f"{n_iter} iterations)")
    print(f"  F:  median={df['F'].median():.2f}  "
          f"95% CI=[{df['F'].quantile(0.025):.2f}, {df['F'].quantile(0.975):.2f}]")
    print(f"  R²: median={df['R2'].median():.3f}  "
          f"95% CI=[{df['R2'].quantile(0.025):.3f}, {df['R2'].quantile(0.975):.3f}]")
    print(f"  p<0.05 fraction:          {(df['p'] < 0.05).mean():.1%}")
    return df


# ---------------------------------------------------------------------------
# (C) One-sample-per-Q bootstrap
# ---------------------------------------------------------------------------

def one_per_q_bootstrap(rel, meta, n_iter=N_ITER, seed=RAREFY_SEED):
    meta_sp = meta.loc[meta["Study_Group"] == "Space_Flight"]
    q_to_ids = {q: meta_sp.index[meta_sp["location"] == q].tolist()
                for q in sorted(meta_sp["location"].dropna().unique())}
    tr_ids_all = meta.index[meta["Study_Group"] == "Terrestrial_Soil"].tolist()
    n_q = len(q_to_ids)

    rng = np.random.default_rng(seed)
    records = []
    for i in range(n_iter):
        sp_sub = [rng.choice(ids) for _, ids in q_to_ids.items()]
        tr_sub = rng.choice(tr_ids_all, size=n_q, replace=False).tolist()
        combined = rel.loc[sp_sub + tr_sub]
        labels = np.array(["Space"] * n_q + ["Terr"] * n_q)
        res = permanova(combined, labels, n_perm=199, seed=int(rng.integers(1e9)))
        records.append({"iter": i + 1, **res})
    df = pd.DataFrame(records).dropna()
    print("\n(C) One-sample-per-quadrant bootstrap PERMANOVA "
          f"(n={n_q} Space + {n_q} Terr, {n_iter} iterations)")
    print(f"  F:  median={df['F'].median():.2f}  "
          f"95% CI=[{df['F'].quantile(0.025):.2f}, {df['F'].quantile(0.975):.2f}]")
    print(f"  R²: median={df['R2'].median():.3f}  "
          f"95% CI=[{df['R2'].quantile(0.025):.3f}, {df['R2'].quantile(0.975):.3f}]")
    print(f"  p<0.05 fraction:          {(df['p'] < 0.05).mean():.1%}")
    return df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=== 43_nonindependence_check.py ===")
    rel, meta = load_analytical_root_zone()
    print(f"Analytical root zone: {rel.shape[0]} samples × {rel.shape[1]} ASVs")

    wb = within_between_report(rel, meta)
    df_b = quadrant_centroid_permanova(rel, meta)
    df_c = one_per_q_bootstrap(rel, meta)

    rows = [
        {"Diagnostic": "Within-Space mean BC",         "Value": round(wb['within_sp'], 3)},
        {"Diagnostic": "Within-Terrestrial mean BC",   "Value": round(wb['within_tr'], 3)},
        {"Diagnostic": "Between-group mean BC",        "Value": round(wb['between'], 3)},
        {"Diagnostic": "Q-centroid PERMANOVA R² (median)",
         "Value": round(df_b['R2'].median(), 3)},
        {"Diagnostic": "Q-centroid PERMANOVA R² 95% CI low",
         "Value": round(df_b['R2'].quantile(0.025), 3)},
        {"Diagnostic": "Q-centroid PERMANOVA R² 95% CI high",
         "Value": round(df_b['R2'].quantile(0.975), 3)},
        {"Diagnostic": "Q-centroid PERMANOVA p<0.05 fraction",
         "Value": round((df_b['p'] < 0.05).mean(), 3)},
        {"Diagnostic": "1-per-Q bootstrap PERMANOVA R² (median)",
         "Value": round(df_c['R2'].median(), 3)},
        {"Diagnostic": "1-per-Q bootstrap PERMANOVA R² 95% CI low",
         "Value": round(df_c['R2'].quantile(0.025), 3)},
        {"Diagnostic": "1-per-Q bootstrap PERMANOVA R² 95% CI high",
         "Value": round(df_c['R2'].quantile(0.975), 3)},
        {"Diagnostic": "1-per-Q bootstrap PERMANOVA p<0.05 fraction",
         "Value": round((df_c['p'] < 0.05).mean(), 3)},
    ]
    out = config.VERSION_DIR / "SuppTable_S1b_nonindependence.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    wb['within_quadrants'].to_csv(config.VERSION_DIR / "SuppTable_S1b_within_quadrant.csv",
                                  index=False)
    print(f"\nSaved: {out}")
    print("Done.")


if __name__ == "__main__":
    main()
