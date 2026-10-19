"""
45_effect_sizes.py
==================
Compute per-function effect sizes (Cliff's delta, Hedges' g) for the
key FAPROTAX functional categories used in the manuscript, on the
analytical root zone (Space n=20 vs Terrestrial n=20).  Also reports
effect size for the N-cycle completeness metric.

Output: version-2_integrated/SuppTable_S1d_effect_sizes.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config

ROOT_ZONE_TYPES = {"Plant Roots", "Plant_Roots", "Rhizosphere"}

KEY_FUNCTIONS = [
    "nitrification",
    "denitrification",
    "nitrate_reduction",
    "nitrogen_fixation",
    "ureolysis",
    "chemoheterotrophy",
    "fermentation",
    "aromatic_compound_degradation",
]


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
    # small-sample correction
    j = 1 - 3 / (4 * (na + nb) - 9)
    return d * j


def main():
    print("=== 45_effect_sizes.py ===")

    meta = pd.read_csv(config.INTEGRATED_METADATA, sep="\t", index_col=0)
    sp_ids = meta.index[(meta["Study_Group"] == "Space_Flight")
                        & meta["sample_type"].isin(ROOT_ZONE_TYPES)]
    tr_ids = meta.index[meta["Study_Group"] == "Terrestrial_Soil"]

    # FAPROTAX functional table → column-normalized to per-sample relative abundance
    fa = pd.read_csv(config.VERSION_DIR / "FAPROTAX_functional_table.tsv",
                     sep="\t", index_col=0)
    col_sums = fa.sum(axis=0).replace(0, 1)
    fa_ra = fa.div(col_sums, axis=1)

    rows = []
    for func in KEY_FUNCTIONS:
        if func not in fa_ra.index:
            continue
        sp = fa_ra.loc[func, [s for s in sp_ids if s in fa_ra.columns]].dropna().values
        tr = fa_ra.loc[func, [s for s in tr_ids if s in fa_ra.columns]].dropna().values
        if len(sp) == 0 or len(tr) == 0:
            continue
        try:
            _, p = mannwhitneyu(sp, tr, alternative="two-sided")
        except ValueError:
            p = np.nan
        log2fc = np.log2((sp.mean() + 1e-6) / (tr.mean() + 1e-6))
        rows.append({
            "metric":            f"FAPROTAX: {func}",
            "Space_n":           len(sp),
            "Terr_n":            len(tr),
            "Space_mean":        round(sp.mean(), 6),
            "Terr_mean":         round(tr.mean(), 6),
            "Log2FC":            round(log2fc, 2),
            "Cliffs_delta":      round(cliffs_delta(sp, tr), 3),
            "Hedges_g":          round(hedges_g(sp, tr), 3),
            "MWU_p":             f"{p:.1e}" if p is not np.nan else "NA",
        })

    # N-cycle completeness
    nc = pd.read_csv(config.VERSION_DIR / "Fig2B_NCycle_Completeness_Results.csv",
                     index_col=0)
    if "Completeness" in nc.columns:
        sp = nc.loc[[s for s in sp_ids if s in nc.index], "Completeness"].dropna().values
        tr = nc.loc[[s for s in tr_ids if s in nc.index], "Completeness"].dropna().values
        if len(sp) and len(tr):
            try:
                _, p = mannwhitneyu(sp, tr, alternative="two-sided")
            except ValueError:
                p = np.nan
            rows.append({
                "metric":            "N-cycle completeness (5 FAPROTAX categories)",
                "Space_n":           len(sp),
                "Terr_n":            len(tr),
                "Space_mean":        round(sp.mean(), 3),
                "Terr_mean":         round(tr.mean(), 3),
                "Log2FC":            "NA",
                "Cliffs_delta":      round(cliffs_delta(sp, tr), 3),
                "Hedges_g":          round(hedges_g(sp, tr), 3),
                "MWU_p":             f"{p:.1e}" if p is not np.nan else "NA",
            })

    df = pd.DataFrame(rows)
    out = config.VERSION_DIR / "SuppTable_S1d_effect_sizes.csv"
    df.to_csv(out, index=False)
    print(df.to_string(index=False))
    print(f"\nSaved: {out}")


if __name__ == "__main__":
    main()
