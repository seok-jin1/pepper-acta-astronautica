"""
38g_lettuce_effect_sizes.py
===========================
Effect-size summary for lettuce — mirrors pepper's 45_effect_sizes.py.

Computes Cliff's δ + Hedges' g + MWU p for:
  - FAPROTAX nitrification Log2FC
  - Per-sample N-cycle completeness
  - Per-sample nitrifier composite RA
  - α-diversity (Shannon, Observed ASVs)

Both aggregate and mission-stratified are reported.

Produces one unified comparison table that can sit next to pepper's
SuppTable_S1d_effect_sizes.csv for direct parallel reading.

Output: acta_astronautica/lettuce/lettuce_effect_sizes.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config

LETTUCE_DIR = config.VEGGIE_LETTUCE_DIR
OUT_DIR = config.AA_DIR / "lettuce"
OUT_DIR.mkdir(parents=True, exist_ok=True)


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


def effect_row(metric, sp, tr, subset):
    sp = np.asarray(sp, dtype=float); tr = np.asarray(tr, dtype=float)
    if len(sp) == 0 or len(tr) == 0:
        return None
    try:
        _, p = mannwhitneyu(sp, tr, alternative="two-sided")
    except ValueError:
        p = np.nan
    log2fc = (np.log2((np.mean(sp) + 1e-6) / (np.mean(tr) + 1e-6))
              if np.mean(sp) + np.mean(tr) > 0 else np.nan)
    return {
        "Subset":        subset,
        "Metric":        metric,
        "Space_n":       len(sp),
        "Terr_n":        len(tr),
        "Space_mean":    round(float(np.mean(sp)), 4),
        "Terr_mean":     round(float(np.mean(tr)), 4),
        "Log2FC":        round(log2fc, 2) if not np.isnan(log2fc) else "NA",
        "Cliffs_delta":  round(cliffs_delta(sp, tr), 3),
        "Hedges_g":      round(hedges_g(sp, tr), 3),
        "MWU_p":         f"{p:.3e}" if not np.isnan(p) else "NA",
    }


def main():
    print("=" * 60)
    print("38g: Effect-size summary for lettuce")
    print("=" * 60)

    rows = []

    # ── (1) Nitrifier composite RA (from 38c output) ─────────────────────────
    nit_file = OUT_DIR / "lettuce_nitrifier_persample.csv"
    if nit_file.exists():
        nit = pd.read_csv(nit_file)
        for subset_name, subset_mask in [
            ("aggregate_n4x4", nit.index == nit.index),
            ("OSD-267_VEG-01A", nit["OSD"] == "OSD-267"),
            ("OSD-269_VEG-03A", nit["OSD"] == "OSD-269"),
        ]:
            sub = nit[subset_mask]
            sp = sub.loc[sub["Study_Group"] == "Space_Flight",
                         "composite_nitrifier_RA"].values
            tr = sub.loc[sub["Study_Group"] == "Terrestrial_Soil",
                         "composite_nitrifier_RA"].values
            r = effect_row("Composite nitrifier RA (%)", sp, tr, subset_name)
            if r:
                rows.append(r)

    # ── (2) N-cycle completeness (from 38d output) ───────────────────────────
    comp_file = OUT_DIR / "lettuce_ncycle_completeness_persample.csv"
    if comp_file.exists():
        comp = pd.read_csv(comp_file)
        for subset_name, subset_mask in [
            ("aggregate_n4x4", comp.index == comp.index),
            ("OSD-267_VEG-01A", comp["OSD"] == "OSD-267"),
            ("OSD-269_VEG-03A", comp["OSD"] == "OSD-269"),
        ]:
            sub = comp[subset_mask]
            sp = sub.loc[sub["Study_Group"] == "Space_Flight", "completeness"].values
            tr = sub.loc[sub["Study_Group"] == "Terrestrial_Soil", "completeness"].values
            r = effect_row("N-cycle completeness (%)", sp, tr, subset_name)
            if r:
                rows.append(r)

    # ── (3) α-diversity (from 38b output) ────────────────────────────────────
    alpha_file = OUT_DIR / "lettuce_alpha_diversity.csv"
    if alpha_file.exists():
        alpha = pd.read_csv(alpha_file)
        for subset_name in alpha["Subset"].unique():
            sub = alpha[alpha["Subset"] == subset_name]
            for metric in ["Shannon", "Observed_ASVs"]:
                sp = sub.loc[sub["Study_Group"] == "Space_Flight", metric].values
                tr = sub.loc[sub["Study_Group"] == "Terrestrial_Soil", metric].values
                r = effect_row(metric, sp, tr, subset_name)
                if r:
                    rows.append(r)

    # ── (4) Per-sample redundancy (from 38d output) — nitrification focus ────
    redund_file = OUT_DIR / "lettuce_persample_redundancy.csv"
    if redund_file.exists():
        red = pd.read_csv(redund_file, index_col=0)
        for subset_name, subset_mask in [
            ("aggregate_n4x4", red.index == red.index),
            ("OSD-267_VEG-01A", red["OSD"] == "OSD-267"),
            ("OSD-269_VEG-03A", red["OSD"] == "OSD-269"),
        ]:
            sub = red[subset_mask]
            if "nitrification" in sub.columns:
                sp = sub.loc[sub["Study_Group"] == "Space_Flight",
                             "nitrification"].values
                tr = sub.loc[sub["Study_Group"] == "Terrestrial_Soil",
                             "nitrification"].values
                r = effect_row("Redundancy: nitrification (# genera)",
                               sp, tr, subset_name)
                if r:
                    rows.append(r)

    df = pd.DataFrame(rows)
    if df.empty:
        print("\n⚠️  No upstream outputs found yet. Run 38b/38c/38d first.")
        return

    out = OUT_DIR / "lettuce_effect_sizes.csv"
    df.to_csv(out, index=False)
    print(df.to_string(index=False))
    print(f"\nSaved: {out}")

    print("\nDone.")


if __name__ == "__main__":
    main()
