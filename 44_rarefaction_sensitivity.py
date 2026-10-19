"""
44_rarefaction_sensitivity.py
=============================
Rarefaction-depth sensitivity check for the core nitrifier signal.

At three rarefaction depths (1k, 5k, 10k reads/sample), compute on the
analytical root zone (Space n=20 + Terr n=20):

  - Nitrifier composite prevalence (Space %, Terr %) for any of
    {Nitrosomonas, Nitrosospira, Nitrospira, Nitrobacter}
  - Nitrification Log2FC from genus-level relative abundance
    (mean Space / mean Terr with 1e-6 pseudocount)
  - Sample retention at that depth

Output: version-2_integrated/SuppTable_S1c_rarefaction_sensitivity.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config

DEPTHS = [1000, 5000, 10000]
RANDOM_SEED = config.RANDOM_SEED
NITRIFIER_GENERA = ["Nitrosomonas", "Nitrosospira", "Nitrospira", "Nitrobacter"]
ROOT_ZONE_TYPES = {"Plant Roots", "Plant_Roots", "Rhizosphere"}


def rarefy_sample(counts: np.ndarray, depth: int, rng) -> np.ndarray:
    """Subsample without replacement to `depth` reads and return counts."""
    counts = counts.astype(int)
    total = counts.sum()
    if total < depth:
        return None  # drop sample
    pool = np.repeat(np.arange(len(counts)), counts)
    chosen = rng.choice(pool, size=depth, replace=False)
    return np.bincount(chosen, minlength=len(counts)).astype(float)


def parse_genus(taxon_str):
    if pd.isna(taxon_str):
        return "Unassigned"
    for part in reversed(str(taxon_str).split(";")):
        p = part.strip()
        if p.startswith("g__"):
            name = p.split("__", 1)[-1].strip()
            if name and name not in ("uncultured", "metagenome"):
                return name
    return "Unassigned"


def main():
    print("=== 44_rarefaction_sensitivity.py ===")

    ft = pd.read_csv(config.FEATURE_TABLE_CLEAN, sep="\t", skiprows=1, index_col=0)
    tax = pd.read_csv(config.EXPORTED_TAXONOMY / "taxonomy.tsv", sep="\t", index_col=0)
    meta = pd.read_csv(config.INTEGRATED_METADATA, sep="\t", index_col=0)

    sp_ids = meta.index[(meta["Study_Group"] == "Space_Flight")
                        & meta["sample_type"].isin(ROOT_ZONE_TYPES)]
    tr_ids = meta.index[meta["Study_Group"] == "Terrestrial_Soil"]
    keep = [s for s in list(sp_ids) + list(tr_ids) if s in ft.columns]
    ft_keep = ft[keep]

    shared = ft_keep.index.intersection(tax.index)
    ft_keep = ft_keep.loc[shared]
    genus = tax.loc[shared, "Taxon"].apply(parse_genus)

    rows = []
    for depth in DEPTHS:
        rng = np.random.default_rng(RANDOM_SEED)
        rarefied = {}
        for s in ft_keep.columns:
            out = rarefy_sample(ft_keep[s].values, depth, rng)
            if out is not None:
                rarefied[s] = out
        ra_df = pd.DataFrame(rarefied, index=ft_keep.index)
        # convert counts → relative abundance (%)
        ra_df = ra_df.div(ra_df.sum(axis=0), axis=1) * 100

        # aggregate to genus level
        ra_df["Genus"] = genus.values
        genus_ra = ra_df.groupby("Genus").sum().drop(index="Unassigned", errors="ignore")

        # restrict to nitrifier genera actually present
        present = [g for g in NITRIFIER_GENERA if g in genus_ra.index]
        comp_ra = genus_ra.loc[present].sum(axis=0) if present else pd.Series(
            0.0, index=genus_ra.columns)

        sp_in = [s for s in sp_ids if s in comp_ra.index]
        tr_in = [s for s in tr_ids if s in comp_ra.index]

        sp_vals = comp_ra.loc[sp_in].values
        tr_vals = comp_ra.loc[tr_in].values

        prev_sp = (sp_vals > 0).mean() * 100
        prev_tr = (tr_vals > 0).mean() * 100

        log2fc = np.log2((sp_vals.mean() + 1e-6) / (tr_vals.mean() + 1e-6))

        rows.append({
            "Rarefaction_depth": depth,
            "Space_n_retained": len(sp_in),
            "Terr_n_retained":  len(tr_in),
            "Nitrifier_genera_detected": ", ".join(present) if present else "none",
            "Space_prevalence_%":        round(prev_sp, 1),
            "Terr_prevalence_%":         round(prev_tr, 1),
            "Nitrifier_composite_Log2FC": round(log2fc, 2),
        })
        print(f"\n  Depth {depth}: Space retained {len(sp_in)}/20, "
              f"Terr retained {len(tr_in)}/20")
        print(f"    Space nitrifier prevalence: {prev_sp:.0f}% "
              f"({int((sp_vals>0).sum())}/{len(sp_in)})")
        print(f"    Terr  nitrifier prevalence: {prev_tr:.0f}% "
              f"({int((tr_vals>0).sum())}/{len(tr_in)})")
        print(f"    Nitrifier composite Log2FC: {log2fc:.2f}")

    df = pd.DataFrame(rows)
    out = config.VERSION_DIR / "SuppTable_S1c_rarefaction_sensitivity.csv"
    df.to_csv(out, index=False)
    print(f"\nSaved: {out}")


if __name__ == "__main__":
    main()
