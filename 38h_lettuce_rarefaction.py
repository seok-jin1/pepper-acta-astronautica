"""
38h_lettuce_rarefaction.py
==========================
Rarefaction-depth sensitivity for lettuce — mirrors pepper's
44_rarefaction_sensitivity.py.

Evaluates nitrifier composite prevalence + Log2FC + sample retention at
multiple rarefaction depths (100, 500, 1000, 5000). Given the low depth of
OSD-267 ground root samples, the depth that retains all 8 samples may be
quite low; this script documents that explicitly.

Uses the DADA2-generated RAW (un-rarefied) counts as input.

Output: acta_astronautica/lettuce/lettuce_rarefaction_sensitivity.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config

LETTUCE_DIR = config.VEGGIE_LETTUCE_DIR
EXPORTS_DIR = LETTUCE_DIR / "exports"
OUT_DIR = config.AA_DIR / "lettuce"
OUT_DIR.mkdir(parents=True, exist_ok=True)

DEPTHS = [100, 500, 1000, 2000, 5000]
RANDOM_SEED = config.RANDOM_SEED
NITRIFIER_GENERA = ["Nitrosomonas", "Nitrosospira", "Nitrospira", "Nitrobacter"]


def rarefy_sample(counts: np.ndarray, depth: int, rng) -> np.ndarray:
    counts = counts.astype(int)
    total = counts.sum()
    if total < depth:
        return None
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
            if name and name not in ("uncultured", "metagenome", "unidentified"):
                return name
    return "Unassigned"


def main():
    print("=" * 60)
    print("38h: Rarefaction-depth sensitivity for lettuce")
    print("=" * 60)

    # Load RAW (un-rarefied) table — this may be a separate file.
    # If we only have rarefied output, fall back to that with a flag.
    raw_ft_path = EXPORTS_DIR / "feature-table-raw.tsv"
    rare_ft_path = EXPORTS_DIR / "feature-table.tsv"
    if raw_ft_path.exists():
        ft_path = raw_ft_path
        print(f"  Using raw (un-rarefied) table: {ft_path}")
    else:
        ft_path = rare_ft_path
        print(f"  NOTE: raw table not found, using rarefied table: {ft_path}")

    ft = pd.read_csv(ft_path, sep="\t", skiprows=1, index_col=0)
    tax = pd.read_csv(EXPORTS_DIR / "taxonomy.tsv", sep="\t", index_col=0)
    meta = pd.read_csv(LETTUCE_DIR / "metadata_harmonized.tsv", sep="\t", index_col=0)

    sample_ids = [s for s in meta.index if s in ft.columns]
    ft = ft[sample_ids]
    print(f"  Input: {ft.shape[0]} ASVs × {ft.shape[1]} samples")
    print(f"  Per-sample depths: min={ft.sum().min()}, median={ft.sum().median()},"
          f" max={ft.sum().max()}")

    shared = ft.index.intersection(tax.index)
    ft = ft.loc[shared]
    genus = tax.loc[shared, "Taxon"].apply(parse_genus)

    sp_all = meta.index[meta["Study_Group"] == "Space_Flight"]
    tr_all = meta.index[meta["Study_Group"] == "Terrestrial_Soil"]

    rows = []
    for depth in DEPTHS:
        rng = np.random.default_rng(RANDOM_SEED)
        rare = {}
        for s in ft.columns:
            out = rarefy_sample(ft[s].values, depth, rng)
            if out is not None:
                rare[s] = out
        if len(rare) == 0:
            rows.append({
                "Rarefaction_depth":  depth,
                "Samples_retained":   0,
                "Reason":             "all samples below depth",
            })
            continue

        ra_df = pd.DataFrame(rare, index=ft.index)
        ra_df = ra_df.div(ra_df.sum(axis=0), axis=1) * 100
        ra_df["Genus"] = genus.values
        genus_ra = ra_df.groupby("Genus").sum().drop(index="Unassigned", errors="ignore")

        present = [g for g in NITRIFIER_GENERA if g in genus_ra.index]
        comp_ra = (genus_ra.loc[present].sum(axis=0) if present
                   else pd.Series(0.0, index=genus_ra.columns))

        sp_in = [s for s in sp_all if s in comp_ra.index]
        tr_in = [s for s in tr_all if s in comp_ra.index]
        sp_vals = comp_ra.loc[sp_in].values
        tr_vals = comp_ra.loc[tr_in].values

        prev_sp = (sp_vals > 0).mean() * 100 if len(sp_vals) else np.nan
        prev_tr = (tr_vals > 0).mean() * 100 if len(tr_vals) else np.nan
        log2fc = (np.log2((sp_vals.mean() + 1e-6) / (tr_vals.mean() + 1e-6))
                  if len(sp_vals) and len(tr_vals) else np.nan)

        rows.append({
            "Rarefaction_depth":         depth,
            "Space_n_retained":          len(sp_in),
            "Terr_n_retained":           len(tr_in),
            "Nitrifier_genera_detected": ", ".join(present) if present else "none",
            "Space_prevalence_%":        round(prev_sp, 1) if not np.isnan(prev_sp) else "NA",
            "Terr_prevalence_%":         round(prev_tr, 1) if not np.isnan(prev_tr) else "NA",
            "Nitrifier_composite_Log2FC": round(log2fc, 2) if not np.isnan(log2fc) else "NA",
        })
        print(f"\n  Depth {depth}: Space {len(sp_in)}/4 retained, Terr {len(tr_in)}/4 retained")
        if present:
            print(f"    Nitrifier prevalence: Space {prev_sp:.0f}%, Terr {prev_tr:.0f}%")

    df = pd.DataFrame(rows)
    out = OUT_DIR / "lettuce_rarefaction_sensitivity.csv"
    df.to_csv(out, index=False)
    print(f"\nSaved: {out}")

    print("\nDone.")


if __name__ == "__main__":
    main()
