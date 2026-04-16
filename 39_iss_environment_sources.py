"""
Script 39: ISS environmental source-pool cross-reference.

Compares the genus-level microbiota of the ISS Space_Flight pepper rhizosphere
(our study) against the ISS-MO Microbial Tracking-1A/1B surface dataset
(OSD-65, Checinska Sielaff et al. 2019). Extends FEAST source-tracking (script
22-24) by asking whether pepper rhizosphere Space_Flight genera are plausibly
derived from the ISS internal environment.

Outputs:
    config.AA_ISS_ENV_SOURCES_CSV  -- overlap statistics
    acta_astronautica/iss_env_genus_overlap.csv  -- genus-level intersection
    acta_astronautica/iss_env_denit_genera.csv   -- N-cycle genus focus
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import config

ISS_ENV_DIR = config.ISS_ENV_DIR
ISS_COUNTS = ISS_ENV_DIR / "processed" / "GLDS-65_GAmplicon_16S-counts.tsv"
ISS_TAX = ISS_ENV_DIR / "processed" / "GLDS-65_GAmplicon_16S-taxonomy.tsv"

NCYCLE_GENERA_SET = {
    "Rhodanobacter", "Nitrosomonas", "Nitrosospira", "Nitrospira", "Nitrobacter",
    "Pseudomonas", "Bradyrhizobium", "Azospirillum", "Bacillus", "Paenibacillus",
    "Burkholderia", "Variovorax", "Sphingomonas",
}


def _strip_silva_suffix(s: str) -> str:
    return re.sub(r"_\d+$", "", str(s).strip())


def load_iss_genera() -> pd.Series:
    """Return a Series indexed by genus; value = total read count across all ISS samples."""
    ct = pd.read_csv(ISS_COUNTS, sep="\t", index_col=0)
    tx = pd.read_csv(ISS_TAX, sep="\t", index_col=0)
    shared = ct.index.intersection(tx.index)
    ct = ct.loc[shared]
    genera = tx.loc[shared, "genus"].fillna("Unassigned").apply(_strip_silva_suffix)
    totals = ct.sum(axis=1)
    return totals.groupby(genera).sum().sort_values(ascending=False)


def load_pepper_space_genera() -> pd.Series:
    """Return Space_Flight pepper rhizosphere genus total counts."""
    tbl = pd.read_csv(config.FEATURE_TABLE_CLEAN, sep="\t", skiprows=1, index_col=0)
    tax = pd.read_csv(config.TAXONOMY_FILE, sep="\t", index_col=0)
    meta = pd.read_csv(config.INTEGRATED_METADATA, sep="\t", index_col=0)

    space_samples = meta.index[meta["Study_Group"] == "Space_Flight"].tolist()
    space_cols = [c for c in tbl.columns if c in space_samples]
    tbl_s = tbl[space_cols]

    shared = tbl_s.index.intersection(tax.index)
    tbl_s = tbl_s.loc[shared]
    genus_strings = tax.loc[shared, "Taxon"].apply(
        lambda t: next(
            (p.strip().replace("g__", "") for p in str(t).split(";") if "g__" in p),
            "Unassigned",
        )
    ).apply(_strip_silva_suffix)
    totals = tbl_s.sum(axis=1)
    return totals.groupby(genus_strings).sum().sort_values(ascending=False)


def compute_overlap(iss: pd.Series, pepper_space: pd.Series) -> dict:
    iss_present = set(iss[iss > 0].index) - {"Unassigned", ""}
    pepper_present = set(pepper_space[pepper_space > 0].index) - {"Unassigned", ""}
    overlap = iss_present & pepper_present
    pepper_only = pepper_present - iss_present
    iss_only = iss_present - pepper_present
    return {
        "n_iss_genera": len(iss_present),
        "n_pepper_space_genera": len(pepper_present),
        "n_overlap": len(overlap),
        "n_pepper_only": len(pepper_only),
        "n_iss_only": len(iss_only),
        "pct_pepper_recoverable_from_iss_env": round(
            100 * len(overlap) / max(1, len(pepper_present)), 2
        ),
        "overlap_genera_sample": sorted(overlap)[:20],
        "pepper_only_sample": sorted(pepper_only)[:20],
    }


def main() -> int:
    config.AA_DIR.mkdir(parents=True, exist_ok=True)
    print("[39] Loading ISS env genera (OSD-65)...", flush=True)
    iss = load_iss_genera()
    print(f"[39] ISS env: {len(iss[iss > 0])} genera present, top 5:")
    print(iss.head(5).to_string())

    print("\n[39] Loading pepper Space rhizosphere genera...", flush=True)
    pepper_space = load_pepper_space_genera()
    print(f"[39] Pepper Space: {len(pepper_space[pepper_space > 0])} genera present, top 5:")
    print(pepper_space.head(5).to_string())

    print("\n[39] Computing overlap...")
    overlap = compute_overlap(iss, pepper_space)
    for k, v in overlap.items():
        if isinstance(v, list):
            print(f"  {k}: {v[:5]}... (n={len(v)})")
        else:
            print(f"  {k}: {v}")

    summary = pd.Series({k: v for k, v in overlap.items() if not isinstance(v, list)})
    summary.to_csv(config.AA_ISS_ENV_SOURCES_CSV)
    print(f"\n[39] Wrote {config.AA_ISS_ENV_SOURCES_CSV}")

    # Detailed intersection table
    both = iss.index.intersection(pepper_space.index)
    detail = pd.DataFrame({
        "iss_env_count": iss.reindex(both, fill_value=0),
        "pepper_space_count": pepper_space.reindex(both, fill_value=0),
    }).sort_values("pepper_space_count", ascending=False)
    detail.to_csv(config.AA_DIR / "iss_env_genus_overlap.csv")

    # Focus on N-cycle genera
    ncycle = pd.DataFrame(index=sorted(NCYCLE_GENERA_SET))
    ncycle["iss_env_count"] = [int(iss.get(g, 0)) for g in ncycle.index]
    ncycle["pepper_space_count"] = [int(pepper_space.get(g, 0)) for g in ncycle.index]
    ncycle["in_iss_env"] = ncycle["iss_env_count"] > 0
    ncycle["in_pepper_space"] = ncycle["pepper_space_count"] > 0
    ncycle.to_csv(config.AA_DIR / "iss_env_ncycle_genera.csv")
    print(f"[39] Wrote N-cycle genus table:")
    print(ncycle.to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
