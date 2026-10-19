"""
Script 38: Veggie lettuce downstream analysis for Acta Astronautica Prediction P3.

Runs FAPROTAX + N-cycle completeness + functional redundancy on Veggie lettuce
data (OSD-267 + OSD-269), then merges with pepper results for cross-mission
comparison.

Depends on:
  - Script 37 QIIME2 export (feature-table.tsv, taxonomy.tsv)
  - FAPROTAX 1.2.12 installation (same as script 12)
  - Pepper results: FAPROTAX_group_comparison.csv, NCycle, FunctionalRedundancy

Outputs:
  - config.AA_VEGGIE_RESULTS_CSV: merged pepper+lettuce FAPROTAX Log2FC table
  - config.AA_DIR / "veggie_ncycle_completeness.csv"
  - config.AA_DIR / "veggie_functional_redundancy.csv"
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import config

FAPROTAX_DIR = config.FAPROTAX_DIR
COLLAPSE_SCRIPT = config.FAPROTAX_COLLAPSE_SCRIPT
FAPROTAX_DB = config.FAPROTAX_DB

LETTUCE_DIR = config.VEGGIE_LETTUCE_DIR

NCYCLE_STEPS = [
    "nitrogen_fixation",
    "nitrification",
    "denitrification",
    "nitrate_reduction",
    "ureolysis",
]

# SILVA v138.1+ renamed phyla; FAPROTAX 1.2.12 uses old names.
# Apply these substitutions to taxonomy strings before FAPROTAX input.
SILVA_PHYLUM_MAP = {
    "Pseudomonadota": "Proteobacteria",
    "Actinomycetota": "Actinobacteria",
    "Bacillota": "Firmicutes",
    "Bacteroidota": "Bacteroidetes",
    "Cyanobacteriota": "Cyanobacteria",
    "Verrucomicrobiota": "Verrucomicrobia",
    "Planctomycetota": "Planctomycetes",
    "Bdellovibrionota": "Bdellovibrionota",  # no change
    "Chloroflexota": "Chloroflexi",
    "Acidobacteriota": "Acidobacteria",
    "Myxococcota": "Deltaproteobacteria",
    "Desulfobacterota": "Deltaproteobacteria",
}


def load_lettuce_data():
    """Load and merge GeneLab-processed counts + taxonomy from OSD-267 and OSD-269.

    GeneLab tables use split taxonomy columns (domain, phylum, ..., species)
    rather than QIIME2-style semicolon-delimited 'Taxon' strings. We
    reconstruct the Taxon string for FAPROTAX compatibility.
    """
    meta = pd.read_csv(LETTUCE_DIR / "metadata_harmonized.tsv", sep="\t", index_col=0)
    root_samples = set(meta.index)

    dfs, taxs = [], []
    for osd, count_file, tax_file in [
        ("OSD-267",
         LETTUCE_DIR / "OSD-267" / "processed" / "GLDS-267_GAmplicon_counts_16S_GLAmpSeq.tsv",
         LETTUCE_DIR / "OSD-267" / "processed" / "GLDS-267_GAmplicon_taxonomy_16S_GLAmpSeq.tsv"),
        ("OSD-269",
         LETTUCE_DIR / "OSD-269" / "processed" / "GLDS-269_GAmplicon_counts.tsv",
         LETTUCE_DIR / "OSD-269" / "processed" / "GLDS-269_GAmplicon_taxonomy.tsv"),
    ]:
        ct = pd.read_csv(count_file, sep="\t", index_col=0)
        # Keep only root samples
        ct = ct[[c for c in ct.columns if c in root_samples]]
        # Prefix ASV IDs to avoid collision between datasets
        ct.index = osd + "_" + ct.index.astype(str)
        dfs.append(ct)

        tx = pd.read_csv(tax_file, sep="\t", index_col=0)
        # Reconstruct SILVA-style Taxon string from split columns
        levels = ["domain", "phylum", "class", "order", "family", "genus", "species"]
        prefixes = ["d__", "p__", "c__", "o__", "f__", "g__", "s__"]
        import re
        def build_taxon(row):
            parts = []
            for lev, pfx in zip(levels, prefixes):
                val = str(row.get(lev, "")).strip()
                if val in ("", "nan", "NA", "Unassigned"):
                    val = ""
                else:
                    val = re.sub(r"_\d+$", "", val)
                    # Translate new SILVA phylum names to old names for FAPROTAX 1.2.12
                    if lev == "phylum" and val in SILVA_PHYLUM_MAP:
                        val = SILVA_PHYLUM_MAP[val]
                parts.append(pfx + val)
            return "; ".join(parts)
        tx["Taxon"] = tx.apply(build_taxon, axis=1)
        tx.index = osd + "_" + tx.index.astype(str)
        taxs.append(tx[["Taxon"]])

    tbl = pd.concat(dfs, axis=0).fillna(0).astype(int)
    tax = pd.concat(taxs, axis=0)
    return tbl, tax, meta


def run_faprotax_on_lettuce(tbl, tax, meta):
    keep = meta.index[meta["Study_Group"].isin(["Space_Flight", "Terrestrial_Soil"])]
    tbl = tbl[[c for c in tbl.columns if c in keep]]

    shared = tbl.index.intersection(tax.index)
    tbl = tbl.loc[shared]

    def clean_taxon(t):
        if pd.isna(t):
            return "Unassigned"
        return "; ".join(p.strip() for p in t.split(";"))

    tax_clean = tax.loc[shared, "Taxon"].apply(clean_taxon)
    tbl_copy = tbl.copy()
    tbl_copy["taxon"] = tax_clean.values
    collapsed = tbl_copy.groupby("taxon").sum()
    collapsed.index.name = "OTU_ID"

    inp = config.AA_DIR / "lettuce_faprotax_input.tsv"
    collapsed.to_csv(inp, sep="\t")

    func_out = config.AA_DIR / "lettuce_FAPROTAX_functional_table.tsv"
    report = config.AA_DIR / "lettuce_FAPROTAX_report.txt"
    cmd = [
        sys.executable,
        str(COLLAPSE_SCRIPT),
        "-i", str(inp),
        "-o", str(func_out),
        "-g", str(FAPROTAX_DB),
        "-r", str(report),
        "--row_names_are_in_column", "OTU_ID",
        "-f", "-v",
    ]
    print("[38] Running FAPROTAX on lettuce data...", flush=True)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("STDERR:", r.stderr[-2000:])
        raise RuntimeError("FAPROTAX failed on lettuce")
    return func_out


def compute_lettuce_log2fc(func_path, meta):
    func = pd.read_csv(func_path, sep="\t", index_col=0, comment="#")
    func_rel = func.div(func.sum(axis=0).replace(0, 1), axis=1) * 100

    sp = [c for c in func_rel.columns if c in meta.index[meta["Study_Group"] == "Space_Flight"]]
    tr = [c for c in func_rel.columns if c in meta.index[meta["Study_Group"] == "Terrestrial_Soil"]]

    eps = 1e-6
    result = pd.DataFrame({
        "Space_Flight_mean": func_rel[sp].mean(axis=1) if sp else 0,
        "Terrestrial_Soil_mean": func_rel[tr].mean(axis=1) if tr else 0,
    })
    result["Log2FC"] = np.log2(
        (result["Space_Flight_mean"] + eps) / (result["Terrestrial_Soil_mean"] + eps)
    )
    result.index.name = "group"
    return result


def compute_ncycle_completeness(func_path, meta):
    func = pd.read_csv(func_path, sep="\t", index_col=0, comment="#")
    func_binary = (func > 0).astype(int)

    rows = []
    for sid in func.columns:
        if sid not in meta.index:
            continue
        present = [s for s in NCYCLE_STEPS if s in func_binary.index and func_binary.loc[s, sid] > 0]
        rows.append({
            "sample_id": sid,
            "Study_Group": meta.loc[sid, "Study_Group"],
            "n_steps": len(present),
            "completeness": len(present) / len(NCYCLE_STEPS) * 100,
            "present_steps": ";".join(present),
        })
    return pd.DataFrame(rows)


def compute_functional_redundancy(tbl, tax, meta):
    keep = meta.index[meta["Study_Group"].isin(["Space_Flight", "Terrestrial_Soil"])]
    tbl = tbl[[c for c in tbl.columns if c in keep]]
    shared = tbl.index.intersection(tax.index)
    tbl = tbl.loc[shared]

    genus_col = tax.loc[shared, "Taxon"].apply(
        lambda t: next((p.strip().replace("g__", "") for p in str(t).split(";") if "g__" in p), "Unknown")
    )
    tbl_g = tbl.copy()
    tbl_g["genus"] = genus_col.values
    genus_agg = tbl_g.groupby("genus").sum()

    target_functions = ["nitrification", "denitrification", "nitrogen_fixation", "chemoheterotrophy"]
    nitrifier_genera = {"Nitrosomonas", "Nitrosospira", "Nitrospira", "Nitrobacter"}
    denitrifier_genera = {"Rhodanobacter", "Pseudomonas", "Azospirillum", "Burkholderia", "Paenibacillus"}
    nfixer_genera = {"Bradyrhizobium", "Azospirillum", "Rhizobium", "Frankia", "Beijerinckia"}
    chemohet_genera = set(genus_agg.index)

    func_genera = {
        "nitrification": nitrifier_genera,
        "denitrification": denitrifier_genera,
        "nitrogen_fixation": nfixer_genera,
        "chemoheterotrophy": chemohet_genera,
    }

    rows = []
    for func, genera_set in func_genera.items():
        present = genera_set.intersection(genus_agg.index)
        for grp in ["Space_Flight", "Terrestrial_Soil"]:
            grp_samples = [c for c in genus_agg.columns if c in meta.index[meta["Study_Group"] == grp]]
            n_contributing = 0
            for g in present:
                if genus_agg.loc[g, grp_samples].sum() > 0:
                    n_contributing += 1
            rows.append({
                "function": func,
                "Study_Group": grp,
                "n_contributing_genera": n_contributing,
            })
    return pd.DataFrame(rows)


def merge_pepper_lettuce_long(pepper_df, lettuce_df):
    pepper_df = pepper_df.assign(crop="pepper")
    lettuce_df = lettuce_df.assign(crop="lettuce")
    return pd.concat([pepper_df, lettuce_df], ignore_index=True)


def main():
    config.AA_DIR.mkdir(parents=True, exist_ok=True)
    print("[38] Loading lettuce data...", flush=True)
    tbl, tax, meta = load_lettuce_data()
    print(f"[38] Samples: {tbl.shape[1]}, ASVs: {tbl.shape[0]}", flush=True)

    # FAPROTAX
    func_path = run_faprotax_on_lettuce(tbl, tax, meta)
    lettuce_log2fc = compute_lettuce_log2fc(func_path, meta)
    print("[38] Lettuce FAPROTAX Log2FC (key functions):")
    for fn in ["nitrification", "denitrification", "nitrogen_fixation", "chemoheterotrophy"]:
        if fn in lettuce_log2fc.index:
            print(f"  {fn}: {lettuce_log2fc.loc[fn, 'Log2FC']:.2f}")

    # Merge with pepper
    pepper_log2fc = pd.read_csv(
        config.VERSION_DIR / "FAPROTAX_group_comparison.csv", index_col=0
    )
    pepper_slim = pepper_log2fc[["Log2FC"]].rename(columns={"Log2FC": "log2fc"})
    pepper_slim.index.name = "function"
    lettuce_slim = lettuce_log2fc[["Log2FC"]].rename(columns={"Log2FC": "log2fc"})
    lettuce_slim.index.name = "function"
    merged = merge_pepper_lettuce_long(pepper_slim.reset_index(), lettuce_slim.reset_index())
    merged.to_csv(config.AA_VEGGIE_RESULTS_CSV, index=False)
    print(f"[38] Wrote {config.AA_VEGGIE_RESULTS_CSV}", flush=True)

    # N-cycle completeness
    ncycle = compute_ncycle_completeness(func_path, meta)
    ncycle.to_csv(config.AA_DIR / "veggie_ncycle_completeness.csv", index=False)
    print("[38] N-cycle completeness (lettuce):")
    print(ncycle.groupby("Study_Group")["completeness"].mean().to_string())

    # Functional redundancy
    redund = compute_functional_redundancy(tbl, tax, meta)
    redund.to_csv(config.AA_DIR / "veggie_functional_redundancy.csv", index=False)
    print("[38] Functional redundancy (lettuce):")
    print(redund.to_string())


if __name__ == "__main__":
    main()
