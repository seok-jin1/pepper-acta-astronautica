"""
65_v4_osd65.py — ISS interior surfaces (ISS-MO, OSD-65): nitrifier lineages with detection limits

OSD-65 is used only as presence/absence context: GeneLab-processed V4 ASV table (A519F/802R, NextSeq
single-end, DADA2 + DECIPHER/SILVA r138). Its 46 libraries are split by name, following the study protocol:
  surface   : wipe samples from 8 ISS locations x 2 sampling sessions, with and without PMA
              (names starting with the location number, e.g. 1-F4_S2, 1P-F4P_S1)
  control   : environmental-control wipes (C*, CP*), Maxwell extraction controls (MC*), GFC / FGFC2
  unlabelled: F1-F4 (no location in the name; reported separately, not used)
Nitrifier lineages are the same pre-specified set as in 63/64 (family/genus/class level).

Inputs : $ISS_ENV_DIR/processed/GLDS-65_GAmplicon_16S-{counts,taxonomy}.tsv
Outputs: $V4_WORK/osd65/{libraries.tsv, nitrifier_by_library.tsv, summary.tsv}
"""
import os
import re
from pathlib import Path

import pandas as pd
from scipy import stats

ENV = Path(os.environ.get("ISS_ENV_DIR", Path.home() / "pepper-aa-work/version-2_integrated/external/iss_environment"))
OUT = Path(os.environ.get("V4_WORK", Path.home() / "pepper-v4" / "work")) / "osd65"
OUT.mkdir(parents=True, exist_ok=True)

ct = pd.read_csv(ENV / "processed" / "GLDS-65_GAmplicon_16S-counts.tsv", sep="\t", index_col=0)
tx = pd.read_csv(ENV / "processed" / "GLDS-65_GAmplicon_16S-taxonomy.tsv", sep="\t", index_col=0).reindex(ct.index)
tx = tx.fillna("NA")

organelle = (tx["family"] == "Mitochondria") | (tx["order"] == "Chloroplast") | (tx["class"] == "Chloroplast")
bacteria_archaea = tx["domain"].isin(["Bacteria", "Archaea"])
use = bacteria_archaea & ~organelle

LINEAGES = {
    "AOB Nitrosomonadaceae": tx["family"] == "Nitrosomonadaceae",
    "AOB Nitrosococcus": tx["genus"] == "Nitrosococcus",
    "NOB Nitrospira": tx["genus"] == "Nitrospira",
    "NOB Nitrobacter": tx["genus"] == "Nitrobacter",
    "NOB Nitrococcus": tx["genus"] == "Nitrococcus",
    "NOB Nitrotoga": tx["genus"].str.contains("Nitrotoga"),
    "NOB Nitrolancea": tx["genus"] == "Nitrolancea",
    "AOA Nitrososphaeria": tx["class"] == "Nitrososphaeria",
}
CANON = tx["genus"].isin(["Nitrosomonas", "Nitrosospira", "Nitrospira", "Nitrobacter"])


def kind(name):
    if re.match(r"^\d", name):
        return "surface"
    if re.match(r"^(C|CP|MC|GFC|FGFC)", name):
        return "control"
    return "unlabelled"


lib = pd.DataFrame(index=ct.columns)
lib["kind"] = [kind(s) for s in lib.index]
lib["location"] = [re.match(r"^(\d+)", s).group(1) if kind(s) == "surface" else "" for s in lib.index]
lib["PMA_in_name"] = [bool(re.match(r"^\d+a?P", s)) for s in lib.index]
lib["total_reads"] = ct.sum()
lib["bacterial_reads"] = ct.loc[use].sum()
lib["organelle_fraction"] = 1 - lib["bacterial_reads"] / lib["total_reads"].where(lib["total_reads"] > 0)
lib["genus_assigned_fraction"] = ct.loc[use & (tx["genus"] != "NA")].sum() / lib["bacterial_reads"].where(lib["bacterial_reads"] > 0)
for name, m in LINEAGES.items():
    lib[name] = ct.loc[use & m].sum()
lib["any nitrifier lineage"] = lib[list(LINEAGES)].sum(axis=1)
lib["canonical genera"] = ct.loc[use & CANON].sum()
lib.to_csv(OUT / "nitrifier_by_library.tsv", sep="\t")


def cp(k, n):
    lo = stats.beta.ppf(0.025, k, n - k + 1) if k > 0 else 0.0
    hi = stats.beta.ppf(0.975, k + 1, n - k) if k < n else 1.0
    return f"{lo:.1%}-{hi:.1%}"


rows = []
for grp, sub in [("surface (all)", lib[lib.kind == "surface"]),
                 ("surface, PMA-treated", lib[(lib.kind == "surface") & lib.PMA_in_name]),
                 ("surface, untreated", lib[(lib.kind == "surface") & ~lib.PMA_in_name]),
                 ("surface, >=1,000 bacterial reads", lib[(lib.kind == "surface") & (lib.bacterial_reads >= 1000)]),
                 ("controls", lib[lib.kind == "control"]),
                 ("unlabelled F1-F4", lib[lib.kind == "unlabelled"])]:
    pooled = int(sub.bacterial_reads.sum())
    for setname in ["any nitrifier lineage", "canonical genera"]:
        k = int((sub[setname] > 0).sum()); reads = int(sub[setname].sum())
        upper = stats.chi2.ppf(0.975, 2 * (reads + 1)) / 2 / pooled if pooled else float("nan")
        rows.append({"group": grp, "set": setname, "libraries": len(sub),
                     "bacterial reads (pooled)": pooled, "median bacterial reads": int(sub.bacterial_reads.median()) if len(sub) else 0,
                     "positive": f"{k}/{len(sub)}", "prev_CI95": cp(k, len(sub)), "nitrifier reads": reads,
                     "pooled RA upper 95%": upper})
summ = pd.DataFrame(rows)
summ.to_csv(OUT / "summary.tsv", sep="\t", index=False)

det = tx.loc[use & pd.concat(LINEAGES, axis=1).any(axis=1), ["phylum", "class", "family", "genus"]].copy()
det["reads_surface"] = ct.loc[det.index, lib.index[lib.kind == "surface"]].sum(axis=1)
det["reads_control"] = ct.loc[det.index, lib.index[lib.kind == "control"]].sum(axis=1)
det["reads_unlabelled"] = ct.loc[det.index, lib.index[lib.kind == "unlabelled"]].sum(axis=1)
det["libraries"] = ct.loc[det.index].apply(lambda r: ", ".join(r.index[r > 0]), axis=1)
det.to_csv(OUT / "nitrifier_asvs.tsv", sep="\t")

with pd.option_context("display.width", 230, "display.max_columns", 20, "display.max_colwidth", 80):
    print(lib.groupby("kind").agg(n=("total_reads", "size"), median_total=("total_reads", "median"),
                                  median_bacterial=("bacterial_reads", "median"),
                                  median_organelle=("organelle_fraction", "median"),
                                  median_genus_assigned=("genus_assigned_fraction", "median")).to_string(), "\n")
    print("surface libraries < 1,000 bacterial reads:", int(((lib.kind == "surface") & (lib.bacterial_reads < 1000)).sum()), "\n")
    print(summ.to_string(index=False), "\n")
    print("nitrifier-lineage ASVs anywhere in OSD-65:\n", det.to_string() if len(det) else "  (none)")
