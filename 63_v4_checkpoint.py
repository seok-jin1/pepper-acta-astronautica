"""
63_v4_checkpoint.py — checkpoint 1 after the V4 re-processing (60-62)

Answers two questions before any re-writing:
  1. Do the re-processed studies now share ASVs (same V4 window)?
  2. Does the core observation survive: canonical nitrifier lineages not detected in ISS root-zone samples?

Nitrifier lineages are fixed in advance (genus/family level, SILVA 138 names):
  AOB  Nitrosomonadaceae (incl. Nitrosomonas, Nitrosospira, uncultured MND1/Ellin6067/IS-44/mle1-7/DSSD61),
       Nitrosococcus
  NOB  Nitrospira, Nitrobacter, Nitrococcus, Candidatus Nitrotoga, Nitrolancea
  AOA  class Nitrososphaeria
Outputs: $V4_WORK/checkpoint/{sample_groups.tsv, nitrifier_by_sample.tsv, nitrifier_prevalence.tsv, shared_asvs.tsv}
"""
import os
import re
from pathlib import Path

import pandas as pd

WORK = Path(os.environ.get("V4_WORK", Path.home() / "pepper-v4" / "work"))
REPO = Path(__file__).resolve().parent
OUT = WORK / "checkpoint"
OUT.mkdir(parents=True, exist_ok=True)

LINEAGES = {
    "AOB Nitrosomonadaceae": r"f__Nitrosomonadaceae",
    "AOB Nitrosococcus": r"g__Nitrosococcus(;|$)",
    "NOB Nitrospira": r"g__Nitrospira(;|$)",
    "NOB Nitrobacter": r"g__Nitrobacter(;|$)",
    "NOB Nitrococcus": r"g__Nitrococcus(;|$)",
    "NOB Nitrotoga": r"g__Candidatus_Nitrotoga(;|$)",
    "NOB Nitrolancea": r"g__Nitrolancea(;|$)",
    "AOA Nitrososphaeria": r"c__Nitrososphaeria",
}

ft = pd.read_csv(WORK / "dada2" / "feature-table.tsv", sep="\t", index_col=0)
tax = pd.read_csv(WORK / "taxonomy" / "taxonomy.tsv", sep="\t", index_col=0)["Taxon"].reindex(ft.index).fillna("Unassigned")
organelle = tax.str.contains("Mitochondria|Chloroplast", regex=True)
ft = ft.loc[~organelle]
tax = tax.loc[ft.index]

# ---------------------------------------------------------------- sample groups
osd = pd.read_csv(REPO / "osd772_metadata.tsv", sep="\t")
norm = lambda s: re.sub(r"[-_]", "", s.lower())
osd_type = {norm(r["sample-id"]): r["sample_type"] for _, r in osd.iterrows()}
osd_loc = {norm(r["sample-id"]): r["location"] for _, r in osd.iterrows()}
srr = pd.read_csv(REPO / "external_mapping.txt", sep="\t", header=None, index_col=0)[1].to_dict()

rows = []
for s in ft.columns:
    if s.startswith("SRR"):
        rows.append((s, "PRJNA1145089", "Terrestrial field", "Rhizosphere", srr.get(s, "?").split("_")[0]))
    elif re.match(r"(Flight|Ground)Root", s) or s.startswith("V03A"):
        study = "OSD-267" if re.match(r"(Flight|Ground)Root", s) else "OSD-269"
        cond = "Flight" if re.search(r"Flight|Flt", s) else "Ground"
        rows.append((s, study, f"Lettuce {cond}", "Root", study))
    else:
        key = norm(re.sub(r"_S\d+$", "", s))
        st = osd_type.get(key)
        if st is None:      # sample-id spelling differs between file names and metadata
            cand = [k for k in osd_type if k in key or key in k]
            st = osd_type[cand[0]] if len(cand) == 1 else "unknown"
            key = cand[0] if len(cand) == 1 else key
        grp = "ISS pepper" if osd_loc.get(key, "") in ("Q1", "Q2", "Q3", "Q4") else "Ground seed"
        rows.append((s, "OSD-772", grp, st, osd_loc.get(key, "NA")))
groups = pd.DataFrame(rows, columns=["sample", "study", "group", "compartment", "cluster"]).set_index("sample")
groups["reads"] = ft.sum().reindex(groups.index)
groups.to_csv(OUT / "sample_groups.tsv", sep="\t")
print(groups.groupby(["study", "group", "compartment"]).size().to_string(), "\n")

# ---------------------------------------------------------------- 1. shared ASVs
pres = ft > 0
by_study = {st: pres.loc[:, groups.index[groups.study == st]].any(axis=1) for st in groups.study.unique()}
iss = pres.loc[:, groups.index[groups.group == "ISS pepper"]].any(axis=1)
ter = by_study["PRJNA1145089"]
shared = ft.index[iss & ter]
print(f"ASVs: ISS pepper {int(iss.sum())}, Terrestrial {int(ter.sum())}, shared {len(shared)}")
share_reads_ter = ft.loc[shared, groups.index[groups.study == 'PRJNA1145089']].values.sum() / ft.loc[:, groups.index[groups.study == 'PRJNA1145089']].values.sum()
share_reads_iss = ft.loc[shared, groups.index[groups.group == 'ISS pepper']].values.sum() / ft.loc[:, groups.index[groups.group == 'ISS pepper']].values.sum()
print(f"reads in shared ASVs: Terrestrial {share_reads_ter:.1%}, ISS {share_reads_iss:.1%}")
pd.DataFrame({"taxon": tax.loc[shared], "reads_total": ft.loc[shared].sum(axis=1)}).sort_values(
    "reads_total", ascending=False).to_csv(OUT / "shared_asvs.tsv", sep="\t")

# ---------------------------------------------------------------- 2. nitrifier screen
reads = {}
for name, pat in LINEAGES.items():
    hit = tax.str.contains(pat, regex=True)
    reads[name] = ft.loc[hit].sum()
nit = pd.DataFrame(reads)
nit["any nitrifier"] = nit.sum(axis=1)
nit = groups.join(nit)
nit.to_csv(OUT / "nitrifier_by_sample.tsv", sep="\t")

cols = list(LINEAGES) + ["any nitrifier"]
prev = nit.groupby(["group", "compartment"]).apply(
    lambda d: pd.Series({**{"n": len(d), "min_reads": int(d.reads.min())},
                         **{c: f"{int((d[c] > 0).sum())}/{len(d)}" for c in cols}}))
prev.to_csv(OUT / "nitrifier_prevalence.tsv", sep="\t")
with pd.option_context("display.width", 250, "display.max_columns", 20):
    print(prev.to_string())

rz = nit[(nit.group == "ISS pepper") & nit.compartment.isin(["Plant Roots", "Rhizosphere"])]
print(f"\nCORE CHECK — ISS root zone (n={len(rz)}): samples with any nitrifier read = {(rz['any nitrifier'] > 0).sum()}; "
      f"total nitrifier reads = {int(rz['any nitrifier'].sum())} of {int(rz.reads.sum())}")
print(rz.loc[rz["any nitrifier"] > 0, cols].to_string() if (rz["any nitrifier"] > 0).any() else "  (none)")
