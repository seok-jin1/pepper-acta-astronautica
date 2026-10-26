"""
66_v4_figures.py — manuscript figures and table data for the V4 re-analysis (outputs of 61-65)

Fig. 1  amplicon coverage of each dataset on the 16S rRNA gene and the shared 240-nt V4 window
Fig. 2  nitrifier-lineage relative abundance per sample, by compartment, with prevalence (k/n)
Fig. 3  within-mission comparison: (a) substrate vs root zone per growth quadrant, (b) lineage composition
Fig. 4  detection power: P(>=1 nitrifier read) vs assumed relative abundance at observed sequencing depths
Tables  table1_datasets.tsv, table2_detection.tsv
"""
import os
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

WORK = Path(os.environ.get("V4_WORK", Path.home() / "pepper-v4" / "work"))
OUT = WORK / "figures"
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.spines.top": False, "axes.spines.right": False})
C_ISS, C_TER, C_LET, C_SURF = "#c0392b", "#1f5f99", "#7f7f7f", "#d68910"

ft = pd.read_csv(WORK / "dada2" / "feature-table.tsv", sep="\t", index_col=0)
tax = pd.read_csv(WORK / "taxonomy" / "taxonomy.tsv", sep="\t", index_col=0)["Taxon"].reindex(ft.index).fillna("Unassigned")
keep = tax.str.contains(r"^d__(?:Bacteria|Archaea)") & ~tax.str.contains("Mitochondria|Chloroplast")  # as in 65 and 67
ft, tax = ft.loc[keep], tax.loc[keep]
groups = pd.read_csv(WORK / "checkpoint" / "sample_groups.tsv", sep="\t", index_col=0)
groups.loc[groups.index.str.startswith("UnsanSeed"), ["group", "compartment"]] = ["Ground seed", "seed"]
groups["reads"] = ft.sum().reindex(groups.index)
groups["quadrant"] = groups.index.str.extract(r"^(Q\d)", expand=False).fillna("").to_numpy()
ANY = (r"f__Nitrosomonadaceae|g__Nitrosococcus(?:;|$)|g__Nitrospira(?:;|$)|g__Nitrobacter(?:;|$)|g__Nitrococcus(?:;|$)|"
       r"g__Candidatus_Nitrotoga(?:;|$)|g__Nitrolancea(?:;|$)|c__Nitrososphaeria")
CANON = r"g__(?:Nitrosomonas|Nitrosospira|Nitrospira)(?:;|$)"
nit_mask = tax.str.contains(ANY, regex=True)
nit = ft.loc[nit_mask].sum()
can = ft.loc[tax.str.contains(CANON, regex=True)].sum()


def unit(row):
    if row.group == "ISS pepper":
        return "ISS root zone" if row.compartment in ("Plant Roots", "Rhizosphere") else f"ISS {row.compartment}"
    return row.group


groups["unit"] = groups.apply(unit, axis=1)
groups["RA"] = nit.reindex(groups.index) / groups["reads"]
groups["canRA"] = can.reindex(groups.index) / groups["reads"]

osd65 = pd.read_csv(WORK / "osd65" / "nitrifier_by_library.tsv", sep="\t", index_col=0)
osd65 = osd65[osd65.kind == "surface"]
osd65_ra = osd65["any nitrifier lineage"] / osd65["bacterial_reads"]
osd65_can = osd65["canonical genera"] / osd65["bacterial_reads"]
ECOLI = r"$\it{E.\ coli}$"

# ------------------------------------------------------------------ Fig. 1
fig, ax = plt.subplots(figsize=(7.0, 2.6))
spans = [("Field pepper rhizosphere (PRJNA1145089)\nV3–V4, 341F–805R, merged reads 407–449 nt", 358, 784, C_TER),
         ("ISS pepper, APH PH-04 (OSD-772)\nV4, 515F–806R, 2 × 250 nt", 534, 786, C_ISS),
         ("Veggie lettuce (OSD-267, OSD-269)\nV4, 515F–806R, 2 × 250 nt", 534, 786, C_LET),
         ("ISS interior surfaces (OSD-65)\nV4, A519F–802R, single reads", 534, 784, C_SURF)]
for i, (lab, a, b, c) in enumerate(spans):
    y = len(spans) - i
    ax.plot([a, b], [y, y], color=c, lw=6, solid_capstyle="butt", alpha=0.85)
    ax.text(330, y, lab, ha="right", va="center", fontsize=7)
ax.axvspan(534, 773, color="#2e8b57", alpha=0.12)
ax.text(653, len(spans) + 0.75, f"shared 240-nt V4 window ({ECOLI} 534–773)\npepper and lettuce reads trimmed to it, denoised per study",
        ha="center", fontsize=7, color="#1e6b44")
for pos, lab in [(341, "341F"), (515, "515F"), (806, "806R")]:
    ax.axvline(pos, color="#999999", lw=0.6, ls=":")
    ax.text(pos, 0.35, lab, ha="center", fontsize=7, color="#666666")
ax.set_xlim(320, 840); ax.set_ylim(0.1, len(spans) + 1.3)
ax.set_yticks([]); ax.spines["left"].set_visible(False)
ax.set_xlabel(f"Position on the 16S rRNA gene ({ECOLI} numbering)")
fig.subplots_adjust(left=0.42, right=0.98, bottom=0.2, top=0.97)
fig.savefig(OUT / "Fig1_study_design_v4_window.png", dpi=300); fig.savefig(OUT / "Fig1_study_design_v4_window.pdf")
plt.close(fig)

# ------------------------------------------------------------------ Fig. 2
order = [("ISS root zone", "ISS root zone", C_ISS), ("ISS Arcillite", "ISS arcillite", C_ISS),
         ("ISS Water", "ISS water", C_ISS), ("ISS Wick", "ISS wick", C_ISS), ("ISS Foam", "ISS foam", C_ISS),
         ("ISS Swab", "ISS carrier swabs", C_ISS), ("OSD-65", "ISS interior surfaces", C_SURF),
         ("Lettuce Flight", "Lettuce roots, flight", C_LET), ("Lettuce Ground", "Lettuce roots, ground", C_LET),
         ("Terrestrial field", "Field pepper rhizosphere", C_TER)]
floor = 3e-7
rng = np.random.default_rng(7)
fig, ax = plt.subplots(figsize=(7.0, 3.6))
for i, (u, lab, c) in enumerate(order):
    if u == "OSD-65":
        ra, cra = osd65_ra, osd65_can
    else:
        sel = groups.unit == u
        if u == "Lettuce Ground":
            sel &= groups.reads >= 100000
        ra, cra = groups.loc[sel, "RA"], groups.loc[sel, "canRA"]
    k, n, kc = int((ra > 0).sum()), len(ra), int((cra > 0).sum())
    x = i + rng.uniform(-0.2, 0.2, n)
    y = ra.clip(lower=floor).to_numpy()
    ax.scatter(x, y, s=16, color=c, alpha=0.85, edgecolor="white", linewidth=0.3, zorder=3)
    hit = (cra > 0).to_numpy()
    ax.scatter(x[hit], y[hit], s=34, facecolor="none", edgecolor="black", linewidth=0.7, zorder=4)
    ax.text(i, 1.13, f"{k}/{n}", ha="center", fontsize=7, fontweight="bold" if k == 0 else "normal", transform=ax.get_xaxis_transform())
    ax.text(i, 1.04, f"{kc}/{n}", ha="center", fontsize=7, fontweight="bold" if kc == 0 else "normal", transform=ax.get_xaxis_transform())
ax.set_yscale("log"); ax.set_ylim(1.5e-7, 0.2)
ax.axhline(floor, color="#bbbbbb", lw=0.6, ls=":")
ax.text(len(order) - 0.45, floor * 1.6, "no nitrifier read (plotted at an arbitrary floor)", ha="right", fontsize=7, color="#888888")
ax.set_xticks(range(len(order)))
ax.set_xticklabels([o[1] for o in order], fontsize=7, rotation=35, ha="right", rotation_mode="anchor")
ax.set_ylabel("Relative abundance\n(nitrifier lineages)")
ax.text(-0.6, 1.13, "positive, nitrifier lineages:", ha="right", fontsize=7, transform=ax.get_xaxis_transform())
ax.text(-0.6, 1.04, "positive, canonical genera:", ha="right", fontsize=7, transform=ax.get_xaxis_transform())
fig.subplots_adjust(left=0.27, right=0.99, bottom=0.3, top=0.86)
fig.savefig(OUT / "Fig2_nitrifiers_by_compartment.png", dpi=300); fig.savefig(OUT / "Fig2_nitrifiers_by_compartment.pdf")
plt.close(fig)

# ------------------------------------------------------------------ Fig. 3
from matplotlib.lines import Line2D
fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0), gridspec_kw={"width_ratios": [1, 1.45]})
ax = axes[0]
for j, q in enumerate(["Q1", "Q2", "Q3", "Q4"]):
    pts = []
    for u, mk, dx in [("ISS Arcillite", "s", -0.12), ("ISS root zone", "o", 0.12)]:
        idx = groups.index[(groups.unit == u) & (groups.quadrant == q)]
        x_, n_ = nit[idx].sum(), groups.loc[idx, "reads"].sum()
        if x_ > 0:
            ax.scatter(j + dx, x_ / n_, marker=mk, s=28, color=C_ISS, zorder=3)
            pts.append(x_ / n_)
        else:
            up = 3.689 / n_
            ax.scatter(j + dx, up, marker="v", s=30, facecolor="white", edgecolor=C_ISS, linewidth=0.9, zorder=3)
            pts.append(up)
    ax.plot([j - 0.12, j + 0.12], pts, color=C_ISS, lw=0.8, alpha=0.6)
ax.set_yscale("log"); ax.set_xticks(range(4)); ax.set_xticklabels(["Q1", "Q2", "Q3", "Q4"])
ax.set_xlim(-0.5, 3.5); ax.set_ylim(2e-6, 2e-2); ax.set_ylabel("Pooled relative abundance of\nnitrifier lineages per quadrant")
ax.legend(handles=[Line2D([], [], marker="s", ls="", color=C_ISS, label="arcillite"),
                   Line2D([], [], marker="o", ls="", color=C_ISS, label="root zone"),
                   Line2D([], [], marker="v", ls="", markerfacecolor="white", markeredgecolor=C_ISS, label="not detected (95% bound)")],
          fontsize=7, frameon=False, loc="upper center", ncol=1, handletextpad=0.3, borderaxespad=0.2)
ax.set_title("(a) Per growth quadrant", fontsize=8, loc="left")
ax = axes[1]
lin = tax[nit_mask].str.replace(r"^.*?f__([^;]+)(?:; g__([^;]+))?.*$", lambda m: (m.group(2) or m.group(1)), regex=True)
lin = lin.replace({"Candidatus_Nitrocosmicus": "Ca. Nitrosocosmicus"})
cols = {"Terrestrial field": "Field pepper\nrhizosphere", "ISS Water": "ISS\nwater", "ISS Arcillite": "ISS\narcillite", "ISS root zone": "ISS\nroot zone"}
comp = pd.DataFrame({lab: ft.loc[nit_mask, groups.index[groups.unit == u]].sum(axis=1).groupby(lin).sum() for u, lab in cols.items()}).fillna(0)
share = comp / comp.sum()
top = share.max(axis=1).sort_values(ascending=False).index[:9]
share = pd.concat([share.loc[top], share.drop(top).sum().to_frame("other lineages").T])
COLORS = {"DSSD61": "#8c564b", "Ellin6067": "#9edae5", "Nitrospira": "#ff7f0e", "MND1": "#ffbb78", "966-1": "#2ca02c",
          "Ca. Nitrosocosmicus": "#98df8a", "mle1-7": "#9467bd", "Nitrolancea": "#c5b0d5", "Nitrosospira": "#e377c2",
          "other lineages": "#c7c7c7"}
ITALIC = {"Nitrospira": r"$\it{Nitrospira}$", "Nitrosospira": r"$\it{Nitrosospira}$", "Nitrolancea": r"$\it{Nitrolancea}$",
          "Ca. Nitrosocosmicus": r"$\it{Ca.}$ Nitrosocosmicus"}
bottom = np.zeros(share.shape[1])
for lname, row in share.iterrows():
    ax.bar(range(share.shape[1]), row.values, bottom=bottom, color=COLORS.get(lname, "#dddddd"), label=ITALIC.get(lname, lname), width=0.65)
    bottom += row.values
for i, c in enumerate(share.columns):
    ax.text(i, 1.02, f"{int(comp[c].sum()):,} reads", ha="center", fontsize=7)
ax.set_xticks(range(share.shape[1])); ax.set_xticklabels(share.columns, fontsize=7)
ax.set_ylabel("Share of nitrifier-lineage reads"); ax.set_ylim(0, 1.1)
ax.legend(fontsize=7, frameon=False, bbox_to_anchor=(1.0, 1.0), loc="upper left")
ax.set_title("(b) Lineage composition", fontsize=8, loc="left")
fig.tight_layout()
fig.savefig(OUT / "Fig3_within_mission.png", dpi=300); fig.savefig(OUT / "Fig3_within_mission.pdf")
plt.close(fig)

# ------------------------------------------------------------------ Fig. 4
rz = groups.loc[groups.unit == "ISS root zone"]
ter = groups.loc[groups.unit == "Terrestrial field"]
ra_grid = np.logspace(-7, -1, 300)
fig, ax = plt.subplots(figsize=(5.2, 3.3))
for depth, lab, ls in [(rz.reads.min(), f"shallowest ISS root-zone sample ({int(rz.reads.min()):,} reads)", "-"),
                       (rz.reads.median(), f"median ISS root-zone sample ({int(rz.reads.median()):,} reads)", "--")]:
    ax.plot(ra_grid, 1 - (1 - ra_grid) ** depth, ls, color=C_ISS, lw=1.2, label=lab)
pooled_upper = 3.689 / rz.reads.sum()
ax.axvspan(1e-7, pooled_upper, color="#999999", alpha=0.15)
ax.text(1.15e-7, 0.97, "pooled 95% bound,\ncanonical genera,\n20 samples together", fontsize=7, color="#555555", va="top")
obs = nit[rz.index].sum() / rz.reads.sum()
ax.axvline(obs, color=C_ISS, lw=0.8, ls=":")
ax.text(obs * 1.12, 0.42, "ISS root zone, nitrifier lineages (pooled)", fontsize=7, color=C_ISS, rotation=90, va="center")
for val, lab, ls, dx, ha in [(ter.canRA.mean(), "field mean, canonical genera", "--", 0.86, "right"),
                             (ter.RA.mean(), "field mean, nitrifier lineages", "-", 1.14, "left")]:
    ax.axvline(val, color=C_TER, lw=1, ls=ls)
    ax.text(val * dx, 0.42, lab, ha=ha, va="center", fontsize=7, color=C_TER, rotation=90)
ax.set_xscale("log"); ax.set_xlabel("Assumed relative abundance in a sample")
ax.set_ylabel("P(≥1 read in that sample)")
ax.legend(fontsize=7, frameon=False, loc="lower right")
fig.tight_layout()
fig.savefig(OUT / "Fig4_detection_power.png", dpi=300); fig.savefig(OUT / "Fig4_detection_power.pdf")
plt.close(fig)

# ------------------------------------------------------------------ tables
track = pd.read_csv(WORK / "dada2" / "track.tsv", sep="\t")
t1 = []
for study, sub in groups.groupby("study"):
    tr = track[track.study == study]
    asvs = int((ft[sub.index] > 0).any(axis=1).sum())
    for (grp, cmpt), s2 in sub.groupby(["group", "compartment"]):
        t1.append({"study": study, "group": grp, "compartment": cmpt, "n": len(s2),
                   "median reads (bacterial)": int(s2.reads.median()), "min reads": int(s2.reads.min()),
                   "study ASVs": asvs, "study read retention": f"{tr.nonchim.sum() / tr.input.sum():.1%}"})
pd.DataFrame(t1).to_csv(OUT / "table1_datasets.tsv", sep="\t", index=False)
print(pd.DataFrame(t1).to_string(index=False))
print("figures:", sorted(p.name for p in OUT.glob("*.png")))
