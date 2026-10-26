"""
64_v4_detection_stats.py — steps 3-4 of the V4 re-analysis

  A. Detection limits: per-sample depth, rule-of-three bounds, pooled abundance with exact Poisson CI,
     prevalence with Clopper-Pearson CI, repeated rarefaction (100x) at common depths.
  B. Cluster-level tests: ISS growth quadrants (Q1-Q4) vs terrestrial cultivar plots (XY6/21/39/40/42).
  C. Within-ISS compartments (same mission, no cross-study confound): substrate/water vs root zone,
     overall and paired by quadrant; identity of nitrifier ASVs shared between compartments.
  D. Veggie lettuce flight vs ground roots (same hardware family, different mission).

Inputs : $V4_WORK/dada2/feature-table.tsv, $V4_WORK/taxonomy/taxonomy.tsv, $V4_WORK/checkpoint/sample_groups.tsv
Outputs: $V4_WORK/stats/*.tsv and fig_nitrifier_by_compartment.png
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

WORK = Path(os.environ.get("V4_WORK", Path.home() / "pepper-v4" / "work"))
OUT = WORK / "stats"
OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(42)
N_ITER = 100

ft = pd.read_csv(WORK / "dada2" / "feature-table.tsv", sep="\t", index_col=0)
tax = pd.read_csv(WORK / "taxonomy" / "taxonomy.tsv", sep="\t", index_col=0)["Taxon"].reindex(ft.index).fillna("Unassigned")
keep = ~tax.str.contains("Mitochondria|Chloroplast")
ft, tax = ft.loc[keep], tax.loc[keep]
groups = pd.read_csv(WORK / "checkpoint" / "sample_groups.tsv", sep="\t", index_col=0)
unsan = groups.index.str.startswith("UnsanSeed")
groups.loc[unsan, ["group", "compartment", "cluster"]] = ["Ground seed", "seed (unsanitized)", "NA"]
groups["reads"] = ft.sum().reindex(groups.index)

ANY = (r"f__Nitrosomonadaceae|g__Nitrosococcus(?:;|$)|g__Nitrospira(?:;|$)|g__Nitrobacter(?:;|$)|g__Nitrococcus(?:;|$)|"
       r"g__Candidatus_Nitrotoga(?:;|$)|g__Nitrolancea(?:;|$)|c__Nitrososphaeria")
CANON = r"g__Nitrosomonas(?:;|$)|g__Nitrosospira(?:;|$)|g__Nitrospira(?:;|$)|g__Nitrobacter(?:;|$)"
sets = {"any nitrifier lineage": tax.str.contains(ANY, regex=True), "canonical genera": tax.str.contains(CANON, regex=True)}


def label(row):
    if row.group == "ISS pepper" and row.compartment in ("Plant Roots", "Rhizosphere"):
        return "ISS root zone"
    if row.group == "ISS pepper":
        return f"ISS {row.compartment}"
    return row.group


groups["unit"] = groups.apply(label, axis=1)
UNITS = ["ISS root zone", "ISS Arcillite", "ISS Water", "ISS Wick", "ISS Foam", "ISS Swab",
         "Terrestrial field", "Lettuce Flight", "Lettuce Ground"]
counts = {k: ft.loc[m].sum() for k, m in sets.items()}


def cp(k, n):
    lo = stats.beta.ppf(0.025, k, n - k + 1) if k > 0 else 0.0
    hi = stats.beta.ppf(0.975, k + 1, n - k) if k < n else 1.0
    return lo, hi


def poisson_ci(k):
    lo = stats.chi2.ppf(0.025, 2 * k) / 2 if k > 0 else 0.0
    return lo, stats.chi2.ppf(0.975, 2 * (k + 1)) / 2


# ---------------------------------------------------------------- A. detection limits
rows = []
for unit in UNITS:
    idx = groups.index[groups.unit == unit]
    depth = groups.loc[idx, "reads"]
    for setname, c in counts.items():
        k = int((c[idx] > 0).sum()); n = len(idx)
        tot = int(c[idx].sum()); pooled = int(depth.sum())
        plo, phi = poisson_ci(tot)
        lo, hi = cp(k, n)
        rows.append({"unit": unit, "set": setname, "n": n, "depth_min": int(depth.min()), "depth_median": int(depth.median()),
                     "positive": k, "prevalence": f"{k}/{n}", "prev_CI95": f"{lo:.1%}-{hi:.1%}",
                     "reads": tot, "pooled_depth": pooled, "pooled_RA": tot / pooled,
                     "pooled_RA_CI95": f"{plo / pooled:.2e}-{phi / pooled:.2e}",
                     "mean_sample_RA": float((c[idx] / depth).mean()),
                     "rule_of_three_RA_upper_per_sample_if_zero": 3 / depth.median()})
det = pd.DataFrame(rows)
det.to_csv(OUT / "detection_limits.tsv", sep="\t", index=False)


def rarefied_prevalence(depth):
    out = []
    for unit in UNITS:
        idx = [s for s in groups.index[groups.unit == unit] if groups.loc[s, "reads"] >= depth]
        if not idx:
            continue
        hits = {k: [] for k in sets}
        for _ in range(N_ITER):
            pos = {k: 0 for k in sets}
            for s in idx:
                sub = rng.multivariate_hypergeometric(ft[s].values.astype(np.int64), depth)
                for k, m in sets.items():
                    pos[k] += int(sub[m.values].sum() > 0)
            for k in sets:
                hits[k].append(pos[k])
        for k in sets:
            h = np.array(hits[k])
            out.append({"depth": depth, "unit": unit, "set": k, "n": len(idx),
                        "positive_median": float(np.median(h)), "positive_min": int(h.min()), "positive_max": int(h.max())})
    return out


rare = pd.DataFrame(rarefied_prevalence(20000) + rarefied_prevalence(1000))
rare.to_csv(OUT / "rarefied_prevalence.tsv", sep="\t", index=False)

# power: chance of >=1 nitrifier read at the ISS root-zone depth if abundance were a fraction of terrestrial
ter_ra = det.query("unit == 'Terrestrial field' and set == 'any nitrifier lineage'").mean_sample_RA.iloc[0]
rz_depth = groups.loc[groups.unit == "ISS root zone", "reads"]
power = pd.DataFrame([{"assumed_RA": ter_ra * f, "fraction_of_terrestrial": f,
                       "P(>=1 read) at min root-zone depth": 1 - (1 - ter_ra * f) ** rz_depth.min(),
                       "expected reads in pooled root zone": ter_ra * f * rz_depth.sum()} for f in (1, 0.1, 0.01, 0.001)])
power.to_csv(OUT / "detection_power.tsv", sep="\t", index=False)


# ---------------------------------------------------------------- B. cluster-level tests
def cluster_table(unit):
    idx = groups.index[groups.unit == unit]
    d = pd.DataFrame({"cluster": groups.loc[idx, "cluster"].astype(str), "depth": groups.loc[idx, "reads"]})
    for k, c in counts.items():
        d[k] = c[idx]
    g = d.groupby("cluster").sum()
    for k in counts:
        g[f"{k} RA"] = g[k] / g["depth"]
    return g


iss_cl = cluster_table("ISS root zone")
ter_cl = cluster_table("Terrestrial field")
crow = []
for k in counts:
    a_pos, b_pos = int((iss_cl[k] > 0).sum()), int((ter_cl[k] > 0).sum())
    f = stats.fisher_exact([[a_pos, len(iss_cl) - a_pos], [b_pos, len(ter_cl) - b_pos]])
    u = stats.mannwhitneyu(iss_cl[f"{k} RA"], ter_cl[f"{k} RA"], alternative="two-sided", method="exact")
    s_iss = counts[k][groups.index[groups.unit == "ISS root zone"]]
    s_ter = counts[k][groups.index[groups.unit == "Terrestrial field"]]
    fs = stats.fisher_exact([[int((s_iss > 0).sum()), int((s_iss == 0).sum())], [int((s_ter > 0).sum()), int((s_ter == 0).sum())]])
    crow.append({"set": k, "ISS quadrants positive": f"{a_pos}/{len(iss_cl)}", "terrestrial plots positive": f"{b_pos}/{len(ter_cl)}",
                 "Fisher p (clusters)": f.pvalue, "MWU p (cluster RA, exact)": u.pvalue, "Fisher p (samples, secondary)": fs.pvalue,
                 "ISS cluster RA": "; ".join(f"{i}:{v:.1e}" for i, v in iss_cl[f"{k} RA"].items()),
                 "terrestrial cluster RA": "; ".join(f"{i}:{v:.1e}" for i, v in ter_cl[f"{k} RA"].items())})
cl = pd.DataFrame(crow)
cl.to_csv(OUT / "cluster_tests.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- C. within-ISS compartments
k = "any nitrifier lineage"
comp_rows = []
r = counts[k][groups.index[groups.unit == "ISS root zone"]]
ra_r = r / groups.loc[r.index, "reads"]
for unit in ["ISS Arcillite", "ISS Water", "ISS Wick", "ISS Foam", "ISS Swab"]:
    a = counts[k][groups.index[groups.unit == unit]]
    fa = stats.fisher_exact([[int((a > 0).sum()), int((a == 0).sum())], [int((r > 0).sum()), int((r == 0).sum())]])
    ra_a = a / groups.loc[a.index, "reads"]
    u = stats.mannwhitneyu(ra_a, ra_r, alternative="two-sided")
    comp_rows.append({"compartment": unit, "positive": f"{int((a > 0).sum())}/{len(a)}", "root zone positive": f"{int((r > 0).sum())}/{len(r)}",
                      "Fisher p": fa.pvalue, "mean RA": ra_a.mean(), "root-zone mean RA": ra_r.mean(),
                      "RA ratio vs root zone": ra_a.mean() / ra_r.mean() if ra_r.mean() > 0 else np.inf, "MWU p": u.pvalue})
comp = pd.DataFrame(comp_rows)
comp.to_csv(OUT / "iss_compartments.tsv", sep="\t", index=False)

pq = []
for q in ["Q1", "Q2", "Q3", "Q4"]:
    def qra(unit):
        idx = groups.index[(groups.unit == unit) & (groups.cluster.astype(str) == q)]
        return counts[k][idx].sum() / groups.loc[idx, "reads"].sum(), len(idx)
    (ra_arc, n_arc), (ra_rz, n_rz) = qra("ISS Arcillite"), qra("ISS root zone")
    pq.append({"quadrant": q, "Arcillite RA": ra_arc, "n_arc": n_arc, "root zone RA": ra_rz, "n_rz": n_rz, "Arcillite > root zone": ra_arc > ra_rz})
pq = pd.DataFrame(pq)
sign_p = stats.binomtest(int(pq["Arcillite > root zone"].sum()), len(pq), 0.5, alternative="greater").pvalue
pq.to_csv(OUT / "iss_quadrant_pairs.tsv", sep="\t", index=False)

asv_rows = []
for asv in ft.index[sets[k]]:
    parts = tax[asv].split("; ")
    row = {"asv": asv[:10], "taxon": "; ".join(parts[-2:])}
    for unit in ["ISS root zone", "ISS Arcillite", "ISS Water", "ISS Foam", "Terrestrial field", "Lettuce Flight", "Lettuce Ground"]:
        row[unit] = int(ft.loc[asv, groups.index[groups.unit == unit]].sum())
    if any(v for kk, v in row.items() if kk.startswith(("ISS", "Lettuce"))):
        asv_rows.append(row)
asv_tab = pd.DataFrame(asv_rows).sort_values("ISS Arcillite", ascending=False)
asv_tab.to_csv(OUT / "iss_nitrifier_asvs.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- D. lettuce
lt = groups[groups.study.isin(["OSD-267", "OSD-269"])].copy()
for kk, c in counts.items():
    lt[kk] = c[lt.index]
lt["upper RA if zero (3/depth)"] = 3 / lt["reads"]
lt.to_csv(OUT / "lettuce_controls.tsv", sep="\t")

# ---------------------------------------------------------------- draft figure
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(9, 4.5))
floor = 1e-7
for i, unit in enumerate(UNITS):
    idx = groups.index[groups.unit == unit]
    ra = (counts[k][idx] / groups.loc[idx, "reads"]).clip(lower=floor)
    x = i + rng.uniform(-0.18, 0.18, len(ra))
    color = "#c0392b" if unit.startswith("ISS") else ("#2471a3" if unit.startswith("Terr") else "#7d7d7d")
    ax.scatter(x, ra, s=18, alpha=0.75, color=color)
ax.set_yscale("log")
ax.axhline(floor, color="#bbbbbb", lw=0.8, ls=":")
ax.text(len(UNITS) - 0.5, floor * 1.4, "not detected", ha="right", fontsize=8, color="#888888")
ax.set_xticks(range(len(UNITS)))
ax.set_xticklabels([u.replace("ISS ", "ISS\n").replace("Lettuce ", "Lettuce\n").replace(" field", "\nfield") for u in UNITS], fontsize=8)
ax.set_ylabel("Nitrifier-lineage relative abundance (per sample)")
ax.set_title("DRAFT - nitrifier lineages by compartment (V4-harmonised, organelles removed)", fontsize=9)
fig.tight_layout()
fig.savefig(OUT / "fig_nitrifier_by_compartment.png", dpi=200)

with pd.option_context("display.width", 220, "display.max_columns", 30, "display.max_colwidth", 70):
    print("== A. detection limits\n", det[["unit", "set", "n", "depth_min", "prevalence", "prev_CI95", "reads", "pooled_RA", "pooled_RA_CI95", "mean_sample_RA"]].to_string(index=False))
    print("\n== A. rarefied prevalence (median [min-max] positive samples over 100 draws)")
    for _, rr in rare.iterrows():
        print(f"   depth {rr.depth:>6} {rr.unit:18s} {rr.set:22s} {rr.positive_median:>5.1f} [{rr.positive_min}-{rr.positive_max}] / {rr.n}")
    print("\n== A. detection power at ISS root-zone depth\n", power.to_string(index=False))
    print("\n== B. cluster-level tests\n", cl.drop(columns=["ISS cluster RA", "terrestrial cluster RA"]).to_string(index=False))
    print("   ISS quadrants:", cl["ISS cluster RA"].iloc[0])
    print("   terrestrial plots:", cl["terrestrial cluster RA"].iloc[0])
    print("\n== C. within-ISS compartments vs root zone\n", comp.to_string(index=False))
    print(pq.to_string(index=False))
    print(f"   sign test (Arcillite > root zone in quadrants), one-sided p = {sign_p:.3f}")
    print("\n== C. nitrifier ASVs seen in ISS/lettuce samples\n", asv_tab.head(12).to_string(index=False))
    print("\n== D. lettuce roots\n", lt[["group", "reads", "any nitrifier lineage", "canonical genera", "upper RA if zero (3/depth)"]].to_string())
