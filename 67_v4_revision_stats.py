"""
67_v4_revision_stats.py — analyses added after internal review round 2, and every number used in manuscript v5

Adds to 64/65:
  - nitrifier definitions: canonical bacterial genera resolvable in this V4 window (Nitrosomonas, Nitrosospira,
    Nitrospira; Nitrobacter is not resolvable from Bradyrhizobium/Afipia and is checked separately), a sensitivity
    set that adds cultured AOA genera, the broad set of 63/64, and the broad set without Group 1.1c and Nitrolancea
  - denominators restricted to Bacteria + Archaea (organelles, Eukaryota and Unassigned removed), as in 65
  - positivity thresholds of >=1, >=2 and >=5 reads
  - two-stage cluster bootstrap (quadrant/cultivar, then sample) for pooled relative abundance and fold differences
  - within-mission comparisons at quadrant level (paired), with sample-level values as descriptive only
  - power against the canonical field level, a cross-talk check, near-identical field variants, seed screen
  - Nitrobacter resolvability against the SILVA 138 V4-240 reference used to train the classifier
  - an illustrative nitrogen budget (PH-04 fertiliser, crew urine N, alkalinity, reactor size, free ammonia, start-up)
  - manuscript tables (Table 1, 2, 3, 5) and a per-sample supplementary table, all written from code

Inputs : $V4_WORK/{dada2,taxonomy,checkpoint,osd65}, $V4_REF/silva-138-99-{v4-240-seqs,tax}.qza
Outputs: $V4_WORK/revision/*.tsv, nitrifier_asvs.fasta, numbers.tsv (one row per number quoted in the text)
"""
import io
import math
import os
import re
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

WORK = Path(os.environ.get("V4_WORK", Path.home() / "pepper-v4" / "work"))
REF = Path(os.environ.get("V4_REF", Path.home() / "pepper-v4" / "ref"))
OUT = WORK / "revision"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20261021
NUM = []  # (key, value, note)


def note(key, value, text=""):
    NUM.append({"key": key, "value": value, "note": text})


# ------------------------------------------------------------------ data
ft = pd.read_csv(WORK / "dada2" / "feature-table.tsv", sep="\t", index_col=0)
tax = pd.read_csv(WORK / "taxonomy" / "taxonomy.tsv", sep="\t", index_col=0)["Taxon"].reindex(ft.index).fillna("Unassigned")
keep = tax.str.contains(r"^d__(?:Bacteria|Archaea)") & ~tax.str.contains("Mitochondria|Chloroplast")
note("asvs_before_filter", int(len(ft)))
ft, tax = ft.loc[keep], tax.loc[keep]
note("asvs_after_filter", int(len(ft)), "Bacteria + Archaea, organelles removed")
fasta, name = {}, None
for line in open(WORK / "dada2" / "asv.fasta"):
    line = line.strip()
    if line.startswith(">"):
        name = line[1:]
    else:
        fasta[name] = line

groups = pd.read_csv(WORK / "checkpoint" / "sample_groups.tsv", sep="\t", index_col=0)
groups.loc[groups.index.str.startswith("UnsanSeed"), ["group", "compartment"]] = ["Ground seed", "seed"]
groups["depth"] = ft.sum().reindex(groups.index)


def unit(row):
    if row.group == "ISS pepper":
        if row.compartment == "Plant Roots":
            return "ISS root sections"
        if row.compartment == "Rhizosphere":
            return "ISS root-substrate"
        return f"ISS {row.compartment}"
    return row.group


groups["unit"] = groups.apply(unit, axis=1)
groups["quadrant"] = np.where(groups.study == "OSD-772", groups.index.str.extract(r"^(Q\d)", expand=False), "")
groups.loc[groups.study == "PRJNA1145089", "quadrant"] = groups.loc[groups.study == "PRJNA1145089", "cluster"]

SETS = {
    "canonical": r"g__(?:Nitrosomonas|Nitrosospira|Nitrospira)(?:;|$)",
    "canonical+AOA": r"g__(?:Nitrosomonas|Nitrosospira|Nitrospira|Candidatus_Nitrocosmicus|Candidatus_Nitrosotalea|"
                     r"Nitrososphaera|Candidatus_Nitrososphaera|Nitrosopumilus|Candidatus_Nitrosopumilus|"
                     r"Candidatus_Nitrosoarchaeum|Nitrosotenuis|Candidatus_Nitrosotenuis|Candidatus_Nitrosocaldus)(?:;|$)",
    "broad": r"f__Nitrosomonadaceae|g__Nitrosococcus(?:;|$)|g__Nitrospira(?:;|$)|g__Nitrobacter(?:;|$)|g__Nitrococcus(?:;|$)|"
             r"g__Candidatus_Nitrotoga(?:;|$)|g__Nitrolancea(?:;|$)|c__Nitrososphaeria",
}
masks = {k: tax.str.contains(v, regex=True) for k, v in SETS.items()}
masks["broad-refined"] = masks["broad"] & ~tax.str.contains(r"f__Group_1\.1c|g__Nitrolancea")
masks["Nitrobacter"] = tax.str.contains(r"g__Nitrobacter(?:;|$)")
counts = {k: ft.loc[m].sum() for k, m in masks.items()}
for k in masks:
    groups[k] = counts[k].reindex(groups.index)

UNITS = {
    "ISS root zone": groups.unit.isin(["ISS root sections", "ISS root-substrate"]),
    "ISS root sections": groups.unit == "ISS root sections",
    "ISS root-substrate": groups.unit == "ISS root-substrate",
    "ISS arcillite": groups.unit == "ISS Arcillite",
    "ISS water": groups.unit == "ISS Water",
    "ISS wick": groups.unit == "ISS Wick",
    "ISS foam": groups.unit == "ISS Foam",
    "ISS carrier swabs": groups.unit == "ISS Swab",
    "ISS carrier, all 75": groups.unit.isin(["ISS root sections", "ISS root-substrate", "ISS Arcillite", "ISS Water",
                                              "ISS Wick", "ISS Foam", "ISS Swab"]),
    "Lettuce roots, flight": groups.group == "Lettuce Flight",
    "Lettuce roots, ground": groups.group == "Lettuce Ground",
    "Lettuce roots, ground (>=100,000 reads)": (groups.group == "Lettuce Ground") & (groups.depth >= 100000),
    "Field pepper rhizosphere": groups.study == "PRJNA1145089",
}


def cp(k, n):
    lo = stats.beta.ppf(0.025, k, n - k + 1) if k > 0 else 0.0
    hi = stats.beta.ppf(0.975, k + 1, n - k) if k < n else 1.0
    return lo, hi


def pois(x, n):
    lo = stats.chi2.ppf(0.025, 2 * x) / 2 / n if x > 0 else 0.0
    hi = stats.chi2.ppf(0.975, 2 * (x + 1)) / 2 / n
    return lo, hi


# ------------------------------------------------------------------ detection by unit and definition
rows = []
for uname, sel in UNITS.items():
    sub = groups[sel]
    for d in ["canonical", "canonical+AOA", "broad", "broad-refined", "Nitrobacter"]:
        x, n = int(sub[d].sum()), int(sub.depth.sum())
        k1, k2, k5 = int((sub[d] >= 1).sum()), int((sub[d] >= 2).sum()), int((sub[d] >= 5).sum())
        lo, hi = cp(k1, len(sub))
        plo, phi = pois(x, n)
        rows.append({"unit": uname, "definition": d, "n": len(sub), "depth_median": float(sub.depth.median()),
                     "depth_min": int(sub.depth.min()), "depth_max": int(sub.depth.max()), "pooled_reads": n,
                     "pos1": k1, "pos2": k2, "pos5": k5, "prev_lo": lo, "prev_hi": hi, "nit_reads": x,
                     "pooled_RA": x / n, "RA_lo": plo, "RA_hi": phi,
                     "mean_sample_RA": float((sub[d] / sub.depth).mean()),
                     "per_sample_bound_min": 3 / sub.depth.max(), "per_sample_bound_max": 3 / sub.depth.min()})
det = pd.DataFrame(rows)

o65 = pd.read_csv(WORK / "osd65" / "nitrifier_by_library.tsv", sep="\t", index_col=0)
s65 = o65[o65.kind == "surface"]
for d, col in [("canonical", "canonical genera"), ("broad", "any nitrifier lineage")]:
    x, n = int(s65[col].sum()), int(s65.bacterial_reads.sum())
    lo, hi = cp(int((s65[col] > 0).sum()), len(s65))
    plo, phi = pois(x, n)
    det.loc[len(det)] = {"unit": "ISS interior surfaces (OSD-65)", "definition": d, "n": len(s65),
                         "depth_median": float(s65.bacterial_reads.median()), "depth_min": int(s65.bacterial_reads.min()),
                         "depth_max": int(s65.bacterial_reads.max()), "pooled_reads": n,
                         "pos1": int((s65[col] >= 1).sum()), "pos2": int((s65[col] >= 2).sum()), "pos5": int((s65[col] >= 5).sum()),
                         "prev_lo": lo, "prev_hi": hi, "nit_reads": x, "pooled_RA": x / n, "RA_lo": plo, "RA_hi": phi,
                         "mean_sample_RA": float((s65[col] / s65.bacterial_reads).mean()),
                         "per_sample_bound_min": 3 / s65.bacterial_reads.max(), "per_sample_bound_max": 3 / s65.bacterial_reads.min()}
det.to_csv(OUT / "detection_by_definition.tsv", sep="\t", index=False)


def D(u, d):
    return det[(det.unit == u) & (det.definition == d)].iloc[0]


rz, fd, arc = D("ISS root zone", "canonical"), D("Field pepper rhizosphere", "canonical"), D("ISS arcillite", "broad")
note("rz_n", int(rz.n)); note("rz_pooled_reads", int(rz.pooled_reads))
note("rz_depth_min", int(rz.depth_min)); note("rz_depth_median", rz.depth_median); note("rz_depth_max", int(rz.depth_max))
note("rz_canonical_upper", rz.RA_hi, "3.69/N, exact two-sided 95% Poisson limit for zero")
note("rz_canonical_prev_hi", rz.prev_hi)
note("rz_per_sample_bound_range", f"{rz.per_sample_bound_min:.2e}-{rz.per_sample_bound_max:.2e}", "3/N per sample")
note("field_canonical_pooled_RA", fd.pooled_RA); note("field_canonical_RA_CI", f"{fd.RA_lo:.3e}-{fd.RA_hi:.3e}")
note("field_canonical_reads", int(fd.nit_reads)); note("field_canonical_mean_RA", fd.mean_sample_RA)
note("fold_field_vs_rz_canonical_upper", fd.pooled_RA / rz.RA_hi)
carrier = D("ISS carrier, all 75", "canonical")
note("carrier_n", int(carrier.n)); note("carrier_canonical_pos", int(carrier.pos1)); note("carrier_pooled_reads", int(carrier.pooled_reads))
note("carrier_canonical_upper", carrier.RA_hi)
for u in ["ISS root zone", "ISS arcillite", "ISS water", "Field pepper rhizosphere"]:
    b = D(u, "broad")
    note(f"{u}|broad", f"{b.pos1}/{b.n}; >=2: {b.pos2}; >=5: {b.pos5}; reads {b.nit_reads}/{b.pooled_reads}; RA {b.pooled_RA:.3e} ({b.RA_lo:.2e}-{b.RA_hi:.2e}); mean {b.mean_sample_RA:.3e}")
aoa = D("ISS root zone", "canonical+AOA")
note("rz_canonical+AOA", f"{aoa.pos1}/{aoa.n} ({aoa.prev_lo:.1%}-{aoa.prev_hi:.1%}); reads {aoa.nit_reads}; RA {aoa.pooled_RA:.2e} ({aoa.RA_lo:.1e}-{aoa.RA_hi:.1e})")
note("Nitrobacter_reads_all_samples", int(ft.loc[masks["Nitrobacter"]].sum().sum()))

# ------------------------------------------------------------------ cluster-level tests between environments
rz_s = groups[UNITS["ISS root zone"]]
fd_s = groups[UNITS["Field pepper rhizosphere"]]
arc_s = groups[UNITS["ISS arcillite"]]


def cluster_ra(sub, d):
    g = sub.groupby("quadrant")
    return g[d].sum() / g.depth.sum()


tests = []
for d in ["canonical", "canonical+AOA", "broad", "broad-refined"]:
    a, b = cluster_ra(rz_s, d), cluster_ra(fd_s, d)
    ka, kb = int((a > 0).sum()), int((b > 0).sum())
    fp = stats.fisher_exact([[ka, len(a) - ka], [kb, len(b) - kb]])[1]
    mw = stats.mannwhitneyu(a, b, alternative="two-sided", method="exact").pvalue
    tests.append({"comparison": "ISS root-zone quadrants vs field cultivars", "definition": d,
                  "positive": f"{ka}/{len(a)} vs {kb}/{len(b)}", "Fisher_p": fp, "MWU_exact_p": mw,
                  "ISS_cluster_RA": "; ".join(f"{q}:{v:.1e}" for q, v in a.items()),
                  "field_cluster_RA": "; ".join(f"{q}:{v:.1e}" for q, v in b.items())})
a, b = cluster_ra(arc_s, "broad"), cluster_ra(fd_s, "broad")
tests.append({"comparison": "ISS arcillite quadrants vs field cultivars (substrate vs soil)", "definition": "broad",
              "positive": f"{int((a > 0).sum())}/4 vs {int((b > 0).sum())}/5",
              "Fisher_p": stats.fisher_exact([[int((a > 0).sum()), 4 - int((a > 0).sum())], [int((b > 0).sum()), 5 - int((b > 0).sum())]])[1],
              "MWU_exact_p": stats.mannwhitneyu(a, b, alternative="two-sided", method="exact").pvalue,
              "ISS_cluster_RA": "; ".join(f"{q}:{v:.1e}" for q, v in a.items()),
              "field_cluster_RA": "; ".join(f"{q}:{v:.1e}" for q, v in b.items())})
tests = pd.DataFrame(tests)
tests.to_csv(OUT / "cluster_tests_v5.tsv", sep="\t", index=False)
for _, r in tests.iterrows():
    note(f"test|{r.comparison}|{r.definition}", f"{r.positive}; Fisher {r.Fisher_p:.4f}; MWU {r.MWU_exact_p:.4f}")

# ------------------------------------------------------------------ two-stage cluster bootstrap
rng = np.random.default_rng(SEED)
B = 10000


def boot_pooled(sub, d, rng, B):
    clusters = [(g[d].to_numpy(float), g.depth.to_numpy(float)) for _, g in sub.groupby("quadrant")]
    out = np.empty(B)
    for i in range(B):
        x = n = 0.0
        for c in rng.integers(0, len(clusters), len(clusters)):
            cx, cn = clusters[c]
            j = rng.integers(0, len(cx), len(cx))
            x += cx[j].sum(); n += cn[j].sum()
        out[i] = x / n
    return out


boot_rows = []
for d in ["broad", "canonical"]:
    bi, bf = boot_pooled(rz_s, d, rng, B), boot_pooled(fd_s, d, rng, B)
    ratio = np.where(bi > 0, bf / np.where(bi > 0, bi, 1), np.inf)
    boot_rows.append({"definition": d, "ISS_root_zone_pooled_RA_2.5": np.percentile(bi, 2.5),
                      "ISS_root_zone_pooled_RA_97.5": np.percentile(bi, 97.5),
                      "field_pooled_RA_2.5": np.percentile(bf, 2.5), "field_pooled_RA_97.5": np.percentile(bf, 97.5),
                      "fold_field_over_ISS_2.5": np.percentile(ratio, 2.5), "share_boot_ISS_zero": float((bi == 0).mean())})
    if d == "canonical":  # zero everywhere: bound the ISS side by the per-quadrant 3.69/N limits
        q = rz_s.groupby("quadrant").depth.sum()
        boot_rows[-1]["per_quadrant_upper_min"] = float((3.689 / q).min())
        boot_rows[-1]["per_quadrant_upper_max"] = float((3.689 / q).max())
boot = pd.DataFrame(boot_rows)
boot.to_csv(OUT / "cluster_bootstrap.tsv", sep="\t", index=False)
bb = boot.set_index("definition")
note("boot_rz_broad_CI", f"{bb.loc['broad', 'ISS_root_zone_pooled_RA_2.5']:.2e}-{bb.loc['broad', 'ISS_root_zone_pooled_RA_97.5']:.2e}")
note("boot_fold_broad_lower", bb.loc["broad", "fold_field_over_ISS_2.5"], "2.5th percentile of field/ISS pooled RA ratio")
note("rz_canonical_per_quadrant_upper", f"{bb.loc['canonical', 'per_quadrant_upper_min']:.2e}-{bb.loc['canonical', 'per_quadrant_upper_max']:.2e}")
note("fold_field_canonical_vs_quadrant_upper", f"{fd.pooled_RA / bb.loc['canonical', 'per_quadrant_upper_max']:.0f}-{fd.pooled_RA / bb.loc['canonical', 'per_quadrant_upper_min']:.0f}")
b_rz, b_fd = D("ISS root zone", "broad"), D("Field pepper rhizosphere", "broad")
note("fold_field_vs_rz_broad_point", b_fd.pooled_RA / b_rz.pooled_RA)
note("fold_field_vs_rz_broad_refined_point", D("Field pepper rhizosphere", "broad-refined").pooled_RA / D("ISS root zone", "broad-refined").pooled_RA)
note("fold_field_vs_arcillite_broad", b_fd.pooled_RA / arc.pooled_RA, "matched substrate vs soil")
rs = D("ISS root-substrate", "broad")
note("rootsubstrate_broad", f"{rs.pos1}/{rs.n}; reads {rs.nit_reads}/{rs.pooled_reads}; RA {rs.pooled_RA:.2e}; fold field {b_fd.pooled_RA / rs.pooled_RA:.0f}")
rsec = D("ISS root sections", "broad")
note("rootsections_broad", f"{rsec.pos1}/{rsec.n}; reads {rsec.nit_reads}/{rsec.pooled_reads}; RA {rsec.pooled_RA:.2e}")

# ------------------------------------------------------------------ within mission, quadrant level
wq = []
rz_q = rz_s.groupby("quadrant")
for comp, cname in [("ISS arcillite", "Arcillite"), ("ISS water", "Water"), ("ISS wick", "Wick"), ("ISS foam", "Foam"), ("ISS carrier swabs", "Swab")]:
    cs = groups[UNITS[comp]]
    cq = cs.groupby("quadrant")
    pairs = []
    for q in sorted(cs.quadrant.unique()):
        c_ra = cq.get_group(q)["broad"].sum() / cq.get_group(q).depth.sum()
        r_ra = rz_q.get_group(q)["broad"].sum() / rz_q.get_group(q).depth.sum()
        pairs.append((q, c_ra, r_ra))
    pq = pd.DataFrame(pairs, columns=["quadrant", "comp_RA", "rz_RA"])
    kc, kr = int((pq.comp_RA > 0).sum()), int((pq.rz_RA > 0).sum())
    diff = pq.comp_RA - pq.rz_RA
    npos, nneg = int((diff > 0).sum()), int((diff < 0).sum())
    sign_one = stats.binomtest(npos, npos + nneg, alternative="greater").pvalue if npos + nneg else np.nan
    sign_two = stats.binomtest(npos, npos + nneg, alternative="two-sided").pvalue if npos + nneg else np.nan
    pooled_ratio = (cs["broad"].sum() / cs.depth.sum()) / (rz_s.loc[rz_s.quadrant.isin(pq.quadrant), "broad"].sum() /
                                                            rz_s.loc[rz_s.quadrant.isin(pq.quadrant)].depth.sum())
    # sample level, descriptive only
    x = cs["broad"] / cs.depth; y = rz_s["broad"] / rz_s.depth
    ks, kz = int((cs["broad"] > 0).sum()), int((rz_s["broad"] > 0).sum())
    perm = stats.permutation_test((x.values, y.values), lambda a, b: stats.mannwhitneyu(a, b).statistic,
                                  permutation_type="independent", n_resamples=20000, alternative="two-sided",
                                  random_state=SEED).pvalue
    wq.append({"comparison": f"{cname} vs root zone", "quadrants": ",".join(pq.quadrant),
               "positive_quadrants": f"{kc}/{len(pq)} vs {kr}/{len(pq)}",
               "Fisher_quadrants_p": stats.fisher_exact([[kc, len(pq) - kc], [kr, len(pq) - kr]])[1],
               "pairs_comp_higher": f"{npos}/{npos + nneg}", "sign_one_sided_p": sign_one, "sign_two_sided_p": sign_two,
               "pooled_ratio": pooled_ratio,
               "quadrant_ratios": "; ".join(f"{q}:{(c / r if r > 0 else np.inf):.1f}" for q, c, r in pairs),
               "quadrant_comp_RA": "; ".join(f"{q}:{c:.1e}" for q, c, r in pairs),
               "quadrant_rz_RA": "; ".join(f"{q}:{r:.1e}" for q, c, r in pairs),
               "samples_positive (descriptive)": f"{ks}/{len(cs)} vs {kz}/{len(rz_s)}",
               "sample_Fisher_p (descriptive)": stats.fisher_exact([[ks, len(cs) - ks], [kz, len(rz_s) - kz]])[1],
               "sample_MWU_perm_p (descriptive)": perm})
wq = pd.DataFrame(wq)
# quadrant bootstrap of the arcillite/root-zone pooled ratio
qs = ["Q1", "Q2", "Q3", "Q4"]
ar_q = arc_s.groupby("quadrant")
rr = []
for _ in range(B):
    pick = rng.choice(qs, 4)
    xa = sum(ar_q.get_group(q)["broad"].sum() for q in pick); na = sum(ar_q.get_group(q).depth.sum() for q in pick)
    xr = sum(rz_q.get_group(q)["broad"].sum() for q in pick); nr = sum(rz_q.get_group(q).depth.sum() for q in pick)
    rr.append((xa / na) / (xr / nr) if xr > 0 else np.inf)
rr = np.array(rr)
wq.loc[wq.comparison == "Arcillite vs root zone", "pooled_ratio_boot_2.5"] = np.percentile(rr, 2.5)
wq.loc[wq.comparison == "Arcillite vs root zone", "boot_share_rz_zero"] = float(np.isinf(rr).mean())
wq.to_csv(OUT / "within_mission_quadrants.tsv", sep="\t", index=False)
for _, r in wq.iterrows():
    note(f"within|{r.comparison}", f"quadrants {r.positive_quadrants} (Fisher {r.Fisher_quadrants_p:.3f}); pairs {r.pairs_comp_higher} "
         f"(sign one {r.sign_one_sided_p:.4f}, two {r.sign_two_sided_p:.4f}); pooled ratio {r.pooled_ratio:.1f}; {r.quadrant_ratios}; "
         f"samples {r['samples_positive (descriptive)']} (Fisher {r['sample_Fisher_p (descriptive)']:.2e}, perm MWU {r['sample_MWU_perm_p (descriptive)']:.2e})")
note("arc_rz_ratio_boot_lower", np.percentile(rr, 2.5))

# ------------------------------------------------------------------ power against canonical and broad field levels
pw = []
for d in ["canonical", "broad"]:
    fmean = D("Field pepper rhizosphere", d).mean_sample_RA
    for frac in [1, 0.1, 0.01, 0.001]:
        p = fmean * frac
        pw.append({"definition": d, "fraction_of_field_mean": frac, "assumed_RA": p,
                   "P_shallowest": 1 - (1 - p) ** rz.depth_min, "P_median": 1 - (1 - p) ** rz.depth_median,
                   "expected_pooled_reads": p * rz.pooled_reads,
                   "P_le_observed_pooled": stats.poisson.cdf(int(D("ISS root zone", d).nit_reads), p * rz.pooled_reads)})
pw = pd.DataFrame(pw)
pw.to_csv(OUT / "power_v5.tsv", sep="\t", index=False)
for _, r in pw.iterrows():
    if r.fraction_of_field_mean in (0.01, 0.001):
        note(f"power|{r.definition}|{r.fraction_of_field_mean}", f"RA {r.assumed_RA:.2e}; P shallow {r.P_shallowest:.3f}; P median {r.P_median:.3f}; expected {r.expected_pooled_reads:.0f}")
for f in [1e-4, 1e-5]:
    note(f"reads_for_95pct_detection|{f}", math.log(0.05) / math.log(1 - f))
note("reads_needed_ge3_95pct_factor", stats.gamma.ppf(0.95, 3), "N*f such that P(X>=3)=0.95 under Poisson")

# ------------------------------------------------------------------ cross-talk check and ASV-level origin
dssd = [a for a in ft.index if a.startswith("ae1fae98aa")][0]
nc = [a for a in ft.index if a.startswith("8dac088009")][0]
rate = ft.loc[dssd, rz_s.index].sum() / rz.pooled_reads
for u in ["ISS wick", "ISS carrier swabs"]:
    cs = groups[UNITS[u]]
    lam = rate * cs.depth.sum()
    note(f"crosstalk|{u}", f"observed {int(ft.loc[dssd, cs.index].sum())}; expected {lam:.1f} at root-zone rate; P(0) {math.exp(-lam):.1e}")


def mism(a, b):
    return sum(x != y for x, y in zip(a, b))


field_asvs = ft.index[ft[fd_s.index].sum(axis=1) > 0]
iss_nit = ft.index[masks["broad"] & (ft[groups.index[UNITS["ISS carrier, all 75"]]].sum(axis=1) > 0)]
near = []
for a in iss_nit:
    best = min(field_asvs, key=lambda b: mism(fasta[a], fasta[b]))
    near.append({"ISS_asv": a[:10], "ISS_taxon": tax[a].split(";")[-2].strip(), "ISS_carrier_reads": int(ft.loc[a, groups.index[UNITS["ISS carrier, all 75"]]].sum()),
                 "nearest_field_asv": best[:10], "mismatches": mism(fasta[a], fasta[best]), "field_reads": int(ft.loc[best, fd_s.index].sum()),
                 "field_taxon": tax[best].split(";")[-2].strip()})
pd.DataFrame(near).to_csv(OUT / "iss_nitrifier_nearest_field.tsv", sep="\t", index=False)
for r in near:
    note(f"nearest_field|{r['ISS_asv']}", f"{r['ISS_taxon']}; {r['mismatches']} mismatches; field reads {r['field_reads']}")
seeds = groups[groups.group.isin(["Ground seed"]) | (groups.compartment == "seed")]
sx = ft.loc[masks["broad"], seeds.index]
for s in seeds.index[sx.sum() > 0]:
    hit = sx[s][sx[s] > 0]
    note(f"seed|{s}", "; ".join(f"{a[:10]} {tax[a].split(';')[-2].strip()} {int(v)}" for a, v in hit.items()))
note("seed_libraries_screened", len(seeds))

# ------------------------------------------------------------------ Nitrobacter resolvability (SILVA 138 V4-240)
def read_qza(path, suffix):
    with zipfile.ZipFile(path) as z:
        nm = [n for n in z.namelist() if n.endswith(suffix)][0]
        return z.read(nm).decode()


try:
    ref_txt = read_qza(REF / "silva-138-99-v4-240-seqs.qza", "dna-sequences.fasta")
    ref_tax = pd.read_csv(io.StringIO(read_qza(REF / "silva-138-99-tax.qza", "taxonomy.tsv")), sep="\t", index_col=0)["Taxon"]
    ref, nm = {}, None
    for line in ref_txt.splitlines():
        if line.startswith(">"):
            nm = line[1:].split()[0]
        else:
            ref[nm] = ref.get(nm, "") + line.strip()
    rgen = ref_tax.str.extract(r"g__([^;]+)")[0]
    seq2gen = {}
    for k, s in ref.items():
        seq2gen.setdefault(s, set()).add(rgen.get(k, "NA"))
    nb = {s for k, s in ref.items() if rgen.get(k) == "Nitrobacter"}
    shared = [s for s in nb if seq2gen[s] - {"Nitrobacter"}]
    note("Nitrobacter_ref_distinct_v4", len(nb)); note("Nitrobacter_ref_identical_to_other_genera", len(shared))
    hits = [a for a in ft.index if fasta[a] in nb]
    rzr = int(ft.loc[hits, rz_s.index].sum().sum()) if hits else 0
    note("ISS_asvs_identical_to_Nitrobacter_ref", f"{len(hits)} ASVs; root-zone reads {rzr} ({rzr / rz.pooled_reads:.1%}); classified as "
         + ", ".join(tax[a].split(';')[-1].strip() or tax[a].split(';')[-2].strip() for a in hits))
    canon_ref = {s for k, s in ref.items() if rgen.get(k) in ("Nitrosomonas", "Nitrosospira", "Nitrospira")}
    rz_asvs = ft.index[ft[rz_s.index].sum(axis=1) > 0]
    enc = lambda seqs: np.array([np.frombuffer(q[:240].ljust(240, "N").encode(), dtype=np.uint8) for q in seqs])
    A, R = enc([fasta[a] for a in rz_asvs]), enc(sorted(canon_ref))
    mins = np.array([(A[i] != R).sum(axis=1).min() for i in range(len(A))])
    note("rz_asvs_min_mismatch_to_canonical_refs", int(mins.min()), f"{len(canon_ref)} distinct canonical reference sequences; {len(rz_asvs)} root-zone ASVs")
except FileNotFoundError:
    note("Nitrobacter_check", "reference qza not found; skipped")

# ------------------------------------------------------------------ shared ASVs and community context
iss_all = groups.index[groups.group == "ISS pepper"]
let = groups.index[groups.study.isin(["OSD-267", "OSD-269"])]
pres = lambda cols: ft[cols].sum(axis=1) > 0
note("shared_ISSpepper_field", int((pres(iss_all) & pres(fd_s.index)).sum()), "all OSD-772 flight samples incl. plant tissues")
note("shared_ISSpepper_field_ISSread_share", ft.loc[pres(iss_all) & pres(fd_s.index), iss_all].sum().sum() / ft[iss_all].sum().sum())
note("shared_ISSpepper_field_fieldread_share", ft.loc[pres(iss_all) & pres(fd_s.index), fd_s.index].sum().sum() / ft[fd_s.index].sum().sum())
note("shared_ISSpepper_lettuce", int((pres(iss_all) & pres(let)).sum()))
note("shared_rootzone_field", int((pres(rz_s.index) & pres(fd_s.index)).sum()))
gen = tax.str.extract(r"g__([^;]+)")[0].fillna("unassigned")
gen = gen.where(~gen.str.startswith("uncultured"), "uncultured")
for u in ["ISS root zone", "ISS arcillite", "Field pepper rhizosphere"]:
    cols = groups.index[UNITS[u]]
    rich = (ft[cols] > 0).sum()
    rar = []
    rr_rng = np.random.default_rng(SEED)
    for c in cols:
        v = ft[c].values
        if v.sum() >= 20000:
            rar.append(int((rr_rng.multivariate_hypergeometric(v, 20000) > 0).sum()))
    pooled = ft[cols].sum(axis=1).groupby(gen).sum() / ft[cols].sum().sum()
    meanra = (ft[cols] / ft[cols].sum()).groupby(gen).sum().mean(axis=1)
    top = pooled.sort_values(ascending=False).head(6)
    note(f"richness|{u}", f"median {rich.median():.0f} (range {rich.min()}-{rich.max()}); total {int((ft[cols].sum(axis=1) > 0).sum())}; rarefied 20k median {np.median(rar):.0f}")
    note(f"topgenera_pooled|{u}", "; ".join(f"{g} {v:.1%} (mean {meanra[g]:.1%})" for g, v in top.items()))

for st in ["OSD-772", "PRJNA1145089", "OSD-267", "OSD-269"]:
    note(f"asvs_study|{st}", int(pres(groups.index[groups.study == st]).sum()))

# ------------------------------------------------------------------ repeated rarefaction (100x) on the same denominator
rar_rows = []
rr_rng = np.random.default_rng(SEED)
for depth_r in [20000, 1000]:
    for u in ["ISS root zone", "ISS arcillite", "Field pepper rhizosphere"]:
        cols = [c for c in groups.index[UNITS[u]] if groups.loc[c, "depth"] >= depth_r]
        res = {d: [] for d in ["canonical", "broad"]}
        mats = {c: ft[c].to_numpy() for c in cols}
        for _ in range(100):
            hits = {d: 0 for d in res}
            for c in cols:
                sub = rr_rng.multivariate_hypergeometric(mats[c], depth_r)
                for d in res:
                    hits[d] += int(sub[masks[d].to_numpy()].sum() > 0)
            for d in res:
                res[d].append(hits[d])
        for d in res:
            v = np.array(res[d])
            rar_rows.append({"depth": depth_r, "unit": u, "definition": d, "n": len(cols),
                             "median": float(np.median(v)), "min": int(v.min()), "max": int(v.max())})
rar = pd.DataFrame(rar_rows)
rar.to_csv(OUT / "rarefaction_v5.tsv", sep="\t", index=False)
for _, r in rar.iterrows():
    note(f"rarefy|{r.depth}|{r.unit}|{r.definition}", f"median {r['median']:.0f}/{r.n} (range {r['min']}-{r['max']})")

# ------------------------------------------------------------------ nitrogen budget (illustrative)
fert = pd.DataFrame([  # Khodadad et al. 2026, Table 1; N split from current product labels
    ("Coated Ca(NO3)2 12-0-0 (T70)", 1.6, 0.12, 0.0, 0.12),
    ("Nutricote 14-4-14 (T100)", 4.0, 0.14, 0.075, 0.065),
    ("Nutricote 0-0-19 + Mg (T180)", 1.6, 0.0, 0.0, 0.0),
    ("Nutricote 14-4-14 (T180)", 5.6, 0.14, 0.075, 0.065)],
    columns=["product", "g_per_quadrant", "N_frac", "NH4N_frac", "NO3N_frac"])
nq = (fert.g_per_quadrant * fert.N_frac).sum(); nh4 = (fert.g_per_quadrant * fert.NH4N_frac).sum(); no3 = (fert.g_per_quadrant * fert.NO3N_frac).sum()
area, days, quadrants = 0.2, 137, 4  # 48 seeds at 240 plants m-2; 137 days; 4 quadrants
n_carrier = nq * quadrants
urine = (7, 16)  # g N per crew member per day (Verbeelen et al. 2021)
mol = tuple(u / 14.007 for u in urine)
koh = tuple(m * 56.11 for m in mol)
c3 = tuple(u / 1350 * 1000 for u in urine)  # L bed at 1.35 kg N m-3 d-1 loading (MELiSSA C-III)
mabr = tuple(u * 0.5 / 1.8 for u in urine)  # m2 membrane, half the N nitrified at 1.8 g N m-2 d-1 (Udert & Waechter 2012)
def tan_for_fa(fa, ph, t=25):  # Anthonisen et al. 1976: FA (mg NH3 L-1) = 17/14 * TAN * 10^pH / (exp(6344/(273+T)) + 10^pH)
    return fa * 14 / 17 * (math.exp(6344 / (273 + t)) + 10 ** ph) / 10 ** ph
doublings = math.log2(b_fd.pooled_RA / b_rz.pooled_RA)
budget = pd.DataFrame([
    ("PH-04 fertiliser N per quadrant", f"{nq:.2f} g N ({nh4:.2f} g NH4+-N, {nh4 / nq:.0%}; {no3:.2f} g NO3−-N)", "Khodadad et al. 2026 Table 1; product labels"),
    ("PH-04 fertiliser N per science carrier", f"{n_carrier:.1f} g N on {area} m² = {n_carrier / area:.0f} g N m⁻² per crop cycle ({n_carrier / area / days:.2f} g N m⁻² d⁻¹ over {days} d)", "4 quadrants; 48 seeds at 240 plants m⁻²"),
    ("Crew urine N", f"{urine[0]}–{urine[1]} g N per crew member per day = {urine[0] / n_carrier:.1f}–{urine[1] / n_carrier:.1f} PH-04 carrier loads per crew-day", "Verbeelen et al. 2021"),
    ("Crop area that would take up one crew member's urine N at the PH-04 supply rate", f"{urine[0] / (n_carrier / area / days):.0f}–{urine[1] / (n_carrier / area / days):.0f} m²", "derived"),
    ("Alkalinity", f"Ureolysis releases ≈1 eq base per mol N; nitrification consumes 2. Without added base about half of the ammonium is oxidised (NH4+:NO3− ≈ 1:1); full nitrification needs ≈{mol[0]:.1f}–{mol[1]:.1f} mol base per crew-day (≈{koh[0]:.0f}–{koh[1]:.0f} g KOH)", "stoichiometry; Udert & Wächter 2012"),
    ("Upstream reactor size per crew member", f"MELiSSA C-III fixed bed at ≤1.35 kg N m⁻³ d⁻¹: ≈{c3[0]:.0f}–{c3[1]:.0f} L of bed; urine MABR at 1.8 g N m⁻² d⁻¹ (half nitrified): ≈{mabr[0]:.0f}–{mabr[1]:.0f} m² of membrane", "Verbeelen et al. 2021; Gòdia et al. 2002; Udert & Wächter 2012"),
    ("NOB inhibition by free ammonia (0.1–1.0 mg NH3 L⁻¹)", f"total ammonium ≈{tan_for_fa(0.1, 7.0):.0f}–{tan_for_fa(1.0, 7.0):.0f} mg N L⁻¹ at pH 7.0; ≈{tan_for_fa(0.1, 7.5):.0f}–{tan_for_fa(1.0, 7.5):.0f} at pH 7.5 (25 °C)", "Anthonisen et al. 1976"),
    ("Start-up from the observed root-zone level to the field level", f"{b_fd.pooled_RA / b_rz.pooled_RA:,.0f}-fold ≈ {doublings:.0f} doublings; {doublings * math.log(2) / 0.5:.0f}–{doublings * math.log(2) / 0.1:.0f} d at net growth 0.5–0.1 d⁻¹ (cycle: {days} d)", "relative abundance as proxy; Prosser 1989"),
], columns=["Quantity", "Value", "Basis"])
budget.to_csv(OUT / "table5_nitrogen_budget.tsv", sep="\t", index=False)
note("fert_N_quadrant", nq); note("fert_NH4_share", nh4 / nq); note("doublings", doublings)

# ------------------------------------------------------------------ manuscript tables
track = pd.read_csv(WORK / "dada2" / "track.tsv", sep="\t")
ret = {s: f"{t.nonchim.sum() / t.input.sum():.1%}" for s, t in track.groupby("study")}
dm = lambda sel: f"{groups[sel].depth.median():,.0f} ({groups[sel].depth.min():,.0f})"
t1 = pd.DataFrame([
    ["OSD-772", "ISS, Advanced Plant Habitat (PH-04, 2021); chile pepper; one plant per quadrant after day 28",
     "Arcillite with controlled-release fertiliser. Root sections 16, root–substrate samples 4, arcillite 16, wick 16, foam 12, water 3, carrier swabs 8 ^b^",
     "V4, 515F/806R; MiSeq 2 × 250", f"Root zone {dm(UNITS['ISS root zone'])}; arcillite {dm(UNITS['ISS arcillite'])}", ret["OSD-772"]],
    ["PRJNA1145089", "Field, Hunan, China (2021); five pepper cultivars", "Soil. Rhizosphere by root shaking, 4 plants × 5 cultivars = 20",
     "V3–V4, 341F/805R; NovaSeq (V4 part extracted)", dm(UNITS["Field pepper rhizosphere"]), ret["PRJNA1145089"]],
    ["OSD-267", "ISS Veggie, VEG-01A; lettuce", "Arcillite pillows with Nutricote. Flight roots 2, ground-control roots 2", "V4, 515F/806R; MiSeq",
     f"Flight {dm((groups.study == 'OSD-267') & (groups.group == 'Lettuce Flight'))}; ground {dm((groups.study == 'OSD-267') & (groups.group == 'Lettuce Ground'))}", ret["OSD-267"]],
    ["OSD-269", "ISS Veggie, VEG-03A (2016); lettuce", "Arcillite pillows with Nutricote. Flight roots 2, ground-control roots 2", "V4, 515F/806R; MiSeq",
     f"Flight {dm((groups.study == 'OSD-269') & (groups.group == 'Lettuce Flight'))}; ground {dm((groups.study == 'OSD-269') & (groups.group == 'Lettuce Ground'))}", ret["OSD-269"]],
    ["OSD-65", "ISS interior surfaces, ISS Microbial Observatory (2015)", "Surface wipes, eight locations in two sessions plus one additional location, with and without PMA: 34 libraries ^c^",
     "V4, A519F/802R; GeneLab-processed table", f"{s65.bacterial_reads.median():,.0f} ({s65.bacterial_reads.min():,.0f})", "Not re-processed"],
], columns=["Dataset", "System and crop", "Substrate and samples used (n)", "Amplicon, platform", "Bacterial and archaeal reads per sample, median (min)", "Reads retained ^a^"])
t1.to_csv(OUT / "table1.tsv", sep="\t", index=False)


def sci(v, digits=2):
    if v == 0:
        return "0"
    e = int(math.floor(math.log10(abs(v))))
    m = v / 10 ** e
    return f"{m:.{digits - 1}f} × 10^−{-e}^" if e < 0 else f"{v:.{digits}g}"


def ra_cell(r):
    if r.nit_reads == 0:
        return f"≤{sci(r.RA_hi, 3)}"
    e = int(math.floor(math.log10(r.pooled_RA)))
    f = 10 ** e
    lo = r.RA_lo / f
    return f"{r.pooled_RA / f:.1f} ({lo:.1f}–{r.RA_hi / f:.1f}) × 10^−{-e}^" if lo >= 0.1 else f"{r.pooled_RA / f:.1f} ({lo:.2f}–{r.RA_hi / f:.1f}) × 10^−{-e}^"


def prev_cell(r):
    lo = "0" if r.prev_lo == 0 else f"{r.prev_lo * 100:.1f}"
    hi = "100" if r.prev_hi >= 1 else f"{r.prev_hi * 100:.1f}"
    return f"{int(r.pos1)}/{int(r.n)} ({lo}–{hi}%)"


t2rows = []
for u, lab in [("ISS root zone", "ISS root zone"), ("ISS root sections", "  root sections"), ("ISS root-substrate", "  root–substrate samples"),
               ("ISS arcillite", "ISS arcillite"), ("ISS water", "ISS water"), ("ISS wick", "ISS wick"), ("ISS foam", "ISS foam"),
               ("ISS carrier swabs", "ISS carrier swabs"), ("ISS carrier, all 75", "ISS carrier, all compartments"),
               ("ISS interior surfaces (OSD-65)", "ISS interior surfaces (OSD-65)"), ("Lettuce roots, flight", "Lettuce roots, flight"),
               ("Lettuce roots, ground", "Lettuce roots, ground control ^a^"), ("Field pepper rhizosphere", "Field pepper rhizosphere")]:
    c, b = D(u, "canonical"), D(u, "broad")
    t2rows.append([lab, f"{int(c.n)}", f"{c.depth_median:,.0f} ({int(c.depth_min):,})",
                   prev_cell(c), ra_cell(c), prev_cell(b), ra_cell(b)])
t2 = pd.DataFrame(t2rows, columns=["Group", "n", "Reads per sample, median (min)", "Canonical genera: positive (95% CI)",
                                   "Canonical genera: pooled RA (95% CI)", "Nitrifier lineages: positive (95% CI)",
                                   "Nitrifier lineages: pooled RA (95% CI)"])
t2.to_csv(OUT / "table2.tsv", sep="\t", index=False)

t3 = []
tt = tests.set_index(["comparison", "definition"])
r0 = tt.loc[("ISS root-zone quadrants vs field cultivars", "canonical")]
t3.append(["Primary", "Canonical genera, ISS quadrants vs field cultivars", r0.positive, f"Fisher p = {r0.Fisher_p:.3f}"])
t3.append(["Secondary", "Canonical genera, cluster relative abundance", "", f"exact MWU p = {r0.MWU_exact_p:.3f} ^b^"])
r1 = tt.loc[("ISS root-zone quadrants vs field cultivars", "canonical+AOA")]
t3.append(["Sensitivity", "Canonical genera plus cultured AOA genera", r1.positive, f"Fisher p = {r1.Fisher_p:.3f}; MWU p = {r1.MWU_exact_p:.3f}"])
r2 = tt.loc[("ISS root-zone quadrants vs field cultivars", "broad")]
t3.append(["Secondary", "Nitrifier lineages, ISS quadrants vs field cultivars", r2.positive, f"Fisher p = {r2.Fisher_p:.2f}; MWU p = {r2.MWU_exact_p:.3f}"])
r3 = tt.loc[("ISS arcillite quadrants vs field cultivars (substrate vs soil)", "broad")]
t3.append(["Secondary", "Nitrifier lineages, ISS arcillite quadrants vs field cultivars", r3.positive, f"Fisher p = {r3.Fisher_p:.2f}; MWU p = {r3.MWU_exact_p:.3f}"])
wa = wq.set_index("comparison").loc["Arcillite vs root zone"]
t3.append(["Secondary", "Nitrifier lineages, arcillite vs root zone within quadrants", f"{wa.pairs_comp_higher} quadrants higher in arcillite",
           f"sign test p = {wa.sign_one_sided_p:.4f} one-sided ({wa.sign_two_sided_p:.3f} two-sided); positive quadrants {wa.positive_quadrants}, Fisher p = {wa.Fisher_quadrants_p:.2f}"])
t3 = pd.DataFrame(t3, columns=["Role", "Comparison (unit: growth quadrant or cultivar)", "Positive units", "Test"])
t3.to_csv(OUT / "table3_tests.tsv", sep="\t", index=False)

# ------------------------------------------------------------------ supplement: per-sample table and nitrifier ASVs
sup = groups[groups.study.isin(["OSD-772", "PRJNA1145089", "OSD-267", "OSD-269"])][["study", "group", "compartment", "quadrant", "depth",
                                                                                       "canonical", "canonical+AOA", "broad", "broad-refined"]].copy()
lin_lab = tax.str.replace(r"^.*?f__([^;]+)(?:; g__([^;]+))?.*$", lambda m: (m.group(2) or m.group(1)), regex=True)
for lname, sub in ft.loc[masks["broad"]].groupby(lin_lab[masks["broad"]]):
    sup[f"reads: {lname}"] = sub.sum().reindex(sup.index).fillna(0).astype(int)
sup.rename(columns={"depth": "bacterial+archaeal reads"}).to_csv(OUT / "supplementary_per_sample.tsv", sep="\t")
with open(OUT / "nitrifier_asvs.fasta", "w") as fh:
    for a in ft.index[masks["broad"] & (ft[groups.index[UNITS["ISS carrier, all 75"] | UNITS["Lettuce roots, flight"]]].sum(axis=1) > 0)]:
        fh.write(f">{a} {tax[a]}\n{fasta[a]}\n")

pd.DataFrame(NUM).to_csv(OUT / "numbers.tsv", sep="\t", index=False)
with pd.option_context("display.width", 250, "display.max_colwidth", 230):
    print(pd.DataFrame(NUM).to_string(index=False))
