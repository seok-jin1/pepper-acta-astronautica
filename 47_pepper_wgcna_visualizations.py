"""
47_pepper_wgcna_visualizations.py
==================================
Publication-quality supplementary figures from pepper WGCNA (power=12).

Produces three figures (idea 1, 2, 3 from brainstorm):

  (A) Module × Function heatmap, Terrestrial vs Space split
  (B) Module eigen-taxon trajectory across samples
  (C) Meta-network of module eigengenes (6-node graph)

Inputs (from 46_pepper_wgcna.R):
  - acta_astronautica/pepper_wgcna/pepper_wgcna_module_assignments.csv
  - acta_astronautica/pepper_wgcna/pepper_wgcna_module_eigengenes.csv
  - acta_astronautica/pepper_wgcna/pepper_wgcna_module_trait_corr.csv
  - version-2_integrated/exported_table_clean/feature-table.tsv
  - version-2_integrated/exported_taxonomy/taxonomy.tsv
  - version-2_integrated/integrated_metadata.tsv
  - version-2_integrated/FAPROTAX_report.txt

Outputs (acta_astronautica/pepper_wgcna/):
  - SuppS6a_module_function_heatmap.{pdf,png}
  - SuppS6b_eigentaxon_trajectory.{pdf,png}
  - SuppS6c_module_metanetwork.{pdf,png}
  - SuppS6_WGCNA_combined.{pdf,png}  (3-panel combined)
"""

from __future__ import annotations

import re
import sys
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config

OUT_DIR = config.AA_DIR / "pepper_wgcna"
OUT_DIR.mkdir(parents=True, exist_ok=True)

ROOT_ZONE = {"Plant Roots", "Plant_Roots", "Rhizosphere"}

# N-cycle focus ordered
FOCUS_FUNCTIONS = [
    "nitrification",
    "nitrogen_fixation",
    "denitrification",
    "nitrate_reduction",
    "ureolysis",
    "fermentation",
    "aerobic_chemoheterotrophy",
]
FUNCTION_LABELS = {
    "nitrification":              "Nitrification",
    "nitrogen_fixation":          "N₂ fixation",
    "denitrification":            "Denitrification",
    "nitrate_reduction":          "Nitrate reduction",
    "ureolysis":                  "Ureolysis",
    "fermentation":               "Fermentation",
    "aerobic_chemoheterotrophy":  "Aerobic chemo-\nheterotrophy",
}

# Consistent module ordering (larger Terr-correlated first)
MODULE_ORDER_PREF = ["turquoise", "green", "grey", "blue", "yellow", "brown"]


def parse_genus(tx):
    if pd.isna(tx):
        return "Unassigned"
    for part in reversed(str(tx).split(";")):
        p = part.strip()
        if p.startswith("g__"):
            name = p.split("__", 1)[-1].strip()
            if name and name not in ("uncultured", "metagenome", "unidentified"):
                return name
    return "Unassigned"


def parse_faprotax_report(path):
    func2genera = {}
    current_func = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            m = re.match(r"^# (\w+) \(\d+ records\):", line)
            if m:
                current_func = m.group(1)
                func2genera[current_func] = set()
                continue
            if current_func and line.startswith("    "):
                for part in line.split(";"):
                    part = part.strip()
                    if part.startswith("g__"):
                        name = part.split("__", 1)[-1].strip()
                        if name and name not in ("uncultured", "metagenome",
                                                 "unidentified", ""):
                            func2genera[current_func].add(name)
    return func2genera


def load_genus_ra():
    ft = pd.read_csv(config.FEATURE_TABLE_CLEAN, sep="\t", skiprows=1, index_col=0)
    tax = pd.read_csv(config.TAXONOMY_FILE, sep="\t", index_col=0)
    meta = pd.read_csv(config.INTEGRATED_METADATA, sep="\t", index_col=0)

    shared = ft.index.intersection(tax.index)
    ft = ft.loc[shared]
    genus = tax.loc[shared, "Taxon"].apply(parse_genus)
    ft = ft.copy()
    ft["Genus"] = genus.values
    g = ft.groupby("Genus").sum()
    g = g.drop(index="Unassigned", errors="ignore")
    col_sums = g.sum(axis=0).replace(0, 1)
    g = g.div(col_sums, axis=1) * 100
    return g, meta


def filter_root_samples(genus_ra, meta):
    sp = meta.index[(meta["Study_Group"] == "Space_Flight")
                    & meta["sample_type"].isin(ROOT_ZONE)]
    tr = meta.index[meta["Study_Group"] == "Terrestrial_Soil"]
    keep = [s for s in list(sp) + list(tr) if s in genus_ra.columns]
    return genus_ra[keep], list(sp), list(tr)


# ─────────────────────────────────────────────────────────────────────────────
# Helper: contrast-aware heatmap annotation
# ─────────────────────────────────────────────────────────────────────────────

def _annotate_heatmap_contrast(ax, data, vmax, threshold_ratio=0.55,
                                fontsize=9, dark_color="#1a1a1a",
                                light_color="#ffffff"):
    """Write integer values inside every cell of a seaborn heatmap with
    automatic contrast: dark text on light cells, light text on dark cells.
    Call AFTER sns.heatmap(..., annot=False)."""
    threshold = vmax * threshold_ratio
    nrows, ncols = data.shape
    values = data.values
    for i in range(nrows):
        for j in range(ncols):
            v = values[i, j]
            color = light_color if v >= threshold else dark_color
            ax.text(j + 0.5, i + 0.5, f"{int(v)}",
                    ha="center", va="center",
                    color=color, fontsize=fontsize,
                    fontweight="bold" if v >= threshold else "normal")


# ─────────────────────────────────────────────────────────────────────────────
# FIG A: Module × Function heatmap (Terr vs Space)
# ─────────────────────────────────────────────────────────────────────────────

def fig_module_function(module_df, genus_ra, sp_ids, tr_ids, func2genera,
                        trait_df):
    """For each (module, function) cell:
         cell value = number of module-member genera assigned to the function
                      AND detected (>0 RA) in that subset.
    Two side-by-side heatmaps: Terr (left), Space (right).
    """
    mods = [m for m in MODULE_ORDER_PREF if m in module_df["module_color"].unique()]
    module_members = {m: set(module_df.loc[module_df["module_color"] == m,
                                            "genus"]) for m in mods}

    # For presence check, require >0 RA in ≥1 sample of subset
    sp_cols = [s for s in sp_ids if s in genus_ra.columns]
    tr_cols = [s for s in tr_ids if s in genus_ra.columns]

    def cell_count(members, function, subset_cols):
        assigned = members & func2genera.get(function, set())
        if not assigned:
            return 0
        present = [g for g in assigned if g in genus_ra.index]
        if not present:
            return 0
        subset = genus_ra.loc[present, subset_cols]
        detected = (subset > 0).any(axis=1)
        return int(detected.sum())

    heat_terr = pd.DataFrame(
        {f: [cell_count(module_members[m], f, tr_cols) for m in mods]
         for f in FOCUS_FUNCTIONS},
        index=mods,
    )
    heat_space = pd.DataFrame(
        {f: [cell_count(module_members[m], f, sp_cols) for m in mods]
         for f in FOCUS_FUNCTIONS},
        index=mods,
    )

    vmax = max(heat_terr.values.max(), heat_space.values.max(), 1)

    fig, axes = plt.subplots(1, 2, figsize=(13, max(4, 0.5 * len(mods) + 2)),
                             constrained_layout=True)

    common_kwargs = dict(cmap="YlOrRd", vmin=0, vmax=vmax,
                         annot=False,
                         linewidths=0.6, linecolor="white",
                         cbar_kws={"label": "# module genera\ndetected"},
                         xticklabels=[FUNCTION_LABELS[f] for f in FOCUS_FUNCTIONS])

    def trait_tag(m):
        row = trait_df[trait_df["module"] == f"ME{m}"]
        if len(row) == 0:
            return ""
        r = float(row.iloc[0]["cor_Space"])
        return f" (r={r:+.2f})"

    ylabels = [f"{m}" + trait_tag(m) for m in mods]

    sns.heatmap(heat_terr, ax=axes[0], yticklabels=ylabels, **common_kwargs)
    _annotate_heatmap_contrast(axes[0], heat_terr, vmax)
    axes[0].set_title("Terrestrial ($n=20$)", fontsize=12, fontweight="bold")
    axes[0].set_ylabel("WGCNA module")
    axes[0].set_xlabel("")
    axes[0].tick_params(axis='x', labelrotation=30)

    sns.heatmap(heat_space, ax=axes[1], yticklabels=ylabels, **common_kwargs)
    _annotate_heatmap_contrast(axes[1], heat_space, vmax)
    axes[1].set_title("Space ($n=20$ root-zone)", fontsize=12, fontweight="bold")
    axes[1].set_ylabel("")
    axes[1].set_xlabel("")
    axes[1].tick_params(axis='x', labelrotation=30)

    for ax in axes:
        plt.setp(ax.get_xticklabels(), ha="right")

    fig.suptitle("WGCNA module × FAPROTAX function: active-genera count per subset",
                 fontsize=12, fontweight="bold")

    for ext in ("pdf", "png"):
        fig.savefig(OUT_DIR / f"SuppS6a_module_function_heatmap.{ext}",
                    dpi=300 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: SuppS6a_module_function_heatmap.pdf/.png")
    return heat_terr, heat_space


# ─────────────────────────────────────────────────────────────────────────────
# FIG B: Eigen-taxon trajectory across samples
# ─────────────────────────────────────────────────────────────────────────────

def fig_eigentaxon_trajectory(eigen_df, meta, trait_df):
    # Join metadata (Study_Group, location/Q if present)
    join = eigen_df.join(meta[["Study_Group", "location", "sample_type"]],
                         how="inner")

    # Order samples: Terr first (grouped), then Space grouped by Q1-Q4
    def sort_key(row):
        grp_rank = 0 if row["Study_Group"] == "Terrestrial_Soil" else 1
        loc = row.get("location", "")
        return (grp_rank, str(loc), row.name)
    join = join.loc[sorted(join.index, key=lambda i: sort_key(join.loc[i]))]

    mods = [c for c in eigen_df.columns if c.startswith("ME")]
    # Re-order by preference
    mod_order = [f"ME{m}" for m in MODULE_ORDER_PREF if f"ME{m}" in mods]
    mod_order += [m for m in mods if m not in mod_order]

    n_mod = len(mod_order)
    fig, axes = plt.subplots(n_mod, 1, figsize=(12, 1.6 * n_mod),
                             sharex=True, constrained_layout=True)
    if n_mod == 1:
        axes = [axes]

    colors = {"Space_Flight": "#3498db", "Terrestrial_Soil": "#2ecc71"}
    x = np.arange(len(join))

    for ax, me in zip(axes, mod_order):
        vals = join[me].values
        cvec = [colors[g] for g in join["Study_Group"]]
        ax.axhline(0, color="gray", lw=0.5)
        ax.bar(x, vals, color=cvec, edgecolor="black", linewidth=0.3, alpha=0.85)
        # Vertical line between Terr and Space groups
        boundary = sum(join["Study_Group"] == "Terrestrial_Soil") - 0.5
        ax.axvline(boundary, color="black", lw=1.0, ls="--", alpha=0.7)

        # Annotate module name + trait correlation
        row = trait_df[trait_df["module"] == me]
        if len(row):
            r = float(row.iloc[0]["cor_Space"])
            p = float(row.iloc[0]["p_Space"])
            cor_txt = f"r$_{{\\mathrm{{Space}}}}$={r:+.2f}, p={p:.1g}"
        else:
            cor_txt = ""
        ax.set_ylabel(me.replace("ME", ""), rotation=0, ha="right", va="center",
                      fontweight="bold")
        ax.text(0.01, 0.95, cor_txt, transform=ax.transAxes,
                ha="left", va="top", fontsize=9,
                bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                          edgecolor="gray", alpha=0.8))

    axes[0].set_title("Module eigen-taxon trajectory across samples "
                      "($n=20$ Terrestrial | $n=20$ Space)",
                      fontsize=12, fontweight="bold")
    axes[-1].set_xlabel("Sample (ordered: Terrestrial then Space)")
    axes[-1].set_xticks([])

    # Legend
    handles = [mpatches.Patch(color=c, label=g.replace("_", " "))
               for g, c in colors.items()]
    axes[0].legend(handles=handles, loc="upper right", fontsize=9, frameon=True)

    for ext in ("pdf", "png"):
        fig.savefig(OUT_DIR / f"SuppS6b_eigentaxon_trajectory.{ext}",
                    dpi=300 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: SuppS6b_eigentaxon_trajectory.pdf/.png")


# ─────────────────────────────────────────────────────────────────────────────
# FIG C: Meta-network of module eigengenes
# ─────────────────────────────────────────────────────────────────────────────

def fig_metanetwork(eigen_df, module_df, trait_df, r_threshold=0.3):
    """Module eigengene correlations → 6-node network."""
    mods = [c for c in eigen_df.columns if c.startswith("ME")]
    mod_order = [f"ME{m}" for m in MODULE_ORDER_PREF if f"ME{m}" in mods]
    mod_order += [m for m in mods if m not in mod_order]
    eig = eigen_df[mod_order]

    # ME-ME correlations
    corr = eig.corr(method="pearson")

    # Module sizes
    size_map = module_df["module_color"].value_counts().to_dict()

    # Trait correlation (Space)
    trait_map = {r["module"]: float(r["cor_Space"])
                 for _, r in trait_df.iterrows()}

    G = nx.Graph()
    for m in mod_order:
        mod = m.replace("ME", "")
        G.add_node(m, size=size_map.get(mod, 1),
                   trait_cor=trait_map.get(m, 0.0),
                   color=mod)

    for m1, m2 in combinations(mod_order, 2):
        r = corr.loc[m1, m2]
        if abs(r) >= r_threshold:
            G.add_edge(m1, m2, weight=r)

    fig, ax = plt.subplots(figsize=(9, 8), constrained_layout=True)
    pos = nx.spring_layout(G, k=1.6, seed=42, iterations=100)

    # Node colors = module color name
    node_colors = [G.nodes[n]["color"] for n in G.nodes]
    node_sizes = [400 + 80 * G.nodes[n]["size"] for n in G.nodes]

    # Node edge (border) color = Space enrichment (red = Space-enriched,
    # blue = Terr-enriched)
    node_edge_colors = []
    for n in G.nodes:
        tc = G.nodes[n]["trait_cor"]
        if tc > 0.3:
            node_edge_colors.append("#d62728")  # red = Space
        elif tc < -0.3:
            node_edge_colors.append("#1f77b4")  # blue = Terr
        else:
            node_edge_colors.append("gray")

    # Edges
    edge_colors = []
    edge_widths = []
    for u, v, d in G.edges(data=True):
        edge_colors.append("#2ca02c" if d["weight"] > 0 else "#e377c2")
        edge_widths.append(1 + 4 * abs(d["weight"]))

    nx.draw_networkx_edges(G, pos, ax=ax, width=edge_widths,
                            edge_color=edge_colors, alpha=0.7)
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors,
                            node_size=node_sizes, edgecolors=node_edge_colors,
                            linewidths=3, alpha=0.95)

    # Node labels = "module\n(size; r=±.xx)"
    labels = {n: f"{n.replace('ME','')}\n(n={G.nodes[n]['size']}; "
                 f"r={G.nodes[n]['trait_cor']:+.2f})"
              for n in G.nodes}
    nx.draw_networkx_labels(G, pos, labels=labels, ax=ax, font_size=9,
                             font_weight="bold")

    # Edge labels = correlation value
    edge_labels = {(u, v): f"{d['weight']:+.2f}"
                   for u, v, d in G.edges(data=True)}
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, ax=ax,
                                  font_size=8, bbox=dict(facecolor="white",
                                                          edgecolor="none",
                                                          pad=1))

    ax.set_title("Module-eigengene meta-network\n"
                 "(node color = module; border = trait direction; "
                 f"edges if $|r| \\geq {r_threshold}$)",
                 fontsize=12, fontweight="bold")
    ax.axis("off")

    # Legend
    legend_elems = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="w",
                   markeredgecolor="#d62728", markersize=15, linewidth=3,
                   label="Space-enriched ($r_{\\mathrm{Space}}>0.3$)"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="w",
                   markeredgecolor="#1f77b4", markersize=15, linewidth=3,
                   label="Terrestrial-enriched ($r_{\\mathrm{Space}}<-0.3$)"),
        plt.Line2D([0], [0], color="#2ca02c", lw=3, label="Positive ME correlation"),
        plt.Line2D([0], [0], color="#e377c2", lw=3, label="Negative ME correlation"),
    ]
    ax.legend(handles=legend_elems, loc="lower left", fontsize=8, frameon=True)

    for ext in ("pdf", "png"):
        fig.savefig(OUT_DIR / f"SuppS6c_module_metanetwork.{ext}",
                    dpi=300 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: SuppS6c_module_metanetwork.pdf/.png")
    return corr


# ─────────────────────────────────────────────────────────────────────────────
# Combined 3-panel figure
# ─────────────────────────────────────────────────────────────────────────────

def fig_combined(module_df, genus_ra, sp_ids, tr_ids, func2genera,
                  eigen_df, meta, trait_df, r_threshold=0.3):
    fig = plt.figure(figsize=(16, 17), constrained_layout=True)
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1.5, 1.4], hspace=0.25)

    # ── Panel A: Module × Function (2 subplots in row 0) ────────────────────
    axA_terr = fig.add_subplot(gs[0, 0])
    axA_space = fig.add_subplot(gs[0, 1])

    mods = [m for m in MODULE_ORDER_PREF if m in module_df["module_color"].unique()]
    module_members = {m: set(module_df.loc[module_df["module_color"] == m,
                                            "genus"]) for m in mods}
    sp_cols = [s for s in sp_ids if s in genus_ra.columns]
    tr_cols = [s for s in tr_ids if s in genus_ra.columns]

    def cell_count(members, function, subset_cols):
        assigned = members & func2genera.get(function, set())
        if not assigned:
            return 0
        present = [g for g in assigned if g in genus_ra.index]
        if not present:
            return 0
        subset = genus_ra.loc[present, subset_cols]
        return int((subset > 0).any(axis=1).sum())

    heat_terr = pd.DataFrame(
        {f: [cell_count(module_members[m], f, tr_cols) for m in mods]
         for f in FOCUS_FUNCTIONS}, index=mods)
    heat_space = pd.DataFrame(
        {f: [cell_count(module_members[m], f, sp_cols) for m in mods]
         for f in FOCUS_FUNCTIONS}, index=mods)
    vmax = max(heat_terr.values.max(), heat_space.values.max(), 1)

    def trait_tag(m):
        row = trait_df[trait_df["module"] == f"ME{m}"]
        if len(row) == 0:
            return ""
        return f" (r={float(row.iloc[0]['cor_Space']):+.2f})"

    ylabels = [f"{m}" + trait_tag(m) for m in mods]
    xlabels = [FUNCTION_LABELS[f] for f in FOCUS_FUNCTIONS]

    sns.heatmap(heat_terr, ax=axA_terr, cmap="YlOrRd", vmin=0, vmax=vmax,
                annot=False, linewidths=0.6, linecolor="white",
                yticklabels=ylabels, xticklabels=xlabels,
                cbar=False)
    _annotate_heatmap_contrast(axA_terr, heat_terr, vmax)
    axA_terr.set_title("(a) Terrestrial ($n=20$)",
                       fontsize=11, fontweight="bold", loc="left")
    axA_terr.set_ylabel("WGCNA module")
    plt.setp(axA_terr.get_xticklabels(), rotation=30, ha="right")

    sns.heatmap(heat_space, ax=axA_space, cmap="YlOrRd", vmin=0, vmax=vmax,
                annot=False, linewidths=0.6, linecolor="white",
                yticklabels=ylabels, xticklabels=xlabels,
                cbar_kws={"label": "# module genera\ndetected"})
    _annotate_heatmap_contrast(axA_space, heat_space, vmax)
    axA_space.set_title("(b) Space ($n=20$ root-zone)",
                        fontsize=11, fontweight="bold", loc="left")
    axA_space.set_ylabel("")
    plt.setp(axA_space.get_xticklabels(), rotation=30, ha="right")

    # ── Panel B: Eigen-taxon trajectory (row 1, full width) ─────────────────
    # use gs[1, :] → one big subplot area, then subdivide internally using
    # a sub-GridSpec
    gsB = gs[1, :].subgridspec(1, 1)
    axB_host = fig.add_subplot(gsB[0, 0])
    axB_host.set_title("(c) Module eigen-taxon trajectory across samples "
                       "($n=20$ Terrestrial | $n=20$ Space)",
                       fontsize=11, fontweight="bold", loc="left")
    axB_host.axis("off")

    inner_gs = gsB[0, 0].subgridspec(6, 1, hspace=0.12)

    join = eigen_df.join(meta[["Study_Group", "location", "sample_type"]],
                         how="inner")

    def sort_key(row):
        grp_rank = 0 if row["Study_Group"] == "Terrestrial_Soil" else 1
        loc = row.get("location", "")
        return (grp_rank, str(loc), row.name)
    join = join.loc[sorted(join.index, key=lambda i: sort_key(join.loc[i]))]
    x = np.arange(len(join))
    colors = {"Space_Flight": "#3498db", "Terrestrial_Soil": "#2ecc71"}

    mod_cols = [c for c in eigen_df.columns if c.startswith("ME")]
    mod_order = [f"ME{m}" for m in MODULE_ORDER_PREF if f"ME{m}" in mod_cols]
    mod_order += [m for m in mod_cols if m not in mod_order]

    for i, me in enumerate(mod_order):
        ax = fig.add_subplot(inner_gs[i])
        vals = join[me].values
        cvec = [colors[g] for g in join["Study_Group"]]
        ax.axhline(0, color="gray", lw=0.5)
        ax.bar(x, vals, color=cvec, edgecolor="black", linewidth=0.3, alpha=0.85)
        boundary = sum(join["Study_Group"] == "Terrestrial_Soil") - 0.5
        ax.axvline(boundary, color="black", lw=1.0, ls="--", alpha=0.7)
        row = trait_df[trait_df["module"] == me]
        if len(row):
            cor_txt = f"{me.replace('ME','')}  (r={float(row.iloc[0]['cor_Space']):+.2f})"
        else:
            cor_txt = me
        ax.set_ylabel(cor_txt, rotation=0, ha="right", va="center", fontsize=9,
                      fontweight="bold")
        if i < len(mod_order) - 1:
            ax.set_xticks([])
        else:
            ax.set_xticks([])
            ax.set_xlabel("Sample (Terrestrial → Space)")

    # Legend in first panel
    handles = [mpatches.Patch(color=c, label=g.replace("_", " "))
               for g, c in colors.items()]
    fig.axes[2].legend(handles=handles, loc="upper right", fontsize=8)

    # ── Panel C: Meta-network (row 2, full width) ──────────────────────────
    axC = fig.add_subplot(gs[2, :])
    eig = eigen_df[mod_order]
    corr = eig.corr(method="pearson")
    size_map = module_df["module_color"].value_counts().to_dict()
    trait_map = {r["module"]: float(r["cor_Space"])
                 for _, r in trait_df.iterrows()}
    G = nx.Graph()
    for m in mod_order:
        mod = m.replace("ME", "")
        G.add_node(m, size=size_map.get(mod, 1),
                   trait_cor=trait_map.get(m, 0.0), color=mod)
    for m1, m2 in combinations(mod_order, 2):
        r = corr.loc[m1, m2]
        if abs(r) >= r_threshold:
            G.add_edge(m1, m2, weight=r)

    pos = nx.spring_layout(G, k=1.4, seed=42, iterations=100)
    node_colors = [G.nodes[n]["color"] for n in G.nodes]
    node_sizes = [300 + 60 * G.nodes[n]["size"] for n in G.nodes]
    node_edge_colors = []
    for n in G.nodes:
        tc = G.nodes[n]["trait_cor"]
        if tc > 0.3:
            node_edge_colors.append("#d62728")
        elif tc < -0.3:
            node_edge_colors.append("#1f77b4")
        else:
            node_edge_colors.append("gray")
    edge_colors = ["#2ca02c" if d["weight"] > 0 else "#e377c2"
                   for _, _, d in G.edges(data=True)]
    edge_widths = [1 + 4 * abs(d["weight"]) for _, _, d in G.edges(data=True)]

    nx.draw_networkx_edges(G, pos, ax=axC, width=edge_widths,
                            edge_color=edge_colors, alpha=0.7)
    nx.draw_networkx_nodes(G, pos, ax=axC, node_color=node_colors,
                            node_size=node_sizes, edgecolors=node_edge_colors,
                            linewidths=2.5, alpha=0.95)
    labels = {n: f"{n.replace('ME','')}\nn={G.nodes[n]['size']}\nr={G.nodes[n]['trait_cor']:+.2f}"
              for n in G.nodes}
    nx.draw_networkx_labels(G, pos, labels=labels, ax=axC, font_size=8,
                             font_weight="bold")
    edge_labels = {(u, v): f"{d['weight']:+.2f}"
                   for u, v, d in G.edges(data=True)}
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, ax=axC,
                                  font_size=7,
                                  bbox=dict(facecolor="white", edgecolor="none",
                                            pad=1))

    axC.set_title("(d) Module-eigengene meta-network\n"
                  f"(borders: red = Space-enriched, blue = Terrestrial-enriched; "
                  f"edges: $|r| \\geq {r_threshold}$; green = positive, pink = negative)",
                  fontsize=11, fontweight="bold", loc="left")
    axC.axis("off")

    for ext in ("pdf", "png"):
        fig.savefig(OUT_DIR / f"SuppS6_WGCNA_combined.{ext}",
                    dpi=300 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: SuppS6_WGCNA_combined.pdf/.png")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    print("=" * 60)
    print("47: WGCNA publication-quality supplementary figures")
    print("=" * 60)

    module_df = pd.read_csv(OUT_DIR / "pepper_wgcna_module_assignments.csv")
    eigen_df = pd.read_csv(OUT_DIR / "pepper_wgcna_module_eigengenes.csv",
                            index_col=0)
    trait_df = pd.read_csv(OUT_DIR / "pepper_wgcna_module_trait_corr.csv")

    print(f"\n  Modules: {list(module_df['module_color'].unique())}")
    print(f"  Samples × ME: {eigen_df.shape}")

    func2genera = parse_faprotax_report(config.VERSION_DIR / "FAPROTAX_report.txt")
    print(f"  FAPROTAX functions parsed: {len(func2genera)}")

    genus_ra, meta = load_genus_ra()
    genus_ra_rz, sp_ids, tr_ids = filter_root_samples(genus_ra, meta)
    print(f"  Root-zone samples: Space={len(sp_ids)}, Terr={len(tr_ids)}")

    print("\n[A] Module × Function heatmap (Terr vs Space)...")
    fig_module_function(module_df, genus_ra_rz, sp_ids, tr_ids, func2genera,
                         trait_df)

    print("\n[B] Eigen-taxon trajectory...")
    fig_eigentaxon_trajectory(eigen_df, meta, trait_df)

    print("\n[C] Module meta-network...")
    fig_metanetwork(eigen_df, module_df, trait_df)

    print("\n[D] Combined 4-panel figure...")
    fig_combined(module_df, genus_ra_rz, sp_ids, tr_ids, func2genera,
                  eigen_df, meta, trait_df)

    print("\nAll 4 figures saved to acta_astronautica/pepper_wgcna/")


if __name__ == "__main__":
    main()
