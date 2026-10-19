"""
Utilities for the robustness and source-tracking expansion analyses.
"""

from __future__ import annotations

from collections import OrderedDict
from pathlib import Path
import re

import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import fisher_exact, mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests

import config


SPACE_GROUP = "Space_Flight"
TERR_GROUP = "Terrestrial_Soil"
GROUND_SEED_GROUP = "Ground_Seed"

ROBUSTNESS_DIR = config.VERSION_DIR / "robustness"
SOURCE_TRACKING_DIR = config.VERSION_DIR / "source_tracking"

FILTER_DEFINITIONS = OrderedDict(
    [
        (
            "ALL",
            {
                "label": "ALL",
                "space_sample_types": None,
                "description": "All ISS samples",
            },
        ),
        (
            "ROOT_ZONE",
            {
                "label": "ROOT_ZONE",
                "space_sample_types": ["Rhizosphere", "Plant Roots", "Arcillite"],
                "description": "Root-adjacent substrate and sink compartments",
            },
        ),
        (
            "ROOTS_RHIZO",
            {
                "label": "ROOTS_RHIZO",
                "space_sample_types": ["Rhizosphere", "Plant Roots"],
                "description": "Closest practical root-associated match",
            },
        ),
        (
            "STRICT_RHIZO",
            {
                "label": "STRICT_RHIZO",
                "space_sample_types": ["Rhizosphere"],
                "description": "Strict rhizosphere-only sensitivity analysis",
            },
        ),
    ]
)

PRIMARY_METRICS = OrderedDict(
    [
        (
            "nitrification_abundance",
            {
                "label": "Nitrification abundance",
                "kind": "continuous",
            },
        ),
        (
            "any_nitrifier_present",
            {
                "label": "Any nitrifier present",
                "kind": "binary",
            },
        ),
        (
            "ncycle_completeness",
            {
                "label": "N-cycle completeness",
                "kind": "continuous",
            },
        ),
        (
            "nitrification_redundancy",
            {
                "label": "Nitrification redundancy",
                "kind": "continuous",
            },
        ),
    ]
)

DEFAULT_NITRIFIER_GENERA = [
    "966-1",
    "Nitrosospira",
    "MND1",
    "Ellin6067",
    "IS-44",
    "mle1-7",
    "Nitrococcus",
    "GOUTA6",
    "DSSD61",
    "Candidatus_Nitrocosmicus",
]

PRIMARY_SOURCE_TYPES = ["Arcillite", "Wick", "Water", "seed", "Swab", "Foam"]
AUXILIARY_SOURCE_TYPES = ["fruit", "leaf", "stem"]
SINK_TYPES = ["Plant Roots", "Rhizosphere"]

NETWORK_R_THRESHOLD = 0.4
NETWORK_P_THRESHOLD = 0.05
NETWORK_TOP_N = 15


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def load_metadata() -> pd.DataFrame:
    return pd.read_csv(config.INTEGRATED_METADATA, sep="\t", index_col=0)


def load_feature_table() -> pd.DataFrame:
    return pd.read_csv(config.FEATURE_TABLE_CLEAN, sep="\t", skiprows=1, index_col=0)


def load_taxonomy() -> pd.DataFrame:
    return pd.read_csv(config.TAXONOMY_FILE, sep="\t", index_col=0)


def extract_genus(taxon_str: str) -> str:
    if pd.isna(taxon_str):
        return "Unassigned"
    for part in reversed(str(taxon_str).split(";")):
        part = part.strip()
        if part.startswith("g__"):
            name = part.split("__", 1)[-1].strip()
            return name if name else "Unassigned"
    return "Unassigned"


def aggregate_to_genus_counts(
    feature_table: pd.DataFrame | None = None,
    taxonomy: pd.DataFrame | None = None,
) -> pd.DataFrame:
    if feature_table is None:
        feature_table = load_feature_table()
    if taxonomy is None:
        taxonomy = load_taxonomy()

    shared = feature_table.index.intersection(taxonomy.index)
    table = feature_table.loc[shared].copy()
    taxonomy = taxonomy.loc[shared].copy()
    taxonomy["Genus"] = taxonomy["Taxon"].apply(extract_genus)
    table["Genus"] = taxonomy["Genus"]
    genus_counts = table.groupby("Genus").sum()
    return genus_counts.drop(index="Unassigned", errors="ignore")


def to_relative_abundance(count_table: pd.DataFrame) -> pd.DataFrame:
    totals = count_table.sum(axis=0).replace(0, np.nan)
    return count_table.div(totals, axis=1).fillna(0) * 100.0


def load_genus_relative() -> pd.DataFrame:
    return to_relative_abundance(aggregate_to_genus_counts())


def load_faprotax_relative() -> pd.DataFrame:
    func = pd.read_csv(
        config.VERSION_DIR / "FAPROTAX_functional_table.tsv",
        sep="\t",
        index_col=0,
        comment="#",
    )
    return func.div(func.sum(axis=0), axis=1).fillna(0) * 100.0


def parse_faprotax_report(report_path: Path | None = None) -> dict[str, set[str]]:
    if report_path is None:
        report_path = config.VERSION_DIR / "FAPROTAX_report.txt"
    if not report_path.exists():
        return {}

    func2genera: dict[str, set[str]] = {}
    current_func: str | None = None

    with open(report_path) as handle:
        for line in handle:
            line = line.rstrip()
            match = re.match(r"^# (\w+) \(\d+ records\):", line)
            if match:
                current_func = match.group(1)
                func2genera[current_func] = set()
                continue
            if current_func and line.startswith("    "):
                genus = None
                for part in line.split(";"):
                    part = part.strip()
                    if part.startswith("g__"):
                        name = part.split("__", 1)[-1].strip()
                        if name and name not in {"uncultured", "metagenome", "unidentified", ""}:
                            genus = name
                if genus:
                    func2genera[current_func].add(genus)

    return func2genera


def get_function_genera(
    function_name: str,
    fallback: list[str] | None = None,
    report_path: Path | None = None,
) -> list[str]:
    parsed = parse_faprotax_report(report_path=report_path)
    genera = sorted(parsed.get(function_name, set()))
    if genera:
        return genera
    return list(fallback or [])


NITRIFIER_GENERA = get_function_genera(
    "nitrification",
    fallback=DEFAULT_NITRIFIER_GENERA,
)


def load_ncycle_results() -> pd.DataFrame:
    candidate_paths = [
        config.VERSION_DIR / "Fig1c_NCycle_Completeness_Results.csv",
        config.VERSION_DIR / "Fig2B_NCycle_Completeness_Results.csv",
    ]
    for path in candidate_paths:
        if path.exists():
            return pd.read_csv(path, index_col=0)
    raise FileNotFoundError(
        "Could not find N-cycle completeness results in either Fig1c_* or Fig2B_* naming."
    )


def load_redundancy_results() -> pd.DataFrame:
    path = config.VERSION_DIR / "Fig1d_Functional_Redundancy_Results.csv"
    return pd.read_csv(path, index_col=0)


def load_alpha_diversity() -> pd.DataFrame:
    path = config.VERSION_DIR / "exported_diversity" / "alpha-diversity.tsv"
    alpha = pd.read_csv(path, sep="\t", index_col=0)
    alpha.columns = ["faith_pd"]
    return alpha


def get_filter_sample_ids(meta: pd.DataFrame, filter_name: str) -> tuple[list[str], list[str]]:
    definition = FILTER_DEFINITIONS[filter_name]
    terr_ids = meta.index[meta["Study_Group"] == TERR_GROUP].tolist()
    if definition["space_sample_types"] is None:
        space_mask = meta["Study_Group"] == SPACE_GROUP
    else:
        space_mask = (
            (meta["Study_Group"] == SPACE_GROUP)
            & meta["sample_type"].isin(definition["space_sample_types"])
        )
    space_ids = meta.index[space_mask].tolist()
    return space_ids, terr_ids


def get_filter_summary(meta: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for filter_name, definition in FILTER_DEFINITIONS.items():
        space_ids, terr_ids = get_filter_sample_ids(meta, filter_name)
        rows.append(
            {
                "filter": filter_name,
                "description": definition["description"],
                "space_sample_types": (
                    "ALL"
                    if definition["space_sample_types"] is None
                    else ", ".join(definition["space_sample_types"])
                ),
                "n_space": len(space_ids),
                "n_terr": len(terr_ids),
            }
        )
    return pd.DataFrame(rows)


def per_sample_nitrifier_presence(
    genus_rel: pd.DataFrame, sample_ids: list[str], nitrifier_genera: list[str] | None = None
) -> pd.Series:
    if nitrifier_genera is None:
        nitrifier_genera = NITRIFIER_GENERA
    present_genera = genus_rel.index.intersection(nitrifier_genera)
    sample_ids = [sample for sample in sample_ids if sample in genus_rel.columns]
    if len(present_genera) == 0:
        return pd.Series(0, index=sample_ids, dtype=int)
    values = (genus_rel.loc[present_genera, sample_ids] > 0).any(axis=0).astype(int)
    return pd.Series(values, index=sample_ids, dtype=int)


def prevalence_by_genus(
    genus_rel: pd.DataFrame, sample_ids: list[str], target_genera: list[str]
) -> pd.DataFrame:
    sample_ids = [sample for sample in sample_ids if sample in genus_rel.columns]
    rows = []
    for genus in target_genera:
        if genus in genus_rel.index:
            values = (genus_rel.loc[genus, sample_ids] > 0).astype(int)
            prevalence = float(values.mean()) if len(values) else np.nan
            mean_abundance = float(genus_rel.loc[genus, sample_ids].mean()) if len(values) else np.nan
        else:
            prevalence = 0.0
            mean_abundance = 0.0
        rows.append(
            {
                "genus": genus,
                "prevalence": prevalence,
                "mean_abundance": mean_abundance,
            }
        )
    return pd.DataFrame(rows)


def mean_difference(values_x: np.ndarray, values_y: np.ndarray) -> float:
    return float(np.mean(values_x) - np.mean(values_y))


def cliffs_delta(values_x: np.ndarray, values_y: np.ndarray) -> float:
    x = np.asarray(values_x, dtype=float)
    y = np.asarray(values_y, dtype=float)
    if len(x) == 0 or len(y) == 0:
        return np.nan
    greater = 0
    less = 0
    for value in x:
        greater += np.sum(value > y)
        less += np.sum(value < y)
    return float((greater - less) / (len(x) * len(y)))


def bootstrap_metric(
    values_x: np.ndarray,
    values_y: np.ndarray,
    stat_fn,
    n_boot: int = 2000,
    seed: int = config.RANDOM_SEED,
) -> tuple[float, float, float]:
    x = np.asarray(values_x, dtype=float)
    y = np.asarray(values_y, dtype=float)
    if len(x) == 0 or len(y) == 0:
        return np.nan, np.nan, np.nan
    observed = float(stat_fn(x, y))
    rng = np.random.default_rng(seed)
    estimates = []
    for _ in range(n_boot):
        x_boot = rng.choice(x, size=len(x), replace=True)
        y_boot = rng.choice(y, size=len(y), replace=True)
        estimates.append(stat_fn(x_boot, y_boot))
    lo, hi = np.quantile(estimates, [0.025, 0.975])
    return observed, float(lo), float(hi)


def fisher_p_from_binary(values_x: np.ndarray, values_y: np.ndarray) -> tuple[float, float]:
    x = np.asarray(values_x, dtype=int)
    y = np.asarray(values_y, dtype=int)
    contingency = np.array(
        [
            [int(np.sum(x == 1)), int(np.sum(x == 0))],
            [int(np.sum(y == 1)), int(np.sum(y == 0))],
        ]
    )
    odds_ratio, p_value = fisher_exact(contingency)
    return float(odds_ratio), float(p_value)


def mann_whitney_p(values_x: np.ndarray, values_y: np.ndarray) -> float:
    x = np.asarray(values_x, dtype=float)
    y = np.asarray(values_y, dtype=float)
    if len(x) == 0 or len(y) == 0:
        return np.nan
    return float(mannwhitneyu(x, y, alternative="two-sided").pvalue)


def select_top_shared_genera(
    space_rel: pd.DataFrame,
    terr_rel: pd.DataFrame,
    n: int = NETWORK_TOP_N,
    force_include: str | None = None,
) -> list[str]:
    shared = space_rel.index.intersection(terr_rel.index)
    combined_mean = (space_rel.loc[shared].mean(axis=1) + terr_rel.loc[shared].mean(axis=1)) / 2.0
    top_genera = combined_mean.sort_values(ascending=False).head(n).index.tolist()
    if force_include and force_include in shared and force_include not in top_genera:
        if len(top_genera) >= n:
            top_genera[-1] = force_include
        else:
            top_genera.append(force_include)
    return top_genera


def build_network(
    genus_rel: pd.DataFrame,
    top_genera: list[str],
    r_threshold: float = NETWORK_R_THRESHOLD,
    p_threshold: float = NETWORK_P_THRESHOLD,
    use_fdr: bool = False,
) -> nx.Graph:
    """Build co-occurrence network using matrix-based Spearman (matches script 09)."""
    raw_data = genus_rel.loc[top_genera]

    graph = nx.Graph()
    graph.add_nodes_from(top_genera)

    # Filter out constant-variance genera before matrix Spearman
    variances = raw_data.var(axis=1)
    variable_genera = [g for g in top_genera if variances.get(g, 0) > 0]
    n_genera = len(variable_genera)

    if n_genera < 2:
        return graph

    data = raw_data.loc[variable_genera].values.T  # (samples, genera)

    if n_genera == 2:
        corr_val, p_val = spearmanr(data[:, 0], data[:, 1])
        corr_matrix = np.array([[1.0, corr_val], [corr_val, 1.0]])
        pval_matrix = np.array([[0.0, p_val], [p_val, 0.0]])
    else:
        corr_matrix, pval_matrix = spearmanr(data)
        if np.ndim(corr_matrix) == 0:
            return graph

    triu_idx = np.triu_indices(n_genera, k=1)
    pvals = pval_matrix[triu_idx].copy()

    if use_fdr and len(pvals) > 0:
        valid = ~np.isnan(pvals)
        if valid.any():
            _, corrected, _, _ = multipletests(pvals[valid], method="fdr_bh")
            pvals[valid] = corrected
        pval_matrix = np.zeros_like(pval_matrix)
        pval_matrix[triu_idx] = pvals
        pval_matrix = pval_matrix + pval_matrix.T

    for idx, (i, j) in enumerate(zip(*triu_idx)):
        r = corr_matrix[i, j]
        p = pval_matrix[i, j]
        if np.isnan(r) or np.isnan(p):
            continue
        if abs(r) > r_threshold and p < p_threshold:
            graph.add_edge(variable_genera[i], variable_genera[j], weight=float(r), pvalue=float(p))
    return graph


def network_summary(
    graph: nx.Graph,
    group: str,
    filter_name: str,
    top_genera: list[str],
    n_samples: int,
) -> dict[str, object]:
    degrees = dict(graph.degree())
    top_hub = "none"
    hub_degree = 0
    if graph.number_of_nodes() > 0:
        top_hub, hub_degree = max(degrees.items(), key=lambda item: (item[1], item[0]))
    negative_edges = sum(1 for _, _, edge in graph.edges(data=True) if edge["weight"] < 0)
    avg_degree = float(np.mean([degree for _, degree in graph.degree()])) if graph.number_of_nodes() else 0.0
    return {
        "filter": filter_name,
        "group": group,
        "n_samples": n_samples,
        "n_nodes": graph.number_of_nodes(),
        "n_edges": graph.number_of_edges(),
        "density": float(nx.density(graph)) if graph.number_of_nodes() > 1 else 0.0,
        "avg_degree": avg_degree,
        "negative_edges": negative_edges,
        "top_hub": top_hub,
        "hub_degree": hub_degree,
        "rhodanobacter_degree": int(degrees.get("Rhodanobacter", 0)),
        "top_genera": ", ".join(top_genera),
    }


def bootstrap_columns(df: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    columns = df.columns.to_numpy()
    sampled = rng.choice(columns, size=len(columns), replace=True)
    return df.loc[:, sampled]


def collapse_feast_output(feast_matrix: pd.DataFrame, sink_meta: pd.DataFrame) -> pd.DataFrame:
    long_df = feast_matrix.copy()
    long_df["SampleID"] = long_df.index
    long_df = long_df.melt(id_vars="SampleID", var_name="Source", value_name="Contribution")
    return long_df.merge(sink_meta, on="SampleID", how="left")
