"""
Script 35: PICRUSt2 / FAPROTAX Consensus Harmonization

Harmonizes MetaCyc pathway-level Log2FC predictions from PICRUSt2 with
FAPROTAX function-level Log2FC into a consensus table for the Acta
Astronautica manuscript (Prediction P1).

Coverage caveat (documented 2026-04-15)
---------------------------------------
The default PICRUSt2 MetaCyc reference pathways used here cover
denitrification and broad heterotrophic/fermentative metabolism well,
but do NOT include canonical nitrification (ammonia/nitrite oxidation)
or nitrogen-fixation pathways. Those functions are therefore validated
via FAPROTAX plus the GTDB genome-level N-cycle gene evidence produced
by script 36 (Fig 3b), rather than by pathway-level consensus here.

Outputs:
    - acta_astronautica/consensus_function_log2fc.csv
    - figures/aa/Fig2d_consensus.{png,pdf}
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import config


# ---------------------------------------------------------------------------
# Function -> MetaCyc pathway mapping
#
# Curated 2026-04-15 against the actual pathway IDs present in the PICRUSt2
# output for this study. Each listed ID has been verified present in
# path_abun_unstrat.tsv.gz; functions left empty (nitrification, nitrogen
# fixation) have no PICRUSt2 MetaCyc pathway coverage and are handled by
# FAPROTAX + GTDB genome evidence elsewhere in the pipeline.
# ---------------------------------------------------------------------------
_FUNCTION_MAPPING = {
    # Canonical denitrification pathways (three hits).
    "denitrification": [
        "DENITRIFICATION-PWY",
        "PWY-5384",
        "PWY-7084",
    ],
    # Broad aerobic/anaerobic heterotrophic metabolism.
    "chemoheterotrophy": [
        "GLYCOLYSIS",
        "TCA",
        "PWY-6142",
    ],
    # Mixed-acid and sugar fermentation pathways.
    "fermentation": [
        "FERMENTATION-PWY",
        "ANAEROFRUCAT-PWY",
        "P122-PWY",
        "P461-PWY",
        "P161-PWY",
        "P241-PWY",
        "P108-PWY",
    ],
    # Aromatic compound degradation (PGPR-relevant trait).
    "aromatic_compound_degradation": [
        "P221-PWY",
        "P23-PWY",
        "P42-PWY",
    ],
    # Nitrification and nitrogen fixation: retained as keys with empty
    # ID lists so the coverage gap is explicit in diagnostics; both will
    # be dropped by aggregate_by_function.
    "nitrification": [],
    "nitrogen_fixation": [],
}


def build_function_mapping() -> dict:
    """Return the curated function -> [MetaCyc pathway IDs] mapping."""
    # Return a copy to avoid mutation of the module-level seed.
    return {k: list(v) for k, v in _FUNCTION_MAPPING.items()}


def aggregate_by_function(pathway_series: pd.Series, mapping: dict) -> pd.Series:
    """Aggregate a pathway-indexed Series into a function-indexed Series.

    For each function, sum the values of any pathway IDs that are present
    in ``pathway_series``. Pathways absent from the series are silently
    skipped. Functions with zero matched pathways are dropped.
    """
    out = {}
    for function, pathway_ids in mapping.items():
        matched = [p for p in pathway_ids if p in pathway_series.index]
        if not matched:
            continue
        out[function] = float(pathway_series.loc[matched].sum())
    return pd.Series(out, name="aggregated")


def consensus(faprotax_log2fc: pd.Series, picrust_log2fc: pd.Series) -> pd.DataFrame:
    """Build a consensus table comparing FAPROTAX and PICRUSt2 Log2FCs.

    Joins on shared function names (inner join). Columns:
        faprotax_log2fc, picrust_log2fc, agreement, consensus_log2fc

    ``agreement`` is True when both values share the same sign (both
    negative, both positive, or both zero). ``consensus_log2fc`` is the
    arithmetic mean of the two Log2FC values.
    """
    df = pd.DataFrame({
        "faprotax_log2fc": faprotax_log2fc,
        "picrust_log2fc": picrust_log2fc,
    }).dropna()
    df["agreement"] = np.sign(df["faprotax_log2fc"]) == np.sign(df["picrust_log2fc"])
    df["consensus_log2fc"] = df[["faprotax_log2fc", "picrust_log2fc"]].mean(axis=1)
    return df


# ---------------------------------------------------------------------------
# Pipeline helpers
# ---------------------------------------------------------------------------
def compute_picrust_function_log2fc(epsilon: float = 1e-6):
    """Compute function-level Log2FC (Space/Terrestrial) from PICRUSt2 pathways.

    Returns
    -------
    log2fc : pd.Series
        Function-indexed Series of log2((space+eps)/(terr+eps)).
    matched : dict
        Diagnostic dict {function: [matched_pathway_ids]}.
    """
    df = pd.read_csv(config.PICRUST2_PATHWAYS, sep="\t",
                     index_col=0, compression="gzip")
    meta = pd.read_csv(config.INTEGRATED_METADATA, sep="\t", index_col=0)

    common = df.columns.intersection(meta.index)
    df = df[common]
    meta = meta.loc[common]

    # Convert to relative abundance per sample (column-wise) as in script 07.
    df_rel = df.div(df.sum(axis=0), axis=1) * 100

    space_samples = meta.index[meta["Study_Group"] == "Space_Flight"]
    terr_samples = meta.index[meta["Study_Group"] == "Terrestrial_Soil"]

    space_mean = df_rel[space_samples].mean(axis=1)
    terr_mean = df_rel[terr_samples].mean(axis=1)

    mapping = build_function_mapping()

    # Diagnostic: which pathway IDs actually matched
    matched = {
        function: [p for p in pathway_ids if p in df_rel.index]
        for function, pathway_ids in mapping.items()
    }

    space_agg = aggregate_by_function(space_mean, mapping)
    terr_agg = aggregate_by_function(terr_mean, mapping)

    # Functions present in both aggregated series
    shared = space_agg.index.intersection(terr_agg.index)
    log2fc = np.log2((space_agg.loc[shared] + epsilon) /
                     (terr_agg.loc[shared] + epsilon))
    log2fc.name = "picrust_log2fc"
    return log2fc, matched


def load_faprotax_function_log2fc() -> pd.Series:
    """Load FAPROTAX function-level Log2FC from the comparison CSV."""
    path = config.VERSION_DIR / "FAPROTAX_group_comparison.csv"
    df = pd.read_csv(path)
    s = df.set_index("group")["Log2FC"]
    s.name = "faprotax_log2fc"
    return s


def _plot_consensus(consensus_df: pd.DataFrame, out_png: Path, out_pdf: Path):
    """Side-by-side bar plot of FAPROTAX vs PICRUSt2 log2fc per function."""
    sort_col = ("consensus_log2fc"
                if "consensus_log2fc" in consensus_df.columns
                else "faprotax_log2fc")
    plot_df = consensus_df.sort_values(sort_col)[["faprotax_log2fc",
                                                  "picrust_log2fc"]]

    functions = plot_df.index.tolist()
    x = np.arange(len(functions))
    width = 0.4

    fig, ax = plt.subplots(figsize=(max(6, 1.0 * len(functions) + 2), 5))
    ax.bar(x - width / 2, plot_df["faprotax_log2fc"], width,
           label="FAPROTAX", color="#2b8cbe")
    ax.bar(x + width / 2, plot_df["picrust_log2fc"], width,
           label="PICRUSt2", color="#e34a33")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(functions, rotation=30, ha="right")
    ax.set_ylabel("Log2 Fold Change (Space / Terrestrial)")
    ax.set_title("FAPROTAX vs PICRUSt2 consensus (Acta Astronautica P1)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf, dpi=300)
    plt.close(fig)


def main():
    config.AA_DIR.mkdir(parents=True, exist_ok=True)
    config.AA_FIG_DIR.mkdir(parents=True, exist_ok=True)

    print("[35] Computing PICRUSt2 function-level Log2FC...")
    picrust_log2fc, matched = compute_picrust_function_log2fc()
    print(f"[35] matched pathways: {matched}")

    # Drop functions where zero pathways matched (already excluded by aggregator,
    # but report which ones were dropped explicitly).
    dropped = [fn for fn, ids in matched.items() if not ids]
    if dropped:
        print(f"[35] Dropped (no MetaCyc IDs matched): {dropped}")

    print("[35] Loading FAPROTAX Log2FC...")
    faprotax_log2fc = load_faprotax_function_log2fc()

    print("[35] Building consensus table...")
    consensus_df = consensus(faprotax_log2fc, picrust_log2fc)
    consensus_df = consensus_df.sort_values("consensus_log2fc")

    out_csv = config.AA_CONSENSUS_CSV
    consensus_df.to_csv(out_csv)
    print(f"[35] Wrote consensus CSV -> {out_csv} ({len(consensus_df)} rows)")
    print(consensus_df.to_string())

    out_png = config.AA_FIG_DIR / "Fig2d_consensus.png"
    out_pdf = config.AA_FIG_DIR / "Fig2d_consensus.pdf"
    _plot_consensus(consensus_df, out_png, out_pdf)
    print(f"[35] Wrote figures -> {out_png}, {out_pdf}")


if __name__ == "__main__":
    main()
