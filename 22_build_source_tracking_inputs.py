"""
Prepare genus- and ASV-level source-tracking inputs for FEAST.
"""

from __future__ import annotations

import pandas as pd

import config
import analysis_expansion_utils as aeu


def build_source_sink_metadata(space_meta: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    records = []
    sink_counter = 0
    for sample_id, row in space_meta.iterrows():
        sample_type = row["sample_type"]
        if sample_type in aeu.SINK_TYPES:
            role = "Sink"
            sink_counter += 1
            feast_id = sink_counter
            analysis_set = "primary"
        elif sample_type in aeu.PRIMARY_SOURCE_TYPES:
            role = "Source"
            feast_id = pd.NA
            analysis_set = "primary"
        elif sample_type in aeu.AUXILIARY_SOURCE_TYPES:
            role = "Source"
            feast_id = pd.NA
            analysis_set = "auxiliary"
        else:
            role = "Exclude"
            feast_id = pd.NA
            analysis_set = "excluded"
        records.append(
            {
                "SampleID": sample_id,
                "Env": sample_type,
                "SourceSink": role,
                "id": feast_id,
                "sample_type": sample_type,
                "location": row["location"],
                "sanitization": row["sanitization"],
                "analysis_set": analysis_set,
            }
        )
    metadata = pd.DataFrame(records)
    primary = metadata[metadata["analysis_set"] == "primary"].copy()
    primary = primary[primary["SourceSink"].isin(["Source", "Sink"])].copy()
    return metadata, primary


def build_function_linkage(space_meta: pd.DataFrame) -> pd.DataFrame:
    func_rel = aeu.load_faprotax_relative()
    ncycle_df = aeu.load_ncycle_results()
    redundancy_df = aeu.load_redundancy_results()

    rows = []
    for sample_id, row in space_meta.iterrows():
        rows.append(
            {
                "SampleID": sample_id,
                "sample_type": row["sample_type"],
                "location": row["location"],
                "nitrification_abundance": float(func_rel.loc["nitrification", sample_id])
                if sample_id in func_rel.columns
                else pd.NA,
                "nitrogen_fixation_abundance": float(func_rel.loc["nitrogen_fixation", sample_id])
                if sample_id in func_rel.columns
                else pd.NA,
                "ureolysis_abundance": float(func_rel.loc["ureolysis", sample_id])
                if sample_id in func_rel.columns
                else pd.NA,
                "ncycle_completeness": float(ncycle_df.loc[sample_id, "Completeness"])
                if sample_id in ncycle_df.index
                else pd.NA,
                "nitrification_redundancy": float(redundancy_df.loc[sample_id, "nitrification"])
                if sample_id in redundancy_df.index
                else pd.NA,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    print("=" * 70)
    print("Build Source-Tracking Inputs")
    print("=" * 70)

    outdir = aeu.ensure_dir(aeu.SOURCE_TRACKING_DIR)
    meta = aeu.load_metadata()
    space_meta = meta[meta["Study_Group"] == aeu.SPACE_GROUP].copy()

    feature_table = aeu.load_feature_table()
    taxonomy = aeu.load_taxonomy()
    genus_counts = aeu.aggregate_to_genus_counts(feature_table=feature_table, taxonomy=taxonomy)

    asv_counts_space = feature_table.loc[:, space_meta.index].T
    genus_counts_space = genus_counts.loc[:, space_meta.index].T

    metadata_all, metadata_primary = build_source_sink_metadata(space_meta)
    linkage_df = build_function_linkage(space_meta)

    asv_counts_space.to_csv(outdir / "asv_counts_space.tsv", sep="\t", index=True, index_label="SampleID")
    genus_counts_space.to_csv(outdir / "genus_counts_space.tsv", sep="\t", index=True, index_label="SampleID")
    metadata_all.to_csv(outdir / "source_sink_metadata.tsv", sep="\t", index=False)
    metadata_primary.to_csv(outdir / "primary_source_sink_metadata.tsv", sep="\t", index=False)
    linkage_df.to_csv(outdir / "function_linkage_space.tsv", sep="\t", index=False)

    print("Saved:")
    print(f"  {outdir / 'asv_counts_space.tsv'}")
    print(f"  {outdir / 'genus_counts_space.tsv'}")
    print(f"  {outdir / 'source_sink_metadata.tsv'}")
    print(f"  {outdir / 'primary_source_sink_metadata.tsv'}")
    print(f"  {outdir / 'function_linkage_space.tsv'}")


if __name__ == "__main__":
    main()
