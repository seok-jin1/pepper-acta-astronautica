import pandas as pd
from importlib import import_module

m = import_module("35_picrust2_faprotax_consensus")


def test_build_function_mapping_has_denitrification_coverage():
    mapping = m.build_function_mapping()
    assert "denitrification" in mapping
    assert "DENITRIFICATION-PWY" in mapping["denitrification"]
    # Nitrification and nitrogen_fixation are retained as keys but have no
    # PICRUSt2 pathway coverage (documented gap); aggregator drops them.
    assert mapping.get("nitrification") == []
    assert mapping.get("nitrogen_fixation") == []


def test_aggregate_by_function_routes_pathways_correctly():
    mapping = {
        "denitrification": ["DENITRIFICATION-PWY", "PWY-7084"],
        "chemoheterotrophy": ["GLYCOLYSIS", "TCA"],
        "nitrification": [],  # no coverage
    }
    pathways = pd.Series({
        "DENITRIFICATION-PWY": 1.0,
        "PWY-7084": 0.5,
        "GLYCOLYSIS": 2.0,
        "TCA": 3.0,
        "UNRELATED-PWY": 99.0,
    })
    agg = m.aggregate_by_function(pathways, mapping)
    assert agg["denitrification"] == 1.5
    assert agg["chemoheterotrophy"] == 5.0
    assert "nitrification" not in agg.index  # dropped (empty mapping)
    assert "UNRELATED-PWY" not in agg.index


def test_consensus_direction_agreement():
    faprotax_log2fc = pd.Series(
        {"denitrification": -3.0, "chemoheterotrophy": 0.5, "fermentation": -1.0}
    )
    picrust_log2fc = pd.Series(
        {"denitrification": -2.5, "chemoheterotrophy": 0.8, "fermentation": -0.9}
    )
    result = m.consensus(faprotax_log2fc, picrust_log2fc)
    assert (result["agreement"] == True).all()
    assert result.loc["denitrification", "consensus_log2fc"] < 0
    assert result.loc["chemoheterotrophy", "consensus_log2fc"] > 0


def test_consensus_disagreement_flag():
    faprotax = pd.Series({"f1": -1.0, "f2": 2.0})
    picrust = pd.Series({"f1": 0.5, "f2": 2.5})
    result = m.consensus(faprotax, picrust)
    assert result.loc["f1", "agreement"] == False
    assert result.loc["f2", "agreement"] == True
