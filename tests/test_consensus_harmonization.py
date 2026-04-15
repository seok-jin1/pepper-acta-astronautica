import pandas as pd
from importlib import import_module
m = import_module("35_picrust2_faprotax_consensus")

def test_map_picrust_to_function_groups_nitrification():
    # PICRUSt2 MetaCyc pathways relevant to nitrification
    pathways = pd.Series({
        "PWY-6523": 2.0,      # ammonia oxidation I
        "PWY-5674": 1.5,      # nitrate reduction IV (denitrification)
        "UNRELATED-PWY": 0.1,
    })
    mapping = m.build_function_mapping()
    assert "nitrification" in mapping
    assert "PWY-6523" in mapping["nitrification"]
    agg = m.aggregate_by_function(pathways, mapping)
    assert agg["nitrification"] == 2.0
    assert agg["denitrification"] == 1.5
    assert "UNRELATED-PWY" not in agg.index

def test_consensus_direction_agreement():
    faprotax_log2fc = pd.Series({"nitrification": -10.0, "denitrification": -3.0, "chemoheterotrophy": 0.5})
    picrust_log2fc  = pd.Series({"nitrification":  -6.0, "denitrification": -2.5, "chemoheterotrophy": 0.8})
    result = m.consensus(faprotax_log2fc, picrust_log2fc)
    assert (result["agreement"] == True).all()
    assert result.loc["nitrification", "consensus_log2fc"] < 0
