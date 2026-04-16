import pandas as pd
from importlib import import_module

m = import_module("38_veggie_lettuce_figures")


def test_merge_pepper_lettuce_long():
    pepper = pd.DataFrame({"function": ["nitrification", "denitrification"], "log2fc": [-10.7, -3.0]})
    lettuce = pd.DataFrame({"function": ["nitrification", "chemoheterotrophy"], "log2fc": [-8.1, 0.4]})
    merged = m.merge_pepper_lettuce_long(pepper, lettuce)
    assert set(merged["crop"].unique()) == {"pepper", "lettuce"}
    assert len(merged) == 4


def test_ncycle_completeness_format():
    ncycle_steps = m.NCYCLE_STEPS
    assert "nitrification" in ncycle_steps
    assert "denitrification" in ncycle_steps
    assert len(ncycle_steps) == 5
