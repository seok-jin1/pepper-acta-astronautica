from importlib import import_module
import pandas as pd

m = import_module("36_gtdb_ncycle_genes")


def test_detect_genes_from_gff_text_amoa():
    gff = (
        "##gff-version 3\n"
        "ctg1\tRefSeq\tCDS\t1\t800\t.\t+\t0\t"
        "ID=cds1;product=methane monooxygenase/ammonia monooxygenase subunit A\n"
        "ctg1\tRefSeq\tCDS\t900\t1600\t.\t+\t0\t"
        "ID=cds2;product=hypothetical protein\n"
    )
    hits = m.detect_genes_from_gff_text(gff)
    assert "amoA" in hits
    assert "nxrA" not in hits


def test_detect_genes_from_gff_text_nosz():
    gff = (
        "c1\tRefSeq\tCDS\t1\t900\t.\t+\t0\t"
        "ID=cds1;product=nitrous-oxide reductase\n"
    )
    hits = m.detect_genes_from_gff_text(gff)
    assert "nosZ" in hits


def test_build_matrix():
    genera = ["GenusA", "GenusB"]
    genes = ["amoA", "nxrA", "nirK", "nirS", "nosZ"]
    hits = {"GenusA": {"amoA", "nxrA"}, "GenusB": {"nirK", "nosZ"}}
    mat = m.build_matrix(genera, genes, hits)
    assert mat.shape == (2, 5)
    assert mat.loc["GenusA", "amoA"] == 1
    assert mat.loc["GenusA", "nirK"] == 0
    assert mat.loc["GenusB", "nosZ"] == 1
