"""
50_nitrospira_phylogeny_reconstructed.py
========================================

*** RECONSTRUCTED SCRIPT (written 2026-10-19, after the 2026-04-23 snapshot) ***

The Nitrospira phylogeny reported in the manuscript (Fig. 5) was
produced with ad-hoc commands that were never saved to this repository. This script
re-implements that analysis from the Methods description and from the intermediate
files that were kept:

    acta_astronautica/nitrospira_asvs.fasta                     (19 ASV sequences)
    acta_astronautica/nitrospira_tree/combined_unaligned.fasta  (ASVs + references)
    acta_astronautica/nitrospira_tree/combined_aligned.fasta    (MAFFT output)
    acta_astronautica/nitrospira_tree/nitrospira_tree.nwk       (FastTree output)
    acta_astronautica/nitrospira_tree/nitrospira_asv_clade.csv  (clade calls)
    figures/aa/SuppS1e_nitrospira_tree.{pdf,png}

It is NOT the original code. Run with --compare to check its output against the
saved files above.

Steps
-----
1. Select ASVs whose genus is g__Nitrospira (config.TAXONOMY_FILE) and that carry
   reads in config.FEATURE_TABLE_CLEAN.
2. Sequences come from version-2_integrated/exported_sequences/dna-sequences.fasta
   (exported by 05_qiime2_pipeline.sh). If that file is absent, the saved
   acta_astronautica/nitrospira_asvs.fasta is used. ASV IDs are 'ASV_' + the first
   8 hex digits of the DADA2 feature ID (MD5 of the sequence); this is checked.
3. Reference set: acta_astronautica/nitrospira_tree/reference_set.fasta (13 records:
   SILVA 138 canonical NOB references, NCBI comammox 16S sequences and two
   Leptospirillum ferriphilum outgroup copies), copied from the saved
   combined_unaligned.fasta because the original download step was not saved.
4. Alignment with MAFFT v7 (--auto); tree with FastTree 2 (-nt -gtr, GTR+CAT).
5. The tree is rooted on the Leptospirillum outgroup and each ASV takes the clade
   of the reference at minimum patristic distance.
6. Tree figure with clade-coloured reference labels.

Requirements: biopython, pandas, matplotlib, plus `mafft` and `FastTree` on PATH
(both ship with the qiime2-amplicon conda env), or point the MAFFT / FASTTREE
environment variables at the binaries.

Usage
-----
    python 50_nitrospira_phylogeny_reconstructed.py [--outdir DIR] [--compare]
"""

from __future__ import annotations

import argparse
import hashlib
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

import config

ROOT = Path(__file__).resolve().parent
SAVED_DIR = ROOT / "acta_astronautica" / "nitrospira_tree"
SAVED_ASVS = ROOT / "acta_astronautica" / "nitrospira_asvs.fasta"
REFERENCE_FASTA = SAVED_DIR / "reference_set.fasta"
REP_SEQS_FASTA = config.VERSION_DIR / "exported_sequences" / "dna-sequences.fasta"
DEFAULT_OUTDIR = SAVED_DIR / "reconstructed"

# reference-name substring -> clade label (labels as in nitrospira_asv_clade.csv)
REFERENCE_CLADES = [
    ("Leptospirillum", "outgroup"),
    ("kreftii", "comammox-clade-B"),
    ("inopinata", "comammox-clade-A"),
    ("nitrosa", "comammox-clade-A"),
    ("nitrificans", "comammox-clade-A"),
    ("defluvii", "lineage-I-NOB"),
    ("moscoviensis", "lineage-II-NOB"),
    ("japonica", "lineage-II-NOB"),
    ("marina", "lineage-IV-NOB"),
]
CLADE_COLOURS = {
    "lineage-I-NOB": ("#2e86c1", "lineage I (N. defluvii)"),
    "lineage-II-NOB": ("#1e8449", "lineage II (N. moscoviensis / N. japonica)"),
    "lineage-IV-NOB": ("#7d3c98", "lineage IV (N. marina)"),
    "comammox-clade-A": ("#d35400", "comammox clade A (N. inopinata / nitrosa / nitrificans)"),
    "comammox-clade-B": ("#943126", "comammox clade B (N. kreftii)"),
    "outgroup": ("#999999", "outgroup (Leptospirillum ferriphilum)"),
}


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            records[name] = ""
        elif name is not None:
            records[name] += line.strip()
    return records


def write_fasta(records: dict[str, str], path: Path) -> None:
    with path.open("w") as fh:
        for name, seq in records.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


def asv_id(feature_id: str) -> str:
    return "ASV_" + feature_id[:8]


def reference_clade(name: str) -> str:
    for key, clade in REFERENCE_CLADES:
        if key in name:
            return clade
    raise ValueError(f"reference with no clade rule: {name}")


def tool(env_var: str, *candidates: str) -> str:
    path = os.environ.get(env_var) or next((shutil.which(c) for c in candidates if shutil.which(c)), None)
    if not path:
        sys.exit(f"{candidates[0]} not found: put it on PATH or set {env_var}")
    return path


def select_nitrospira_asvs() -> list[str]:
    tax = pd.read_csv(config.TAXONOMY_FILE, sep="\t", index_col=0)
    genus = tax["Taxon"].str.extract(r"g__([^;]+)")[0].str.strip()
    nitrospira = set(tax.index[genus == "Nitrospira"])
    table = pd.read_csv(config.FEATURE_TABLE_CLEAN, sep="\t", skiprows=1, index_col=0)
    with_reads = set(table.index[table.sum(axis=1) > 0])
    return sorted(nitrospira & with_reads)


def asv_sequences(feature_ids: list[str]) -> tuple[dict[str, str], str]:
    wanted = {asv_id(f): f for f in feature_ids}
    if REP_SEQS_FASTA.exists():
        source = str(REP_SEQS_FASTA.relative_to(ROOT))
        rep = read_fasta(REP_SEQS_FASTA)
        seqs = {asv_id(f): rep[f] for f in rep if f in set(feature_ids)}  # keeps rep-seqs order
    else:
        source = f"{SAVED_ASVS.relative_to(ROOT)} (rep-seqs export not found)"
        saved = read_fasta(SAVED_ASVS)
        seqs = {k: v for k, v in saved.items() if k in wanted}  # keeps saved order
    missing = set(wanted) - set(seqs)
    if missing:
        sys.exit(f"no sequence for: {sorted(missing)}")
    for name, seq in seqs.items():
        if hashlib.md5(seq.upper().encode()).hexdigest()[:8] != name[4:]:
            sys.exit(f"{name}: sequence does not match its MD5-derived ID")
    return seqs, source


def build_tree(unaligned: Path, aligned: Path, newick: Path) -> None:
    mafft = tool("MAFFT", "mafft")
    fasttree = tool("FASTTREE", "FastTree", "fasttree")
    with aligned.open("w") as out:
        subprocess.run([mafft, "--auto", "--quiet", str(unaligned)], stdout=out, check=True)
    with newick.open("w") as out:
        subprocess.run([fasttree, "-nt", "-gtr", "-quiet", str(aligned)], stdout=out,
                       stderr=subprocess.DEVNULL, check=True)


def assign_clades(newick: Path, references: list[str], asvs: list[str]):
    from Bio import Phylo

    tree = Phylo.read(str(newick), "newick")
    outgroup = [r for r in references if reference_clade(r) == "outgroup"]
    tree.root_with_outgroup(*outgroup)
    labelled = [r for r in references if reference_clade(r) != "outgroup"]
    rows = []
    for asv in asvs:
        dists = {ref: tree.distance(asv, ref) for ref in labelled}
        nearest = min(dists, key=dists.get)
        rows.append({"ASV": asv, "nearest_reference": nearest,
                     "branch_distance": dists[nearest], "clade": reference_clade(nearest)})
    return tree, pd.DataFrame(rows).sort_values("ASV").reset_index(drop=True)


def plot_tree(tree, clades: pd.DataFrame, references: list[str], stem: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from Bio import Phylo

    def short(name: str) -> str:
        if name.startswith("ASV_"):
            return name
        return "Leptospirillum" if name.startswith("Leptospirillum") else name.split("_")[0]

    colours = {short(r): CLADE_COLOURS[reference_clade(r)][0] for r in references}
    counts = clades["clade"].value_counts()
    fig, ax = plt.subplots(figsize=(12, 11))
    Phylo.draw(tree, axes=ax, do_show=False,
               label_func=lambda c: short(c.name) if c.is_terminal() and c.name else None,
               label_colors=lambda label: colours.get(label, "#555555"),
               branch_labels=lambda c: None)
    summary = ", ".join(f"{k} (n={counts[k]})" for k in CLADE_COLOURS if k in counts)
    ax.set_title(f"RECONSTRUCTED - phylogenetic placement of {len(clades)} Nitrospira ASVs\n{summary}",
                 fontsize=11)
    ax.set_xlabel("Substitutions per site (FastTree GTR+CAT)")
    ax.legend(handles=[Patch(color=c, label=l) for c, l in CLADE_COLOURS.values()],
              loc="lower right", fontsize=8)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def splits(tree, taxa: set[str]) -> set[frozenset]:
    anchor = sorted(taxa)[0]
    out = set()
    for clade in tree.get_nonterminals():
        side = frozenset(t.name for t in clade.get_terminals())
        if anchor in side:
            side = frozenset(taxa - side)
        if 1 < len(side) < len(taxa) - 1:
            out.add(side)
    return out


def compare(outdir: Path, seqs: dict[str, str], clades: pd.DataFrame, tree) -> None:
    from Bio import Phylo

    print("\n=== comparison with the saved outputs ===")
    saved_seqs = read_fasta(SAVED_ASVS)
    same_seqs = saved_seqs.keys() == seqs.keys() and all(saved_seqs[k].upper() == seqs[k].upper() for k in seqs)
    print(f"ASV set and sequences identical to nitrospira_asvs.fasta: {same_seqs}")

    saved_aln = {k: v.upper() for k, v in read_fasta(SAVED_DIR / "combined_aligned.fasta").items()}
    new_aln = {k: v.upper() for k, v in read_fasta(outdir / "combined_aligned.fasta").items()}
    print(f"alignment identical to combined_aligned.fasta: {saved_aln == new_aln}")

    saved_tree = Phylo.read(str(SAVED_DIR / "nitrospira_tree.nwk"), "newick")
    taxa = {t.name for t in tree.get_terminals()}
    a, b = splits(saved_tree, taxa), splits(tree, taxa)
    print(f"tree topology (unrooted splits): {len(a & b)}/{len(a | b)} shared, "
          f"Robinson-Foulds distance = {len(a ^ b)}")

    saved = pd.read_csv(SAVED_DIR / "nitrospira_asv_clade.csv").sort_values("ASV").reset_index(drop=True)
    merged = saved.merge(clades, on="ASV", suffixes=("_saved", "_new"))
    same_clade = (merged["clade_saved"] == merged["clade_new"]).sum()
    same_ref = (merged["nearest_reference_saved"] == merged["nearest_reference_new"]).sum()
    print(f"clade calls matching nitrospira_asv_clade.csv: {same_clade}/{len(merged)}; "
          f"nearest reference matching: {same_ref}/{len(merged)}")
    diff = merged[merged["clade_saved"] != merged["clade_new"]]
    if len(diff):
        print(diff[["ASV", "clade_saved", "clade_new"]].to_string(index=False))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--compare", action="store_true", help="compare with the saved outputs")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    feature_ids = select_nitrospira_asvs()
    seqs, source = asv_sequences(feature_ids)
    references = read_fasta(REFERENCE_FASTA)
    print(f"{len(seqs)} Nitrospira ASVs with reads; sequences from {source}")
    print(f"{len(references)} reference sequences from {REFERENCE_FASTA.relative_to(ROOT)}")

    write_fasta(seqs, args.outdir / "nitrospira_asvs.fasta")
    unaligned = args.outdir / "combined_unaligned.fasta"
    write_fasta({**references, **seqs}, unaligned)
    aligned = args.outdir / "combined_aligned.fasta"
    newick = args.outdir / "nitrospira_tree.nwk"
    build_tree(unaligned, aligned, newick)

    tree, clades = assign_clades(newick, list(references), list(seqs))
    clades.to_csv(args.outdir / "nitrospira_asv_clade.csv", index=False)
    plot_tree(tree, clades, list(references), args.outdir / "SuppS1e_nitrospira_tree_reconstructed")
    print(clades["clade"].value_counts().to_string())
    print(f"outputs written to {args.outdir}")

    if args.compare:
        compare(args.outdir, seqs, clades, tree)


if __name__ == "__main__":
    main()
