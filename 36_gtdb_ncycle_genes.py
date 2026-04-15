"""
Script 36: GTDB N-cycle gene matrix via NCBI Datasets CLI (Strategy A)

Purpose
-------
Provide an independent, genome-level validation layer for Acta Astronautica
Prediction P2 (nitrogen-cycle functional redundancy in the ISS rhizosphere).
Task 2 showed that PICRUSt2 pathway-level prediction lacks coverage for
nitrification and canonical N-fixation, so this script is load-bearing for
the nitrification evidence.

Strategy
--------
For each target rhizosphere genus, fetch one reference RefSeq genome using
the NCBI Datasets CLI, parse the PGAP-annotated ``genomic.gff`` file, and
mark N-cycle gene presence by **substring regex matching on the ``product=``
field** (the ``gene=`` attribute is frequently missing in PGAP annotations
for N-cycle enzymes, but ``product=`` is reliable).

Limitations
-----------
- For each genus the script fetches either the type-strain reference or
  (via species override) a canonical species; ``fetch_genome_gffs``
  unions gene hits across all genomic.gff files in the bundle, so the
  matrix reflects "at least one canonical strain in this genus carries
  the gene" rather than presence in every strain.
- The ``nxrA`` pattern intentionally also matches membrane-bound nitrate
  reductase alpha subunit (``narG``): these enzymes are homologous and
  PGAP annotations often do not distinguish the nitrifying NO2->NO3
  (nxrA) from the denitrifying NO3->NO2 (narG) function. Positive hits
  in this column should therefore be read as "has a narG/nxrA-type
  alpha subunit".
- Nitrospira reference genomes annotate their nitrite oxidoreductase
  alpha subunit as the generic product ``molybdopterin oxidoreductase
  family protein``, which the regex does not match. The matrix may
  therefore under-report nxrA in Nitrospira; the Methods section
  should cite Daims et al. 2015 / Lucker et al. 2010 for literature
  confirmation of nxr presence in this lineage.
- Absence can reflect incomplete assembly or annotation rather than a
  true gene-content gap.

Outputs
-------
- ``config.AA_GTDB_MATRIX_CSV`` -- 13 genera x 5 gene presence/absence
- ``config.AA_DIR / "gtdb_ncycle_matrix_sources.csv"`` -- provenance
  (GCF accession used or failure reason) per genus.

No plotting is done here; figure composition is Task 9.
"""

from __future__ import annotations

import logging
import re
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path
from typing import Iterable

import pandas as pd

import config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
log = logging.getLogger("gtdb_ncycle")


TARGET_GENERA = [
    "Rhodanobacter",
    "Nitrosomonas",
    "Nitrosospira",
    "Nitrospira",
    "Nitrobacter",
    "Pseudomonas",
    "Bradyrhizobium",
    "Azospirillum",
    "Bacillus",
    "Paenibacillus",
    "Burkholderia",
    "Variovorax",
    "Sphingomonas",
]

# Species-level overrides for genera whose genus-level RefSeq fetch fails
# (either no --reference assembly or fallback download is too large to
# complete reliably). The scientifically informative representative is
# specified explicitly so sanity checks remain meaningful.
SPECIES_OVERRIDES: dict[str, str] = {
    "Nitrospira": "Nitrospira moscoviensis",
    "Bradyrhizobium": "Bradyrhizobium japonicum",
    "Bacillus": "Bacillus subtilis",
}

NCYCLE_GENES = ["amoA", "nxrA", "nirK", "nirS", "nosZ"]

NCYCLE_PRODUCT_PATTERNS: dict[str, list[str]] = {
    "amoA": [r"ammonia monooxygenase.*subunit a"],
    "nxrA": [
        r"nitrite oxidoreductase.*(subunit a|alpha)",
        r"nitrate reductase.*(narG|subunit alpha)",
    ],
    "nirK": [r"copper.containing nitrite reductase", r"\bnirk\b"],
    "nirS": [r"cytochrome cd1.*nitrite reductase", r"\bnirs\b"],
    "nosZ": [r"nitrous.oxide reductase", r"\bnosz\b"],
}

CACHE_DIR = Path("/tmp/gtdb_cache")


# ---------------------------------------------------------------------------
# Core parsing (unit-tested)
# ---------------------------------------------------------------------------


def _compile_patterns() -> dict[str, list[re.Pattern]]:
    return {
        gene: [re.compile(p, re.IGNORECASE) for p in pats]
        for gene, pats in NCYCLE_PRODUCT_PATTERNS.items()
    }


_COMPILED = _compile_patterns()


def detect_genes_from_gff_text(gff_text: str) -> set[str]:
    """Scan a GFF3 text blob for N-cycle product strings and return hit gene set."""
    hits: set[str] = set()
    for line in gff_text.splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 9:
            continue
        attrs = parts[8]
        # Extract product= value (up to ';' or end of line)
        m = re.search(r"product=([^;]+)", attrs)
        if not m:
            continue
        product = m.group(1)
        for gene, patterns in _COMPILED.items():
            if gene in hits:
                continue
            for pat in patterns:
                if pat.search(product):
                    hits.add(gene)
                    break
    return hits


def build_matrix(
    genera: list[str],
    genes: list[str],
    hits_by_genus: dict[str, Iterable[str]],
) -> pd.DataFrame:
    """Assemble a genera x genes presence/absence DataFrame of 0/1 ints."""
    rows = []
    for g in genera:
        ghits = set(hits_by_genus.get(g, set()))
        rows.append({gene: int(gene in ghits) for gene in genes})
    return pd.DataFrame(rows, index=genera, columns=genes)


# ---------------------------------------------------------------------------
# Genome fetching
# ---------------------------------------------------------------------------


def _run_datasets_download(taxon: str, zip_path: Path, reference_only: bool) -> bool:
    """Invoke `datasets` CLI; return True if zip was produced.

    ``taxon`` may be a genus or a species name; the command is the same.
    """
    cmd = [
        "datasets",
        "download",
        "genome",
        "taxon",
        taxon,
        "--assembly-source",
        "RefSeq",
        "--include",
        "gff3",
        "--filename",
        str(zip_path),
    ]
    if reference_only:
        cmd.insert(5, "--reference")
    log.info("datasets: %s", " ".join(cmd))
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    except Exception as e:  # noqa: BLE001
        log.warning("datasets invocation failed for %s: %s", genus, e)
        return False
    if r.returncode != 0:
        log.warning(
            "datasets nonzero exit for %s (ref=%s): %s",
            taxon,
            reference_only,
            r.stderr.strip()[:300],
        )
        return False
    return zip_path.exists() and zip_path.stat().st_size > 0


def fetch_genome_gffs(genus: str, cache_dir: Path) -> tuple[list[Path], list[str]]:
    """Fetch RefSeq genome GFF(s) for a genus. Returns (gff_paths, accessions).

    Strategy:
      1. If a species-level override is defined for the genus, query that
         species directly with ``--reference``.
      2. Else try ``--reference`` at the genus level.
      3. If that returns nothing, try the species override (if defined).
      4. Never fall back to the unconstrained genus-level (all assemblies)
         download — those bundles can be hundreds of MB and truncate.

    All ``genomic.gff`` files inside the downloaded bundle are returned; the
    caller is expected to union hits across them (represents "the genus has
    capacity to perform the function in at least one reference strain").
    """
    gdir = cache_dir / genus
    gdir.mkdir(parents=True, exist_ok=True)

    # Reuse cached GFFs if present
    existing = list(gdir.rglob("genomic.gff"))
    if existing:
        accs = [a for a in (_accession_from_path(p) for p in existing) if a]
        log.info("cache hit for %s -> %d GFF(s) (%s)", genus, len(existing), ",".join(accs[:3]))
        return existing, accs

    species = SPECIES_OVERRIDES.get(genus)
    zip_path = gdir / "dataset.zip"
    attempts: list[tuple[str, bool]] = []
    if species:
        attempts.append((species, True))  # species + --reference first
    attempts.append((genus, True))        # genus + --reference
    if species and not any(a[0] == species for a in attempts[:1]):
        attempts.append((species, True))

    got = False
    for taxon, ref_only in attempts:
        if zip_path.exists():
            zip_path.unlink()
        if _run_datasets_download(taxon, zip_path, reference_only=ref_only):
            got = True
            log.info("fetched %s via taxon=%s ref=%s", genus, taxon, ref_only)
            break

    if not got:
        log.warning("no genome recoverable for %s (tried %s)", genus, attempts)
        return [], []

    extract_dir = gdir / "extracted"
    if extract_dir.exists():
        shutil.rmtree(extract_dir)
    extract_dir.mkdir()
    try:
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(extract_dir)
    except zipfile.BadZipFile:
        log.warning("bad zip for %s", genus)
        return [], []

    gffs = sorted(extract_dir.rglob("genomic.gff"))
    if not gffs:
        log.warning("no genomic.gff inside zip for %s", genus)
        return [], []
    accs = [a for a in (_accession_from_path(p) for p in gffs) if a]
    return gffs, accs


def _accession_from_path(gff: Path) -> str | None:
    for part in gff.parts[::-1]:
        if part.startswith("GCF_") or part.startswith("GCA_"):
            return part
    return None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> int:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    config.AA_DIR.mkdir(parents=True, exist_ok=True)

    hits_by_genus: dict[str, set[str]] = {}
    source_rows: list[dict] = []

    for genus in TARGET_GENERA:
        log.info("=== %s ===", genus)
        gff_paths, accessions = fetch_genome_gffs(genus, CACHE_DIR)
        if not gff_paths:
            hits_by_genus[genus] = set()
            source_rows.append(
                {
                    "genus": genus,
                    "accession": "",
                    "n_gffs": 0,
                    "status": "no_genome_recovered",
                    "hits": "",
                }
            )
            continue
        union_hits: set[str] = set()
        for gff in gff_paths:
            try:
                text = gff.read_text(errors="replace")
            except Exception as e:  # noqa: BLE001
                log.warning("cannot read %s: %s", gff, e)
                continue
            union_hits |= detect_genes_from_gff_text(text)
        hits_by_genus[genus] = union_hits
        log.info(
            "%s hits (union over %d GFF): %s",
            genus, len(gff_paths), sorted(union_hits),
        )
        source_rows.append(
            {
                "genus": genus,
                "accession": ";".join(accessions[:5]),
                "n_gffs": len(gff_paths),
                "status": "ok",
                "hits": ";".join(sorted(union_hits)),
            }
        )

    matrix = build_matrix(TARGET_GENERA, NCYCLE_GENES, hits_by_genus)
    matrix.index.name = "genus"
    matrix.to_csv(config.AA_GTDB_MATRIX_CSV)
    log.info("wrote %s", config.AA_GTDB_MATRIX_CSV)

    sources = pd.DataFrame(source_rows)
    sources_path = config.AA_DIR / "gtdb_ncycle_matrix_sources.csv"
    sources.to_csv(sources_path, index=False)
    log.info("wrote %s", sources_path)

    log.info("\n%s", matrix.to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
