"""
Configuration file for ISS Chile Pepper Microbiome Analysis

This module provides centralized configuration for all analysis scripts,
ensuring portability and reproducibility across different computing environments.

Usage:
    import config
    table = pd.read_csv(config.FEATURE_TABLE_CLEAN, sep='\t')

Last Updated: 2025-02-09 (Peer Review Phase 2)
"""

from pathlib import Path
import sys

# ============================================================================
# Directory Structure
# ============================================================================

# Base directory (repository root)
BASE_DIR = Path(__file__).parent.resolve()

# Version-specific data directory
VERSION_DIR = BASE_DIR / "version-2_integrated"

# Output directories
ARCHIVE_DIR = BASE_DIR / "archive"
TESTS_DIR = BASE_DIR / "tests"

# Current Biology Correspondence output directory
CURRENT_BIOLOGY_DIR = BASE_DIR / "current_biology"

# ============================================================================
# Data Files - Feature Tables
# ============================================================================

# Combined dataset (ISS + Terrestrial)
EXPORTED_TABLE_CLEAN = VERSION_DIR / "exported_table_clean"
FEATURE_TABLE_CLEAN = EXPORTED_TABLE_CLEAN / "feature-table.tsv"

# ISS samples only
EXPORTED_TABLE_SPACE = VERSION_DIR / "exported_table_space"
FEATURE_TABLE_SPACE = EXPORTED_TABLE_SPACE / "feature-table.tsv"

# Terrestrial samples only
EXPORTED_TABLE_TERR = VERSION_DIR / "exported_table_terrestrial"
FEATURE_TABLE_TERR = EXPORTED_TABLE_TERR / "feature-table.tsv"

# ============================================================================
# Data Files - Taxonomy
# ============================================================================

EXPORTED_TAXONOMY = VERSION_DIR / "exported_taxonomy"
TAXONOMY_FILE = EXPORTED_TAXONOMY / "taxonomy.tsv"

# ============================================================================
# Data Files - Metadata
# ============================================================================

INTEGRATED_METADATA = VERSION_DIR / "integrated_metadata.tsv"

# ============================================================================
# Data Files - Diversity Metrics
# ============================================================================

# Exported diversity (TSV format, ready for analysis)
EXPORTED_DIVERSITY = VERSION_DIR / "exported_diversity"
SHANNON_VECTOR = EXPORTED_DIVERSITY / "shannon.tsv"
OBSERVED_FEATURES = EXPORTED_DIVERSITY / "observed_features.tsv"

# QIIME2 artifacts (for reference, not typically used directly in Python)
CORE_METRICS = VERSION_DIR / "core-metrics-results"

# ============================================================================
# Data Files - PICRUSt2 Functional Predictions
# ============================================================================

PICRUST2_DIR = VERSION_DIR / "picrust2_out"

# Pathway abundance (unstratified)
PICRUST2_PATHWAYS = PICRUST2_DIR / "pathways_out" / "path_abun_unstrat.tsv.gz"

# Enzyme predictions
PICRUST2_EC = PICRUST2_DIR / "EC_metagenome_out" / "pred_metagenome_unstrat.tsv.gz"

# KEGG ortholog predictions
PICRUST2_KO = PICRUST2_DIR / "KO_metagenome_out" / "pred_metagenome_unstrat.tsv.gz"

# Pathway descriptions
PICRUST2_PATH_DESC = PICRUST2_DIR / "pathways_out" / "path_abun_descriptions.tsv"

# ============================================================================
# Output Files - Figures
# ============================================================================

# Main figures
MAIN_FIG_PREFIX = VERSION_DIR / "Fig"

# Supplementary figures
SUPP_FIG_PREFIX = VERSION_DIR / "SuppS"

# Network visualizations
NETWORK_PREFIX = VERSION_DIR / "Network"

# ============================================================================
# Analysis Parameters
# ============================================================================

# Random seed for reproducibility (CRITICAL for peer review)
RANDOM_SEED = 42

# Network analysis parameters
NETWORK_PARAMS = {
    "r_threshold": 0.4,        # Spearman |r| > 0.4, raw p < 0.05 (genus-level; Barberán et al., 2012)
    "p_threshold": 0.05,       # Raw p-value threshold (FDR not applied at genus level; see 10_network_analysis.py)
    "top_n_asvs": 100,        # Top N most abundant ASVs for network construction
    "min_prevalence": 0.05,   # Minimum proportion of samples ASV must appear in
}

# Diversity analysis parameters
DIVERSITY_PARAMS = {
    "rarefaction_depth": 1000,   # From 06_plot_depth.py analysis; retains ~86% of samples
    "permanova_permutations": 999,
}

# PICRUSt2 filtering parameters
PICRUST2_PARAMS = {
    "min_pathway_abundance": 0.01,  # Minimum relative abundance for pathways
    "top_pathways": 30,            # Number of top pathways to display
}

# Statistical testing parameters
STATS_PARAMS = {
    "alpha": 0.05,                 # Significance level
    "fdr_method": "fdr_bh",       # Benjamini-Hochberg FDR correction
    "mannwhitneyu_alternative": "two-sided",
}

# ============================================================================
# Validation Functions
# ============================================================================

def validate_paths(verbose=True):
    """
    Validate that all required input files exist.

    Args:
        verbose (bool): If True, print detailed validation results

    Returns:
        bool: True if all critical files exist, False otherwise

    Raises:
        SystemExit: If critical files are missing (when called directly)
    """
    critical_files = {
        "Feature Table (Clean)": FEATURE_TABLE_CLEAN,
        "Feature Table (Space)": FEATURE_TABLE_SPACE,
        "Feature Table (Terrestrial)": FEATURE_TABLE_TERR,
        "Taxonomy": TAXONOMY_FILE,
        "Metadata": INTEGRATED_METADATA,
        "Shannon Diversity": SHANNON_VECTOR,
    }

    optional_files = {
        "PICRUSt2 Pathways": PICRUST2_PATHWAYS,
        "PICRUSt2 EC": PICRUST2_EC,
        "PICRUSt2 KO": PICRUST2_KO,
    }

    missing_critical = []
    missing_optional = []

    # Check critical files
    for name, path in critical_files.items():
        if not path.exists():
            missing_critical.append(f"  - {name}: {path}")
        elif verbose:
            print(f"✅ {name}")

    # Check optional files
    for name, path in optional_files.items():
        if not path.exists():
            missing_optional.append(f"  - {name}: {path}")
        elif verbose:
            print(f"✅ {name}")

    # Report results
    if missing_critical:
        print("\n❌ ERROR: Missing required files:")
        for item in missing_critical:
            print(item)
        print("\n💡 Hint: Run previous analysis steps first (see README.md workflow)")
        if not verbose:  # If called as validation in a script
            sys.exit(1)
        return False

    if missing_optional and verbose:
        print("\n⚠️  Warning: Missing optional files (PICRUSt2 analysis):")
        for item in missing_optional:
            print(item)
        print("   → These are only needed for functional analysis (Main Fig 4, 6)")

    if verbose:
        print("\n✅ All critical input files validated successfully!")

    return True


def validate_output_dirs(create=True):
    """
    Validate that output directories exist, optionally creating them.

    Args:
        create (bool): If True, create missing directories

    Returns:
        bool: True if all directories exist or were created
    """
    required_dirs = [
        VERSION_DIR,
        EXPORTED_TABLE_CLEAN,
        EXPORTED_TABLE_SPACE,
        EXPORTED_TABLE_TERR,
        EXPORTED_TAXONOMY,
        CORE_METRICS,
    ]

    for dir_path in required_dirs:
        if not dir_path.exists():
            if create:
                dir_path.mkdir(parents=True, exist_ok=True)
                print(f"✅ Created directory: {dir_path}")
            else:
                print(f"❌ Missing directory: {dir_path}")
                return False

    return True


def get_figure_path(figure_type, figure_number, description="", ext="png"):
    """
    Generate standardized figure file paths.

    Args:
        figure_type (str): "Main" or "Supp"
        figure_number (int): Figure number
        description (str): Optional description to append to filename
        ext (str): File extension (default: png)

    Returns:
        Path: Full path to figure file

    Example:
        >>> get_figure_path("Main", 1, "AlphaDiversity")
        PosixPath('.../version-2_integrated/Main_Fig1_AlphaDiversity.png')
    """
    prefix = MAIN_FIG_PREFIX if figure_type.lower() == "main" else SUPP_FIG_PREFIX

    if description:
        filename = f"{prefix.name}{figure_number}_{description}.{ext}"
    else:
        filename = f"{prefix.name}{figure_number}.{ext}"

    return VERSION_DIR / filename


# ============================================================================
# Session Information
# ============================================================================

def get_session_info():
    """
    Get current Python and package versions for reproducibility.

    Returns:
        dict: Dictionary containing version information
    """
    import pandas as pd
    import numpy as np
    import scipy
    import matplotlib
    import seaborn as sns
    from datetime import datetime

    info = {
        "timestamp": datetime.now().isoformat(),
        "python_version": sys.version,
        "python_executable": sys.executable,
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "matplotlib": matplotlib.__version__,
        "seaborn": sns.__version__,
    }

    return info


def log_session_info(output_file="SESSION_INFO.txt"):
    """
    Write session information to file for reproducibility documentation.

    Args:
        output_file (str): Path to output file
    """
    info = get_session_info()

    with open(BASE_DIR / output_file, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("ISS Chile Pepper Microbiome Analysis - Session Information\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"Analysis Timestamp: {info['timestamp']}\n\n")
        f.write(f"Python Version:\n{info['python_version']}\n\n")
        f.write(f"Python Executable: {info['python_executable']}\n\n")
        f.write("Package Versions:\n")
        f.write(f"  - pandas: {info['pandas']}\n")
        f.write(f"  - numpy: {info['numpy']}\n")
        f.write(f"  - scipy: {info['scipy']}\n")
        f.write(f"  - matplotlib: {info['matplotlib']}\n")
        f.write(f"  - seaborn: {info['seaborn']}\n")
        f.write("\n" + "=" * 70 + "\n")

    print(f"✅ Session information written to: {output_file}")


# ============================================================================
# Self-Test (run when executed directly)
# ============================================================================

# ============================================================================
# Acta Astronautica redesign (2026-04-15) — external datasets and outputs
# ============================================================================

# Veggie lettuce cross-mission validation (Prediction P3)
VEGGIE_LETTUCE_OSDS = ["OSD-267", "OSD-269"]
VEGGIE_LETTUCE_DIR = VERSION_DIR / "external" / "veggie_lettuce"

# ISS environmental source pool (Fig 4d)
ISS_ENV_OSD = "OSD-65"
ISS_ENV_DIR = VERSION_DIR / "external" / "iss_environment"

# Outputs for the Acta Astronautica manuscript
AA_DIR = BASE_DIR / "acta_astronautica"
AA_CONSENSUS_CSV = AA_DIR / "consensus_function_log2fc.csv"
AA_GTDB_MATRIX_CSV = AA_DIR / "gtdb_ncycle_matrix.csv"
AA_VEGGIE_RESULTS_CSV = AA_DIR / "veggie_lettuce_results.csv"
AA_ISS_ENV_SOURCES_CSV = AA_DIR / "iss_env_sources.csv"

AA_FIG_DIR = BASE_DIR / "figures" / "aa"


if __name__ == "__main__":
    print("=" * 70)
    print("ISS Chile Pepper Analysis - Configuration Validation")
    print("=" * 70)
    print(f"\nBase Directory: {BASE_DIR}")
    print(f"Version Directory: {VERSION_DIR}")
    print(f"Random Seed: {RANDOM_SEED}")
    print("\n" + "=" * 70)
    print("Validating Input Files:")
    print("=" * 70 + "\n")

    validate_paths(verbose=True)

    print("\n" + "=" * 70)
    print("Generating Session Info:")
    print("=" * 70 + "\n")

    log_session_info()
