#!/bin/bash
# Rename R script outputs for Current Biology figure scheme
mkdir -p current_biology
cp version-2_integrated/SuppS2_PERMANOVA.png current_biology/SuppS1c_PERMANOVA.png 2>/dev/null
cp version-2_integrated/SuppS2_PERMANOVA.pdf current_biology/SuppS1c_PERMANOVA.pdf 2>/dev/null
cp version-2_integrated/SuppS4_BNTI_RCbray.png current_biology/SuppS2c_BNTI_RCbray.png 2>/dev/null
cp version-2_integrated/SuppS4_BNTI_RCbray.pdf current_biology/SuppS2c_BNTI_RCbray.pdf 2>/dev/null
echo "R script outputs copied to current_biology/"
