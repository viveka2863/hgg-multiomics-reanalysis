#!/bin/bash
# C3 control arm: authors' WGCNA, paper settings, run in a scratch copy. Run from repo root.
# WGCNA/dynamicTreeCut are deterministic (no RNG), so no seed is needed.
set -euo pipefail
rm -rf scratch/JUMPnc && mkdir -p scratch && cp -R JUMPn/JUMPnc scratch/JUMPnc
cd scratch/JUMPnc
mv hgg_whl_cluster_output.txt shipped_hgg_whl_cluster_output.txt
mv hgg_whl_eigeneME_output.txt shipped_hgg_whl_eigeneME_output.txt
# minimal edits: working dir (line 8) and no graphics windows
sed -i '' -e '8s|.*|workingDir = ".";|' -e 's|^sizeGrWindow|#sizeGrWindow|' wgcna_v1.0.R
Rscript wgcna_v1.0.R > run.log 2>&1
