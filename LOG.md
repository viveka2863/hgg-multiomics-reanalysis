# LOG

## 2026-10-18 C3 control arm (authors' WGCNA, paper settings)
- Installed BiocManager, WGCNA, dplyr, limma, impute, preprocessCore, GO.db, AnnotationDbi; library(WGCNA) loads. sessionInfo in env/sessionInfo.txt.
- Copied JUMPn/JUMPnc to scratch/JUMPnc; edits in the copy only: line 8 workingDir ".", sizeGrWindow lines commented, shipped outputs renamed shipped_*. Ran with Rscript (run.log). JUMPn/ and JUMPp/ untouched.
- No seed set: WGCNA TOM/hclust/cutreeDynamic have no RNG. Not tuned.
- Metric: ARI (CLAUDE.md default); Reproduced = ARI >= 0.95 vs shipped output. Rows aligned by order (shipped files have no accession column; checked intensities identical to the input, TRUE).
- Result: shipped vs rerun ARI 0.9993 (2 proteins differ); paper vs rerun 0.9507.
- Surprise/finding: paper labels (WP-C1..5) only map approximately onto WGCNA colours: WP-C1/C2 are not cleanly separated by blue/yellow (60 WP-C1 in yellow, 22 WP-C2 in blue). Paper-vs-shipped ARI is 0.9507, barely above the 0.95 threshold, so the paper's final labels differ somewhat from this script's shipped output (the script in the repo may not be exactly what produced the figure). Reported as found, not adjusted.
- Unexpected: shipped output has no protein accessions; script writes no accession column either.
