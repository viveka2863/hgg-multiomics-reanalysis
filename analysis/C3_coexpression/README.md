# C3: co-expression clusters (WGCNA)

- **Kept:** Authors' WGCNA settings (signed, power 16, deepSplit 2, merge 0.20, minModuleSize 100) as the control arm. It reproduces their shipped labels (ARI 0.9993) and matches the paper's WP-C labels for 96.9% of proteins (ARI 0.951).
- **Compared:** merge height {0.15, 0.20, 0.25} x deepSplit {1, 2, 3}; log2 input; unscaled input (authors' per-group rescaling removed); hierarchical (Ward.D2) and k-means at k = 5; WGCNA on the top 8,000 proteins; leave-one-sample-out. Pre-registered in LOG.md before the run.
- **Added:** Check of the undocumented per-group rescaling in HGG_whl_input (PD x0.886, TRK x0.843, Cortex x1); silhouette analysis showing k = 5 is not data-driven.
- **Result:** No variant matches the paper's labels better than the control. A cluster rising Cortex < PDGFRA < NTRK (WP-C1) is found in 14/14 variants (Jaccard 0.80-0.92) and all 10 leave-one-out runs (ARI 0.87-0.98). The number of clusters is not stable: merge 0.25 gives 3, log2 and unscaled give 4; five clusters are two pairs (WP-C1+C2, WP-C3+C4).
- **Verdict:** C3a (ordered tumour-up cluster, WP-C1): robust. C3b (exactly five clusters): method-dependent. Limits: n = 10, features pre-selected as differential, ARI vs paper measures similarity to their method.

Details: original/ (control arm), improved/ (script, tables, figures), results/metrics.csv.
