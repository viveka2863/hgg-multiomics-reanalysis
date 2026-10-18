# C7: PI3K-AKT pathway activity score a(P)

- **Kept:** The paper's formula a(P) = sum(Ci*Fi)/sqrt(k), k = 22, observed value 1.45 (taken from the paper), and the authors' bootstrap null (22 sites drawn from all 28,527 phosphosites, 10,000 draws).
- **Compared:** The same bootstrap with the pool restricted to sites with |log2FC NTRK vs PDGFRA| >= tau, tau = 0 to 1.5. Pre-registered in LOG.md before the run.
- **Added:** Restatement of a(P) = 1.45 as a mean signed change of 0.31 log2 units (about 1.24-fold) per protein; critical a(P) needed for p < 0.05 at each tau.
- **Result:** The authors' p-value reproduces (p = 0.0019 vs 0.003 reported). It rises above 0.05 once the null is drawn from sites changing by >= 0.4 log2 (p = 0.057) and reaches 0.17 at tau = 1.0. The observed score is above the null mean at every tau.
- **Verdict:** method-dependent (pre-registered rule). The direction (NTRK > PDGFRA) is consistent across nulls, but its statistical significance depends on how the comparison sites are selected. Limits: the 22 proteins and their values are not in the repo, so a selection-matched null and a sample-label permutation (35 splits, minimum one-sided p = 1/35) await Supplementary Data 4h and 2 (part B).

Details: improved/ (script, table, figure), results/metrics.csv.
