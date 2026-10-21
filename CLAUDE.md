# Project context
Robustness audit of Wang et al. 2019 (Nat Commun 10:3718): which conclusions about PDGFRA- vs NTRK1-driven
high-grade glioma survive when the data mining method changes? Each claim (C0-C12) is handled as
Keep (paper's method, credited) / Compare (original vs alternatives) / Extend (new analysis),
and ends with a verdict: robust | method-dependent | not supported.
Primary contrast: NTRK vs PDGFRA (tumour vs cortex is confounded by cell composition).

# Rules
- NEVER modify JUMPn/ or JUMPp/ (authors' code). Run their scripts only in copies under scratch/ (gitignored).
- Our work goes in analysis/<claim>/{original,improved,extensions,compare}/ plus a README.md with five lines:
  Kept / Compared / Added / Result / Verdict.
- Control arm = the authors' method with minimal edits needed to run on macOS. Do NOT tune it.
- Every comparison includes the paper-settings arm, so differences can be attributed to the method.
- Set random seeds. Save sessionInfo() to env/. Keep committed outputs small (tables, figures; no big intermediates).
- Write metrics and thresholds in LOG.md BEFORE running a comparison. Log every decision and surprise in LOG.md.
- Commit after each claim with a clear message. Do not push unless asked. Do not change git config.
- If a result differs from the paper, report it as a finding. Do not adjust things until it matches.

# Metrics (defaults; edit before first comparison)
Feature lists: Jaccard. Rankings/scores: Spearman. Clusterings: adjusted Rand index (ARI).
Stability: leave-one-sample-out. Agree = Jaccard/ARI/Spearman >= 0.8. Stable = result kept in >= 9/10 runs.
Reproduced = ARI >= 0.95 against shipped output.

# Environment
macOS arm64, R 4.6.1, Java 11, Perl 5.34. Use Rscript. Commit email is already set (noreply).

# Data (all inside the repo)
- JUMPn/JUMPnc/HGG_whl_input.txt: 4,703 DE proteins; Cluster = paper's WP-C1..5; columns Cortex1-3, PD1-4, TRK1-3 (raw TMT intensities)
- JUMPn/JUMPnk/Integration/: whl_hgg_v5 (~13.4k proteins), pho_hgg_uniq_v5 (~32k sites, incl. contaminants), ikap_hgg_v5 (90 kinases)
- JUMPn/JUMPnk/PhosphoWholeNormalization/: Mouse_Kinase_Substrate.txt (~12.7k), normalised matrix (984 kinase-assigned sites)
- JUMPn/JUMPna/.RData: z = log2FC of 28,527 phosphosites (for the a(P) score)
- RNA-seq is NOT in the repo (GEO GSE114331 / Supplementary Data 3).
- Sample names differ across files (Cortex/PD/TRK, PDGFRA_1a.., Data1-10, sig126..). Build analysis/00_sample_key first.

# Known issues to audit (already found)
1. WGCNA script: merge height 0.2 (paper says 0.15), raw not log intensities, run only on pre-selected DE proteins, softPower 16.
2. Kinase limma sample lists use PDGFRA1,2,2,4 and NTRK2,3 (3 PDGFRA incl. duplicate, 2 NTRK, 3 control), not 4/3/3.
3. a(P) bootstrap null draws from all 28,527 sites (p ~0.002); nulls from comparably large changes give ~0.09-0.17.
   Permutation test with 4 vs 3 samples cannot give p below 1/35.
4. Enrichment script default FDR cutoff 0.2 (paper: 0.05). Only 984 sites (~3%) have a known kinase.
5. Methods say ANOVA + permutation for DE; the repo's statistics are limma-based.
6. C3 control arm: the authors' script reproduces its shipped output (ARI 0.9993) but matches the paper's WP-C labels
   only 96.9% (ARI 0.951). Reassigning proteins by highest eigengene correlation (kME) made it worse (ARI 0.904): rejected.
7. Shipped/rerun WGCNA outputs have no protein accessions; rows align to HGG_whl_input by order.
   Always attach Protein_Accession from the input and assert the intensity columns are identical.

# Results convention
Every comparison script (1) writes tidy CSVs with Protein_Accession (or site ID) into its claim folder,
(2) saves figures as PNG, (3) appends rows to results/metrics.csv with columns: claim, variant, metric, value, n.
Scripts are the source of truth; no notebooks for analysis.

# Evidence rule
Nothing is reported unless a script in this repo reproduces it. Checks done outside the repo (notes, sandbox, ad hoc)
are leads only. Findings read from the authors' code are labelled "read from code" and cite file and line.
