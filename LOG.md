# LOG

## 2026-10-18 C3 control arm (authors' WGCNA, paper settings)
- Installed BiocManager, WGCNA, dplyr, limma, impute, preprocessCore, GO.db, AnnotationDbi; library(WGCNA) loads. sessionInfo in env/sessionInfo.txt.
- Copied JUMPn/JUMPnc to scratch/JUMPnc; edits in the copy only: line 8 workingDir ".", sizeGrWindow lines commented, shipped outputs renamed shipped_*. Ran with Rscript (run.log). JUMPn/ and JUMPp/ untouched.
- No seed set: WGCNA TOM/hclust/cutreeDynamic have no RNG. Not tuned.
- Metric: ARI (CLAUDE.md default); Reproduced = ARI >= 0.95 vs shipped output. Rows aligned by order (shipped files have no accession column; checked intensities identical to the input, TRUE).
- Result: shipped vs rerun ARI 0.9993 (2 proteins differ); paper vs rerun 0.9507.
- Surprise/finding: paper labels (WP-C1..5) only map approximately onto WGCNA colours: WP-C1/C2 are not cleanly separated by blue/yellow (60 WP-C1 in yellow, 22 WP-C2 in blue). Paper-vs-shipped ARI is 0.9507, barely above the 0.95 threshold, so the paper's final labels differ somewhat from this script's shipped output (the script in the repo may not be exactly what produced the figure). Reported as found, not adjusted.
- Unexpected: shipped output has no protein accessions; script writes no accession column either.

## 2026-10-18 00_sample_key (sample/channel mapping across all files)
Script: analysis/00_sample_key/make_sample_key.R (no RNG). Outputs: sample_key.csv, mapping_evidence.csv (112 rows), input_vs_channel_spearman.csv, transgene_check.csv. sessionInfo: env/sessionInfo_00_sample_key.txt.
Metrics set before running: exact value identity (rounded 2 dp), constant-ratio fraction (|ratio/median-1|<1e-4), Spearman as fallback; accept a mapping only if best match is one-to-one and the runner-up is clearly lower.

VERIFIED BY DATA
- whl_hgg_v5 columns == hgg_whole_v1 channels sig126..sig131 in order: 13,419/13,419 values exactly identical per column, runner-up 0. (v5 has no accession, only Gene; v1 is the bridge.)
- pho_hgg_uniq_v5 columns == raw site file sig126..131 in order, exact; only ~50% of v5 values appear verbatim in the raw file (rest are presumably collapsed/averaged duplicates), runner-up 0, so the column mapping is unambiguous. Not checked: how rows were collapsed.
- jpeng_hgg_phosphosite_matrix Data1..10 == raw site sig126..131, all 69,808 rows exactly identical. Data_k = channel k.
- HGG_whl_input -> channels: input intensities are not raw v1 values. For the Cortex columns they are exactly identical (98.5% of 4,702 proteins; rest rounding). For PD and TRK columns, input = channel x a constant per group (PD 0.8857, TRK 0.8433, Cortex 1.0000) for 99.9% of proteins; runner-up <=0.002. Spearman agrees (best rho >=0.9997, runner-up <=0.9984). Result: PD1-4 = sig126,127N,127C,128N; TRK1-3 = sig128C,129N,129C; Cortex1-3 = sig130N,130C,131.
- ikap_hgg_v5 vs IKAP_result: all 90 kinases identical per column. The repeated header "PDGFRA2" (3rd column of ikap_hgg_v5) is IKAP_result PDGFRA3. IKAP_result is ordered CNTRL, PDGFRA, NTRK (v5 files: PDGFRA, NTRK, CNTRL); map by name, never by position.
- Normalised matrices: Data8-10 have column mean ~0 and Data1-7 ~ +1 (log2 vs controls), so Data8-10 = controls.
- Transgenes: whole-proteome GFP high in ch1-4, mCherry high in ch5-7, PDGFRA_D842V protein highest in ch1-4 and lowest in ch8-10; input file shows the same by label (PDGFRA_D842V ~5x in PD, TPM3_NTRK1 ~5x in TRK).
  Site level (raw file): PDGFRA Y849 3.5x in ch1-4 vs rest; NTRK1 S501 22x in ch5-7 vs rest. Other transgene sites are weak or non-specific: PDGFRA S767 peaks in ch5-6 (NTRK-like), TPM3-part sites S180/Y179 peak in ch8-10 (probably shared with endogenous Tpm3 peptides). Only 2/10 sites discriminate cleanly.

FINDING: the channel layout in the task text ("channels 1-3 cortex, 4-7 PDGFRA, 8-10 NTRK", from Fig 1c/d) does NOT match the data. All files agree on PDGFRA = ch1-4, NTRK = ch5-7, control cortex = ch8-10. HGG_whl_input only lists the columns in a different order (Cortex, PD, TRK). Whether the paper's figure uses its own ordering is not checked (paper not in repo); reported, not adjusted.
FINDING: HGG_whl_input is hgg_whole_v1 rescaled per group (PD x0.886, TRK x0.843), which shifts tumour-vs-cortex log2FC by -0.18 and -0.25. Does not change within-tumour comparisons (NTRK vs PDGFRA shifts by 0.07 log2) but means input and v1 are not interchangeable for FC.
Surprise (my own bug, fixed): first Spearman matrix was transposed (rows channels, columns input), which briefly suggested a different mapping; caught because the input file's own transgene rows contradicted it.

ASSUMPTIONS (not proven by data)
- ikap_hgg_v5 / IKAP_result tumour columns follow channel order (PDGFRA1-4 = ch1-4, NTRK1-3 = ch5-7). Only control columns are supported by data (IKAP vs substrate means: CNTRL1-3 best match Data8-10, rho 0.79/0.61/0.84); tumour columns give weak rho 0.2-0.4, not discriminative. Support is the header text plus the same order as whl_hgg_v5.
- hgg_phosphosite_matrix_normalized_whole (step 4 output) keeps Data1..10 = channels 1..10. Correlation with its input matrix is weak for Data1-7 (rho 0.16-0.23, not diagnostic because step 4 subtracts protein-level change); strong for Data8-10 (0.68-0.83).
- hgg_whole_normalized Data1..10 = channels 1..10 by name and control centring only.
- Sample naming only: PDGFRA_1a / 1b are mapped to PD1 / PD2 by position and by the exact data match in ikap/v5. The names suggest ch1 and ch2 might be two channels of one tumour (technical pair) while PDGFRA_2, _3 are separate tumours; the repo does not say. If true, PDGFRA has 3 independent tumours, not 4. Unresolved.
- Cortex1-3 naming order within the group follows the data match; their biological identity (independent cortex samples) is assumed.

AMBIGUITIES FLAGGED
- ikap_hgg_v5 repeats header "PDGFRA2" (3rd column = PDGFRA3). KinaseActivityHeatmap/PDGFRA.txt lists PDGFRA1,2,2,4 and NTRK_Filtered.txt lists NTRK2,3 (known issue 2): by name that selects PDGFRA2 twice and omits PDGFRA3 and NTRK1.
- Column order differs between files: IKAP_result (CNTRL, PDGFRA, NTRK) and HGG_whl_input (Cortex, PD, TRK) vs v5 files / raw (PDGFRA, NTRK, CNTRL).
- hgg_whole_v1 transgene protein "cu|TMP3_NTRK1" is flat across channels (~46-55k), while the input row "cu|TPM3_NTRK1|TMP3_NTRK1" (different accession spelling) is 5x in TRK; accession mismatch is why 1 of 4,703 input proteins did not match v1. Not investigated further.
- NTRK-transgene site names are spelled TMP3_NTRK1 in the site file but TPM3_NTRK1 elsewhere.

## 2026-10-18 PRE-REGISTRATION: C3 compare arm (written BEFORE running c3_compare.R)
- Q-A: Which settings best recover the paper's WP-C labels? Metrics: ARI and % identical (each cluster mapped to its majority paper label). Ranking only, no pass/fail.
- Q-B: Does "a cluster rising cortex < PDGFRA < NTRK" survive? Per variant, per-cluster mean of per-protein z-scored intensities by group. Claim holds if some cluster has Cortex < PD < TRK AND Jaccard >= 0.8 with the paper's WP-C1 members (all-features variant: members restricted to the 4,703 DE proteins).
- Q-C: Stability. Leave-one-sample-out (10 runs, control settings): ARI of each run's labels vs full-data control labels. Stable if ARI >= 0.8 in >= 9 of 10 runs.
- Variants fixed in advance: control (power 16, deepSplit 2, merge 0.20, minModuleSize 100); grid merge {0.15,0.20,0.25} x deepSplit {1,2,3}; log2 input; unscaled input (per-group rescaling removed); hierarchical Ward.D2 on 1-Pearson, k=5; k-means k=5 (seed 1, nstart 25) with silhouette k=2..8; all-features (top 8,000 proteins by MAD of log2, per-sample scale factors estimated from the data); leave-one-out.
- Control check: the rerun must reproduce the shipped labels (ARI >= 0.999) or the script stops.

## 2026-10-18 C3 compare arm RESULTS (c3_compare.R; summary_variants.csv, results/metrics.csv)
Verified by data: control reproduces shipped labels (ARI 0.9993); channel-mapping Spearman >= 0.9997; rescaling factors PD 0.8857, TRK 0.8433, Cortex 1.0.
- Q-A: no variant recovers the paper's labels better than the control (ARI 0.9507, 96.9% identical). Merge 0.15 = control (ARI 1.0000), so it does not explain the 147 differing proteins; kME reassignment also rejected. Residual unexplained.
- Q-B: holds in 14/14 variants (Jaccard with WP-C1 0.80-0.92; lowest k-means 0.801, log2 0.810, unscaled 0.828).
- Q-C: 10/10 leave-one-out runs ARI >= 0.8 (range 0.869-0.975; lowest dropping PD2, PD3, TRK2).
- Structure (cross-tabs): five clusters are two pairs (WP-C1+C2, WP-C3+C4); merge 0.25 gives 3 clusters; unscaled input merges 306/307 WP-C2 proteins into WP-C1; silhouette peaks at k=2.
- Exploratory, not pre-registered: 463 proteins (control green+black) have PDGFRA >= NTRK.
- Limits: n=10; DE pre-selection (circularity); ARI vs paper = similarity to their method; LOO only for control settings; all-features includes 4,327/4,703 DE proteins.

## 2026-10-18 PRE-REGISTRATION: C7 part A (written BEFORE running c7_null_sensitivity.R)
Disclosure: informal Python checks earlier in this project (thresholds 0.5 and 1.0) already suggested p > 0.05 under selection-matched nulls. This entry formalises them; the tau grid and decision rule below were fixed before this run.
Data: JUMPn/JUMPna/.RData, z = log2FC NTRK vs PDGFRA of 28,527 phosphosites (no site names; the 22 pathway proteins and their Fi are NOT in the repo). Observed a(P) = 1.45 is taken from the paper. Paper Methods say Fi is relative to cortex; the repo says NTRK vs PDGFRA (logged discrepancy).
- Q-A: reproduce the authors' null (22 values drawn without replacement from all sites; positions 4, 8, 21 sign-flipped; 10,000 draws; seed 1). Expect p about 0.002-0.003 (paper 0.003; saved workspace 0.0022).
- Q-B: sensitivity of p to the null: draw from sites with |log2FC| >= tau, tau = 0, 0.1, ..., 1.5. Report p(tau), the 95th-percentile critical a(P), and the first tau where p >= 0.05. Note: the paper's phosphosite DE rule (>= 2-fold, |log2FC| >= 1, in at least one of three pairwise comparisons) is not a threshold on NTRK-vs-PDGFRA alone, so tau = 1.0 overstates the selection and tau = 0 understates it; the true selection lies between.
- Q-C: restate a(P) = 1.45 as mean signed log2FC per protein = 1.45 / sqrt(22).
- Decision rule: robust if p < 0.05 for all tau in [0, 1.0]; method-dependent if p crosses 0.05 within [0, 1.0]; not supported if p >= 0.05 at tau = 0.
- Part B (needs Supplementary Data 4h for the 22 proteins' Fi and Ci, and Supplementary Data 2 for DE flags): reproduce a(P), sample-label permutation null (35 splits; minimum one-sided p = 1/35), bootstrap CI over samples. To be pre-registered when the data are downloaded.

## 2026-10-18 PRE-REGISTRATION: C5 kinase activity (written BEFORE running any C5 script)
Prior knowledge: the coded sample lists (PDGFRA.txt, NTRK_Filtered.txt in KinaseActivityHeatmap) were read from the repo before this plan was written, so Q-A is a confirmatory check, not a blind test. The rules below were fixed before the C5 repo run.
Data: IKAP_result_v2_20160708.txt (90 kinases x 10 samples), ALL_LIMMA_ANOVA_HGG_Result_20160708.txt (published F, P, adj.P and contrasts group2-group1 = NTRK-PDGFRA, group3-group2 = CNTRL-NTRK, group3-group1 = CNTRL-PDGFRA), the normalised phospho matrix (984 sites, KINASE column with one or more kinases; Data1-4 PDGFRA, 5-7 NTRK, 8-10 cortex; Data8-10 verified by data, Data1-7 assumed). IKAP itself is not in the repo.
- Q-A Reproduction: limma (lmFit + contrasts.fit + eBayes, three groups) on IKAP_result_v2 with the coded lists (PDGFRA1, PDGFRA2, PDGFRA2, PDGFRA4; NTRK2, NTRK3; CNTRL1-3) reproduces the published contrasts, F, P and adj.P. Pass if max abs difference < 1e-6; otherwise report the maximum differences.
- Q-B Sample effect: same analysis with all 10 samples (PDGFRA1-4, NTRK1-3, CNTRL1-3). Report (i) kinases with F-test adj.P < 0.05 under each, (ii) kinases with NTRK-vs-PDGFRA contrast adj.P < 0.05 (BH over 90 kinases) under each, and the Jaccard overlap of each pair of sets. Rule: robust if Jaccard >= 0.8 for both; method-dependent otherwise; not supported if the NTRK-vs-PDGFRA set is empty with all 10 samples.
- Q-C: record counts of published adj.P < 0.05, 0.01, 0.001 (62, 48, 30) against the paper's 41. No resolution without Supplementary Data 4e (deferred).
- Q-D Alternative score: substrate average = per sample, mean of Data values over a kinase's substrate sites (a site counts for every kinase listed), kinases with >= 3 sites. Compare with IKAP for kinases in both: (i) per-kinase Pearson correlation across the 10 samples (median, fraction >= 0.8); (ii) Spearman across kinases of the NTRK-minus-PDGFRA mean difference; (iii) Jaccard of kinases with NTRK-vs-PDGFRA contrast adj.P < 0.05 (limma, all 10 samples). Agree = >= 0.8 (CLAUDE.md).
- Q-E Stability: leave-one-sample-out (10 runs, all-10-sample analysis): Jaccard of the NTRK-vs-PDGFRA set vs the full set; stable if Jaccard >= 0.8 in >= 9 of 10 runs.
- Q-F Descriptive only: ranks of PRKAA1, PRKAA2, EEF2K (CRISPR-essential) among the kinases by NTRK-vs-PDGFRA effect and by F under each method; 3 positives, no statistics.
- Out of scope: AGC/CMGC superfamily test (not fully specified in the paper); IKAP rerun.
