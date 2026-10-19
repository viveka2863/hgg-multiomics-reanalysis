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

## 2026-10-18 C7 part A RESULTS (c7_null_sensitivity.R; c7_null_sensitivity.csv, results/metrics.csv)
- Q-A reproduced: authors' null gives p = 0.0019 (MC SE 0.0004; paper 0.003, saved workspace 0.0022); observed a(P) = 1.45 is 2.87 null SDs above the null mean (0.086).
- Q-B: p crosses 0.05 between tau = 0.3 (p = 0.030) and 0.4 (p = 0.057, MC SE 0.002); p = 0.084 at 0.5 and 0.170 at 1.0. Critical a(P) for p < 0.05 rises from 0.87 (tau 0) to 1.49 (tau 0.4) and 2.43 (tau 1.0); observed 1.45. Pre-registered decision rule -> method-dependent. p is non-monotonic above tau ~1.1 (pools < 1,000 sites): not interpreted.
- Q-C: a(P) = 1.45 = mean signed log2FC 0.309 per protein (about 1.24-fold), not a 1.45-fold change.
- Direction: observed a(P) exceeds the null mean at every tau (z_obs 0.88-2.87).
- Limits: the 22 proteins/Fi/Ci are not in the repo, so selection could not be matched to the paper's DE rule (|log2FC| >= 1 in any of three comparisons); tau = 0 understates and tau = 1 overstates the selection. Null draws are single sites whereas observed Fi are per-protein averages (would narrow the null; untested). Part B needed for a matched null and a sample-label permutation.

## 2026-10-18 C5 RESULTS (analysis/C5_kinase_activity/c5_kinase_activity.R; CSVs in original/ and improved/, results/metrics.csv, env/sessionInfo_C5.txt)
Verified by data (script output):
- Q-A: limma on the 9 coded columns reproduces the published table exactly. Max abs diff: three contrasts <= 4.4e-15, AveExpr 4.7e-15, F 2.2e-12, P 4.9e-15, adj.P 4.9e-15 (max relative <= 1.0e-14). The alternative-coding fallback was not needed.
- Q-B: F adj.P<0.05: 62 coded, 63 all-10 (Jaccard 0.953). NTRK-vs-PDGFRA adj.P<0.05 (BH over 90): 9 coded, 11 all-10 (Jaccard 0.818). Pre-registered rule -> robust (both >= 0.8). Note the NTRK-vs-PDGFRA margin is thin (0.818 vs 0.8).
- Q-C: published adj.P<0.05/0.01/0.001 = 62/48/30 (paper: 41). Counts only.
- Q-D: 101 kinases with >=3 substrate sites (187 with any site); 59 in common with IKAP. Per-kinase Pearson median 0.957, 81.4% >= 0.8. Spearman of NTRK-minus-PDGFRA difference 0.668. Substrate-score limma: 45 kinases adj.P<0.05; IKAP all-10 set on the common kinases has 6; Jaccard 0.133 (IKAP BH computed over 90, then restricted; substrate BH over 101).
- Q-E: LOO Jaccard vs full set: 0.667, 0.833, 0.833, 0.846, 0.636, 0.545, 0.379, 0.917, 0.909, 0.917 (drop PDGFRA1-4, NTRK1-3, CNTRL1-3). >= 0.8 in 6/10 runs -> not stable by the pre-registered criterion. Dropping NTRK3 gives 29 kinases (n_set), i.e. the set is sensitive to which NTRK samples are present.
- Q-F (descriptive, in c5_descriptive_ranks.csv): PRKAA2 ranks 2nd by signed NTRK-PDGFRA effect in all three scores; PRKAA1 and EEF2K rank mid/low by effect (37-47th, 26-47th); EEF2K has only 1 substrate site so no substrate score.
Assumed / not verified: Data1-7 -> PDGFRA1-4, NTRK1-3 mapping (Data8-10 = cortex verified earlier); substrate-score column order mapped to IKAP sample names on that assumption; IKAP not rerun; kinase names matched exactly between IKAP and the KINASE column (no alias mapping); a site counts for every listed kinase (duplicates removed).
Choices within the pre-registration: Q-F ranks given both by signed and absolute effect (rule did not specify); LOO Jaccard defined as 1 if both sets empty (never occurred); NTRK-vs-PDGFRA adj.P is BH of limma's contrast p over the kinases in each analysis.
Fixes to own code: removed an unused alternative-design branch and a placeholder column after Q-A passed exactly; script re-run is idempotent for C5 rows in results/metrics.csv.
Not committed (per instruction). git status --short JUMPn JUMPp: empty.

## 2026-10-18 C5 correction and verdict
- Correction to the C5 RESULTS entry: Q-F ranks by signed NTRK-PDGFRA effect: PRKAA2 2nd (IKAP coded), 2nd (IKAP all-10), 6th (substrate average); PRKAA1 45th, 37th, 34th; EEF2K 47th (coded), 26th (all-10), no substrate score. By F (tumour vs cortex): PRKAA1 9th/3rd/5th, PRKAA2 8th/4th/3rd, EEF2K 20th/17th (coded/all-10; no substrate score).
- Verdict: C5a robust to sample choice; C5b method-dependent and fragile (Q-B passes narrowly at Jaccard 0.818; Q-E fails 6/10; Q-D Jaccard 0.133). The Q-B rule alone gives "robust"; the label also reflects the pre-registered stability and agreement criteria.
- Not verified: Data1-7 mapping; IKAP scale (Fig. 1 axis label assumes log2); the 62-vs-41 gap (needs Supplementary Data 4e).

## 2026-10-18 PRE-REGISTRATION: C1 sample exploration (written BEFORE running c1_exploration.R)
Prior knowledge: earlier sample-key work flagged that PDGFRA_1a/1b (PD1/PD2) might be technical replicates; this was not tested. No C1 result was computed before this entry.
Data: whl_hgg_v5 (whole proteome, 13,419 rows) and pho_hgg_uniq_v5 (phospho, 31,957 rows; 44 contaminant rows with CON_ removed), all values > 0, no missing values. Columns mapped to PD1-4, TRK1-3, Cortex1-3 via the sample key. RNA-seq is not in the repo (C11).
- Q-A Paper setup: PCA on log2 of all features, and Ward clustering of z-scored top-1000-MAD log2 features (Euclidean), recover the three groups (ARI = 1 at k = 3) in both layers. Report PC1/PC2 % variance against the paper (PC1 76% proteome, 66% phospho).
- Q-B Method choice: PCA variants (log2 vs raw; centred vs scaled; all vs top-1000 MAD) with k-means (k = 3, PC1-2, seed 1, nstart 25); hierarchical variants (log2 vs raw; all vs top-1000; Euclidean on z-scored features vs 1 - Pearson; linkage ward.D, ward.D2, average, complete; cut at k = 3). Metric: ARI against the true genotype labels. Rule: robust if ARI = 1 for >= 90% of variants within each analysis type and layer; otherwise method-dependent, listing the failing variants.
- Q-C Genotype separation within tumours only (7 samples, features re-selected within tumours, k = 2): same variants and rule, ARI against PDGFRA vs NTRK.
- Q-D Stability: leave-one-sample-out (10 runs per layer) for (P1) log2, top-1000 MAD re-selected, z-scored Euclidean, ward.D, k = 3, and (P2) PCA log2 all features centred + k-means k = 3. ARI against true labels of the remaining samples. Stable if ARI >= 0.8 in >= 9 of 10 runs.
- Q-E Replicate check: median absolute log2 difference between every sample pair (log2, all features). PD1-PD2 flagged as possible technical replicates if their value is lower than every other within-group pair; then n = 3 independent PDGFRA tumours is flagged for C2/C7.
- Verdict rule: C1 robust if Q-B (all samples and tumours only) is robust and Q-D is stable in both layers; otherwise method-dependent with the failing variants named.

## 2026-10-18 C1 RESULTS (c1_exploration.R; c1_*.csv, results/metrics.csv)
Pre-registered results: Q-A proteome reproduced (PC1/PC2 76.1/12.0; ARI 1.0), phospho PC1/PC2 67.7/13.4 with ARI 0.676 (all samples) and tumour-only PCA ARI 0.432, tumour-only hierarchical ARI 1.0. Q-B: fraction of variants with ARI = 1: proteome 0.34-0.38, phospho 0.22-0.25 (< 0.9) -> method-dependent. Q-D: hierarchical LOO >= 0.8 in 6/10 runs in both layers; PCA+k-means 9/10 proteome, 3/10 phospho. Q-E: PD1-PD2 median abs log2 diff 0.063 (proteome, rank 4/12) and 0.159 (phospho, rank 5/12); lowest pairs are cortex pairs; not flagged.
Exploratory (not pre-registered): successful variants are essentially log2 + top-1000 MAD features; raw intensities and z-scoring all features fail.
Observation for part B: the repeated ARI values 0.676 and 0.432 are the values expected when one NTRK sample is grouped with PDGFRA. Which sample is not established by this run; see the C1 part B pre-registration for provenance.
Implication to check: NTRK-vs-PDGFRA comparisons (C2, C5, C7) may rest on two typical NTRK samples.

## 2026-10-18 PRE-REGISTRATION: C1 part B (written BEFORE running c1b_trk3.R)
Prior knowledge and its source: (1) the repeated ARI values in the C1 run.log (0.6763, 0.4324, 0.643) are the values obtained when one NTRK sample is clustered with the PDGFRA group (worked example reproducible in analysis/C1_exploration/improved/ari_signature_check.R: 0.6763 all samples, 0.4324 tumours only; a PDGFRA sample moving to NTRK would give 0.6591 instead); (2) C5 leave-one-out: dropping NTRK3 raised the NTRK-vs-PDGFRA kinase set from 11 to 29 (c5_loo_ntrk_vs_pdgfra.csv), dropping NTRK1 or NTRK2 gave 7 and 6; (3) the mCherry row in HGG_whl_input.txt is about half as high in TRK3 (143k) as in TRK1/TRK2 (297k, 307k); (4) the paper's Fig. 2a and 2c show one NTRK sample apart. TRK3 as the specific sample was also suggested by an exploratory re-implementation of the PCA and clustering outside the repository (Python, by the assistant); that check is not reproduced in the repo and is not part of the evidence. Because the hypothesis came from the same ten samples, part B is a reproduction of an exploratory observation, not an independent test. All questions below are exploratory; the criteria were fixed before the repo run.
- Q-F1 (which sample): among failing log2 variants (all samples, k = 3), report the misgrouped sample(s). TRK3 is "the" misgrouped sample if it is the only one in >= 80% of failing log2 variants.
- Q-F2 (tumour vs cortex): ARI for k = 2 (Tumour vs Cortex) across all variants; fraction with ARI = 1 per layer and transform.
- Q-F3 (TRK3 profile): per-sample transgene/marker rows (GFP, mCherry, PDGFRA_D842V, TPM3-NTRK1) from HGG_whl_input, and TRK3's position along the PDGFRA-to-TRK1/2 axis in tumour-only PCA (0 = PDGFRA mean, 1 = TRK1/2 mean).
- Q-F4 (specificity): fraction of log2 variants (all samples, k = 3) with ARI = 1 after dropping each single sample. TRK3 is "atypical" if dropping it raises the fraction by >= 0.25 over full data in both layers and gives the highest fraction of all ten drops.
- Q-F5: paper-setup linkage (log2, top-1000, Euclidean z): ARI for ward.D, ward.D2, average, complete, with and without TRK3.

## Figure fixes for the final pass (cosmetic; numbers unaffected)
- C1 part B c1b_tumour_positions: overlapping sample labels (stagger label heights).
- C3 group_profiles: legend sits over some points in the lower-left corner (move legend to "bottom").
- C5 fig1_coded_vs_all10_effect: axis label says "log2-scale IKAP"; the repo does not state the IKAP scale (relabel "IKAP activity units").
- C3 pairwise_ari_heatmap: colour scale is clipped at 0.7-1; say so in the caption.

## 2026-10-18 PRE-REGISTRATION: C2 differential analysis (written BEFORE running any C2 script)
Prior knowledge (all reproducible in the repo): C3 estimated the per-sample rescaling of HGG_whl_input relative to hgg_whole_v1 (PD x0.8857, TRK x0.8433, Cortex x1.0); C1 part B showed TRK3 is atypical (exploratory). No C2 result was computed before this entry.
Data: protein layer = hgg_whole_v1_20160420.txt (accession level, channels sig126..sig131 -> PD1-4, TRK1-3, Cortex1-3 per the sample key; co| contaminants removed); phospho layer = pho_hgg_uniq_v5_20160513.txt (CON_ rows removed; row IDs made unique). Paper counts: 4,703 DE proteins, 6,768 DE phosphosites. Benchmark list (proteins only): the 4,703 accessions in HGG_whl_input.txt. The phospho benchmark needs Supplementary Data 2 (deferred).
Paper-style pipeline (control): one-way ANOVA across Cortex/PDGFRA/NTRK, p from 1000 label permutations with a pooled null, BH adjustment, adj.P < 0.05, and |log2 fold change| >= log2(1.5) (proteins) or log2(2) (phospho) in at least one of the three pairwise comparisons.
- Q-A Reproduction (proteins): paper-style pipeline on the unscaled matrix and on the per-sample rescaled matrix (factors estimated from HGG_whl_input as in C3). Report counts and Jaccard, precision and recall against the 4,703 list. Reproduced if Jaccard >= 0.8.
- Q-B Method choice (both layers): parametric ANOVA, permutation ANOVA (paper-style), limma moderated F; FC cutoff none / 1.5 / 2; FDR 0.05 / 0.01. Rule: robust if every alternative at the paper's thresholds has Jaccard >= 0.8 with the paper-style list in both layers; otherwise method-dependent with the offending variants named. Also report the same lists with and without TRK3 (n = 9) for the limma F pipeline.
- Q-C Rescaling (proteins): Jaccard of paper-style lists with vs without the rescaling, and the number of features that cross the FC cutoff only because of it.
- Q-D Genotype level (primary contrast, both layers): limma moderated t, NTRK vs PDGFRA, adj.P < 0.05 and |log2FC| >= the layer's cutoff; counts, direction (share higher in NTRK), and Jaccard of all-10 vs without-TRK3 lists. Robust if Jaccard >= 0.8.
- Q-E Stability: leave-one-sample-out (10 runs) for the limma F pipeline (3-group list) and the NTRK-vs-PDGFRA list in both layers: Jaccard vs full-data list; stable if >= 0.8 in >= 9 of 10 runs.
- Q-F Phospho normalised to protein: phospho log2 minus matched protein log2 (gene-level sum of whl_hgg_v5 intensities; sites with gene token NA or no match excluded); limma F pipeline on both versions. Report Jaccard of DE lists among matched sites, and the counts of sites DE only unnormalised (abundance-driven) and only normalised. The paper states there was no significant difference; agree if Jaccard >= 0.8.
- Verdict rule: C2a (DE lists) robust if Q-B robust and Q-E stable in both layers; C2b (genotype-level lists) robust if Q-D robust and Q-E stable; C2c (normalisation statement) supported if Q-F agrees. Otherwise method-dependent / not supported with reasons.
