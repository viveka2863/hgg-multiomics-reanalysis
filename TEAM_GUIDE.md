# Team guide: HGG multi-omics robustness audit

Read this first. It explains how the repository is organised, what the original paper did, what we did
instead, and where to find every result. Numbers are copied from the claim READMEs, `LOG.md` and
`results/metrics.csv`. If this guide and a README or CSV disagree, the README or CSV wins; tell the group.

Contents: 1 Project in brief | 2 Repository map | 3 How the work is organised | 4 The original paper |
5 Claim by claim (paper vs ours) | 6 Cross-cutting findings | 7 What we did not cover | 8 How to run things |
9 Working together | 10 Open tasks | 11 Glossary

---

## 1. Project in brief

**Paper:** Wang H. et al., *Deep multiomics profiling of brain tumors identifies signaling networks downstream
of cancer driver genes*, Nat Commun 10:3718 (2019), doi:10.1038/s41467-019-11661-4.
Two mouse high-grade glioma (HGG) models, one driven by PDGFRA D842V and one by the TPM3-NTRK1 fusion, plus
normal cortex, profiled by whole proteome, phosphoproteome (TMT10-plex mass spectrometry) and RNA-seq.

**Our project:** not a reproduction for its own sake. We ask: *which conclusions of the paper survive when
the data mining method changes?* Each paper claim is a numbered "claim" (C1 to C11). For each one we
(1) keep the paper's method, credited, (2) compare it with alternatives, (3) extend it with new analyses, and
(4) end with a verdict. Everything is pre-registered in `LOG.md` before it is run.

**What the course asks for** (group project brief): propose a topic, identify datasets, develop an analysis
strategy using the data mining techniques from the lectures, discuss the findings ("does the use of different
methods strongly impact outcome? is using several methods needed?"), and hand in a 15-20 page report plus a
20-30 slide PowerPoint with voice recording. Rubric (each criterion 0-25): Innovation, Understanding and
application, Teamwork (peer-to-peer exchange), Critical thinking (links to lectures and exercises).
Peer review at the end can reduce your grade. Check the deadline on NTU-Learn.

**The one-sentence result:** the paper's headline numbers can be reproduced only after steps it does not
describe (a per-group rescaling and leaving one NTRK sample out), several claims depend on analysis choices
(reference set, null model, thresholds), and the CRISPR validation cannot validate the rankings because it has
no control gene set.

---

## 2. Repository map

```
hgg-multiomics-reanalysis/
  README.md               short front page (points here)
  TEAM_GUIDE.md           this file
  README_original.md      the authors' README, kept for credit
  CLAUDE.md               rules and context used while working with Claude Code (read it; it holds the
                          "known issues" list and the metric defaults)
  LOG.md                  chronological log: pre-registrations, results, verifications, decisions
  JUMPn/  JUMPp/          the authors' code and data, UNMODIFIED (byte-identical to the original repo)
  analysis/               all our work, one folder per claim (below)
  results/metrics.csv     every quantitative result: claim, variant, metric, value, n
  env/                    sessionInfo per analysis (R version and packages)
  data/external/          NOT in git. The paper's supplementary files (you download them, section 8)
  scratch/                NOT in git. Copies of the authors' files used to run their scripts
```

Inside `analysis/`:

```
analysis/
  00_sample_key/          sample and TMT channel mapping across the five naming schemes
  C1_exploration/         PCA and clustering of the 10 samples
  C2_differential/        differential expression lists and which samples/scaling the paper used
  C3_coexpression/        WGCNA co-expression clusters
  C4_pathways/            pathway enrichment (ORA, FCS, GSEA) and the "phospho-only" claim
  C5_kinase_activity/     IKAP kinase activities
  C6_tf_activity/         transcription factor activity (protein level)
  C7_pi3k_akt/            PI3K-AKT pathway activity score a(P)
  C8_network/             kinase-TF network
  C9_benchmark/           CRISPR benchmark of every ranking
  C10_function_prediction/  guilt-by-association (our extension)
  C11_rna_layer/          RNA vs protein vs phospho, RNA-based TF activity
```

Every claim folder has a `README.md` with five lines: **Kept / Compared / Added / Result / Verdict**.
Start with that file. Inside a folder:
- `original/`: the paper's method (the "paper-settings arm"). Every comparison has one.
- `improved/`: alternatives, controls and follow-ups (part B, part C).
- CSVs: tidy tables. PNGs: figures (PDFs were removed; figures are PNG only).
- `*_lib.R`: shared helper functions for that claim (other claims re-use some, see section 8).

File name prefixes: `cN_` main analysis, `cNb_` follow-up or supplement round ("part B"), `cNc_` second
follow-up ("part C"). In `results/metrics.csv`, variants starting with `partB_`, `partC_`, `suppB_`, `suppC_`
mark those follow-ups.

---

## 3. How the work is organised

**Pre-registration.** Before running a claim we append a "PRE-REGISTRATION" entry to `LOG.md` with the
questions (Q-A, Q-B, ...), metrics and thresholds. Then the analysis runs, then a "RESULTS" entry and a
"verification and verdict" entry follow. We never change a threshold after seeing results. Analyses we thought
of after seeing data are labelled "exploratory" or "post hoc".

**Evidence rule.** Nothing is reported unless a script in the repo reproduces it. Checks done outside the repo
are only leads. Where this guide mentions such a lead, it says so (section 7).

**Default metrics** (from `CLAUDE.md`): lists use Jaccard, rankings use Spearman, clusterings use the adjusted
Rand index (ARI). "Agree" means at least 0.8. "Stable" means the result holds in at least 9 of 10
leave-one-sample-out runs. "Reproduced" for a clustering means ARI at least 0.95 against the shipped output.

**Verdict labels:**
- *robust*: holds across the alternatives and sample removal;
- *method-dependent*: the answer changes with the method or setting;
- *not supported*: fails under the paper's own or fairer settings;
- *reproduced / not reproduced*: used when we tried to recreate a number from the paper.

**Run-everything rule.** Scripts are the source of truth. No notebooks. Seeds are set; `env/` holds sessionInfo.

---

## 4. The original paper

All page numbers are journal pages (the PDF page of the same number).

### Data and design (p. 2-3, Fig. 1)
- Mouse HGG from p53-null astrocytes expressing PDGFRA D842V or TPM3-NTRK1, implanted in nude mice; normal
  cortex as control.
- TMT10-plex: **3 cortex, 4 PDGFRA, 3 NTRK** samples (Fig. 1c shows cortex 1-3, PDGFRA 4-7, NTRK 8-10).
- Whole proteome: 13,860 proteins identified; **13,567** quantified in every sample (11,718 gene products).
- Phosphoproteome: 30,431 sites identified; **28,527** quantified in every sample.
- RNA-seq: 4 cortex, 7 PDGFRA, 7 NTRK samples (two batches; replicates from a second cohort), FPKM.
- Data: GEO GSE114331, ProteomeXchange PXD005360. Code: github.com/hongwang198745/HGG_Source_Code.

### Methods (p. 12-13), step by step
1. **Differential expression:** one-way ANOVA across cortex, NTRK, PDGFRA; p-values by 1,000 permutations
   (Storey), Benjamini-Hochberg adjustment; initial cutoff 0.05; then fold change at least 1.5 (proteins) or
   2 (phosphosites) in at least one comparison; a final FDR estimated by permutation.
2. **Clustering:** WGCNA on the DE features (Pearson, positive correlations only, hybrid dynamic tree cut,
   minimum merge height 0.15), eigengene per cluster, features assigned by correlation.
3. **Pathways:** ClueGO on each cluster; right-sided hypergeometric test, BH q < 0.05; KEGG, WikiPathways and
   Reactome pooled; kappa score 0.5 to group pathways into modules. The reference gene set is not stated.
4. **Kinase activity:** IKAP on phospho data normalised to whole proteome; kinase-substrate links from
   PhosphoSitePlus; filters: 0.2 SD for convergence, substrate size below 3 removed unless supported by an
   upstream kinase, BH p < 0.05.
5. **TF activity:** DE proteins (per cluster) and DE transcripts overlapped with ENCODE TF targets, Fisher exact
   test, BH 0.05; a TF is accepted only if it changes in both transcriptome and whole proteome; then checked
   against TF protein or phosphorylation.
6. **Kinase-TF network:** kinase-substrate and TF-target links; edges **accepted manually** when kinase
   activity, TF phosphosite and target expression all follow NTRK > PDGFRA > Cortex.
7. **Pathway activity a(P):** `a(P) = sum(Ci x Fi) / sqrt(k)` over PI3K-AKT proteins with annotated functional
   phosphosites; bootstrap of 10,000 draws of 22 Fi values from all quantified phosphorylation events;
   `p = sum(a(P) > 1.45) / 10,000`.
8. **CRISPR:** mini library of 15 genes (9 from the omics analysis, 6 from cross-species comparison), three
   dual-gRNA constructs each, triplicates, 15 days; a gene is a hit if 2 of 3 constructs are significantly
   depleted (Student t test).
9. **Normalisation remark (Methods p. 13):** kinase activity used phospho normalised to whole proteome; for the
   pathway analysis they did not normalise and state there is "no significant difference". The Discussion
   (p. 11) mentions "scale normalization" (Supplementary Fig. 11), which the Methods do not describe.

### Claims and where they are (figure, page)
| # | Paper claim | Where | Our claim |
|---|---|---|---|
| P1 | PCA and clustering separate cortex, PDGFRA, NTRK (proteome PC1 76%, PC2 12%; phospho 66%, 15%; RNA 54%, 18%) | Fig. 2a-d, p. 4 | C1, C11 |
| P2 | RNA and protein abundance correlate moderately (R2 = 0.5) | Fig. 2e, p. 4 | C11 |
| P3 | 4,703 DE proteins and 6,768 DE phosphosites (2,301 phosphoproteins) | p. 3, Fig. 3a | C2 |
| P4 | Five whole-proteome (WP-C) and five phospho (PP-C) clusters; WP-C1, PP-C1, PP-C2 follow Cortex < PDGFRA < NTRK | Fig. 3b,c, p. 5 | C3 |
| P5 | 67 modules; cell cycle in WP-C1, PI3K in PP-C2; **most cancer signalling pathways are altered only in the phosphoproteome** | Fig. 3d-f, p. 5 | C4 |
| P6 | Regulatory proteins (TFs, kinases, epigenetic genes) are globally higher in tumours, more in NTRK | p. 3-4, Supp. Fig. 5 | not tested |
| P7 | 41 of 187 kinases reprogrammed; kinase superfamilies; 34 AKT substrates | Fig. 4, p. 6-7 | C5 (41 only) |
| P8 | 47 TFs from RNA, 38 seen by MS, 23 supported by protein or phospho | Fig. 5a,b, p. 8 | C6, C11 |
| P9 | Core kinase-TF network: 5 kinases, 6 TFs | Fig. 5c, p. 8 | C8 |
| P10 | NTRK shows 1.45-fold greater PI3K-AKT activity, P = 0.003; higher Ki67 and shorter survival | Fig. 6, p. 8-9 | C7 (a(P) only) |
| P11 | NTRK upregulates other RTKs (feedback loop) | Fig. 7, p. 9-10 | not tested |
| P12 | CRISPR: 5 of 9 prioritised regulators are needed (Jun, Myc, Prkaa1, Prkaa2, Eef2k yes; Brca1, Ep300, Insr, Jund no) | Fig. 8, p. 10-11 | C9 |

---

## 5. Claim by claim: paper vs ours

Format: what the paper did, what we did, what differs, and where to find it. Verdicts are the current ones
(after all follow-ups). Paths are relative to `analysis/<folder>/`.

### 00_sample_key (setup)
- **Problem:** five naming schemes. PD1-4 = channels sig126, sig127N, sig127C, sig128N = PDGFRA_1a, 1b, 2, 3;
  TRK1-3 = sig128C, sig129N, sig129C = NTRK_1-3; Cortex1-3 = sig130N, sig130C, sig131 = CNTRL_1-3.
- **Known issues:** IKAP tumour columns are only assumed (correlation 0.2-0.4); PDGFRA_1a and 1b might be one
  tumour; NTRK3 is atypical; Fig. 1c channel numbering disagrees with the TMT layout in the repo data.
- **Where:** `00_sample_key/sample_key.csv` (+ `make_sample_key.R`, README).

### Undocumented preprocessing (found in C2, affects everything)
`JUMPn/JUMPnc/HGG_whl_input.txt` equals the raw `hgg_whole_v1` values multiplied by per-group constants
(PD x0.8857, TRK x0.8433, Cortex x1.0). The paper's own Supplementary Data 1 intensities equal these rescaled
values. The Methods do not describe it. Factors: `C2_differential/c2_rescale_factors.csv`.

### C1 Sample exploration (paper P1; Fig. 2a,b)
- **Paper:** PCA on all quantified log-intensities and Ward clustering on the 1,000 most variable features
  separate the three groups.
- **Ours:** reproduced the paper setup (proteome PC1 76.1%, ARI 1.0), then ran 80 variants (16 PCA plus
  k-means, 64 hierarchical: log2 vs raw, all vs top-1,000 features, centring/scaling, distances, linkages),
  at k = 3 and tumours only at k = 2; leave-one-sample-out; replicate check. Part B asked which sample breaks it.
- **Results:** only 22-38% of variants recover three groups (pre-registered bar 90%). **TRK3 is the only
  misgrouped sample in 100% of failing log2 variants** (6/6 proteome, 11/11 phospho); dropping TRK3 gives ARI
  1.0 in all 20 log2 variants. Tumour vs cortex separates in 90-100% of variants. Phospho PC1 is 67.7% (paper 66%).
  PD1 and PD2 are not the closest replicate pair (rank 4-5 of 12).
- **Verdict:** method-dependent. Tumour vs cortex is robust; the PDGFRA vs NTRK split rests on TRK1 and TRK2.
- **Where:** `C1_exploration/improved/` (`c1_exploration.R`, `c1b_trk3.R`, `ari_signature_check.R`; figures
  `c1_paper_setup`, `c1_variant_recovery`, `c1_sample_distance`, `c1b_drop_one`, `c1b_tumour_positions`).

### C2 Differential expression lists (paper P3)
- **Paper:** ANOVA with permutation p, BH, FC 1.5 / 2 gives 4,703 proteins and 6,768 phosphosites.
- **Ours (main):** `original/c2_qa.R` rebuilds the paper pipeline; `improved/c2_main.R` compares parametric
  ANOVA, permutation ANOVA and limma F, FC cutoffs (none, 1.5, 2), FDR 0.05 / 0.01, with and without the
  rescaling, with and without TRK3, leave-one-out, genotype (NTRK vs PDGFRA) lists, and protein-normalised
  phospho.
- **Results (main):** protein list not reproduced (Jaccard 0.549 raw, 0.738 rescaled; 36 variants, none at
  0.8). The test hardly matters at the paper's thresholds because about 11,000 of 13,567 proteins pass FDR before
  the fold-change filter; **the list size is set by the fold-change cutoff**. Rescaling changes the list
  (Jaccard 0.572). Genotype lists: 527 proteins / 751 phosphosites with all samples, 820 / 1,363 without TRK3;
  95% / 97% of the all-sample hits persist, so the low Jaccard (0.59 / 0.53) is mostly a power effect. Phospho
  normalised to protein vs not: Jaccard 0.279 at feature level.
- **Supplement round 1 (part B, the paper's own files):** repo benchmark = the supplement list (Jaccard 1.000);
  the paper's intensities equal our rescaled values; our repo phospho matrix (v5, May 2016) is **not** the
  paper's final matrix (only 71% of keys match; 6,658 matched keys differ), so phospho comparisons use the
  supplement intensities.
- **Part C (which NTRK samples?):** the supplement's NTRK ratio columns are reproduced exactly only with
  **NTRK = TRK1 + TRK2** (median difference 0.0025 vs 0.054 with all three; phospho 0 vs 0.078). ANOVA p-values
  agree best without TRK3 (Spearman 0.959 protein, 0.983 phospho vs 0.814 and 0.900). DE lists reproduce
  without TRK3 (Jaccard 0.824 protein, 0.938 phospho vs 0.741 and 0.757); repo pipeline on the rescaled matrix
  without TRK3 gives 0.825. About 16% of our protein list is still not in the paper's (unexplained). The paper
  describes three NTRK samples and never states the exclusion.
- **Verdict:** C2a lists reproduced only with the rescaling and TRK3 excluded (neither stated); C2b genotype
  lists method-dependent; C2c normalisation statement not supported at the feature level (pathway level: C4).
- **Where:** `C2_differential/` (`c2_*.csv` and `c2_fig1-3` at the folder root; `improved/c2b_*`, `c2c_*`).
  Key tables: `c2_counts_per_variant.csv`, `c2_qa_summary.csv`, `c2_qb_method_jaccard.csv`,
  `improved/c2c_h1_best_subset.csv`, `c2c_h3_de_lists_vs_labelled.csv`. The DE flag files
  (`c2_protein_DE_flags.csv.gz` etc.) are inputs to C4, C6, C9, C10.

### C3 Co-expression clusters (paper P4; Fig. 3b,c)
- **Paper:** WGCNA, merge height 0.15, five WP clusters (1,611 / 268 / 376 / 104 / 2,344 proteins).
- **Ours:** control arm = the authors' own script run unchanged (`original/run_control.sh`, outputs in
  `original/`); 14 variants (merge heights, log vs raw, unscaled) plus 10 leave-one-sample-out runs; ARI against
  the paper's labels.
- **Results:** the script reproduces its shipped output (ARI 0.9993) but matches the paper's labels only 96.9%
  (ARI 0.951). The WP-C1 cluster (Cortex < PDGFRA < NTRK) appears in 14/14 variants (Jaccard 0.80-0.92) and
  10/10 leave-one-out runs (ARI 0.87-0.98). **The number of clusters is setting-dependent** (merge 0.25 gives 3;
  log2 or unscaled gives 4). The script uses merge 0.2 (paper says 0.15) and raw intensities.
- **Verdict:** C3a robust (the main pattern); C3b method-dependent (five clusters).
- **Where:** `C3_coexpression/` (`improved/c3_compare.R`; figures `pairwise_ari_heatmap`, `group_profiles`,
  `loo_stability`, `silhouette`; `original/labels_with_accession.csv` has paper / shipped / rerun labels per
  protein and feeds C4 and C10).

### C4 Pathways and the "phospho-only" claim (paper P5; Fig. 3d-f)
- **Paper:** ClueGO per cluster; cell cycle in WP-C1; most cancer signalling pathways enriched in phospho
  clusters but not whole-proteome clusters.
- **Ours:** enrichment on the **global** DE lists with the 186 KEGG sets in `JUMPn/JUMPnp/kegg_DB`, using
  ORA (Fisher), FCS (Wilcoxon) and GSEA (fgsea); 9 pre-registered signalling sets; variants for background,
  size matching, normalisation, FDR, symbol mapping; leave-one-out. Part B tests a genome-wide reference set.
- **Results:** themes reproduce (cell cycle in WP-C1: 34/98 genes, BH 4.4e-6; 6/6 WP-C5 metabolism themes).
  The phospho-only claim holds in **0 of 8** settings when the background is the set of measured genes (phospho
  DE genes are 43% of detected genes; signalling sets sit at chance). With a genome-wide reference it appears
  (7 of 9 sets enriched, 5 phospho-only). Post hoc: detected phosphoproteins are specifically rich in signalling
  sets (median fold 1.87 vs 1.04 for all sets; proteome 1.33 vs 1.26). Methods disagree on proteome pathway lists
  (ORA vs FCS Jaccard 0.47; vs GSEA 0.32). Normalised vs un-normalised phospho at pathway level, genome-wide
  reference: Jaccard 0.70, Spearman 0.91.
- **Differences from the paper:** global lists instead of per-cluster; KEGG only (the paper pooled
  WikiPathways and Reactome; FoxO, HIF-1, PI3K cascade sets are missing); PP cluster labels not used; the
  paper's reference set is unknown.
- **Verdict:** C4a themes supported; C4b phospho-only claim reference-dependent (not supported against the
  measured universe); C4c normalisation at pathway level partly supported.
- **Where:** `C4_pathways/` (`original/c4_main.R`; `improved/c4_variants.R`, `c4_qe.R`, `c4b_reference.R`,
  `c4b_loo.R`, `c4b_fold_enrichment.R`; `c4_lib.R` is re-used by C6, C10, C11).

### C5 Kinase activity (paper P7; Fig. 4)
- **Paper:** IKAP; 41 of 187 kinases reprogrammed.
- **Ours:** reproduced the shipped table exactly (max difference 4.7e-15); used all 10 samples; leave-one-out;
  a simple substrate-average score vs IKAP; supplement comparison with sheet 4e and Fig. 4a.
- **Results:** the shipped sample lists are coded as PDGFRA 1,2,2,4 and NTRK 2,3 (not 4 / 3 / 3). Using all 10
  samples: F-test kinase sets 62 to 63 (Jaccard 0.953), NTRK-vs-PDGFRA sets 9 to 11 (0.818). Leave-one-out stable
  in 6/10 runs (dropping NTRK3 gives 29 kinases). Substrate average vs IKAP: Pearson 0.957 per kinase but
  contrast-set Jaccard 0.133 (45 vs 6 significant kinases). Supplement: sheet 4e P-values equal the shipped
  coded-list table (Spearman 0.998); 63 kinases at BH < 0.05; keeping those with at least 3 quantified substrates
  gives **41** (the paper's count); the set differs from Source Data 4a (Jaccard 0.783) but matches the printed
  Fig. 4a better (0.822), because the printed figure and the Source Data differ by one kinase (MAPKAPK2 vs MAP3K8).
- **Verdict:** C5a robust (exact reproduction; F-set stable); C5b method-dependent (contrast set, IKAP vs simple
  score); the 41 is reconciled against the printed figure.
- **Where:** `C5_kinase_activity/` (`c5_kinase_activity.R`, `c5b_reconcile_paper_tables.R` at the folder root;
  `improved/c5b_posthoc_printed_fig4a.R` and `c5b_fig4a_labels_from_pdf.csv`).

### C6 TF activity (paper P8; Fig. 5a,b)
- **Paper:** 47 TFs from RNA targets, 23 supported by MS.
- **Ours:** proteome-level ORA on the 115 ENCODE target sets in `JUMPn/JUMPnt/TF_DB`; benchmark = the 23 TFs
  printed in Fig. 5b; three reference sets; genotype lists; ORA vs FCS; leave-one-out; CRISPR TF ranks.
- **Results:** 3 of 23 enriched (E2F4, REST, STAT2); ranking quality AUC 0.726 vs 0.623 for target-set size alone;
  AUC 0.72-0.75 under all three universes but the number of enriched TFs swings 4 / 26 / 6. NTRK-vs-PDGFRA lists
  give AUC 0.51 and 0.66 (small lists). ORA vs FCS Jaccard 0.375. Leave-one-out 8/10.
- **Differences:** global list, protein only, no cluster-wise requirement of both layers (C11 adds RNA).
- **Verdict:** C6a supported with a thin margin; C6b ranking robust, enriched lists universe-dependent;
  C6c genotype-level signal not supported.
- **Where:** `C6_tf_activity/` (`original/c6_qa.R`, `improved/c6_main.R`, `c6_lib.R`).

### C7 PI3K-AKT pathway activity a(P) (paper P10; Fig. 6)
- **Paper:** NTRK has 1.45-fold greater activity, P = 0.003 (bootstrap).
- **Ours, part A:** reproduced the bootstrap (p = 0.002; paper 0.003); the null draws from **all 28,527 sites**;
  restricting the null to comparably large changes pushes p above 0.05 (p = 0.030 at threshold 0.3, 0.057 at 0.4,
  0.17 at 1.0). a(P) = 1.45 equals a mean signed log2 change of 0.309 per protein (about 1.24-fold).
- **Ours, part B (supplement sheet 4h):** a(P) reproduced (PDGFRA 3.2403, NTRK 4.6906). The Methods' "1.45" is
  the **difference** (1.4503); the text's "1.45-fold" is the **ratio** (1.4476). 17 of 22 sites favour NTRK
  (sign test p 0.017, Wilcoxon p 0.021); mean difference 0.309 log2 (bootstrap 95% CI 0.059 to 0.526).
  Bootstrap null from all sites p = 0.0027; from the 6,768 DE sites only p = 0.092. The NTRK values in sheet 4h
  come from TRK1 + TRK2 (C2 part C).
- **Verdict:** method-dependent. The direction is consistent; the significance depends on the null.
- **Where:** `C7_pi3k_akt/improved/` (`c7_null_sensitivity.R`, `c7b_ap_per_site.R`, tables `c7b_c1-c4_*`).

### C8 Kinase-TF network (paper P9; Fig. 5c)
- **Paper:** manual acceptance of edges with NTRK > PDGFRA > Cortex co-activation; 5 kinases, 6 TFs.
- **Ours:** an explicit rule R0 on the candidate edges we can build from the repo (28 edges, 15 kinases, 14 TFs):
  an edge is accepted if kinase activity (IKAP), TF phosphosite and TF-target abundance are all ordered
  NTRK > PDGFRA > Cortex. Variants (2 of 3 layers, without TRK3), leave-one-out, base rates, degree and
  betweenness. Part B: permutation null over all 4,200 label assignments, global offsets, circularity check.
- **Results:** R0 accepts 15 of 28 edges and recovers all 7 recoverable paper nodes (MYC, JUND, EP300, BRCA1,
  AKT1, MAPK3, CDK5) including AKT1 to BRCA1 S686 (JUN, CEBPB, PRKCA, PRKAA1 have no candidate edges). Relaxing
  to 2 of 3 layers accepts 26; leave-one-out holds in 8/10 (fails when TRK1 or TRK2 is dropped). The TF-target
  layer passes every candidate edge because proteome and phospho data have a **global group offset**
  (Cortex < PDGFRA < NTRK, about 0.2 log2). Permutation p = 0.0445 over all assignments, but 0.057 when only the
  genotype labels are permuted. 9 of the 15 accepted edges are circular (IKAP is inferred from the same sites);
  the 6 MYC edges are not and all survive median-centring. CRISPR hit MYC and the three non-hits are all accepted.
- **Verdict:** C8a recovered; C8b not robust; C8c borderline.
- **Where:** `C8_network/` (`c8_lib.R` at the root; `original/c8_original.R`; `improved/c8_improved.R`,
  `c8b_partB.R`, `c8b_posthoc_genotype_only.R`).

### C9 CRISPR benchmark (paper P12; Fig. 8)
- **Paper:** 5 of 9 regulators needed; says this validates the multi-omics approach.
- **Ours:** AUC and exact permutation p for six gene-level scores (protein and phospho limma F and genotype
  contrasts, expression level, phosphosite count) on the 9 genes (B9) and on 15 genes with the cross-species
  ones (B15); type-specific rankings (IKAP kinases, TF enrichment, network acceptance) descriptive only;
  without TRK3.
- **Results:** no score reaches AUC 0.8 with p below 0.0083; best is 0.69 (p = 0.24); expression level is at
  chance; type-specific benchmarks cannot reach p below 0.1 (3 vs 1 kinases, 2 vs 3 TFs). The B15 extra genes
  were selected for a strong NTRK-over-PDGFRA pattern, which biases B15.
- **Verdict:** no ranking validated. The paper's 56% has no comparison set of non-prioritised genes, and the
  Source Data contains only the gRNA library, not the counts, so hit calls cannot be re-derived.
- **Where:** `C9_benchmark/improved/c9_main.R` and `c9_*.csv`.

### C10 Guilt-by-association function prediction (our extension; Week 7 topic)
- **Paper:** does not claim this; it uses clusters to describe function.
- **Ours:** hide a gene's KEGG membership and predict it from its co-expression neighbours (leave-one-gene-out
  neighbour voting) on the 4,699 DE genes; baselines (degree only, cluster vote); a within-group network that
  removes the group shift; Spearman, without TRK3, leave-one-out.
- **Results:** median set AUC 0.657 (47 sets at 0.7 or more, 76 significant); gain over baselines 0.056 and
  0.057; the within-group network scores 0.677, so it is not just the group shift; stable in 9/10 runs. Strong for
  protein complexes (ribosome 0.945, oxidative phosphorylation 0.926, cell cycle 0.796), weak for signalling
  sets (0.54-0.73) and focal adhesion (0.406).
- **Verdict:** informative (narrow margin), adds value over baselines (narrow), not group-driven, stable.
- **Where:** `C10_function_prediction/` (`original/c10a_gba.R`, `improved/c10b_controls.R`, `c10_lib.R`).

### C11 RNA layer (paper P1, P2, P8)
- **Paper:** RNA PCA (PC1 54%, PC2 18%), RNA-protein R2 = 0.5, RNA targets define TF activity.
- **Ours:** group-level log2 ratios of RNA vs protein vs phospho; which RNA samples the supplement ratios use;
  RNA PCA and clustering by group and batch; ORA of RNA DE genes on the ENCODE sets (as C6).
- **Results:** fold changes correlate (Spearman 0.65-0.70 RNA vs protein, 0.60-0.66 phospho vs protein; sign
  agreement 97-100%). **The RNA ratio columns use all 18 samples**, unlike the protein and phospho ones.
  RNA separates genotype (ARI 1.0) and not batch (-0.08). RNA TF AUC 0.715 (protein 0.726); only E2F4 and REST
  are enriched in both layers (2 of 23).
- **Differences:** our RNA PCA uses the top 1,000 variable genes (70.0% / 20.8%), the paper used all genes
  (54% / 18%); the paper's R2 is an abundance correlation, ours is a fold-change correlation.
- **Verdict:** concordant; RNA uses all samples; RNA TF activity is as good as protein, with little overlap.
- **Where:** `C11_rna_layer/` (`original/c11_qa_concordance.R`, `improved/c11_main.R`, `c11_lib.R`,
  `c11_convert_sd3.R`).

---

## 6. Cross-cutting findings (the story for the report)

1. **TRK3 is the atypical NTRK sample.** It breaks clustering (C1), the supplement's ratios, p-values and DE
   lists exclude it (C2 part C), the kinase limma uses coded lists with two NTRK samples (C5), and leave-one-out
   fails when TRK1 or TRK2 is dropped (C5, C6, C8). The NTRK side of a(P) (C7) rests on two samples.
2. **A per-group rescaling is applied before the paper's DE analysis** and is not described in the Methods.
3. **Enrichment results depend on the reference set** (C4, C6): the ranking is robust, the thresholded lists
   are not.
4. **A global group offset (Cortex < PDGFRA < NTRK) makes "ordered across layers" rules easy to satisfy**
   (C8), and the DE list size is set by the fold-change filter, not the test (C2).
5. **Circularity:** IKAP activities are inferred from the same phosphosites used as independent evidence (C8).
6. **Nulls matter:** the a(P) p-value moves from 0.003 to 0.09 with a selection-matched null (C7); the C8
   permutation p moves from 0.044 to 0.057 when only genotype is permuted.
7. **The CRISPR validation cannot test the prioritisation** (C9).
8. **Our repo's phospho matrix is an older version than the paper's supplement** (C2 part B).

**Improved analysis strategy** (the constructive output): report all-sample and TRK3-excluded arms and state
the exclusion; document the normalisation; use measured-universe backgrounds and rank-based metrics next to
thresholded lists; test co-activation rules against group-restricted permutations and base rates; use
selection-matched nulls for activity scores; include non-prioritised control genes in validation screens; use
the supplement as the canonical data version.

---

## 7. What we did not cover (state these as limitations)

- Regulatory protein families (P6), kinase superfamilies and AKT substrates (Fig. 4c,d), Ki67 and survival
  (Fig. 6c,d), RTK feedback loop (Fig. 7), cross-species integration (Supp. Fig. 9).
- Per-cluster pathway and TF enrichment (the paper's design) and PP cluster-level analysis.
- Re-deriving CRISPR hit calls (the Source Data has no gRNA counts).
- Sample-level RNA-protein comparison (no RNA-to-proteome sample map).
- **Leads checked outside the repo (not in any script, do not cite as results):** six of the seven Fig. 7c RTKs
  stay significant if all three NTRK samples are used (only Pdgfra drops; that row may mix human and mouse
  peptides); the Ki67 p = 0.048 (means 0.32 and 0.22, n = 4 each, SD 0.04 and 0.07) becomes about 0.06 with a
  Welch test.

---

## 8. How to run things

**Environment:** macOS arm64, R 4.6.1, run with `Rscript` from the repo root (paths are relative; each script
starts with `stopifnot(file.exists("CLAUDE.md"))` or similar). Packages used: limma; WGCNA (C3); igraph (C8);
fgsea from Bioconductor (C4); readxl (supplement scripts); cluster and matrixStats (C11). See `env/`.

**Supplementary data (needed by `c2b_*`, `c2c_*`, `c5b_*`, `c7b_*`, `c11_*`; not in git):** download from the
article page (Supplementary information) and put in `data/external/`:
`41467_2019_11661_MOESM3_ESM.pdf` (file descriptions), `MOESM4` (Supp Data 1, proteome), `MOESM5` (Supp Data 2,
phospho; use sheet "Supp Data 2a"), `MOESM6` (Supp Data 3, RNA-seq), `MOESM7` (Supp Data 4, sheets 4a-4i),
`MOESM9` (Source Data). `c2b_convert_xlsx.R` and `c11_convert_sd3.R` convert the sheets to CSV under
`data/external/csv/` (also not in git). The Supplementary Information PDF (`MOESM1`, with Supp. Figs.) is not
downloaded yet.

**Who needs what:**
- C4, C6, C8, C9, C10 run from repo files and earlier outputs only.
- C3's control arm runs the authors' scripts on copies in `scratch/` (`original/run_control.sh`).
- Dependencies between claims: C2 DE flags feed C4, C6, C9; C3 labels feed C4, C10; `c4_lib.R` is used by C6,
  C10 and C11; C5, C6, C8 outputs feed C9; C6 outputs feed C11.

**Sanity checks after any change:** `git status --short JUMPn JUMPp` must print nothing;
`results/metrics.csv` gets rows appended, never edited.

---

## 9. Working together

- Never edit `JUMPn/` or `JUMPp/`. Never commit `data/external/` or `scratch/`.
- Work on a branch per person or task and merge through pull requests; keep commits small with clear messages.
- New decisions, surprises and anything exploratory go into `LOG.md` (append only). Do not rewrite old entries.
- If a result differs from the paper it is a finding; do not tune settings until it matches.
- New figures: PNG only, saved by the script that produced them; new numbers: appended to `results/metrics.csv`
  by script.
- Record who did what in `CONTRIBUTIONS.md` (create it) for the teamwork criterion and the peer review.
- Check the module's policy on AI tools and disclose how AI assistance was used in the report.

---

## 10. Open tasks

1. Verdict table and heatmap (one row per claim: paper wording, our verdict, key number, which methods were
   compared); it should answer the brief's two questions: does the method matter, are several methods needed.
2. Figure fixes (cosmetic): C1 part B label overlap in `c1b_tumour_positions`; C3 `group_profiles` legend;
   C5 `fig1_coded_vs_all10_effect` axis label ("IKAP activity units"); C3 ARI heatmap colour scale note;
   C2 fig 1 y-axis label; C8 fig 1 self-loop and label clipping; C9 fig 2 B15 title; C11 fig 1 scatter size.
3. Replace stale text (top-level README status) and add the per-claim package table to `env/`.
4. Decide on optional extras: per-cluster C4 using Supp Data 4b and the PP labels; Ki67 and survival;
   Supplementary Information PDF.
5. Report (15-20 pages) and slides (20-30, with voice recording); assign sections.
6. Link findings to lecture content and exercises (rubric: critical thinking).

---

## 11. Glossary

- **DE:** differentially expressed. **WP-C1..5 / PP-C1..5:** whole-proteome / phosphoproteome co-expression
  clusters from the paper. **TMT:** tandem mass tag (10 samples in one run).
- **ANOVA / BH / FDR:** analysis of variance; Benjamini-Hochberg false discovery rate adjustment.
- **limma:** linear models for omics data (moderated t and F tests).
- **WGCNA:** weighted gene co-expression network analysis. **ARI:** adjusted Rand index (cluster agreement).
- **ORA / FCS / GSEA:** over-representation analysis; functional class scoring; gene set enrichment analysis.
- **IKAP:** algorithm inferring kinase activity from substrate phosphorylation. **a(P):** the paper's PI3K-AKT
  pathway activity score.
- **Jaccard:** overlap of two lists (shared / union). **AUC:** probability a positive outranks a negative.
- **LOO:** leave-one-sample-out. **UNIVERSE / background:** the set of genes the enrichment test is compared against.
- **Pre-registration:** thresholds and questions written in `LOG.md` before running.
