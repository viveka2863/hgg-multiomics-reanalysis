# Data versions used in this audit

Policy (LOG.md, group-review round 1, P3): three data versions are never mixed inside one analysis, and any crossing is stated.
This file defines the versions and lists, from reading the scripts (not from running them), which version each claim and variant uses.
"unclear" means the script does not settle it.

## Definitions

| Version | What it is | Files | Scale |
|---|---|---|---|
| **A** | Repository files, as shipped | proteome: `JUMPn/JUMPnk/PhosphoWholeNormalization/hgg_whole_v1_20160420.txt` (equal to `whl_hgg_v5`; via `load_protein()` in `c2_lib.R`, contaminants `co\|` removed, 13,567 proteins); phospho: `JUMPn/JUMPnk/Integration/pho_hgg_uniq_v5_20160513.txt` (via `load_phospho()`, `CON_` rows removed, 31,913 rows, duplicated ids kept) | TMT intensities per channel, unscaled |
| **B** | Published supplement | Supp Data 1 (MOESM4, 13,567 proteins incl. the 4,703 labelled WP-C1..5), Supp Data 2a (MOESM5, 28,527 sites, 6,768 labelled PP-C), Supp Data 3 (MOESM6, RNA, 18 samples), sheets 4e/4h (MOESM7), Source Data 4a (MOESM9); read from `data/external/csv` | Intensities rescaled per group by the authors (tumour columns x about 0.8857 / 0.8433 for protein, about 0.8848 for phospho; Cortex x 1.0); see `c2b_a_sd1_scaling_check.csv`, `c2b_a_sd2a_scaling_check.csv` |
| **C** | Ours | version A log2 intensities with per-sample median centring (`median_centre()` in `analysis/C1_exploration/improved/c1c_qc.R`: each sample's median is set to the grand median, the median of the 10 sample medians) | log2, per-sample shift only |
| A' (not a pre-registered version) | Repository file that carries rescaled values | `JUMPn/JUMPnc/HGG_whl_input.txt` (4,703 proteins; tumour columns equal A x the C2 rescale factors, Cortex x 1, so its values equal the B protein values; `c2b_a0_benchmark_vs_SD1.csv`, `c2_rescale_factors.csv`) and "A x `c2_rescale_factors.csv`" computed in scripts | rescaled, as B |

Two repository sources are not intensity matrices and have no version: the authors' kinase tables (`IKAP_result_v2`, `ALL_LIMMA_ANOVA_HGG_Result`, `Mouse_Kinase_Substrate`) and the normalised phospho matrix `hgg_phosphosite_matrix_normalized_whole_v1` (phospho normalised to protein; whether it was built from rescaled or unscaled intensities is unclear). The RNA layer (Supp Data 3, "R") is a separate layer, not a version of the MS data.

## Claim / variant to version

| Claim | Variant | Version(s) | Evidence (script) |
|---|---|---|---|
| C0 / 00_sample_key | key built from A files, HGG_whl_input, raw site files | A, A' (ids and channels only) | `make_sample_key.R` |
| C1 | all PCA / hclust variants, paper set-up (part A) | A (`whl_hgg_v5`, `pho_hgg_uniq_v5`) | `c1_exploration.R` |
| C1 | part B (TRK3 profile) | A for clustering and drop-one; marker rows `cu\|...` from `HGG_whl_input` (A') | `c1b_trk3.R` |
| C1 | part C (QC, round 1) | A, B and C each reported separately | `c1c_qc.R` |
| C2 | grid of 18 variants per layer, genotype lists, leave-one-out, phospho normalisation (parts A-C of the original claim) | A (unscaled) | `c2_main.R` |
| C2 | "rescaled" arm and benchmark | A' (A x factors) vs `HGG_whl_input` | `c2_qa.R` |
| C2 | part B (supplement round 1) | B (SD1, SD2a) vs A unscaled and A x factors; phospho: A only (factors not applied) | `c2b_paper_supp_compare.R` |
| C2 | part C (which NTRK samples) | B intensities and ratios for the sheet arm; A x protein factors for the repo arm; sheet 4h (B) | `c2c_ntrk_subset.R` |
| C2 | part D (site mapping, round 1) | A phospho (rescaled by the phospho factors) vs B SD2a | `c2d_site_mapping.R` |
| C3 | control and WGCNA grid (merge heights, deep split), log2, hierarchical, k-means | A' (`HGG_whl_input`, 4,703 DE proteins) | `c3_compare.R`, `original/compare.R` |
| C3 | "unscaled" variant | A (`hgg_whole_v1` on the same 4,703 proteins) | `c3_compare.R` |
| C3 | "all features" variant (top 8,000 by MAD) | A x factors estimated in the script (A') | `c3_compare.R` |
| C4 | original and improved enrichment | DE flags from C2 on A (unscaled); WP-C cluster labels from the paper (`HGG_whl_input.Cluster`, A') | `c4_lib.R`, `original/c4_main.R`, `c4_variants.R` |
| C5 | original (authors' limma, coded 3/2/3) and improved (all 10) | repository normalised phospho matrix and IKAP outputs (no version; scaling unclear) | `c5_kinase_activity.R` |
| C5 | part B and post hoc (reconciliation with Fig. 4a) | B (sheet 4e, Source Data 4a) vs the repository tables | `c5b_reconcile_paper_tables.R`, `c5b_posthoc_printed_fig4a.R` |
| C6 | TF activity | A (C2 flags and `load_protein()`) | `c6_lib.R`, `c6_main.R` |
| C7 | part A (null sensitivity) | authors' `JUMPna/.RData` z (28,527 sites; version unclear, no site names; cannot be aligned to SD2a rows) | `c7_null_sensitivity.R` |
| C7 | part B (per-site a(P)) | B (sheet 4h, SD2a ratios) and the repository z (unclear) | `c7b_ap_per_site.R` |
| C8 | network, original and improved | A (pho_hgg_uniq_v5, hgg_whole_v1) plus the authors' kinase tables | `c8_lib.R` |
| C8 | part B (per-sample median-centred sensitivity) | A with per-sample median centring, close to C but implemented separately in the script; also the repository normalised matrix | `c8_lib.R` (`sweep` on column medians), `c8b_partB.R` |
| C9 | benchmark | A, plus tables from C5, C6, C8 | `c9_main.R` |
| C10 | function prediction (original and controls) | A' (A x factors) restricted to the paper's 4,703 labelled proteins | `c10_lib.R` |
| C11 | RNA concordance (original) | RNA B-layer (SD3) vs protein B (SD1) and vs A x factors (A') over all 10 samples; phospho SD2a (B) | `c11_qa_concordance.R`, `c11_lib.R` |
| C11 | structure and TF activity (improved) | RNA (SD3); protein DE flags from C2 on A | `c11_main.R` |

## Places where versions are mixed

1. **C1 part B**: clustering on A, marker rows from `HGG_whl_input` (A').
2. **C2 part B**: B intensities against A unscaled and A x factors (the comparison itself is the purpose).
3. **C2 part C**: B sheet arm and A x protein factors repo arm in one comparison table; the phospho scaling (0.8848) differs from the protein factors (0.8433 for TRK), so the protein factors must not be applied to phospho (checked in `c2d_*`: 0 sites agree within 0.01 with the protein factors).
4. **C3**: control on A' (rescaled input); the "unscaled" variant on A; the "all features" variant on A' computed from A inside the script; labels from the paper.
5. **C4**: DE flags on A (unscaled) crossed with paper cluster labels from A'. The 4,703 benchmark list is B/A' while the flags are A.
6. **C5**: the repository normalised matrix (version unclear) compared with B tables (sheet 4e, Source Data 4a).
7. **C7**: B sheet 4h and SD2a ratios together with the repository z of unclear version.
8. **C8 part B**: A with a per-sample median centring that is similar to, but not the same code as, version C.
9. **C10**: A' intensities with A-derived protein-to-gene map; labels from the paper.
10. **C11**: RNA (SD3) against protein in B and in A'; TF activity uses A-based flags. RNA and MS samples are not assumed to be paired (`sample_key_extended.csv`).
11. **C9** inherits whatever C5, C6, C8 use.

Decision needed from the group (not made here): C3 and C10 start from the rescaled A' values because the authors' pipeline does; per P4 these are reproduction arms, and any improved arm that should use A or C instead is not yet defined.
