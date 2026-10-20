# C11 Q-B to Q-D (pre-registered in LOG.md): samples behind the SD3 ratios, structure of the RNA data, TF activity from RNA.
# Run from the repo root after original/c11_qa_concordance.R: Rscript analysis/C11_rna_layer/improved/c11_main.R
set.seed(1)
suppressPackageStartupMessages({ library(matrixStats); library(cluster) })
source("analysis/C11_rna_layer/improved/c11_lib.R")
OUT_I <- file.path(C11, "improved")
sd3 <- read_sd3(); F <- sd3$F; R <- sd3$R; meta <- sd3$meta
MET <- list()
add <- function(v, m, x, n) MET[[length(MET) + 1]] <<- m11(v, m, x, n)

# =============================== Q-B: which samples behind the ratios ===============================================
GR <- list(Cortex = 1:4, PDGFRA = 5:11, NTRK = 12:18)
masks <- function(n) lapply(seq_len(2^n - 1), function(m) which(bitwAnd(m, 2^(seq_len(n) - 1)) > 0))   # all non-empty subsets, bit masks
SUB <- lapply(GR, function(ix) masks(length(ix))); stopifnot(lengths(SUB) == c(15, 127, 127))
sub_names <- lapply(names(GR), function(g) vapply(SUB[[g]], function(s) paste(GR[[g]][s], collapse = "+"), "")); names(sub_names) <- names(GR)
full_id <- sapply(SUB, length)                                         # the full set is the last mask (all bits set)
W <- function(g) { n <- length(GR[[g]]); sapply(SUB[[g]], function(s) { w <- numeric(n); w[s] <- 1 / length(s); w }) }
G <- nrow(F)
t0 <- proc.time()
defs <- list(log2_ratio_of_means = function(g) log2(F[, GR[[g]]] %*% W(g)),   # G x nsubsets
             mean_of_log2 = function(g) log2(F[, GR[[g]]]) %*% W(g))
# median |ours - paper| for all (subset of a, subset of b) pairs: rows = subsets of a, cols = subsets of b
med_grid <- function(A, B, rat) { out <- matrix(NA_real_, ncol(A), ncol(B))
  for (i in seq_len(ncol(A))) out[i, ] <- colMedians(abs((A[, i] - rat) - B)); out }
grid <- list()
for (dn in names(defs)) {
  A <- lapply(names(GR), function(g) defs[[dn]](g)); names(A) <- names(GR)
  grid[[dn]] <- list(PC = med_grid(A$PDGFRA, A$Cortex, R[, "PDGFRA_Cortex"]),       # [PDGFRA subset, Cortex subset]
                     NC = med_grid(A$NTRK,   A$Cortex, R[, "NTRK_Cortex"]),         # [NTRK subset, Cortex subset]
                     NP = med_grid(A$NTRK,   A$PDGFRA, R[, "NTRK_PDGFRA"]))         # [NTRK subset, PDGFRA subset]
}
runtime <- (proc.time() - t0)[["elapsed"]]
cat(sprintf("Q-B enumeration (2 definitions; 15 x 127 x 127 subset triples via 3 pairwise grids): %.1f s elapsed\n", runtime))

SEQ <- list(); JOINT <- list(); TOP <- list()
for (dn in names(defs)) {
  gd <- grid[[dn]]
  # step 1: (PDGFRA, Cortex) subsets reproducing PDGFRA/Cortex
  b1 <- which(gd$PC == min(gd$PC), arr.ind = TRUE)[1, ]; ties1 <- sum(gd$PC <= min(gd$PC) + 1e-12)
  p1 <- b1[1]; c1 <- b1[2]
  # step 2: NTRK subset reproducing NTRK/Cortex given the Cortex subset
  n2 <- which.min(gd$NC[, c1]); ties2 <- sum(gd$NC[, c1] <= min(gd$NC[, c1]) + 1e-12)
  # check: NTRK/PDGFRA at (n2, p1)
  allv <- c(PC = gd$PC[full_id["PDGFRA"], full_id["Cortex"]], NC = gd$NC[full_id["NTRK"], full_id["Cortex"]], NP = gd$NP[full_id["NTRK"], full_id["PDGFRA"]])
  bestv <- c(PC = gd$PC[p1, c1], NC = gd$NC[n2, c1], NP = gd$NP[n2, p1])
  full_set <- c(Cortex = c1 == full_id["Cortex"], PDGFRA = p1 == full_id["PDGFRA"], NTRK = n2 == full_id["NTRK"])
  SEQ[[dn]] <- data.frame(definition = dn, ratio_column = c("PDGFRA_Cortex", "NTRK_Cortex", "NTRK_PDGFRA"),
                          role = c("fitted (step 1)", "fitted (step 2)", "check"),
                          cortex_subset = sub_names$Cortex[c1], pdgfra_subset = sub_names$PDGFRA[p1], ntrk_subset = sub_names$NTRK[n2],
                          n_cortex = length(SUB$Cortex[[c1]]), n_pdgfra = length(SUB$PDGFRA[[p1]]), n_ntrk = length(SUB$NTRK[[n2]]),
                          median_abs_diff_best = unname(bestv), median_abs_diff_all_samples = unname(allv),
                          improvement = unname(allv - bestv), n_tied_best_step = c(ties1, ties2, NA),
                          best_equals_full_groups = all(full_set), within_0.005_of_all = unname(allv - bestv) <= 0.005)
  # supplementary joint search (objective = sum of the three medians; NOT part of the pre-registration)
  best_j <- c(Inf, 0, 0, 0)
  for (n in seq_len(127)) { S3 <- gd$PC + matrix(gd$NC[n, ], 127, 15, byrow = TRUE) + matrix(gd$NP[n, ], 127, 15)   # [p, c]
    m <- min(S3); if (m < best_j[1]) { w <- which(S3 == m, arr.ind = TRUE)[1, ]; best_j <- c(m, n, w[1], w[2]) } }
  jn <- best_j[2]; jp <- best_j[3]; jc <- best_j[4]
  JOINT[[dn]] <- data.frame(definition = dn, ratio_column = c("PDGFRA_Cortex", "NTRK_Cortex", "NTRK_PDGFRA"),
                            cortex_subset = sub_names$Cortex[jc], pdgfra_subset = sub_names$PDGFRA[jp], ntrk_subset = sub_names$NTRK[jn],
                            median_abs_diff_joint = c(gd$PC[jp, jc], gd$NC[jn, jc], gd$NP[jn, jp]), median_abs_diff_all_samples = unname(allv),
                            sum_of_three_joint = best_j[1], sum_of_three_all_samples = sum(allv))
  # top 5 per grid for transparency
  for (k in names(gd)) { m <- gd[[k]]; o <- order(m)[1:5]; ij <- arrayInd(o, dim(m))
    rows <- if (k == "PC") list(sub_names$PDGFRA, sub_names$Cortex) else if (k == "NC") list(sub_names$NTRK, sub_names$Cortex) else list(sub_names$NTRK, sub_names$PDGFRA)
    TOP[[paste(dn, k)]] <- data.frame(definition = dn, grid = k, rank = 1:5, subset_row = rows[[1]][ij[, 1]], subset_col = rows[[2]][ij[, 2]], median_abs_diff = m[o]) }
}
SEQ <- do.call(rbind, SEQ); JOINT <- do.call(rbind, JOINT); TOP <- do.call(rbind, TOP); rownames(SEQ) <- rownames(JOINT) <- rownames(TOP) <- NULL
SEQ$runtime_seconds_enumeration <- runtime
write.csv(SEQ, file.path(OUT_I, "c11_qb_best_subsets_sequential.csv"), row.names = FALSE)
write.csv(JOINT, file.path(OUT_I, "c11_qb_best_subsets_joint_supplementary.csv"), row.names = FALSE)
write.csv(TOP, file.path(OUT_I, "c11_qb_top5_per_grid.csv"), row.names = FALSE)
print(SEQ, digits = 4); print(JOINT, digits = 4)
for (i in seq_len(nrow(SEQ))) { v <- paste0("QB_seq_", SEQ$definition[i], "_", SEQ$ratio_column[i])
  add(v, "median_abs_diff_best", SEQ$median_abs_diff_best[i], G); add(v, "median_abs_diff_all_samples", SEQ$median_abs_diff_all_samples[i], G)
  add(v, "within_0.005_of_all_or_full", as.numeric(SEQ$best_equals_full_groups[i] | SEQ$within_0.005_of_all[i]), 1) }
add("QB_runtime", "seconds", runtime, 15 * 127 * 127)

save_fig11(OUT_I, "c11_fig2_best_subset_vs_all", 10, 4.8, function() {
  par(mfrow = c(1, 2), mar = c(7.5, 6.5, 3, 1))
  for (dn in names(defs)) {
    s <- SEQ[SEQ$definition == dn, ]; j <- JOINT[JOINT$definition == dn, ]
    v <- rbind(all = s$median_abs_diff_all_samples, sequential = s$median_abs_diff_best, joint = j$median_abs_diff_joint)
    b <- barplot(v, beside = TRUE, col = c("#999999", "#0072B2", "#E69F00"), names.arg = gsub("_", "/", s$ratio_column), ylab = "", main = gsub("_", " ", dn), cex.main = 0.95, ylim = c(0, max(v) * 1.5), las = 2, cex.names = 0.8)
    mtext("median |ours - paper| log2 ratio", side = 2, line = 5, cex = 0.8)
    legend("topright", c("all 18 samples", "best subset (sequential, pre-registered)", "best subset (joint, supplementary)"), fill = c("#999999", "#0072B2", "#E69F00"), bty = "n", cex = 0.7)
  } })

# =============================== Q-C: structure of the RNA data ======================================================
Lg <- log2(F)
X18 <- top_mad(Lg, 1000); p <- prcomp(t(X18), center = TRUE, scale. = FALSE); ve <- 100 * p$sdev^2 / sum(p$sdev^2)
cat(sprintf("PCA: %d genes; PC1 %.1f%%, PC2 %.1f%%\n", nrow(X18), ve[1], ve[2]))
grp <- setNames(meta$group, meta$sample); bat <- setNames(meta$batch, meta$sample)
tum <- meta$sample[meta$group != "Cortex"]; Xt <- top_mad(Lg[, tum], 1000)
ward3 <- run_ward(X18, 3); ward2 <- run_ward(Xt, 2)
sil <- function(X, lab) { s <- silhouette(as.integer(factor(lab)), zdist(X)); c(overall = mean(s[, 3]), tapply(s[, 3], lab, mean)) }
SIL <- rbind(data.frame(samples = "all 18", labels = "group", level = names(sil(X18, grp[colnames(X18)])), silhouette = sil(X18, grp[colnames(X18)])),
             data.frame(samples = "tumour 14", labels = "group", level = names(sil(Xt, grp[tum])), silhouette = sil(Xt, grp[tum])),
             data.frame(samples = "tumour 14", labels = "batch", level = names(sil(Xt, bat[tum])), silhouette = sil(Xt, bat[tum])))
rownames(SIL) <- NULL
QC <- data.frame(quantity = c("n_genes_top_MAD_all18", "PC1_variance_pct", "PC2_variance_pct", "ARI_ward_k3_vs_group_all18",
                              "ARI_ward_k2_vs_batch_tumour14", "ARI_ward_k2_vs_genotype_tumour14", "batch_driven (ARI batch >= ARI genotype)"),
                 value = c(nrow(X18), ve[1], ve[2], ari(ward3, grp[colnames(X18)]), ari(ward2, bat[tum]), ari(ward2, grp[tum]), NA))
QC$value[7] <- as.numeric(QC$value[5] >= QC$value[6])
write.csv(QC, file.path(OUT_I, "c11_qc_structure.csv"), row.names = FALSE); write.csv(SIL, file.path(OUT_I, "c11_qc_silhouette.csv"), row.names = FALSE)
CL <- rbind(data.frame(analysis = "ward_k3_all18", sample = colnames(X18), cluster = ward3, group = grp[colnames(X18)], batch = bat[colnames(X18)]),
            data.frame(analysis = "ward_k2_tumour14", sample = tum, cluster = ward2, group = grp[tum], batch = bat[tum]))
PC <- data.frame(sample = rownames(p$x), PC1 = p$x[, 1], PC2 = p$x[, 2], group = grp[rownames(p$x)], batch = bat[rownames(p$x)])
write.csv(CL, file.path(OUT_I, "c11_qc_clusters.csv"), row.names = FALSE); write.csv(PC, file.path(OUT_I, "c11_qc_pca_scores.csv"), row.names = FALSE)
print(QC, digits = 4); print(SIL, digits = 3); print(table(CL$analysis, CL$cluster, CL$batch, useNA = "ifany"))
for (i in seq_len(6)) add(paste0("QC_", QC$quantity[i]), "value", QC$value[i], if (i >= 5) 14 else 18)
add("QC_batch_driven", "ARI_batch_ge_ARI_genotype", QC$value[7], 14)
for (i in seq_len(nrow(SIL))) add(paste("QC_silhouette", SIL$samples[i], SIL$labels[i], SIL$level[i]), "silhouette", SIL$silhouette[i], sum(meta$sample %in% (if (SIL$samples[i] == "all 18") meta$sample else tum)))

# =============================== Q-D: TF activity from RNA ===========================================================
cx <- build_c6(); QA6 <- "anova_perm_FC1.5_FDR0.05"
mapped_u <- lapply(cx$mapped, toupper)                                    # case-insensitive matching: upper-cased symbols on both sides
i <- dedup_idx(sd3$gene, rowMeans(F))                                     # one row per gene (highest mean FPKM)
mp <- apply(abs(R[i, ]), 1, max); de_rna <- key(sd3$gene[i])[sd3$p[i] < 0.05 & mp >= log2(1.5)]
univ_rna <- key(sd3$gene[i])[apply(F[i, ], 1, max) > 1]
cat(sprintf("RNA genes (dedup) %d; DE %d; universe (FPKM>1 in >=1 sample) %d; DE in universe %d\n", length(i), length(de_rna), length(univ_rna), length(intersect(de_rna, univ_rna))))
tab_rna <- run_ora(de_rna, univ_rna, mapped_u)
tab_pro <- run_ora(cx$genes[[QA6]], cx$u$own, cx$mapped)
sm_rna <- summarise_tab(tab_rna); sm_pro <- summarise_tab(tab_pro)
ref <- read.csv("analysis/C6_tf_activity/original/c6_qa_summary.csv")
stopifnot(abs(sm_pro$AUC - ref$AUC) < 1e-9, sm_pro$n_enriched == ref$n_enriched)   # protein side reproduces the committed C6 Q-A
enr <- function(t) t$TF[!is.na(t$BH) & t$BH < 0.05]
eR <- enr(tab_rna); eP <- enr(tab_pro); both <- intersect(eR, eP)
cmp <- merge(tab_rna[, c("TF", "benchmark", "size", "overlap", "p", "BH")], tab_pro[, c("TF", "size", "overlap", "p", "BH")], by = "TF", all = TRUE, suffixes = c("_rna", "_protein"))
cmp$benchmark <- cmp$TF %in% BENCH; cmp$enriched_rna <- cmp$TF %in% eR; cmp$enriched_protein <- cmp$TF %in% eP
write.csv(cmp[order(cmp$p_rna), ], file.path(OUT_I, "c11_qd_enrichment_rna_vs_protein.csv"), row.names = FALSE)
QD <- data.frame(layer = c("RNA", "protein (C6 Q-A, own universe)"), n_de_genes = c(length(de_rna), length(cx$genes[[QA6]])),
                 n_universe = c(length(univ_rna), length(cx$u$own)), rbind(sm_rna, sm_pro))
write.csv(QD, file.path(OUT_I, "c11_qd_auc_summary.csv"), row.names = FALSE); print(QD, digits = 3)
OV <- data.frame(n_enriched_rna = length(eR), n_enriched_protein = length(eP), n_both = length(both), jaccard = jaccard(eR, eP),
                 n_bench_in_both = sum(BENCH %in% both), recall_of_23_in_intersection = sum(BENCH %in% both) / 23,
                 n_bench_tested_rna = sm_rna$n_bench_tested, n_bench_tested_protein = sm_pro$n_bench_tested,
                 tfs_both = paste(sort(both), collapse = ";"), bench_in_both = paste(sort(intersect(BENCH, both)), collapse = ";"),
                 bench_enriched_rna = paste(sort(intersect(BENCH, eR)), collapse = ";"), bench_enriched_protein = paste(sort(intersect(BENCH, eP)), collapse = ";"))
write.csv(OV, file.path(OUT_I, "c11_qd_overlap.csv"), row.names = FALSE); print(t(OV))
add("QD_RNA", "AUC", sm_rna$AUC, sm_rna$n_tested); add("QD_RNA", "AUC_size_only", sm_rna$AUC_size_only, sm_rna$n_tested)
add("QD_RNA", "n_enriched", sm_rna$n_enriched, sm_rna$n_tested); add("QD_RNA", "n_bench_enriched", sm_rna$n_bench_enriched, sm_rna$n_bench_tested)
add("QD_RNA", "n_de_genes", length(de_rna), length(univ_rna)); add("QD_protein_C6", "AUC", sm_pro$AUC, sm_pro$n_tested)
add("QD_overlap", "jaccard_enriched_TFs", OV$jaccard, length(union(eR, eP))); add("QD_overlap", "recall_of_23_in_intersection", OV$recall_of_23_in_intersection, 23)
add("QD_overlap", "n_enriched_both", OV$n_both, length(union(eR, eP)))

save_fig11(OUT_I, "c11_fig3_pca_and_tf_auc", 11, 5, function() {
  layout(matrix(1:2, 1), widths = c(2.2, 1)); par(mar = c(4.5, 4.5, 3, 1))
  gc <- c(Cortex = "#009E73", PDGFRA = "#0072B2", NTRK = "#D55E00"); bp <- c(`1` = 16, `2` = 17); bp_c <- 15
  pch <- ifelse(is.na(PC$batch), bp_c, bp[as.character(PC$batch)])
  plot(PC$PC1, PC$PC2, col = gc[PC$group], pch = pch, cex = 1.8, xlab = sprintf("PC1 (%.1f%%)", ve[1]), ylab = sprintf("PC2 (%.1f%%)", ve[2]),
       main = "RNA PCA, log2 FPKM, top 1,000 MAD genes", cex.main = 0.95)
  text(PC$PC1, PC$PC2, sub("S", "", PC$sample), pos = 3, cex = 0.7)
  legend("topleft", c(names(gc), "batch 1", "batch 2", "cortex (no batch)"), col = c(gc, "grey30", "grey30", "grey30"), pch = c(rep(16, 3), 16, 17, 15), bty = "n", cex = 0.8)
  b <- barplot(QD$AUC, names.arg = c("RNA", "protein"), ylim = c(0, 1), col = c("#CC79A7", "#999999"), ylab = "AUC (23 benchmark TFs vs others, ORA p)", main = "TF recovery", cex.main = 0.95)
  abline(h = 0.7, lty = 2, lwd = 1.5); abline(h = 0.5, lty = 3, col = "grey50"); text(b, QD$AUC + 0.04, sprintf("%.2f", QD$AUC), cex = 0.9)
})
put_metrics11(do.call(rbind, MET), c("QB", "QC", "QD"))
writeLines(capture.output(sessionInfo()), "env/sessionInfo_C11.txt")
