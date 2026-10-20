# C10 Q-B..Q-E: baselines, group-structure control, robustness, named sets (pre-registered 2026-10-21)
source("analysis/C10_function_prediction/c10_lib.R")
OUTD <- file.path(C10, "improved")
U <- load_universe(); N <- nrow(U$X); S <- get_sets(rownames(U$X)); sets <- S$main
idxs <- lapply(sets, function(s) match(s, rownames(U$X)))
clus <- as.integer(sub("WP-C", "", U$wp)); stopifnot(!anyNA(clus), all(clus %in% 1:5))
msg <- function(...) cat(sprintf(...), "\n")
W <- net(U$X, "pearson")

# check of the C++ scorer against a naive R implementation (3 sets)
naive <- function(W, idx) {
  nm <- setdiff(seq_len(nrow(W)), idx)
  mean(vapply(idx, function(g) { sj <- rowSums(W[, setdiff(idx, g), drop = FALSE]); sg <- sj[g]; sj <- sj[nm]
    mean(sj < sg) + 0.5 * mean(sj == sg) }, 1))
}
for (i in c(1, 40, 100)) stopifnot(abs(naive(W, idxs[[i]]) - set_auc(function(x) auc_gba(W, x), idxs[[i]])) < 1e-12)
msg("C++ GBA scorer matches the naive R implementation on 3 sets")

qa <- read.csv("analysis/C10_function_prediction/original/c10_qa_per_set.csv", stringsAsFactors = FALSE)
stopifnot(identical(qa$set, names(idxs)))
nulls <- make_null_sets(lengths(idxs), N)
stopifnot(abs(set_auc(function(x) auc_gba(W, x), idxs[[1]]) - qa$auc[1]) < 1e-12)

summ <- function(res) data.frame(median_auc = median(res$auc), n_auc_ge_0.7 = sum(res$auc >= 0.7), n_sig_BH05 = sum(res$BH < 0.05), n_sets = nrow(res))
tabs <- list(); tabs$pearson <- qa[, c("set", "size", "auc", "null_mean", "null_sd", "null_q95", "p", "BH")]

# ---- Q-B baselines ----------------------------------------------------------
deg <- rowSums(W)
tabs$degree <- run_scorer(function(i) auc_deg(deg, i), idxs, nulls)
tabs$cluster_vote <- run_scorer(function(i) auc_clu(clus, i), idxs, nulls)
# ---- Q-C within-group residual network ---------------------------------------
Wg <- net(within_group_X(U$X), "pearson")
tabs$within_group <- run_scorer(function(i) auc_gba(Wg, i), idxs, nulls)
# ---- Q-D robustness -----------------------------------------------------------
Ws <- net(U$X, "spearman")
tabs$spearman <- run_scorer(function(i) auc_gba(Ws, i), idxs, nulls)
W9 <- net(U$X[, setdiff(SAMPLES, "TRK3")], "pearson")
tabs$no_TRK3 <- run_scorer(function(i) auc_gba(W9, i), idxs, nulls)
loo <- do.call(rbind, lapply(SAMPLES, function(s) {
  Wl <- net(U$X[, setdiff(SAMPLES, s)], "pearson")
  a <- vapply(idxs, function(i) set_auc(function(x) auc_gba(Wl, x), i), 1)
  data.frame(left_out = s, median_auc = median(a), n_auc_ge_0.7 = sum(a >= 0.7), n_sets = length(a), ge_0.65 = median(a) >= 0.65)
}))
write.csv(loo, file.path(OUTD, "c10_qd_leave_one_out.csv"), row.names = FALSE)

# ---- outputs ------------------------------------------------------------------
long <- do.call(rbind, lapply(names(tabs), function(n) cbind(network = n, tabs[[n]])))
write.csv(long, file.path(OUTD, "c10_per_set_all_networks.csv"), row.names = FALSE)
nullsum <- aggregate(cbind(null_mean, null_sd, null_q95) ~ network + size, long, mean)
write.csv(nullsum, file.path(OUTD, "c10_null_summary.csv"), row.names = FALSE)
sm <- do.call(rbind, lapply(names(tabs), function(n) cbind(network = n, summ(tabs[[n]]))))
write.csv(sm, file.path(OUTD, "c10_summaries.csv"), row.names = FALSE)

g <- tabs$pearson$auc
gain <- sapply(c("degree", "cluster_vote"), function(b) c(median_of_paired_gain = median(g - tabs[[b]]$auc),
                                                       diff_of_medians = median(g) - median(tabs[[b]]$auc)))
adds_value <- all(gain["median_of_paired_gain", ] >= 0.05)
wg_drop <- median(g) - median(tabs$within_group$auc)
stable <- sum(loo$ge_0.65)
qb <- data.frame(baseline = colnames(gain), t(gain), adds_value_both = adds_value, row.names = NULL)
write.csv(qb, file.path(OUTD, "c10_qb_gain.csv"), row.names = FALSE)
qcd <- data.frame(median_pearson = median(g), median_within_group = median(tabs$within_group$auc), drop = wg_drop,
                  driven_by_group_structure = wg_drop >= 0.10,
                  median_spearman = median(tabs$spearman$auc), median_no_TRK3 = median(tabs$no_TRK3$auc),
                  loo_runs_ge_0.65 = stable, loo_stable = stable >= 9, loo_min = min(loo$median_auc), loo_max = max(loo$median_auc))
write.csv(qcd, file.path(OUTD, "c10_qc_qd_summary.csv"), row.names = FALSE)

# ---- Q-E named sets -------------------------------------------------------------
NAMED <- c("KEGG_CELL_CYCLE", "KEGG_OXIDATIVE_PHOSPHORYLATION", "KEGG_RIBOSOME", "KEGG_SPLICEOSOME",
           "KEGG_CITRATE_CYCLE_TCA_CYCLE", "KEGG_FOCAL_ADHESION", S_WIDE)
ne <- do.call(rbind, lapply(unique(NAMED), function(n) {
  idx <- match(S$all[[n]], rownames(U$X)); k <- length(idx)
  r <- data.frame(set = n, size = k, in_main_ge10 = k >= 10, auc_pearson = NA_real_, auc_within_group = NA_real_,
                  auc_degree = NA_real_, auc_cluster_vote = NA_real_, p_pearson = NA_real_, BH_pearson = NA_real_)
  if (k >= 3) { r$auc_pearson <- set_auc(function(x) auc_gba(W, x), idx); r$auc_within_group <- set_auc(function(x) auc_gba(Wg, x), idx)
    r$auc_degree <- set_auc(function(x) auc_deg(deg, x), idx); r$auc_cluster_vote <- set_auc(function(x) auc_clu(clus, x), idx) }
  if (k >= 10) { j <- match(n, tabs$pearson$set); r$p_pearson <- tabs$pearson$p[j]; r$BH_pearson <- tabs$pearson$BH[j] }
  r }))
write.csv(ne, file.path(OUTD, "c10_qe_named_sets.csv"), row.names = FALSE)

# ---- figures --------------------------------------------------------------------
fig(OUTD, "c10_fig1_auc_vs_baselines", 8, 5, function() {
  o <- order(g); par(mar = c(4.5, 4.5, 1, 1))
  plot(NA, xlim = c(1, length(g)), ylim = range(0.3, 1), xlab = "KEGG sets (sorted by GBA AUC)", ylab = "set AUC")
  points(tabs$degree$auc[o], col = "grey60", pch = 16, cex = .6); points(tabs$cluster_vote$auc[o], col = "#E69F00", pch = 17, cex = .6)
  points(g[o], col = "#0072B2", pch = 16, cex = .8); abline(h = 0.65, lty = 2); abline(h = 0.5, col = "grey85")
  legend("topleft", c("GBA (Pearson)", "degree only", "cluster vote", "0.65"), col = c("#0072B2", "grey60", "#E69F00", "black"),
         pch = c(16, 16, 17, NA), lty = c(NA, NA, NA, 2), bty = "n") })
fig(OUTD, "c10_fig2_gba_vs_within_group", 5, 5, function() {
  par(mar = c(4.5, 4.5, 1, 1)); r <- range(g, tabs$within_group$auc)
  plot(g, tabs$within_group$auc, xlim = r, ylim = r, pch = 16, col = "#0072B2", xlab = "set AUC, Pearson network (10 samples)",
       ylab = "set AUC, within-group residual network"); abline(0, 1, lty = 2) })
fig(OUTD, "c10_fig3_leave_one_out", 6, 4.5, function() {
  par(mar = c(5, 4.5, 1, 1)); b <- barplot(loo$median_auc, names.arg = loo$left_out, las = 2, ylim = c(0, 0.8), ylab = "median set AUC", col = "#0072B2")
  abline(h = 0.65, lty = 2); abline(h = median(g), col = "#D55E00") })

mm <- rbind(
  do.call(rbind, lapply(c("degree", "cluster_vote"), function(b) rbind(
    mrow(paste0("QB_", b), "median_set_auc", median(tabs[[b]]$auc), nrow(qa)),
    mrow(paste0("QB_", b), "median_paired_gain_of_GBA", gain["median_of_paired_gain", b], nrow(qa)),
    mrow(paste0("QB_", b), "n_sets_sig_BH05", sum(tabs[[b]]$BH < 0.05), nrow(qa))))),
  mrow("QC_within_group", "median_set_auc", median(tabs$within_group$auc), nrow(qa)),
  mrow("QC_within_group", "drop_vs_QA", wg_drop, nrow(qa)),
  mrow("QC_within_group", "n_sets_sig_BH05", sum(tabs$within_group$BH < 0.05), nrow(qa)),
  mrow("QD_spearman", "median_set_auc", median(tabs$spearman$auc), nrow(qa)),
  mrow("QD_spearman", "n_sets_sig_BH05", sum(tabs$spearman$BH < 0.05), nrow(qa)),
  mrow("QD_no_TRK3", "median_set_auc", median(tabs$no_TRK3$auc), nrow(qa)),
  mrow("QD_no_TRK3", "n_sets_sig_BH05", sum(tabs$no_TRK3$BH < 0.05), nrow(qa)),
  mrow("QD_leave_one_out", "runs_median_ge_0.65", stable, 10),
  mrow("QD_leave_one_out", "min_median_auc", min(loo$median_auc), 10),
  mrow("QD_leave_one_out", "max_median_auc", max(loo$median_auc), 10))
put_metrics(mm, unique(mm$variant))
writeLines(capture.output(sessionInfo()), "env/sessionInfo_C10.txt")
print(sm); print(qb); print(qcd); print(loo); print(ne)
