# C6 Q-B to Q-F (pre-registered in LOG.md). Run from the repo root after original/c6_qa.R.
set.seed(1)
source("analysis/C6_tf_activity/improved/c6_lib.R")
OUT <- "analysis/C6_tf_activity/improved"
cx <- build_c6(); L <- UNIVERSES; QA <- "anova_perm_FC1.5_FDR0.05"

# ---- ORA for every list x universe (Q-A, Q-B, Q-C) --------------------------------------------------
tabs <- list(); summ <- list()
for (ln in names(cx$genes)) for (un in L) {
  t <- tag(run_ora(cx$genes[[ln]], cx$u[[un]], cx$mapped), list = ln, universe = un, method = "ORA")
  tabs[[paste(ln, un)]] <- t
  summ[[paste(ln, un)]] <- cbind(list = ln, universe = un, method = "ORA", summarise_tab(t))
}
# ---- Q-D: FCS (own universe only; statistic = -log10 limma F adj.P, best feature per gene) -----------
Lp <- log2(load_protein()); st <- gene_stat(Lp, cx$pm)
stopifnot(setequal(names(st), cx$u$own))
# sanity: the all-10 limma F list re-derived here equals the C2 flag
f10 <- de_flag(limma_F_p(Lp, GROUP[colnames(Lp)]), max_pairwise(Lp, GROUP[colnames(Lp)]), 0.05, CUT[["protein"]])
stopifnot(identical(sort(rownames(Lp)[f10]), sort(cx$fl$Feature[cx$fl[["limma_F_FC1.5_FDR0.05"]]])))
tf <- fcs(st, restrict(cx$mapped, cx$u$own)); names(tf)[1] <- "TF"; tf$benchmark <- tf$TF %in% BENCH
tf <- tag(tf, list = "limmaF_stat", universe = "own", method = "FCS")
tabs[["FCS"]] <- tf; summ[["FCS"]] <- cbind(list = "limmaF_stat", universe = "own", method = "FCS", summarise_tab(tf))
# supplementary: ORA on the limma F list (same family as the FCS statistic)
lf <- genes_of(rownames(Lp)[f10], cx$pm)
tl <- tag(run_ora(lf, cx$u$own, cx$mapped), list = "limma_F_FC1.5_FDR0.05", universe = "own", method = "ORA")
tabs[["ORA_limmaF"]] <- tl; summ[["ORA_limmaF"]] <- cbind(list = "limma_F_FC1.5_FDR0.05", universe = "own", method = "ORA", summarise_tab(tl))

alltab <- do.call(rbind, tabs); rownames(alltab) <- NULL
alltab <- alltab[, c("list", "universe", "method", "TF", "benchmark", "size", "overlap", "n_list", "n_bg", "p", "BH")]
write_gz(alltab, file.path(OUT, "c6_enrichment_tables.csv.gz"))
S <- do.call(rbind, summ); rownames(S) <- NULL
write.csv(S, file.path(OUT, "c6_auc_summary.csv"), row.names = FALSE); print(S, digits = 3)
get <- function(ln, un, me = "ORA") alltab[alltab$list == ln & alltab$universe == un & alltab$method == me, ]
enr <- function(t) t$TF[!is.na(t$BH) & t$BH < 0.05]

# ---- Jaccard of enriched TFs ------------------------------------------------------------------------
J <- data.frame(comparison = c("Q-C: DE_all10 vs DE_noTRK3 (ORA, own)", "Q-D: ORA(anova_perm list) vs FCS (own)",
                               "Q-D: ORA(limma_F list) vs FCS (own)"),
                jaccard = c(jaccard(enr(get("DE_all10", "own")), enr(get("DE_noTRK3", "own"))),
                            jaccard(enr(get(QA, "own")), enr(get("limmaF_stat", "own", "FCS"))),
                            jaccard(enr(get("limma_F_FC1.5_FDR0.05", "own")), enr(get("limmaF_stat", "own", "FCS")))),
                n_a = c(length(enr(get("DE_all10", "own"))), length(enr(get(QA, "own"))), length(enr(get("limma_F_FC1.5_FDR0.05", "own")))),
                n_b = c(length(enr(get("DE_noTRK3", "own"))), length(enr(get("limmaF_stat", "own", "FCS"))), length(enr(get("limmaF_stat", "own", "FCS")))))
write.csv(J, file.path(OUT, "c6_jaccard.csv"), row.names = FALSE); print(J)
# Jaccard of enriched TFs across universes (descriptive)
Ju <- do.call(rbind, lapply(names(cx$genes), function(ln) data.frame(list = ln,
        pair = c("own-U1", "own-U3", "U1-U3"),
        jaccard = c(jaccard(enr(get(ln, "own")), enr(get(ln, "U1"))), jaccard(enr(get(ln, "own")), enr(get(ln, "U3"))),
                    jaccard(enr(get(ln, "U1")), enr(get(ln, "U3")))))))
write.csv(Ju, file.path(OUT, "c6_jaccard_across_universes.csv"), row.names = FALSE); print(Ju)

# ---- Q-E: leave-one-sample-out ------------------------------------------------------------------------
loo <- do.call(rbind, lapply(SAMPLES, function(s) {
  keep <- setdiff(SAMPLES, s); g <- GROUP[keep]; x <- Lp[, keep]
  f <- de_flag(limma_F_p(x, g), max_pairwise(x, g), 0.05, CUT[["protein"]])
  ge <- genes_of(rownames(x)[f], cx$pm)
  t <- run_ora(ge, cx$u$own, cx$mapped)
  cbind(left_out = s, n_DE_features = sum(f), n_DE_genes = length(ge), summarise_tab(t))
}))
loo$is_TRK3_run <- loo$left_out == "TRK3"; loo$auc_ge_0.7 <- loo$AUC >= 0.7
write.csv(loo, file.path(OUT, "c6_qe_leave_one_out.csv"), row.names = FALSE); print(loo, digits = 3)

# ---- Q-F: CRISPR ranks (1 = smallest ORA p among tested TFs), main list -------------------------------
cr <- do.call(rbind, lapply(L, function(un) { t <- get(QA, un); t$rank <- rank(t$p, ties.method = "min")
  t <- t[t$TF %in% names(CRISPR), ]
  data.frame(universe = un, TF = t$TF, crispr = CRISPR[t$TF], rank = t$rank, n_tested = nrow(get(QA, un)), p = t$p, BH = t$BH) }))
rownames(cr) <- NULL; write.csv(cr, file.path(OUT, "c6_qf_crispr_ranks.csv"), row.names = FALSE); print(cr)

# ---- figures -----------------------------------------------------------------------------------------
cols <- c(anova_perm_FC1.5_FDR0.05 = "#0072B2", DE_all10 = "#D55E00", DE_noTRK3 = "#009E73")
save_fig("c6_fig1_auc", 8, 5, function() {
  M <- sapply(names(cx$genes), function(ln) sapply(L, function(un) S$AUC[S$list == ln & S$universe == un & S$method == "ORA" & S$list == ln][1]))
  sz <- sapply(L, function(un) S$AUC_size_only[S$list == QA & S$universe == un & S$method == "ORA"])
  par(mar = c(4.5, 4.5, 3, 1)); b <- barplot(t(M), beside = TRUE, ylim = c(0, 1), col = cols, ylab = "AUC (benchmark TFs vs others, ORA p)",
        xlab = "universe", names.arg = L)
  abline(h = 0.7, lty = 2, lwd = 1.5); abline(h = 0.5, lty = 3, col = "grey50")
  mid <- colMeans(b); points(mid, sz, pch = 23, bg = "white", cex = 1.6, lwd = 2)
  legend("topright", c(names(cx$genes), "size-only baseline (main list)", "AUC 0.7"), fill = c(cols, NA, NA), border = NA,
         pch = c(NA, NA, NA, 23, NA), lty = c(NA, NA, NA, NA, 2), bty = "n", cex = 0.75)
  title("Benchmark-TF recovery (AUC) per universe and DE list", cex.main = 0.95)
})
save_fig("c6_fig2_enriched_bench_vs_other", 9, 5, function() {
  d <- S[S$method == "ORA" & S$list %in% names(cx$genes), ]
  v <- rbind(bench = d$n_bench_enriched, other = d$n_enriched - d$n_bench_enriched)
  par(mar = c(7.5, 4.5, 3, 1)); b <- barplot(v, col = c("#CC79A7", "#999999"), ylab = "enriched TFs (BH < 0.05)",
        names.arg = paste(substr(d$list, 1, 11), d$universe, sep = "\n"), las = 2, cex.names = 0.7)
  legend("topleft", c("benchmark (of 23)", "other TFs"), fill = c("#CC79A7", "#999999"), bty = "n")
  title("Enriched TFs: benchmark vs others", cex.main = 0.95)
})
save_fig("c6_fig3_leave_one_out_auc", 8, 5, function() {
  par(mar = c(5, 4.5, 3, 1)); b <- barplot(loo$AUC, names.arg = loo$left_out, ylim = c(0, 1), ylab = "AUC (own universe)",
        xlab = "sample left out", col = ifelse(loo$auc_ge_0.7, "#0072B2", "#BBBBBB"))
  abline(h = 0.7, lty = 2, lwd = 1.5); text(b, loo$AUC + 0.04, sprintf("%.2f", loo$AUC), cex = 0.8)
  title("Leave-one-sample-out: benchmark-TF AUC (Q-A rerun)", cex.main = 0.95)
})

# ---- metrics --------------------------------------------------------------------------------------------
m <- list()
for (i in seq_len(nrow(S))) { v <- paste0("QABCD_", S$method[i], "_", S$list[i], "_", S$universe[i])
  m[[length(m) + 1]] <- m6(v, "AUC", S$AUC[i], S$n_tested[i]); m[[length(m) + 1]] <- m6(v, "AUC_size_only", S$AUC_size_only[i], S$n_tested[i])
  m[[length(m) + 1]] <- m6(v, "n_bench_enriched", S$n_bench_enriched[i], S$n_bench_tested[i])
  m[[length(m) + 1]] <- m6(v, "n_enriched", S$n_enriched[i], S$n_tested[i])
  m[[length(m) + 1]] <- m6(v, "p_overlap_bench_enriched", S$p_overlap_bench_enriched[i], S$n_tested[i]) }
for (i in seq_len(nrow(J))) m[[length(m) + 1]] <- m6(paste0("Jaccard_", i), "jaccard_enriched_TFs", J$jaccard[i], J$n_a[i] + J$n_b[i])
m[[length(m) + 1]] <- m6("QE_LOO", "runs_AUC_ge_0.7", sum(loo$auc_ge_0.7), 10)
m[[length(m) + 1]] <- m6("QE_LOO", "min_AUC", min(loo$AUC), 10)
m[[length(m) + 1]] <- m6("QE_LOO", "AUC_TRK3_dropped", loo$AUC[loo$left_out == "TRK3"], 9)
for (i in seq_len(nrow(cr))) m[[length(m) + 1]] <- m6(paste0("QF_", cr$universe[i], "_", cr$TF[i]), "rank", cr$rank[i], cr$n_tested[i])
c6_put_metrics(do.call(rbind, m))
writeLines(capture.output(sessionInfo()), "env/sessionInfo_C6.txt")
