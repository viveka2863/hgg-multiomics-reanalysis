# ======================================================================
# C5 kinase activity. Run from repo root: Rscript analysis/C5_kinase_activity/c5_kinase_activity.R
# Design and thresholds fixed in LOG.md "PRE-REGISTRATION: C5 kinase activity".
# ======================================================================
options(stringsAsFactors = FALSE)
suppressMessages(library(limma))
set.seed(1)
stopifnot(file.exists("CLAUDE.md"))
D   <- "JUMPn/JUMPnk"
OA  <- "analysis/C5_kinase_activity/original"
OI  <- "analysis/C5_kinase_activity/improved"
dir.create("env", showWarnings = FALSE)

ik  <- read.delim(file.path(D, "KinaseActivityHeatmap/IKAP_result_v2_20160708.txt"), check.names = FALSE, row.names = 1)
pub <- read.delim(file.path(D, "KinaseActivityHeatmap/ALL_LIMMA_ANOVA_HGG_Result_20160708.txt"), check.names = FALSE)
ph  <- read.delim(file.path(D, "PhosphoWholeNormalization/hgg_phosphosite_matrix_normalized_whole_v1_20160420.txt"), check.names = FALSE)
stopifnot(nrow(ik) == 90, ncol(ik) == 10, nrow(ph) == 984, setequal(rownames(pub), rownames(ik)))

metrics <- list()
addm <- function(variant, metric, value, n) metrics[[length(metrics) + 1]] <<-
  data.frame(claim = "C5", variant = variant, metric = metric, value = value, n = n)
jac <- function(a, b) if (length(union(a, b)) == 0) 1 else length(intersect(a, b)) / length(union(a, b))

# groups: group1 = PDGFRA, group2 = NTRK, group3 = CNTRL
run_limma <- function(M, grp) {
  grp <- factor(grp, levels = c("group1", "group2", "group3"))
  design <- model.matrix(~ 0 + grp); colnames(design) <- levels(grp)
  fit <- lmFit(M, design)
  cm  <- makeContrasts("group2...group1" = group2 - group1, "group3...group2" = group3 - group2,
                       "group3...group1" = group3 - group1, levels = design)
  fit <- eBayes(contrasts.fit(fit, cm))
  tt  <- topTable(fit, number = Inf, sort.by = "none")        # F, P.Value, adj.P.Val, AveExpr, contrasts
  # adj.P of the NTRK-vs-PDGFRA contrast (BH over rows)
  tt$p_NTRK_PDGFRA   <- fit$p.value[, "group2...group1"]
  tt$adjp_NTRK_PDGFRA <- p.adjust(tt$p_NTRK_PDGFRA, "BH")
  tt <- cbind(Kinase = rownames(tt), tt); rownames(tt) <- NULL
  tt
}
GRPCOLS <- function(cols) ifelse(grepl("^PDGFRA", cols), "group1", ifelse(grepl("^NTRK", cols), "group2", "group3"))

# ---------------- Q-A -------------------------------------------------
cols9 <- c("PDGFRA1", "PDGFRA2", "PDGFRA2", "PDGFRA4", "NTRK2", "NTRK3", "CNTRL1", "CNTRL2", "CNTRL3")
M9 <- as.matrix(ik[, cols9]); stopifnot(ncol(M9) == 9)
A  <- run_limma(M9, GRPCOLS(cols9))
write.csv(A, file.path(OA, "c5_original_limma_coded_3_2_3.csv"), row.names = FALSE)
cmp_cols <- c("group2...group1", "group3...group2", "group3...group1", "AveExpr", "F", "P.Value", "adj.P.Val")
cmpA <- function(A, label) {
  pb <- pub[A$Kinase, ]
  do.call(rbind, lapply(cmp_cols, function(cn) {
    d <- abs(A[[cn]] - pb[[cn]])
    data.frame(design = label, column = cn, max_abs_diff = max(d), max_rel_diff = max(d / abs(pb[[cn]])))
  }))
}
QA <- cmpA(A, "coded_3_2_3")
print(QA, digits = 4)
exact <- all(QA$max_abs_diff < 1e-6)
cat("Q-A exact (<1e-6) on all columns:", exact, "\n")
if (!exact) stop("Q-A not exact: report differences and stop (pre-registered: no tuning)")
write.csv(QA, file.path(OA, "c5_original_vs_published.csv"), row.names = FALSE)
for (i in seq_len(nrow(QA))) addm(paste0("coded_3_2_3_vs_published_", QA$column[i]), "max_abs_diff", QA$max_abs_diff[i], 90)
for (i in seq_len(nrow(QA))) addm(paste0("coded_3_2_3_vs_published_", QA$column[i]), "max_rel_diff", QA$max_rel_diff[i], 90)

# ---------------- Q-B -------------------------------------------------
cols10 <- c(paste0("PDGFRA", 1:4), paste0("NTRK", 1:3), paste0("CNTRL", 1:3))
M10 <- as.matrix(ik[, cols10]); g10 <- GRPCOLS(cols10)
B <- run_limma(M10, g10)
write.csv(B, file.path(OI, "c5_improved_limma_all10.csv"), row.names = FALSE)
Fset <- function(T) T$Kinase[T$adj.P.Val < 0.05]
Nset <- function(T) T$Kinase[T$adjp_NTRK_PDGFRA < 0.05]
sets <- list(F_coded = Fset(A), F_all10 = Fset(B), N_coded = Nset(A), N_all10 = Nset(B))
cat(sprintf("Q-B (i) F adj.P<0.05: coded %d, all10 %d; (ii) NTRK-vs-PDGFRA adj.P<0.05: coded %d, all10 %d\n",
            length(sets$F_coded), length(sets$F_all10), length(sets$N_coded), length(sets$N_all10)))
jF <- jac(sets$F_coded, sets$F_all10); jN <- jac(sets$N_coded, sets$N_all10)
cat(sprintf("Jaccard F sets %.3f; NTRK-vs-PDGFRA sets %.3f\n", jF, jN))
verdictB <- if (length(sets$N_all10) == 0) "not supported" else if (jF >= 0.8 && jN >= 0.8) "robust" else "method-dependent"
cat("Q-B pre-registered rule ->", verdictB, "\n")
for (v in c("coded_3_2_3", "all_10")) {
  k <- if (v == "all_10") "all10" else "coded"
  addm(v, "n_F_adjP_lt_0.05", length(sets[[paste0("F_", k)]]), 90)
  addm(v, "n_NTRK_vs_PDGFRA_adjP_lt_0.05", length(sets[[paste0("N_", k)]]), 90)
}
addm("coded_vs_all10", "jaccard_F_set", jF, 90); addm("coded_vs_all10", "jaccard_NTRK_vs_PDGFRA_set", jN, 90)
write.csv(data.frame(set = names(sets), n = lengths(sets), kinases = sapply(sets, paste, collapse = ";")),
          file.path(OI, "c5_improved_sets.csv"), row.names = FALSE)

# ---------------- Q-C -------------------------------------------------
QC <- data.frame(threshold = c(0.05, 0.01, 0.001), n_published = sapply(c(0.05, 0.01, 0.001), function(t) sum(pub$adj.P.Val < t)),
                 paper_reported = c(41, NA, NA))
print(QC); write.csv(QC, file.path(OI, "c5_published_counts.csv"), row.names = FALSE)
for (i in 1:3) addm(paste0("published_adjP_lt_", QC$threshold[i]), "n_kinases", QC$n_published[i], 90)

# ---------------- Q-D -------------------------------------------------
dcols <- paste0("Data", 1:10)
stopifnot(all(dcols %in% colnames(ph)))
sub <- lapply(seq_len(nrow(ph)), function(i) trimws(strsplit(ph$KINASE[i], "[,/]")[[1]]))
pairs <- data.frame(site = rep(ph$Phosphosite, lengths(sub)), kinase = unlist(sub))
pairs <- pairs[pairs$kinase != "" & !is.na(pairs$kinase), ]; pairs <- unique(pairs)
nsite <- table(pairs$kinase); kin3 <- names(nsite)[nsite >= 3]
P <- as.matrix(ph[, dcols]); rownames(P) <- ph$Phosphosite
S <- t(sapply(kin3, function(k) colMeans(P[pairs$site[pairs$kinase == k], , drop = FALSE], na.rm = TRUE)))
colnames(S) <- c(paste0("PDGFRA", 1:4), paste0("NTRK", 1:3), paste0("CNTRL", 1:3))   # Data1-4,5-7,8-10 (same order as cols10)
S <- S[, cols10]
cat(sprintf("Substrate score: %d kinases with >=3 sites (of %d with any site)\n", nrow(S), length(nsite)))
write.csv(cbind(Kinase = rownames(S), n_sites = as.integer(nsite[rownames(S)]), as.data.frame(S)),
          file.path(OI, "c5_substrate_average_score.csv"), row.names = FALSE)
common <- intersect(rownames(S), rownames(ik))
cat(sprintf("Kinases in both IKAP and substrate score: %d; IKAP kinases lacking a substrate score: %d\n",
            length(common), 90 - length(common)))
addm("substrate_avg", "n_kinases_ge3_sites", nrow(S), nrow(S)); addm("substrate_avg", "n_common_with_IKAP", length(common), length(common))
r <- sapply(common, function(k) cor(as.numeric(M10[k, ]), as.numeric(S[k, ]), method = "pearson"))
dif <- function(X) rowMeans(X[, g10 == "group2", drop = FALSE]) - rowMeans(X[, g10 == "group1", drop = FALSE])
d_ik <- dif(M10)[common]; d_su <- dif(S)[common]
rho <- suppressWarnings(cor(d_ik, d_su, method = "spearman"))
Cs <- run_limma(S, g10)
write.csv(Cs, file.path(OI, "c5_substrate_limma_all10.csv"), row.names = FALSE)
Nsub <- Cs$Kinase[Cs$adjp_NTRK_PDGFRA < 0.05]
Nik_common <- intersect(sets$N_all10, common)
jD <- jac(Nsub, Nik_common)
cat(sprintf("Q-D (i) Pearson median %.3f, fraction>=0.8 %.3f (n=%d); (ii) Spearman of NTRK-PDGFRA diff %.3f; (iii) substrate set %d, IKAP-all10 set on common %d, Jaccard %.3f\n",
            median(r), mean(r >= 0.8), length(r), rho, length(Nsub), length(Nik_common), jD))
addm("substrate_avg_vs_IKAP", "pearson_per_kinase_median", median(r), length(r))
addm("substrate_avg_vs_IKAP", "pearson_per_kinase_frac_ge_0.8", mean(r >= 0.8), length(r))
addm("substrate_avg_vs_IKAP", "spearman_NTRK_minus_PDGFRA_diff", rho, length(common))
addm("substrate_avg", "n_NTRK_vs_PDGFRA_adjP_lt_0.05", length(Nsub), nrow(S))
addm("substrate_avg_vs_IKAP_all_10", "jaccard_NTRK_vs_PDGFRA_set_common", jD, length(common))
QDtab <- data.frame(Kinase = common, pearson_across_samples = r, diff_IKAP = d_ik, diff_substrate = d_su,
                    n_sites = as.integer(nsite[common]))
write.csv(QDtab, file.path(OI, "c5_ikap_vs_substrate_per_kinase.csv"), row.names = FALSE)

# ---------------- Q-E -------------------------------------------------
loo <- do.call(rbind, lapply(seq_along(cols10), function(i) {
  Tl <- run_limma(M10[, -i], g10[-i]); s <- Nset(Tl)
  data.frame(Kinase = paste0("drop_", cols10[i]), dropped = cols10[i], n_set = length(s),
             jaccard_vs_full = jac(s, sets$N_all10))
}))
names(loo)[1] <- "run"
write.csv(loo, file.path(OI, "c5_loo_ntrk_vs_pdgfra.csv"), row.names = FALSE)
nstable <- sum(loo$jaccard_vs_full >= 0.8)
cat(sprintf("Q-E LOO Jaccard: min %.3f median %.3f; runs >=0.8: %d/10 -> %s\n", min(loo$jaccard_vs_full),
            median(loo$jaccard_vs_full), nstable, if (nstable >= 9) "stable" else "not stable"))
print(loo, digits = 3)
for (i in seq_len(nrow(loo))) {
  addm(paste0("loo_drop_", loo$dropped[i]), "jaccard_vs_full_NTRK_vs_PDGFRA_set", loo$jaccard_vs_full[i], loo$n_set[i])
}
addm("loo", "runs_jaccard_ge_0.8", nstable, 10)

# ---------------- Q-F -------------------------------------------------
targets <- c("PRKAA1", "PRKAA2", "EEF2K")
rk <- function(T, effcol, label, Fcol = "F") {
  e <- T[[effcol]]; names(e) <- T$Kinase; f <- T[[Fcol]]; names(f) <- T$Kinase
  out <- data.frame(Kinase = targets, method = label, n_kinases = nrow(T))
  out$in_method <- targets %in% T$Kinase
  out$effect_NTRK_minus_PDGFRA <- e[targets]
  out$rank_effect_signed_desc <- rank(-e)[targets]; out$rank_effect_abs_desc <- rank(-abs(e))[targets]
  out$F <- f[targets]; out$rank_F_desc <- rank(-f)[targets]
  out
}
QF <- rbind(rk(A, "group2...group1", "IKAP_coded_3_2_3"), rk(B, "group2...group1", "IKAP_all_10"),
            rk(Cs, "group2...group1", "substrate_avg_all_10"))
print(QF, digits = 3)
write.csv(QF, file.path(OI, "c5_descriptive_ranks.csv"), row.names = FALSE)

# ---------------- Figures --------------------------------------------
fig <- function(name, expr, w = 5.5, h = 5) {
  png(file.path(OI, paste0(name, ".png")), width = w, height = h, units = "in", res = 200); expr(); dev.off()
  pdf(file.path(OI, paste0(name, ".pdf")), width = w, height = h); expr(); dev.off()
}
fig("fig1_coded_vs_all10_effect", function() {
  x <- A$`group2...group1`; y <- B$`group2...group1`
  plot(x, y, pch = 16, col = "grey40", xlab = "NTRK - PDGFRA, coded 3/2/3 (log2-scale IKAP)", ylab = "NTRK - PDGFRA, all 10 samples",
       main = sprintf("Spearman %.3f", cor(x, y, method = "spearman")))
  abline(0, 1, lty = 2); abline(h = 0, v = 0, col = "grey80")
  s <- B$adjp_NTRK_PDGFRA < 0.05; points(x[s], y[s], pch = 16, col = "firebrick")
  legend("topleft", legend = "adj.P<0.05 (all 10)", pch = 16, col = "firebrick", bty = "n")
})
fig("fig2_ikap_vs_substrate", function() {
  plot(d_ik, d_su, pch = 16, col = "grey40", xlab = "IKAP: NTRK - PDGFRA mean difference", ylab = "Substrate average: NTRK - PDGFRA",
       main = sprintf("%d kinases; Spearman %.3f", length(common), rho))
  abline(h = 0, v = 0, col = "grey80")
})
fig("fig3_loo_jaccard", function() {
  bp <- barplot(loo$jaccard_vs_full, names.arg = sub("CNTRL", "C", loo$dropped), ylim = c(0, 1.05), las = 2,
                ylab = "Jaccard vs full-data set", main = "Leave-one-sample-out, NTRK vs PDGFRA")
  abline(h = 0.8, lty = 2, col = "firebrick")
})

# ---------------- metrics + session ----------------------------------
mt <- do.call(rbind, metrics)
old <- read.csv("results/metrics.csv"); old <- old[old$claim != "C5", ]       # idempotent re-run
write.csv(rbind(old, mt), "results/metrics.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "env/sessionInfo_C5.txt")
cat("Done.\n")
