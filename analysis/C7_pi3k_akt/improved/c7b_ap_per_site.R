# ======================================================================
# C7 part B (supplement round 1, Q-C): per-site analysis of the a(P) score from Supplementary Data sheet 4h, plus bootstrap nulls.
# Design and thresholds: LOG.md "PRE-REGISTRATION: supplement round 1". Run from the repo root AFTER
# analysis/C2_differential/improved/c2b_convert_xlsx.R:   Rscript analysis/C7_pi3k_akt/improved/c7b_ap_per_site.R
# Null (a) pool: repo JUMPna z (28,527 sites, as in c7_null_sensitivity.R); null (b) pool: the 6,768 labelled (DE) sites of Supp Data 2a.
# ======================================================================
options(stringsAsFactors = FALSE)
set.seed(1)
stopifnot(file.exists("CLAUDE.md"))
OUT <- "analysis/C7_pi3k_akt/improved"
rdx <- function(f) read.csv(file.path("data/external/csv", f), header = FALSE, colClasses = "character", na.strings = "")
msg <- function(...) cat(sprintf(...), "\n")
NB <- 10000

# ---- sheet 4h -------------------------------------------------------------------------------------------------------
h <- rdx("t4h.csv"); stopifnot(h[1, 5] == "Ci", h[1, 8] == "k")
h <- h[-1, ]; stopifnot(nrow(h) == 22)
S <- data.frame(site = h[, 1], F_PDGFRA = as.numeric(h[, 3]), F_NTRK = as.numeric(h[, 4]), Ci = as.numeric(h[, 5]),
                FiCi_PDGFRA_sheet = as.numeric(h[, 6]), FiCi_NTRK_sheet = as.numeric(h[, 7]))
K <- nrow(S); stopifnot(all(S$Ci %in% c(-1, 1)), sum(S$Ci == -1) == 3, as.numeric(h[1, 8]) == 22)
sheet_aP <- c(PDGFRA = as.numeric(h[1, 9]), NTRK = as.numeric(h[1, 10]), ratio = as.numeric(h[1, 11]))
msg("sheet 4h: %d sites, negative-Ci sites at rows %s; sheet a(P): PDGFRA %.4f, NTRK %.4f, ratio %.4f", K, paste(which(S$Ci == -1), collapse = ","), sheet_aP[1], sheet_aP[2], sheet_aP[3])
msg("check: sheet's own Fi*Ci columns equal Ci x Fi (max abs diff): PDGFRA %.2e, NTRK %.2e",
    max(abs(S$FiCi_PDGFRA_sheet - S$Ci * S$F_PDGFRA)), max(abs(S$FiCi_NTRK_sheet - S$Ci * S$F_NTRK)))

# ---- C1: reproduce a(P) -----------------------------------------------------------------------------------------------
aP <- c(PDGFRA = sum(S$Ci * S$F_PDGFRA) / sqrt(K), NTRK = sum(S$Ci * S$F_NTRK) / sqrt(K))
OBS <- aP[["NTRK"]] - aP[["PDGFRA"]]
c1 <- data.frame(quantity = c("aP_PDGFRA", "aP_NTRK", "aP_ratio_NTRK_over_PDGFRA", "aP_difference_NTRK_minus_PDGFRA"),
                 recomputed = c(aP[["PDGFRA"]], aP[["NTRK"]], aP[["NTRK"]] / aP[["PDGFRA"]], OBS),
                 sheet = c(sheet_aP[["PDGFRA"]], sheet_aP[["NTRK"]], sheet_aP[["ratio"]], 4.6906 - 3.2403))
c1$abs_diff <- abs(c1$recomputed - c1$sheet); c1$agree_within_0.01 <- c1$abs_diff <= 0.01
c1$note <- c("", "", "", "sheet value = 4.6906 - 3.2403 (pre-registered observed 1.4503)")
write.csv(c1, file.path(OUT, "c7b_c1_aP_reproduction.csv"), row.names = FALSE); print(c1, row.names = FALSE, digits = 5)
C1_ok <- all(c1$agree_within_0.01[1:2])

# ---- C2: per-site paired analysis ----------------------------------------------------------------------------------------------
S$diff_NTRK_minus_PDGFRA <- S$F_NTRK - S$F_PDGFRA            # = log2 NTRK/PDGFRA at the site
S$d <- S$Ci * S$diff_NTRK_minus_PDGFRA                         # contribution to the a(P) difference (times sqrt(22))
npos <- sum(S$d > 0); nz <- sum(S$d != 0)
sign_p2 <- binom.test(npos, nz, 0.5, alternative = "two.sided")$p.value
sign_p1 <- binom.test(npos, nz, 0.5, alternative = "greater")$p.value
wil <- suppressWarnings(wilcox.test(S$d, mu = 0, alternative = "two.sided", exact = length(unique(abs(S$d))) == length(S$d) && all(S$d != 0)))
wil1 <- suppressWarnings(wilcox.test(S$d, mu = 0, alternative = "greater", exact = length(unique(abs(S$d))) == length(S$d) && all(S$d != 0)))
bs <- replicate(NB, mean(S$d[sample.int(K, K, replace = TRUE)]))
ci <- unname(quantile(bs, c(0.025, 0.975)))
c2 <- data.frame(quantity = c("n_sites", "n_positive", "n_nonzero", "mean_d", "median_d", "boot_CI95_low", "boot_CI95_high", "sign_test_p_two_sided", "sign_test_p_one_sided_greater",
                              "wilcoxon_p_two_sided", "wilcoxon_p_one_sided_greater", "wilcoxon_exact", "check_sum_d_over_sqrtK_equals_observed"),
                 value = c(K, npos, nz, mean(S$d), median(S$d), ci[1], ci[2], sign_p2, sign_p1, wil$p.value, wil1$p.value,
                           as.numeric(length(unique(abs(S$d))) == K), sum(S$d) / sqrt(K) - OBS))
write.csv(c2, file.path(OUT, "c7b_c2_per_site_paired.csv"), row.names = FALSE); print(c2, row.names = FALSE, digits = 5)
write.csv(S[, c("site", "F_PDGFRA", "F_NTRK", "Ci", "diff_NTRK_minus_PDGFRA", "d")], file.path(OUT, "c7b_c2_per_site_values.csv"), row.names = FALSE)

# ---- C3: leave-one-site-out ----------------------------------------------------------------------------------------------------
loo <- data.frame(dropped_site = S$site, diff_div_sqrt21 = sapply(seq_len(K), function(i) sum(S$d[-i]) / sqrt(K - 1)),
                  diff_div_sqrt22 = sapply(seq_len(K), function(i) sum(S$d[-i]) / sqrt(K)))
loo$aP_PDGFRA_sqrt21 <- sapply(seq_len(K), function(i) sum(S$Ci[-i] * S$F_PDGFRA[-i]) / sqrt(K - 1))
loo$aP_NTRK_sqrt21 <- sapply(seq_len(K), function(i) sum(S$Ci[-i] * S$F_NTRK[-i]) / sqrt(K - 1))
write.csv(loo, file.path(OUT, "c7b_c3_leave_one_site_out.csv"), row.names = FALSE)
msg("C3: a(P) difference with one site dropped (formula with sqrt(21)): min %.3f (drop %s), max %.3f (drop %s); all > 0: %s; (sqrt(22) form: %.3f-%.3f)",
    min(loo$diff_div_sqrt21), loo$dropped_site[which.min(loo$diff_div_sqrt21)], max(loo$diff_div_sqrt21), loo$dropped_site[which.max(loo$diff_div_sqrt21)],
    all(loo$diff_div_sqrt21 > 0), min(loo$diff_div_sqrt22), max(loo$diff_div_sqrt22))

# ---- C4: bootstrap nulls ---------------------------------------------------------------------------------------------------------
e <- new.env(); load("JUMPn/JUMPna/.RData", envir = e); zA <- as.numeric(e$z); stopifnot(length(zA) == 28527, !anyNA(zA))
y <- rdx("sd2a.csv"); stopifnot(y[2, 17] == "NTRK/PDGFRA"); y <- y[-(1:2), ]; y <- y[!is.na(y[, 3]), ]; stopifnot(nrow(y) == 28527)
rN <- as.numeric(y[, 17]); lab <- grepl("^PP-C", y[, 1]); stopifnot(sum(lab) == 6768)
zB <- rN[lab]
msg("check: repo z vs Supp Data 2a NTRK/PDGFRA ratio by row position (Pearson / Spearman): %.4f / %.4f; sign agreement %.3f",
    cor(zA, rN, use = "complete.obs"), cor(zA, rN, method = "spearman", use = "complete.obs"), mean(sign(zA) == sign(rN), na.rm = TRUE))
boot <- function(pool) { r <- numeric(NB); for (i in seq_len(NB)) r[i] <- sum(S$Ci * pool[sample.int(length(pool), K)]) / sqrt(K); r }
nullA <- boot(zA); nullB <- boot(zB)
nl <- function(lbl, v, pool) data.frame(null = lbl, pool_n = length(pool), pool_mean = mean(pool), pool_sd = sd(pool), share_pool_positive = mean(pool > 0),
  null_mean = mean(v), null_sd = sd(v), observed = OBS, p_share_ge_observed = mean(v >= OBS), p_plus_one = (sum(v >= OBS) + 1) / (NB + 1), mc_se = sqrt(max(mean(v >= OBS) * (1 - mean(v >= OBS)), 1e-12) / NB),
  crit95 = unname(quantile(v, 0.95)), z_obs = (OBS - mean(v)) / sd(v))
c4 <- rbind(nl("(a) all 28,527 quantified sites (repo JUMPna z)", nullA, zA), nl("(b) 6,768 DE sites of Supp Data 2a (paper's NTRK/PDGFRA ratios)", nullB, zB))
write.csv(c4, file.path(OUT, "c7b_c4_bootstrap_nulls.csv"), row.names = FALSE); print(t(c4), quote = FALSE)
write.csv(data.frame(draw = seq_len(NB), a_P_diff_null_a = nullA, a_P_diff_null_b = nullB), gzfile(file.path(OUT, "c7b_c4_null_draws.csv.gz")), row.names = FALSE)

# ---- verdict (pre-registered rule) --------------------------------------------------------------------------------------------------
ps <- c(sign = sign_p2, wilcoxon = wil$p.value, boot_a = c4$p_share_ge_observed[1], boot_b = c4$p_share_ge_observed[2])
verdict <- if (!C1_ok) "not reproduced (C1 failed)" else if (all(ps < 0.05)) "robust" else "method-dependent"
msg("p-values: %s -> verdict by the pre-registered rule: %s (C1 %s)", paste(names(ps), signif(ps, 3), collapse = "; "), verdict, if (C1_ok) "agrees" else "does not agree")
write.csv(data.frame(test = names(ps), p = ps, lt_0.05 = ps < 0.05, C1_agrees = C1_ok, verdict = verdict), file.path(OUT, "c7b_verdict.csv"), row.names = FALSE)

# ---- figure -----------------------------------------------------------------------------------------------------------------------------
save_fig2 <- function(name, w, h, fn) {
  png(file.path(OUT, paste0(name, ".png")), w, h, units = "in", res = 150); fn(); invisible(dev.off())
}
save_fig2("c7b_fig1_per_site_and_nulls", 14, 5.5, function() {
  layout(matrix(1:3, 1), widths = c(1.4, 1, 1)); par(mar = c(8, 4.5, 3, 1))
  o <- order(S$d); cols <- ifelse(S$Ci == -1, "#E15759", "#4E79A7")
  barplot(S$d[o], names.arg = S$site[o], las = 2, cex.names = 0.7, col = cols[o], ylab = "Ci x (log2FC NTRK - log2FC PDGFRA)",
          main = sprintf("Per-site contributions: %d of %d positive", npos, K))
  abline(h = c(mean(S$d), ci), lty = c(1, 2, 2), col = "black")
  legend("topleft", bty = "n", cex = 0.8, fill = c("#4E79A7", "#E15759"), legend = c("Ci = +1", "Ci = -1 (inhibitory)"))
  mtext(sprintf("mean %.3f, bootstrap 95%% CI %.3f to %.3f", mean(S$d), ci[1], ci[2]), side = 3, line = 0, cex = 0.7)
  par(mar = c(4.5, 4.5, 3, 1))
  for (i in 1:2) { v <- list(nullA, nullB)[[i]]
    hist(v, breaks = 60, freq = FALSE, col = "grey85", border = "grey60", xlim = range(c(v, OBS * 1.05)), xlab = "a(P) difference under the null",
         main = sprintf("Null (%s)\np = %.4f", c("a) all 28,527 sites", "b) 6,768 DE sites")[i], c4$p_share_ge_observed[i]), cex.main = 0.95)
    abline(v = OBS, col = "red", lwd = 2); legend("topright", bty = "n", col = "red", lwd = 2, legend = sprintf("observed %.3f", OBS), cex = 0.85) }
})

# ---- metrics ---------------------------------------------------------------------------------------------------------------------------------
M <- NULL
add <- function(variant, metric, value, n) M <<- rbind(M, data.frame(claim = "C7", variant = variant, metric = metric, value = value, n = n))
add("suppB_C1_aP", "aP_PDGFRA_recomputed", aP[["PDGFRA"]], K); add("suppB_C1_aP", "aP_NTRK_recomputed", aP[["NTRK"]], K)
add("suppB_C1_aP", "max_abs_diff_vs_sheet", max(c1$abs_diff[1:2]), K); add("suppB_C1_aP", "observed_difference", OBS, K)
add("suppB_C2_per_site", "n_positive", npos, K); add("suppB_C2_per_site", "sign_test_p_two_sided", sign_p2, K)
add("suppB_C2_per_site", "wilcoxon_p_two_sided", wil$p.value, K); add("suppB_C2_per_site", "mean_d", mean(S$d), K)
add("suppB_C2_per_site", "boot_CI95_low", ci[1], K); add("suppB_C2_per_site", "boot_CI95_high", ci[2], K)
add("suppB_C3_leave_one_site_out", "min_diff", min(loo$diff_div_sqrt21), K - 1); add("suppB_C3_leave_one_site_out", "max_diff", max(loo$diff_div_sqrt21), K - 1)
for (i in 1:2) { v <- c("suppB_C4_null_all_sites", "suppB_C4_null_DE_sites")[i]
  add(v, "p_share_ge_observed", c4$p_share_ge_observed[i], NB); add(v, "null_mean", c4$null_mean[i], c4$pool_n[i]); add(v, "null_sd", c4$null_sd[i], c4$pool_n[i]); add(v, "crit95", c4$crit95[i], c4$pool_n[i]) }
mf <- "results/metrics.csv"; old <- read.csv(mf, stringsAsFactors = FALSE)
old <- old[!(old$claim == "C7" & grepl("^suppB_", old$variant)), ]
write.csv(rbind(old, M[, names(old)]), mf, row.names = FALSE)
cat("done\n")
