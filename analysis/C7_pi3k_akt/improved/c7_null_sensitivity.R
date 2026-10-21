# ======================================================================
# C7 part A: how sensitive is the a(P) bootstrap p-value to the choice of null?
# Run from the repo root:  Rscript analysis/C7_pi3k_akt/improved/c7_null_sensitivity.R
# Uses only repo data: JUMPn/JUMPna/.RData (z = log2FC NTRK vs PDGFRA, 28,527 sites, no site names).
# Pre-registration is in LOG.md.
# ======================================================================
options(stringsAsFactors = FALSE)
set.seed(1)
stopifnot(file.exists("CLAUDE.md"))
OUT <- "analysis/C7_pi3k_akt/improved"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

e <- new.env(); load("JUMPn/JUMPna/.RData", envir = e)   # separate environment: keeps the saved .Random.seed out of our session
z <- as.numeric(e$z)
stopifnot(length(z) == 28527, !anyNA(z))
K <- 22; NEG <- c(4, 8, 21); AOBS <- 1.45; NB <- 10000
cat(sprintf("Pool: %d sites; mean %.3f, sd %.3f, fraction positive %.3f\n", length(z), mean(z), sd(z), mean(z > 0)))
cat(sprintf("Observed a(P) = %.2f -> mean signed log2FC per protein = %.3f (about %.2f-fold)\n",
            AOBS, AOBS / sqrt(K), 2^(AOBS / sqrt(K))))
cat("Authors' reported p = 0.003; saved workspace p =", e$pvalue, "\n\n")

boot <- function(pool) {
  ap <- numeric(NB)
  for (i in seq_len(NB)) {
    a <- pool[sample.int(length(pool), K)]
    a[NEG] <- -a[NEG]
    ap[i] <- sum(a) / sqrt(K)
  }
  ap
}
summarise <- function(ap, tau, n_pool) {
  p <- mean(ap > AOBS)
  data.frame(tau = tau, pool_n = n_pool, null_mean = mean(ap), null_sd = sd(ap), p_value = p,
             mc_se = sqrt(max(p * (1 - p), 1e-12) / NB), crit95 = unname(quantile(ap, 0.95)),
             z_obs = (AOBS - mean(ap)) / sd(ap))
}

taus <- seq(0, 1.5, by = 0.1)
res <- list(); ap_authors <- NULL
for (tau in taus) {
  pool <- z[abs(z) >= tau]
  stopifnot(length(pool) >= K)
  ap <- boot(pool)
  if (tau == 0) ap_authors <- ap
  s <- summarise(ap, tau, length(pool)); res[[length(res) + 1]] <- s
  cat(sprintf("tau %.1f: pool %5d  p = %.4f  null mean %.3f sd %.3f  95%% critical a(P) = %.2f\n",
              tau, length(pool), s$p_value, s$null_mean, s$null_sd, s$crit95))
}
R <- do.call(rbind, res)
write.csv(R, file.path(OUT, "c7_null_sensitivity.csv"), row.names = FALSE)
cross <- R$tau[which(R$p_value >= 0.05)[1]]

save_fig <- function(name, fun, w, h) {
  png(file.path(OUT, paste0(name, ".png")), width = w, height = h, units = "in", res = 150); fun(); invisible(dev.off())
}
save_fig("c7_null_sensitivity", function() {
  par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))
  hist(ap_authors, breaks = 60, main = "Authors' null (all sites)", xlab = "a(P) from random draws",
       col = "grey85", border = "grey60")
  abline(v = AOBS, col = "red", lwd = 2)
  legend("topright", legend = "observed 1.45", col = "red", lwd = 2, bty = "n", cex = 0.8)
  plot(R$tau, R$p_value, type = "b", pch = 16, ylim = c(0, max(0.2, R$p_value)),
       xlab = "null drawn from |log2FC| >= tau", ylab = "bootstrap p-value", main = "p-value vs null selection")
  abline(h = 0.05, lty = 2, col = "red")
  plot(R$tau, R$crit95, type = "b", pch = 16, xlab = "null drawn from |log2FC| >= tau",
       ylab = "a(P) needed for p < 0.05", main = "Critical value")
  abline(h = AOBS, lty = 2, col = "red")
}, w = 12, h = 4.5)

M <- do.call(rbind, lapply(seq_len(nrow(R)), function(i) data.frame(
  claim = "C7", variant = sprintf("null_tau_%.1f", R$tau[i]),
  metric = c("p_value", "null_mean", "null_sd", "crit95", "z_obs"),
  value = c(R$p_value[i], R$null_mean[i], R$null_sd[i], R$crit95[i], R$z_obs[i]), n = R$pool_n[i])))
M <- rbind(M, data.frame(claim = "C7", variant = "observed",
                         metric = c("a_P", "mean_signed_log2FC"), value = c(AOBS, AOBS / sqrt(K)), n = K))
mf <- "results/metrics.csv"
if (file.exists(mf)) {
  old <- read.csv(mf)
  old <- old[!(old$claim == "C7" & old$variant %in% M$variant), , drop = FALSE]
  M <- rbind(old[, names(M)], M)
}
write.csv(M, mf, row.names = FALSE)

cat("\n================ SUMMARY (paste this back) ================\n")
print(R, row.names = FALSE, digits = 4)
cat(sprintf("\nFirst tau where p >= 0.05: %s\n", ifelse(is.na(cross), "none up to 1.5", sprintf("%.1f", cross))))
cat("Outputs written to", OUT, "and results/metrics.csv\n")