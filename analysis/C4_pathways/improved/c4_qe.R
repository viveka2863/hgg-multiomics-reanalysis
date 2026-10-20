# C4 Q-E: leave-one-sample-out stability of the main-setting Q-B classification. Run from repo root.
set.seed(1)
source("analysis/C4_pathways/improved/c4_lib.R")
ctx <- build_ctx(); sets <- ctx$sets; bg <- ctx$bg
sbP <- sets_in_bg(sets, bg$protein, "table", ctx$ortho); sbH <- sets_in_bg(sets, bg$phospho, "table", ctx$ortho)

one_run <- function(keep, label) {
  g <- GROUP[keep]; de <- list()
  for (ly in c("protein", "phospho")) {
    x <- ctx$L[[ly]][, keep]
    f <- de_flag(limma_F_p(x, g), max_pairwise(x, g), 0.05, CUT[[ly]])
    de[[ly]] <- genes_of(rownames(x)[f], ctx$maps[[ly]])
  }
  tP <- ora(de$protein, bg$protein, sbP); tH <- ora(de$phospho, bg$phospho, sbH)
  cl <- classify_wide(tP, tH)
  list(summ = data.frame(left_out = label, n_samples = length(keep), n_DE_genes_proteome = length(de$protein), n_DE_genes_phospho = length(de$phospho),
                         n_enriched_all_proteome = length(enriched(tP)), n_enriched_all_phospho = length(enriched(tH)), claim_summary(cl)),
       cl = tag(cl, left_out = label))
}
runs <- c(list(one_run(SAMPLES, "none (all 10, limma F)")), lapply(SAMPLES, function(s) one_run(setdiff(SAMPLES, s), s)))
qe <- do.call(rbind, lapply(runs, `[[`, "summ")); cl <- do.call(rbind, lapply(runs, `[[`, "cl"))
qe$is_TRK3_run <- qe$left_out == "TRK3"
write.csv(qe, file.path(OUT, "c4_qe_leave_one_out.csv"), row.names = FALSE)
write.csv(cl, file.path(OUT, "c4_qe_S_wide_classification.csv"), row.names = FALSE)
print(qe)
loo <- qe[qe$n_samples == 9, ]; nh <- sum(loo$claim_holds)
cat(sprintf("\nQ-E: claim holds in %d of 10 LOO runs (stable if >= 9); TRK3-dropped run holds: %s; runs with >=3 S_wide sets enriched: %d\n",
            nh, loo$claim_holds[loo$left_out == "TRK3"], sum(loo$n_enriched_either >= 3)))

save_fig("c4_fig3_leave_one_out", 9, 5, function() {
  par(mar = c(5, 4.5, 3, 1)); v <- ifelse(is.na(loo$phospho_only_share), 0, loo$phospho_only_share)
  b <- barplot(v, names.arg = loo$left_out, ylim = c(0, 1), col = ifelse(loo$claim_holds, "#D55E00", "#BBBBBB"),
               ylab = "phospho-only share of enriched S_wide sets", xlab = "sample left out")
  text(b, v + 0.04, ifelse(is.na(loo$phospho_only_share), "n/a", sprintf("%.2f", loo$phospho_only_share)), cex = 0.8)
  mtext(sprintf("n enriched S_wide: %s", paste(loo$n_enriched_either, collapse = " ")), side = 1, line = 3.8, cex = 0.7)
  abline(h = 0.5, lty = 2, lwd = 1.5)
  title("Leave-one-sample-out: phospho-only share (n/a = no S_wide set enriched)", cex.main = 0.9)
})
m <- rbind(mrow("QE_LOO", "runs_claim_holds", nh, 10),
           mrow("QE_LOO", "runs_ge3_S_wide_enriched", sum(loo$n_enriched_either >= 3), 10),
           mrow("QE_LOO", "max_S_wide_enriched_in_a_run", max(loo$n_enriched_either), 10),
           mrow("QE_all10_limmaF", "S_wide_enriched_either", qe$n_enriched_either[1], 9),
           mrow("QE_all10_limmaF", "claim_holds", as.integer(qe$claim_holds[1]), 9),
           mrow("QE_TRK3_dropped", "claim_holds", as.integer(loo$claim_holds[loo$left_out == "TRK3"]), 9))
put_metrics(m)
writeLines(capture.output(sessionInfo()), "env/sessionInfo_C4.txt")
