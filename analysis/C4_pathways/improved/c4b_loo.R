# C4 part B, Q-G4: leave-one-sample-out stability of the V8 (reference U1) classification. Run from the repo root.
# Same procedure as c4_qe.R (limma F lists, C2 thresholds); only the reference differs.
set.seed(1)
source("analysis/C4_pathways/improved/c4_lib.R")
ctx <- build_ctx(); sets <- ctx$sets
U1 <- sort(unique(ctx$ortho$V2)); sbU1 <- sets_in_bg(sets, U1, "table", ctx$ortho)

one_run <- function(keep, label) {
  g <- GROUP[keep]; de <- list()
  for (ly in c("protein", "phospho")) {
    x <- ctx$L[[ly]][, keep]
    f <- de_flag(limma_F_p(x, g), max_pairwise(x, g), 0.05, CUT[[ly]])
    de[[ly]] <- genes_of(rownames(x)[f], ctx$maps[[ly]])
  }
  tP <- ora(de$protein, U1, sbU1); tH <- ora(de$phospho, U1, sbU1)
  cl <- classify_wide(tP, tH)
  list(summ = data.frame(left_out = label, n_samples = length(keep), n_DE_genes_proteome_in_ref = tP$n_list[1], n_DE_genes_phospho_in_ref = tH$n_list[1],
                         n_enriched_all_proteome = length(enriched(tP)), n_enriched_all_phospho = length(enriched(tH)), claim_summary(cl)),
       cl = tag(cl, left_out = label))
}
runs <- c(list(one_run(SAMPLES, "none (all 10, limma F)")), lapply(SAMPLES, function(s) one_run(setdiff(SAMPLES, s), s)))
qg4 <- do.call(rbind, lapply(runs, `[[`, "summ")); cl <- do.call(rbind, lapply(runs, `[[`, "cl"))
write.csv(qg4, file.path(OUT, "c4b_qg4_leave_one_out.csv"), row.names = FALSE)
write.csv(cl, file.path(OUT, "c4b_qg4_S_wide_classification.csv"), row.names = FALSE)
print(qg4[, c("left_out", "n_DE_genes_proteome_in_ref", "n_DE_genes_phospho_in_ref", "n_enriched_either", "n_phospho_only", "n_proteome_only", "n_both", "phospho_only_share", "claim_holds")])
loo <- qg4[qg4$n_samples == 9, ]; nh <- sum(loo$claim_holds)
cat(sprintf("\nQ-G4: claim holds in %d of 10 LOO runs under V8 (stable if >= 9); all-10 limma-F run holds: %s\n", nh, qg4$claim_holds[1]))

save_fig("c4b_fig3_leave_one_out", 9, 5, function() {
  par(mar = c(5, 4.5, 3, 1)); v <- ifelse(is.na(loo$phospho_only_share), 0, loo$phospho_only_share)
  b <- barplot(v, names.arg = loo$left_out, ylim = c(0, 1), col = ifelse(loo$claim_holds, "#D55E00", "#BBBBBB"),
               ylab = "phospho-only share of enriched S_wide sets", xlab = "sample left out")
  text(b, v + 0.04, ifelse(is.na(loo$phospho_only_share), "n/a", sprintf("%.2f", loo$phospho_only_share)), cex = 0.8)
  mtext(sprintf("n enriched S_wide: %s", paste(loo$n_enriched_either, collapse = " ")), side = 1, line = 3.8, cex = 0.7)
  abline(h = 0.5, lty = 2, lwd = 1.5)
  title("V8 (U1) leave-one-sample-out: phospho-only share (n/a = none enriched)", cex.main = 0.9)
})
m <- rbind(mrow("partB_QG4_LOO_V8", "runs_claim_holds", nh, 10),
           mrow("partB_QG4_LOO_V8", "runs_ge3_S_wide_enriched", sum(loo$n_enriched_either >= 3), 10),
           mrow("partB_QG4_LOO_V8", "min_S_wide_enriched_in_a_run", min(loo$n_enriched_either), 10),
           mrow("partB_QG4_all10_limmaF_V8", "claim_holds", as.integer(qg4$claim_holds[1]), 9))
put_metrics(m)
