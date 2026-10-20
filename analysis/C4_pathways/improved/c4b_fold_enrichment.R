# C4 part B, post hoc (not pre-registered): fold enrichment of the signalling sets among ALL detected genes.
# Run from the repo root: Rscript analysis/C4_pathways/improved/c4b_fold_enrichment.R
options(stringsAsFactors = FALSE)
stopifnot(file.exists("CLAUDE.md"))
OUT <- "analysis/C4_pathways/improved"
d <- read.csv(file.path(OUT, "c4b_qg1_detected_enrichment_tables.csv"))
d$fold <- (d$overlap / d$n_list) / (d$size / d$n_bg)
SW <- c("KEGG_GLIOMA","KEGG_NEUROTROPHIN_SIGNALING_PATHWAY","KEGG_MTOR_SIGNALING_PATHWAY","KEGG_MAPK_SIGNALING_PATHWAY",
        "KEGG_ERBB_SIGNALING_PATHWAY","KEGG_INSULIN_SIGNALING_PATHWAY","KEGG_PATHWAYS_IN_CANCER",
        "KEGG_P53_SIGNALING_PATHWAY","KEGG_VEGF_SIGNALING_PATHWAY")
res <- NULL
for (ref in unique(d$reference)) for (ls in unique(d$list)) {
  x <- d[d$reference == ref & d$list == ls, ]
  res <- rbind(res, data.frame(reference = ref, list = ls,
    median_fold_S_wide = median(x$fold[x$set %in% SW]), median_fold_all_sets = median(x$fold),
    n_S_wide_above_all_sets_median = sum(x$fold[x$set %in% SW] > median(x$fold))))
}
write.csv(res, file.path(OUT, "c4b_fold_enrichment_summary.csv"), row.names = FALSE)
write.csv(d[d$set %in% SW, c("reference","list","set","size","overlap","n_list","n_bg","fold")],
          file.path(OUT, "c4b_fold_enrichment_S_wide.csv"), row.names = FALSE)
cat("Median fold enrichment among ALL detected genes:\n"); print(res, row.names = FALSE, digits = 3)
u <- d[d$reference == "U1" & d$set %in% SW, c("set", "list", "fold")]
w <- reshape(u, idvar = "set", timevar = "list", direction = "wide")
cat("\nS_wide sets with higher fold in detected phospho than in detected proteome (U1):",
    sum(w$fold.detected_phospho > w$fold.detected_proteome), "of", nrow(w), "\n")
m <- rbind(phospho = w$fold.detected_phospho, proteome = w$fold.detected_proteome); colnames(m) <- sub("KEGG_", "", w$set)
fig <- function() {
  par(mar = c(11, 4.5, 3, 1))
  barplot(m, beside = TRUE, las = 2, cex.names = 0.6, col = c("grey30", "grey80"),
          ylab = "fold enrichment vs reference (U1)", main = "Detected genes: enrichment of signalling sets (post hoc)")
  abline(h = 1, lty = 3); legend("topright", legend = rownames(m), fill = c("grey30", "grey80"), bty = "n", cex = 0.8)
}
png(file.path(OUT, "c4b_fig4_fold_enrichment.png"), width = 9, height = 5.5, units = "in", res = 150); fig(); invisible(dev.off())
pdf(file.path(OUT, "c4b_fig4_fold_enrichment.pdf"), width = 9, height = 5.5); fig(); invisible(dev.off())
vn <- paste0("partB_posthoc_fold_", res$reference, "_", res$list)
M <- rbind(data.frame(claim = "C4", variant = vn, metric = "median_fold_S_wide", value = res$median_fold_S_wide, n = 9),
           data.frame(claim = "C4", variant = vn, metric = "median_fold_all_sets", value = res$median_fold_all_sets, n = 186))
mf <- "results/metrics.csv"; old <- read.csv(mf)
old <- old[!(old$claim == "C4" & old$variant %in% M$variant), , drop = FALSE]
write.csv(rbind(old[, names(M)], M), mf, row.names = FALSE)
cat("Outputs written to", OUT, "and results/metrics.csv\n")