# C10 Q-A: guilt-by-association on the Pearson network of the 4,703 DE proteins (pre-registered 2026-10-21)
source("analysis/C10_function_prediction/c10_lib.R")
OUTD <- file.path(C10, "original")
U <- load_universe(); N <- nrow(U$X); sets <- get_sets(rownames(U$X))$main
idxs <- lapply(sets, function(s) match(s, rownames(U$X)))
msg <- function(...) cat(sprintf(...), "\n")
msg("labelled input %d; dropped GN NA %d; dropped duplicate genes %d; universe genes %d; genes with a WP label %d",
    U$n_labelled_input, U$n_dropped_GN_NA, U$n_dropped_dup_gene, N, sum(!is.na(U$wp) & grepl("^WP-C", U$wp)))
msg("sets in KEGG_DB %d; sets with >=10 genes in the universe %d", length(read_sets()), length(sets))
stopifnot(identical(colnames(U$X), SAMPLES))
W <- net(U$X, "pearson")
nulls <- make_null_sets(lengths(idxs), N)
t0 <- Sys.time()
res <- run_scorer(function(i) auc_gba(W, i), idxs, nulls)
msg("Q-A time %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
write.csv(res, file.path(OUTD, "c10_qa_per_set.csv"), row.names = FALSE)
nsum <- aggregate(cbind(null_mean, null_sd, null_q95) ~ size, res, mean)
write.csv(nsum, file.path(OUTD, "c10_qa_null_summary.csv"), row.names = FALSE)
qa <- data.frame(universe = N, genes_with_WP_label = sum(grepl("^WP-C", U$wp)), n_sets_ge10 = nrow(res),
                 median_auc = median(res$auc), n_auc_ge_0.7 = sum(res$auc >= 0.7), n_sig_BH05 = sum(res$BH < 0.05),
                 informative_median_ge_0.65 = median(res$auc) >= 0.65)
write.csv(qa, file.path(OUTD, "c10_qa_summary.csv"), row.names = FALSE)
print(qa)
put_metrics(rbind(mrow("QA_pearson", "universe_genes", N, N), mrow("QA_pearson", "n_sets_ge10", nrow(res), N),
  mrow("QA_pearson", "genes_with_WP_label", qa$genes_with_WP_label, N),
  mrow("QA_pearson", "median_set_auc", qa$median_auc, nrow(res)), mrow("QA_pearson", "n_sets_auc_ge_0.7", qa$n_auc_ge_0.7, nrow(res)),
  mrow("QA_pearson", "n_sets_sig_BH05", qa$n_sig_BH05, nrow(res))), "QA_pearson")
writeLines(capture.output(sessionInfo()), "env/sessionInfo_C10.txt")
