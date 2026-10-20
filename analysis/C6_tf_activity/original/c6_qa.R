# C6 Q-A (original method: ORA of the paper-style protein DE list vs ENCODE TF target sets, own universe). Run from the repo root.
set.seed(1)
source("analysis/C6_tf_activity/improved/c6_lib.R")
OUT <- "analysis/C6_tf_activity/original"
cx <- build_c6()
# mapping coverage of the TF sets and universe sizes
cov <- data.frame(TF = names(cx$sets), n_human = lengths(cx$sets),
                  n_human_with_ortholog = vapply(cx$sets, function(s) sum(s %in% cx$ortho$V1), 1L),
                  n_mouse = lengths(cx$mapped),
                  n_mouse_in_own = lengths(restrict(cx$mapped, cx$u$own)), n_mouse_in_U1 = lengths(restrict(cx$mapped, cx$u$U1)),
                  n_mouse_in_U3 = lengths(restrict(cx$mapped, cx$u$U3)), benchmark = names(cx$sets) %in% BENCH, row.names = NULL)
write.csv(cov, file.path(OUT, "c6_mapping_coverage.csv"), row.names = FALSE)
us <- data.frame(universe = UNIVERSES, n_genes = unname(lengths(cx$u)),
                 list = rep(names(cx$genes), each = 3), n_list_genes = unlist(lapply(cx$genes, function(g) lengths(lapply(cx$u, function(u) intersect(g, u))))))
write.csv(us, file.path(OUT, "c6_universe_sizes.csv"), row.names = FALSE)
cat("TF sets: human symbols with an ortholog:", sum(cov$n_human_with_ortholog), "/", sum(cov$n_human),
    sprintf("(%.1f%%)\n", 100 * sum(cov$n_human_with_ortholog) / sum(cov$n_human)))
print(sapply(cx$u, length)); print(us)
print(summary(cov$n_mouse_in_own)); cat("benchmark with <5 targets in own:", sum(cov$benchmark & cov$n_mouse_in_own < 5), "\n")

tab <- run_ora(cx$genes$anova_perm_FC1.5_FDR0.05, cx$u$own, cx$mapped)
tab <- tag(tab, list = "anova_perm_FC1.5_FDR0.05", universe = "own", method = "ORA")
tab <- tab[order(tab$p), ]
write.csv(tab, file.path(OUT, "c6_qa_enrichment_own_protein.csv"), row.names = FALSE)
s <- summarise_tab(tab); write.csv(s, file.path(OUT, "c6_qa_summary.csv"), row.names = FALSE); print(s)
cat("enriched benchmark TFs:", paste(tab$TF[tab$benchmark & tab$BH < 0.05], collapse = ","), "\n")
