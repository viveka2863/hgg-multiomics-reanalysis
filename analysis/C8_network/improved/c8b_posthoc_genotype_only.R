# C8 part B, post hoc (not pre-registered): permutation null restricted to assignments that keep the true Cortex group.
# Run from the repo root: Rscript analysis/C8_network/improved/c8b_posthoc_genotype_only.R
options(stringsAsFactors = FALSE)
stopifnot(file.exists("CLAUDE.md"))
OUT <- "analysis/C8_network/improved"
d <- read.csv(file.path(OUT, "c8b_h1_permutation_distribution.csv"),
              colClasses = c(PDGFRA = "character", NTRK = "character", Cortex = "character"))
true_cortex <- d$Cortex[d$is_identity][1]
obs  <- d$n_edges[d$is_identity][1]
keep <- d$Cortex == true_cortex
g <- d[keep, ]
res <- data.frame(
  scope = c("all 4,200 assignments", "assignments keeping the true Cortex group", "assignments breaking the Cortex group"),
  n_assignments = c(nrow(d), sum(keep), sum(!keep)),
  n_ge_observed = c(sum(d$n_edges >= obs), sum(g$n_edges >= obs), sum(d$n_edges[!keep] >= obs)))
res$share <- res$n_ge_observed / res$n_assignments
write.csv(res, file.path(OUT, "c8b_posthoc_genotype_only_null.csv"), row.names = FALSE)
write.csv(g[order(-g$n_edges), c("perm", "PDGFRA", "NTRK", "Cortex", "is_identity", "n_edges", "n_paper_nodes")],
          file.path(OUT, "c8b_posthoc_genotype_only_assignments.csv"), row.names = FALSE)
print(res, row.names = FALSE, digits = 3)
cat("\nAccepted edges among the", nrow(g), "assignments keeping the true Cortex group (sorted; observed =", obs, "):\n")
print(sort(g$n_edges, decreasing = TRUE))
M <- data.frame(claim = "C8", variant = "partB_posthoc_genotype_only_null",
                metric = c("n_assignments_ge_observed", "share_ge_observed"),
                value = c(res$n_ge_observed[2], res$share[2]), n = nrow(g))
mf <- "results/metrics.csv"; old <- read.csv(mf)
old <- old[!(old$claim == "C8" & old$variant %in% M$variant), , drop = FALSE]
write.csv(rbind(old[, names(M)], M), mf, row.names = FALSE)
cat("Outputs written to", OUT, "and results/metrics.csv\n")