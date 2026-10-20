# C8 improved arm: Q-B robustness, Q-C informativeness, Q-D network measures, Q-E CRISPR. Run from repo root.
set.seed(1)
source("analysis/C8_network/c8_lib.R")
library(igraph)
out <- "analysis/C8_network/improved/"
D <- load_all(); E <- candidates(D)
lay <- layers(D, S); et <- edge_table(D, E, lay)
N0 <- nodes_of(et, "accept_R0")
mrows <- list(); addm <- function(variant, metric, value, n) mrows[[length(mrows) + 1]] <<- data.frame("C8", variant, metric, value, n)

## ---- Q-B robustness ----
v1 <- nodes_of(et, "accept_V1")
s9 <- setdiff(S, "TRK3")
et9 <- edge_table(D, E, layers(D, s9)); v2 <- nodes_of(et9, "accept_R0")
qb <- data.frame(run = c("V1_2of3", "V2_noTRK3"), n_samples = c(10, 9),
                 n_edges = c(sum(et$accept_V1), sum(et9$accept_R0)), n_nodes = c(length(v1), length(v2)),
                 jaccard_vs_R0 = c(jacc(v1, N0), jacc(v2, N0)))
loo <- do.call(rbind, lapply(S, function(d) {
  e <- edge_table(D, E, layers(D, setdiff(S, d)))
  nd <- nodes_of(e, "accept_R0")
  data.frame(run = paste0("LOO_", d), n_samples = 9, n_edges = sum(e$accept_R0), n_nodes = length(nd), jaccard_vs_R0 = jacc(nd, N0))
}))
qb <- rbind(qb, loo)
write.csv(qb, paste0(out, "c8_qb_robustness.csv"), row.names = FALSE)
n_loo_ok <- sum(loo$jaccard_vs_R0 >= 0.8)
robust <- qb$jaccard_vs_R0[1] >= 0.8 && qb$jaccard_vs_R0[2] >= 0.8 && n_loo_ok >= 9
print(qb)
cat("Q-B: V1 J =", qb$jaccard_vs_R0[1], " V2 J =", qb$jaccard_vs_R0[2], " LOO >=0.8:", n_loo_ok, "/10 -> robust:", robust, "\n")
addm("V1_2of3", "jaccard_nodes_vs_R0", qb$jaccard_vs_R0[1], length(v1))
addm("V2_noTRK3", "jaccard_nodes_vs_R0", qb$jaccard_vs_R0[2], length(v2))
for (i in seq_len(nrow(loo))) addm(loo$run[i], "jaccard_nodes_vs_R0", loo$jaccard_vs_R0[i], loo$n_nodes[i])
addm("LOO", "n_runs_jaccard_ge_0.8", n_loo_ok, 10)
addm("LOO", "min_jaccard", min(loo$jaccard_vs_R0), 10)

## ---- Q-C informativeness ----
ord1 <- ordered_NPC(lay$L1[, "N"], lay$L1[, "P"], lay$L1[, "C"])
ord2 <- ordered_NPC(lay$L2[, "N"], lay$L2[, "P"], lay$L2[, "C"]); ok2 <- !is.na(lay$L2[, "N"] + lay$L2[, "P"] + lay$L2[, "C"])
ev3 <- lay$n_det >= 10
ord3 <- ordered_NPC(lay$L3[ev3, "N"], lay$L3[ev3, "P"], lay$L3[ev3, "C"])
br <- data.frame(layer = c("L1_kinase_activity", "L2_phosphosite", "L3_TF_targets"),
                 n_ordered = c(sum(ord1), sum(ord2[ok2]), sum(ord3)),
                 n_evaluable = c(length(ord1), sum(ok2), sum(ev3)))
br$base_rate <- br$n_ordered / br$n_evaluable
nE <- nrow(et); obs <- sum(et$accept_R0)
expd <- nE * prod(br$base_rate); ratio <- obs / expd
br$candidate_edges_ordered <- c(sum(et$ord_L1), sum(et$ord_L2), sum(et$ord_L3))
write.csv(br, paste0(out, "c8_qc_base_rates.csv"), row.names = FALSE)
qc <- data.frame(n_candidate_edges = nE, observed_accepted = obs, expected_independent = expd, ratio = ratio, informative = ratio >= 2)
write.csv(qc, paste0(out, "c8_qc_observed_vs_expected.csv"), row.names = FALSE)
print(br); print(qc)
for (i in 1:3) addm(br$layer[i], "base_rate_ordered", br$base_rate[i], br$n_evaluable[i])
addm("R0", "expected_accepted_independent", expd, nE); addm("R0", "observed_accepted", obs, nE)
addm("R0", "observed_over_expected", ratio, nE)

## ---- Q-D network measures ----
mk <- function(d) graph_from_data_frame(d, directed = FALSE)
cand <- data.frame(from = paste0("K:", et$Kinase), to = paste0("T:", et$TF))
a <- et[et$accept_R0, ]; accd <- data.frame(from = paste0("K:", a$Kinase), to = paste0("T:", a$TF))
gc <- mk(cand); ga <- mk(accd)    # multigraph: one edge per kinase-TF-site
dc <- degree(gc); da <- setNames(rep(0, length(dc)), names(dc)); da[names(degree(ga))] <- degree(ga)
bw <- betweenness(simplify(gc))
deg <- data.frame(node = names(dc), type = ifelse(startsWith(names(dc), "K:"), "kinase", "TF"),
                  degree_candidate = as.numeric(dc), degree_accepted = as.numeric(da[names(dc)]),
                  betweenness_candidate = as.numeric(bw[names(dc)]), accepted_node = names(dc) %in% N0)
write.csv(deg, paste0(out, "c8_qd_degrees.csv"), row.names = FALSE)
rho_all <- suppressWarnings(cor(deg$degree_candidate, deg$degree_accepted, method = "spearman"))
dd <- deg[deg$accepted_node, ]
rho_acc <- suppressWarnings(cor(dd$degree_candidate, dd$degree_accepted, method = "spearman"))
rho_bw <- suppressWarnings(cor(deg$betweenness_candidate, deg$degree_accepted, method = "spearman"))
cat("Q-D: Spearman degree cand vs acc (all", nrow(deg), "nodes) =", rho_all, "; accepted nodes only (", nrow(dd), ") =", rho_acc,
    "; betweenness(cand) vs accepted degree =", rho_bw, "\n")
addm("R0", "spearman_degree_cand_vs_acc_all_nodes", rho_all, nrow(deg))
addm("R0", "spearman_degree_cand_vs_acc_accepted_nodes", rho_acc, nrow(dd))
addm("R0", "spearman_betweenness_cand_vs_acc_degree", rho_bw, nrow(deg))

## ---- Q-E CRISPR (descriptive) ----
crispr <- data.frame(TF = c("MYC", "BRCA1", "EP300", "JUND"), CRISPR_hit = c(TRUE, FALSE, FALSE, FALSE))
crispr$in_candidate <- paste0("T:", crispr$TF) %in% names(dc)
crispr$accepted_R0 <- paste0("T:", crispr$TF) %in% N0
crispr$n_accepted_edges <- sapply(crispr$TF, function(t) sum(et$accept_R0 & et$TF == t))
crispr$n_candidate_edges <- sapply(crispr$TF, function(t) sum(et$TF == t))
write.csv(crispr, paste0(out, "c8_qe_crispr.csv"), row.names = FALSE)
print(crispr)
addm("R0", "crispr_TFs_accepted", sum(crispr$accepted_R0), 4)
addm("R0", "crispr_hit_MYC_accepted", as.numeric(crispr$accepted_R0[1]), 1)

append_metrics(do.call(rbind, setNames(mrows, NULL)))

## ---- figures ----
f2 <- function() {
  par(mfrow = c(1, 2), mar = c(8, 4, 3, 1))
  b <- barplot(br$base_rate, names.arg = c("kinase\nactivity", "phosphosite", "TF target\nset"), ylim = c(0, max(br$base_rate) * 1.25),
               ylab = "share ordered NTRK > PDGFRA > Cortex", main = "Base rate per layer", col = "grey70")
  text(b, br$base_rate, sprintf("%.3f\n(n=%d)", br$base_rate, br$n_evaluable), pos = 3, cex = 0.7, xpd = NA)
  b2 <- barplot(c(expd, obs), names.arg = c("expected\n(independent)", "observed\n(R0)"), col = c("grey70", "#d95f02"),
                ylim = c(0, max(expd, obs) * 1.25), ylab = "accepted edges", main = sprintf("Ratio obs/exp = %.2f", ratio))
  text(b2, c(expd, obs), sprintf("%.1f", c(expd, obs)), pos = 3)
}
png(paste0(out, "c8_fig2_base_rates.png"), 1400, 750, res = 130); f2(); dev.off()
pdf(paste0(out, "c8_fig2_base_rates.pdf"), 11, 6); f2(); dev.off()
f3 <- function() {
  par(mar = c(8, 4, 3, 1))
  b <- barplot(qb$jaccard_vs_R0, names.arg = qb$run, las = 2, ylim = c(0, 1.05), ylab = "Jaccard of accepted nodes vs R0",
               col = ifelse(qb$jaccard_vs_R0 >= 0.8, "grey60", "#d95f02"), main = "C8 robustness (V1, V2, leave-one-sample-out)")
  abline(h = 0.8, lty = 2, col = "red")
}
png(paste0(out, "c8_fig3_leave_one_out.png"), 1400, 750, res = 130); f3(); dev.off()
pdf(paste0(out, "c8_fig3_leave_one_out.pdf"), 11, 6); f3(); dev.off()
