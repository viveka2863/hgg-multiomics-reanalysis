# C8 part B: Q-H1 label permutation, Q-H2 global offset / centred R0, Q-H3 circularity. Pre-registered in LOG.md. Run from repo root.
set.seed(1)
source("analysis/C8_network/c8_lib.R")
out <- "analysis/C8_network/improved/"
D <- load_all(); E <- candidates(D)
lay <- layers(D, S); et <- edge_table(D, E, lay)
stopifnot(nrow(E) == 28, sum(et$accept_R0) == 15)
paper <- c("T:MYC", "T:JUND", "T:EP300", "T:BRCA1", "K:AKT1", "K:MAPK3", "K:CDK5")
n_paper <- function(acc_idx) { nd <- unique(c(paste0("K:", E$Kinase[acc_idx]), paste0("T:", E$TF[acc_idx]))); sum(paper %in% nd) }
mrows <- list(); addm <- function(variant, metric, value, n) mrows[[length(mrows) + 1]] <<- data.frame("C8", variant, metric, value, n)

## per-sample values of the three layers for the 28 candidate edges (labels only enter through group means)
per_sample <- function(D) {
  X1 <- D$km[E$Kinase, S, drop = FALSE]
  v <- D$pm[D$site_long$row, S, drop = FALSE]
  agg <- rowsum(v, D$site_long$key, na.rm = FALSE) / as.vector(table(D$site_long$key)[rownames(rowsum(v, D$site_long$key))])
  X2 <- agg[E$key, S, drop = FALSE]
  z <- t(scale(t(D$prm[, S]))); z <- z[!is.na(z[, 1]), , drop = FALSE]
  sc <- t(sapply(names(D$tgm), function(t) if (length(D$tgm[[t]]) < 10) rep(NA_real_, 10) else colMeans(z[D$tgm[[t]], , drop = FALSE])))
  colnames(sc) <- S
  list(X1 = X1, X2 = X2, X3 = sc[E$TF, S, drop = FALSE])
}
R0_idx <- function(X, idxN, idxP, idxC) {
  ok <- rep(TRUE, nrow(X$X1))
  for (m in X) ok <- ok & ordered_NPC(rowMeans(m[, idxN, drop = FALSE]), rowMeans(m[, idxP, drop = FALSE]), rowMeans(m[, idxC, drop = FALSE]))
  which(ok)
}
X <- per_sample(D)

## ---- Q-H1 ----
pd <- 1:4; trk <- 5:7; ctx <- 8:10
cc <- combn(10, 3); perm <- list()
for (i in seq_len(ncol(cc))) {
  rest <- setdiff(1:10, cc[, i]); pp <- combn(rest, 4)
  for (j in seq_len(ncol(pp))) perm[[length(perm) + 1]] <- list(C = cc[, i], P = pp[, j], N = setdiff(rest, pp[, j]))
}
stopifnot(length(perm) == 4200)
res <- do.call(rbind, lapply(seq_along(perm), function(k) {
  p <- perm[[k]]; a <- R0_idx(X, p$N, p$P, p$C)
  data.frame(perm = k, PDGFRA = paste(p$P, collapse = ""), NTRK = paste(p$N, collapse = ""), Cortex = paste(p$C, collapse = ""),
             is_identity = setequal(p$P, pd) && setequal(p$N, trk) && setequal(p$C, ctx),
             n_edges = length(a), n_paper_nodes = n_paper(a))
}))
write.csv(res, paste0(out, "c8b_h1_permutation_distribution.csv"), row.names = FALSE)
idrow <- res[res$is_identity, ]; stopifnot(nrow(idrow) == 1, idrow$n_edges == 15, idrow$n_paper_nodes == 7)
no_id <- res[!res$is_identity, ]
p_edges_with <- mean(res$n_edges >= 15); p_edges_wo <- mean(no_id$n_edges >= 15)
p_nodes_with <- mean(res$n_paper_nodes == 7); p_nodes_wo <- mean(no_id$n_paper_nodes == 7)
cat(sprintf("Q-H1: observed 15 edges / 7 nodes; perm edges mean %.2f median %g max %d\n", mean(res$n_edges), median(res$n_edges), max(res$n_edges)))
cat(sprintf("  p(edges>=15) with identity %.5f (%d/4200), without %.5f (%d/4199)\n", p_edges_with, sum(res$n_edges >= 15), p_edges_wo, sum(no_id$n_edges >= 15)))
cat(sprintf("  p(7/7 nodes) with identity %.5f (%d/4200), without %.5f (%d/4199)\n", p_nodes_with, sum(res$n_paper_nodes == 7), p_nodes_wo, sum(no_id$n_paper_nodes == 7)))
print(table(res$n_edges))
addm("partB_H1", "perm_p_edges_ge15_with_identity", p_edges_with, 4200)
addm("partB_H1", "perm_p_edges_ge15_without_identity", p_edges_wo, 4199)
addm("partB_H1", "perm_p_nodes_7of7_with_identity", p_nodes_with, 4200)
addm("partB_H1", "perm_p_nodes_7of7_without_identity", p_nodes_wo, 4199)
addm("partB_H1", "perm_median_edges", median(res$n_edges), 4200)
addm("partB_H1", "perm_max_edges", max(res$n_edges), 4200)

## ---- Q-H2 ----
grpm <- function(M) { cm <- M - rowMeans(M); sapply(c(PDGFRA = "P", NTRK = "N", Cortex = "C"), function(g) mean(cm[, grp_of(S) == g])) }
Mph <- D$pm[rowSums(is.na(D$pm)) == 0, S]; Mpr <- D$prm[, S]
off <- rbind(data.frame(layer = "proteome", t(grpm(Mpr)), n_features = nrow(Mpr)),
             data.frame(layer = "phospho", t(grpm(Mph)), n_features = nrow(Mph)))
write.csv(off, paste0(out, "c8b_h2_group_offsets.csv"), row.names = FALSE); print(off)
for (i in 1:2) for (g in c("PDGFRA", "NTRK", "Cortex")) addm(paste0("partB_H2_", off$layer[i]), paste0("group_offset_", g), off[i, g], off$n_features[i])
# centre each sample by its median over all features of that layer
D2 <- D
D2$pm <- sweep(D$pm, 2, apply(D$pm, 2, median, na.rm = TRUE))
D2$prm <- sweep(D$prm, 2, apply(D$prm, 2, median))
X2c <- per_sample(D2)
lay2 <- layers(D2, S); et2 <- edge_table(D2, E, lay2)
a2 <- which(et2$accept_R0)
stopifnot(identical(unname(a2), unname(R0_idx(X2c, 5:7, 1:4, 8:10))))
write.csv(et2, paste0(out, "c8b_h2_centred_R0_edges.csv"), row.names = FALSE)
nd2 <- nodes_of(et2, "accept_R0"); N0 <- nodes_of(et, "accept_R0")
cat(sprintf("Q-H2: centred R0 accepts %d edges (uncentred 15); paper nodes %d/7; node Jaccard vs R0 %.3f\n", length(a2), n_paper(a2), jacc(nd2, N0)))
cat("  ordered-layer counts centred:", sum(et2$ord_L1), sum(et2$ord_L2), sum(et2$ord_L3), " uncentred:", sum(et$ord_L1), sum(et$ord_L2), sum(et$ord_L3), "\n")
br <- function(lay) {
  o1 <- ordered_NPC(lay$L1[, "N"], lay$L1[, "P"], lay$L1[, "C"])
  o2 <- ordered_NPC(lay$L2[, "N"], lay$L2[, "P"], lay$L2[, "C"]); ok2 <- !is.na(lay$L2[, "N"] + lay$L2[, "P"] + lay$L2[, "C"])
  ev <- lay$n_det >= 10; o3 <- ordered_NPC(lay$L3[ev, "N"], lay$L3[ev, "P"], lay$L3[ev, "C"])
  c(sum(o1), length(o1), sum(o2[ok2]), sum(ok2), sum(o3), sum(ev))
}
b0 <- br(lay); b1 <- br(lay2)
bt <- data.frame(layer = rep(c("L1_kinase_activity", "L2_phosphosite", "L3_TF_targets"), 2), condition = rep(c("uncentred", "centred"), each = 3),
                 n_ordered = c(b0[c(1, 3, 5)], b1[c(1, 3, 5)]), n_evaluable = c(b0[c(2, 4, 6)], b1[c(2, 4, 6)]))
bt$base_rate <- bt$n_ordered / bt$n_evaluable
write.csv(bt, paste0(out, "c8b_h2_base_rates.csv"), row.names = FALSE); print(bt)
exp_c <- 28 * prod(bt$base_rate[4:6])
cat(sprintf("  expected accepted (independent, centred base rates) %.2f; observed centred %d\n", exp_c, length(a2)))
addm("partB_H2_centred", "n_accepted_edges", length(a2), 28)
addm("partB_H2_centred", "paper_nodes_recovered", n_paper(a2), 7)
addm("partB_H2_centred", "jaccard_nodes_vs_R0", jacc(nd2, N0), length(nd2))
addm("partB_H2_centred", "expected_accepted_independent", exp_c, 28)
for (i in 1:6) addm(paste0("partB_H2_", bt$condition[i]), paste0("base_rate_", bt$layer[i]), bt$base_rate[i], bt$n_evaluable[i])
ed <- data.frame(edge = paste0(E$Kinase, ">", E$TF, "_", E$site), R0 = et$accept_R0, centred = et2$accept_R0)
write.csv(ed, paste0(out, "c8b_h2_edges_R0_vs_centred.csv"), row.names = FALSE)

## ---- Q-H3 ----
mx <- read.delim("JUMPn/JUMPnk/PhosphoWholeNormalization/hgg_phosphosite_matrix_normalized_whole_v1_20160420.txt", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(identical(names(mx), c("Phosphosite", paste0("Data", 1:10), "KINASE")), nrow(mx) == 984)
parts <- strsplit(mx$Phosphosite, "_"); stopifnot(all(lengths(parts) == 3))
mkey <- sapply(parts, function(p) paste(p[2], p[3], sep = "_"))
kin <- strsplit(toupper(mx$KINASE), "[,/;[:space:]]+")
pairs <- data.frame(key = rep(mkey, lengths(kin)), K = unlist(kin), stringsAsFactors = FALSE)
epair <- paste(E$key, toupper(E$Kinase)); hit <- epair %in% paste(pairs$key, pairs$K)
site_in <- E$key %in% mkey
ed3 <- data.frame(Kinase = E$Kinase, TF = E$TF, site = E$site, key = E$key, site_in_matrix = site_in, kinase_listed_with_site = hit, accept_R0 = et$accept_R0)
write.csv(ed3, paste0(out, "c8b_h3_circularity.csv"), row.names = FALSE)
cat(sprintf("Q-H3: %d/28 candidate edges have site+kinase together in the IKAP matrix KINASE column (site present for %d/28); among R0-accepted %d/15\n",
            sum(hit), sum(site_in), sum(hit & et$accept_R0)))
addm("partB_H3", "candidate_edges_site_kinase_in_matrix", sum(hit), 28)
addm("partB_H3", "candidate_edges_site_in_matrix", sum(site_in), 28)
addm("partB_H3", "accepted_edges_site_kinase_in_matrix", sum(hit & et$accept_R0), 15)

## ---- figures ----
f1 <- function() {
  tb <- table(factor(res$n_edges, levels = 0:max(c(res$n_edges, 15))))
  b <- barplot(tb, xlab = "accepted edges under R0 (28 candidates)", ylab = "label assignments (of 4,200)", col = "grey70",
               main = sprintf("Label-permutation null; p(>=15) = %.4f with / %.4f without identity", p_edges_with, p_edges_wo))
  abline(v = b[16], col = "#d95f02", lwd = 2); text(b[16], max(tb) * 0.9, "observed = 15", pos = 4, col = "#d95f02")
}
f2 <- function() {
  m <- t(as.matrix(off[, c("Cortex", "PDGFRA", "NTRK")])); colnames(m) <- off$layer
  barplot(m, beside = TRUE, col = c("grey60", "#7570b3", "#d95f02"), ylab = "mean (log2 - feature mean)", main = "Per-group offsets",
          legend.text = rownames(m), args.legend = list(x = "topleft")); abline(h = 0)
}
f3 <- function() {
  tb <- data.frame(R0 = c(sum(et$accept_R0), sum(et$ord_L1), sum(et$ord_L2), sum(et$ord_L3)),
                   centred = c(length(a2), sum(et2$ord_L1), sum(et2$ord_L2), sum(et2$ord_L3)))
  m <- t(as.matrix(tb)); colnames(m) <- c("accepted (all 3)", "L1 ordered", "L2 ordered", "L3 ordered")
  b <- barplot(m, beside = TRUE, col = c("#d95f02", "grey60"), ylim = c(0, 30), ylab = "candidate edges (of 28)", main = "R0 vs centred R0",
               legend.text = rownames(m), args.legend = list(x = "topright"))
  text(b, m, m, pos = 3)
}
for (nm in c("f1", "f2", "f3")) {
  fn <- c(f1 = "c8b_fig1_permutation", f2 = "c8b_fig2_group_offsets", f3 = "c8b_fig3_R0_vs_centred")[[nm]]
  png(paste0(out, fn, ".png"), 1400, 750, res = 130); get(nm)(); dev.off()
  pdf(paste0(out, fn, ".pdf"), 11, 6); get(nm)(); dev.off()
}
append_metrics(do.call(rbind, setNames(mrows, NULL)))
dir.create("env", showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), "env/sessionInfo_C8b.txt")
