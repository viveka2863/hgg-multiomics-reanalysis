# ======================================================================
# C1: do the 10 samples separate by genotype, and does that depend on the method?
# Run from the repo root:  Rscript analysis/C1_exploration/improved/c1_exploration.R
# Layers: whole proteome and phosphoproteome (RNA-seq is not in the repo; see C11).
# Pre-registration is in LOG.md.
# ======================================================================
options(stringsAsFactors = FALSE)
set.seed(1)
stopifnot(file.exists("CLAUDE.md"))
OUT <- "analysis/C1_exploration/improved"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
say <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), sprintf(...))

## ---------------- helpers ----------------
ari <- function(a, b) {
  ok <- !is.na(a) & !is.na(b); a <- a[ok]; b <- b[ok]
  tab <- table(a, b); n <- sum(tab)
  c2 <- function(x) x * (x - 1) / 2
  s <- sum(c2(tab)); sa <- sum(c2(rowSums(tab))); sb <- sum(c2(colSums(tab)))
  e <- sa * sb / c2(n)
  out <- (s - e) / ((sa + sb) / 2 - e)
  if (is.nan(out)) NA_real_ else out
}
save_fig <- function(name, fun, w, h) {
  png(file.path(OUT, paste0(name, ".png")), width = w, height = h, units = "in", res = 150); fun(); invisible(dev.off())
}
METRICS <- list()
add_metric <- function(variant, metric, value, n)
  METRICS[[length(METRICS) + 1]] <<- data.frame(claim = "C1", variant = variant, metric = metric, value = as.numeric(value), n = n)

## ---------------- data ----------------
samples <- c("Cortex1","Cortex2","Cortex3","PD1","PD2","PD3","PD4","TRK1","TRK2","TRK3")
group   <- setNames(sub("[0-9]+$", "", samples), samples)       # Cortex / PD / TRK
tum     <- samples[group[samples] != "Cortex"]
map <- c(PDGFRA_1a = "PD1", PDGFRA_1b = "PD2", PDGFRA_2 = "PD3", PDGFRA_3 = "PD4",
         NTRK_1 = "TRK1", NTRK_2 = "TRK2", NTRK_3 = "TRK3",
         CNTRL_1 = "Cortex1", CNTRL_2 = "Cortex2", CNTRL_3 = "Cortex3")   # from analysis/00_sample_key
read_layer <- function(path, drop_contam) {
  d <- read.table(path, header = TRUE, sep = "\t", quote = "", comment.char = "", check.names = FALSE)
  if (drop_contam) d <- d[!grepl("CON_", d[[1]]), ]
  m <- as.matrix(d[, -1]); storage.mode(m) <- "numeric"
  rownames(m) <- make.unique(as.character(d[[1]]))
  stopifnot(!anyNA(map[colnames(m)]))
  colnames(m) <- map[colnames(m)]
  m[, samples]
}
say("Loading layers")
L <- list(proteome = read_layer("JUMPn/JUMPnk/Integration/whl_hgg_v5_20160513.txt", FALSE),
          phospho  = read_layer("JUMPn/JUMPnk/Integration/pho_hgg_uniq_v5_20160513.txt", TRUE))
stopifnot(all(sapply(L, function(m) all(m > 0) && !anyNA(m))))
cat(sprintf("proteome: %d features; phospho: %d features (contaminants removed)\n", nrow(L$proteome), nrow(L$phospho)))

tf  <- function(m, t) if (t == "log2") log2(m) else m
sel <- function(m, f) if (f == "all") m else m[order(apply(m, 1, mad), decreasing = TRUE)[1:1000], , drop = FALSE]
prep <- function(layer, transform, features, cols) {
  X <- sel(tf(L[[layer]][, cols, drop = FALSE], transform), features)
  X[apply(X, 1, sd) > 0, , drop = FALSE]
}
run_hc <- function(X, dist_type, linkage, k) {
  d <- if (dist_type == "euclid_z") { Z <- scale(t(X)); Z[is.na(Z)] <- 0; dist(Z) } else as.dist(1 - cor(X))
  cutree(hclust(d, method = linkage), k = k)
}

## ---------------- Q-A / Q-B / Q-C: PCA variants ----------------
say("PCA variants")
grid_pca <- expand.grid(layer = names(L), transform = c("log2", "raw"), features = c("all", "top1000"),
                        scaling = c("centred", "scaled"), stringsAsFactors = FALSE)
PCA <- list()
for (i in seq_len(nrow(grid_pca))) {
  g <- grid_pca[i, ]; sc <- g$scaling == "scaled"
  X <- prep(g$layer, g$transform, g$features, samples)
  p <- prcomp(t(X), center = TRUE, scale. = sc); ve <- 100 * p$sdev^2 / sum(p$sdev^2)
  set.seed(1); km <- kmeans(p$x[, 1:2], centers = 3, nstart = 25)
  Xt <- prep(g$layer, g$transform, g$features, tum)
  pt <- prcomp(t(Xt), center = TRUE, scale. = sc)
  set.seed(1); kt <- kmeans(pt$x[, 1:2], centers = 2, nstart = 25)
  nm <- paste("pca", g$layer, g$transform, g$features, g$scaling, sep = "_")
  PCA[[i]] <- data.frame(variant = nm, layer = g$layer, transform = g$transform, features = g$features,
                         scaling = g$scaling, PC1_pct = ve[1], PC2_pct = ve[2],
                         ARI_all3 = ari(km$cluster, group[rownames(p$x)]),
                         ARI_tumour2 = ari(kt$cluster, group[rownames(pt$x)]))
}
PCA <- do.call(rbind, PCA)

## ---------------- hierarchical variants ----------------
say("Hierarchical variants")
grid_hc <- expand.grid(layer = names(L), transform = c("log2", "raw"), features = c("all", "top1000"),
                       dist = c("euclid_z", "cor"), linkage = c("ward.D", "ward.D2", "average", "complete"),
                       stringsAsFactors = FALSE)
HC <- list()
for (i in seq_len(nrow(grid_hc))) {
  g <- grid_hc[i, ]
  X  <- prep(g$layer, g$transform, g$features, samples)
  Xt <- prep(g$layer, g$transform, g$features, tum)
  HC[[i]] <- data.frame(variant = paste("hc", g$layer, g$transform, g$features, g$dist, g$linkage, sep = "_"),
                        layer = g$layer, transform = g$transform, features = g$features, dist = g$dist, linkage = g$linkage,
                        ARI_all3 = ari(run_hc(X, g$dist, g$linkage, 3), group[colnames(X)]),
                        ARI_tumour2 = ari(run_hc(Xt, g$dist, g$linkage, 2), group[colnames(Xt)]))
}
HC <- do.call(rbind, HC)
write.csv(PCA, file.path(OUT, "c1_pca_variants.csv"), row.names = FALSE)
write.csv(HC,  file.path(OUT, "c1_hclust_variants.csv"), row.names = FALSE)
for (i in seq_len(nrow(PCA))) {
  add_metric(PCA$variant[i], "pc1_pct", PCA$PC1_pct[i], 10); add_metric(PCA$variant[i], "pc2_pct", PCA$PC2_pct[i], 10)
  add_metric(PCA$variant[i], "ari_all3", PCA$ARI_all3[i], 10); add_metric(PCA$variant[i], "ari_tumour2", PCA$ARI_tumour2[i], 7)
}
for (i in seq_len(nrow(HC))) {
  add_metric(HC$variant[i], "ari_all3", HC$ARI_all3[i], 10); add_metric(HC$variant[i], "ari_tumour2", HC$ARI_tumour2[i], 7)
}

## ---------------- Q-D: leave-one-sample-out ----------------
say("Leave-one-out")
LOO <- NULL
for (ly in names(L)) {
  for (i in seq_along(samples)) {
    keep <- samples[-i]
    X1 <- prep(ly, "log2", "top1000", keep)
    a1 <- ari(run_hc(X1, "euclid_z", "ward.D", 3), group[keep])
    X2 <- prep(ly, "log2", "all", keep)
    p2 <- prcomp(t(X2), center = TRUE, scale. = FALSE)
    set.seed(1); k2 <- kmeans(p2$x[, 1:2], centers = 3, nstart = 25)
    a2 <- ari(k2$cluster, group[keep])
    LOO <- rbind(LOO, data.frame(layer = ly, dropped = samples[i], ARI_P1_hc = a1, ARI_P2_pca_kmeans = a2))
  }
}
write.csv(LOO, file.path(OUT, "c1_loo_stability.csv"), row.names = FALSE)
for (i in seq_len(nrow(LOO))) {
  add_metric(paste("loo", LOO$layer[i], "P1_hc", "drop", LOO$dropped[i], sep = "_"), "ari", LOO$ARI_P1_hc[i], 9)
  add_metric(paste("loo", LOO$layer[i], "P2_pca", "drop", LOO$dropped[i], sep = "_"), "ari", LOO$ARI_P2_pca_kmeans[i], 9)
}

## ---------------- Q-E: replicate check ----------------
say("Replicate check")
prs <- t(combn(samples, 2)); REP <- NULL; MD <- list()
for (ly in names(L)) {
  X <- log2(L[[ly]])
  md <- apply(prs, 1, function(p) median(abs(X[, p[1]] - X[, p[2]])))
  REP <- rbind(REP, data.frame(layer = ly, sample_i = prs[, 1], sample_j = prs[, 2],
                               same_group = group[prs[, 1]] == group[prs[, 2]], median_abs_log2_diff = md))
  M <- matrix(0, 10, 10, dimnames = list(samples, samples))
  for (k in seq_len(nrow(prs))) { M[prs[k, 1], prs[k, 2]] <- md[k]; M[prs[k, 2], prs[k, 1]] <- md[k] }
  MD[[ly]] <- M
}
write.csv(REP, file.path(OUT, "c1_replicate_check.csv"), row.names = FALSE)
REPSUM <- NULL
for (ly in names(L)) {
  w <- REP[REP$layer == ly & REP$same_group, ]
  w$rank <- rank(w$median_abs_log2_diff)
  r12 <- w[w$sample_i == "PD1" & w$sample_j == "PD2", ]
  low <- w[which.min(w$median_abs_log2_diff), ]
  REPSUM <- rbind(REPSUM, data.frame(layer = ly, PD1_PD2_diff = r12$median_abs_log2_diff, PD1_PD2_rank_of_12 = r12$rank,
                                     lowest_pair = paste(low$sample_i, low$sample_j, sep = "-"),
                                     lowest_diff = low$median_abs_log2_diff, flag_replicates = r12$rank == 1))
  add_metric(paste0("replicate_", ly), "PD1_PD2_median_abs_log2_diff", r12$median_abs_log2_diff, 45)
  add_metric(paste0("replicate_", ly), "PD1_PD2_rank_among_12_within_group_pairs", r12$rank, 12)
}
write.csv(REPSUM, file.path(OUT, "c1_replicate_summary.csv"), row.names = FALSE)

## ---------------- metrics file ----------------
M <- do.call(rbind, METRICS); mf <- "results/metrics.csv"
if (file.exists(mf)) {
  old <- read.csv(mf)
  old <- old[!(old$claim == "C1" & old$variant %in% M$variant), , drop = FALSE]
  M <- rbind(old[, names(M)], M)
}
write.csv(M, mf, row.names = FALSE)

## ---------------- figures ----------------
say("Figures")
COL <- c(Cortex = "purple", PD = "deepskyblue3", TRK = "red3")
PAPER <- lapply(names(L), function(ly) {
  X <- prep(ly, "log2", "all", samples); Xh <- prep(ly, "log2", "top1000", samples)
  Z <- scale(t(Xh)); Z[is.na(Z)] <- 0
  list(pca = prcomp(t(X), center = TRUE, scale. = FALSE), hc = hclust(dist(Z), method = "ward.D"))
}); names(PAPER) <- names(L)
save_fig("c1_paper_setup", function() {
  par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))
  for (ly in names(PAPER)) {
    p <- PAPER[[ly]]$pca; ve <- 100 * p$sdev^2 / sum(p$sdev^2)
    plot(p$x[, 1], p$x[, 2], pch = 16, cex = 1.6, col = COL[group[rownames(p$x)]],
         xlab = sprintf("PC1 (%.0f%%)", ve[1]), ylab = sprintf("PC2 (%.0f%%)", ve[2]), main = paste(ly, "PCA (log2, all features)"))
    text(p$x[, 1], p$x[, 2], rownames(p$x), pos = 3, cex = 0.7)
  }
  for (ly in names(PAPER)) plot(PAPER[[ly]]$hc, main = paste(ly, "Ward.D, top-1000 MAD, z-scored"), xlab = "", sub = "", cex = 0.8)
}, w = 10, h = 9)

frac1 <- function(x) mean(x > 0.999, na.rm = TRUE)
FR <- sapply(names(L), function(ly) c(
  "PCA+k-means, all samples" = frac1(PCA$ARI_all3[PCA$layer == ly]),
  "PCA+k-means, tumours only" = frac1(PCA$ARI_tumour2[PCA$layer == ly]),
  "Hierarchical, all samples" = frac1(HC$ARI_all3[HC$layer == ly]),
  "Hierarchical, tumours only" = frac1(HC$ARI_tumour2[HC$layer == ly])))
save_fig("c1_variant_recovery", function() {
  par(mar = c(5, 4.5, 3, 1))
  barplot(t(FR), beside = TRUE, ylim = c(0, 1.1), las = 2, cex.names = 0.7, col = c("grey40", "grey80"),
          ylab = "fraction of variants with ARI = 1", main = "C1: how often is the genotype recovered?")
  abline(h = 0.9, lty = 2, col = "red"); legend("bottomright", legend = names(L), fill = c("grey40", "grey80"), cex = 0.8, bty = "n")
}, w = 9, h = 5.5)
save_fig("c1_sample_distance", function() {
  par(mfrow = c(1, 2), mar = c(6, 6, 3, 1))
  for (ly in names(MD)) {
    M <- MD[[ly]]; n <- nrow(M)
    image(1:n, 1:n, M[, n:1], col = hcl.colors(50, "viridis", rev = TRUE), axes = FALSE, xlab = "", ylab = "",
          main = paste(ly, ": median |log2 difference|"))
    axis(1, 1:n, colnames(M), las = 2, cex.axis = 0.7); axis(2, 1:n, rev(rownames(M)), las = 2, cex.axis = 0.7)
    for (i in 1:n) for (j in 1:n) text(i, n - j + 1, sprintf("%.2f", M[i, j]), cex = 0.5)
  }
}, w = 12, h = 6)

## ---------------- print summary ----------------
cat("\n================ SUMMARY (paste this back) ================\n")
cat("\nPaper-setup PCA (log2, all features, centred): % variance PC1 / PC2 (paper: proteome 76 / 12, phospho 66 / 15)\n")
print(PCA[PCA$transform == "log2" & PCA$features == "all" & PCA$scaling == "centred", c("layer", "PC1_pct", "PC2_pct", "ARI_all3", "ARI_tumour2")], row.names = FALSE, digits = 4)
cat("\nPaper-setup hierarchical (log2, top1000, euclid_z, ward.D):\n")
print(HC[HC$transform == "log2" & HC$features == "top1000" & HC$dist == "euclid_z" & HC$linkage == "ward.D", c("layer", "ARI_all3", "ARI_tumour2")], row.names = FALSE, digits = 4)
cat("\nFraction of variants with ARI = 1 (rule: >= 0.9):\n"); print(round(FR, 3))
cat("\nVariants with ARI < 1 (all samples or tumours only):\n")
bad <- rbind(PCA[PCA$ARI_all3 < 0.999 | PCA$ARI_tumour2 < 0.999, c("variant", "ARI_all3", "ARI_tumour2")],
             HC[HC$ARI_all3 < 0.999 | HC$ARI_tumour2 < 0.999, c("variant", "ARI_all3", "ARI_tumour2")])
if (nrow(bad)) print(bad, row.names = FALSE, digits = 3) else cat("none\n")
cat("\nLeave-one-out ARI (rule: >= 0.8 in >= 9 of 10):\n"); print(LOO, row.names = FALSE, digits = 3)
for (ly in names(L)) cat(sprintf("  %s: P1 runs >= 0.8: %d/10; P2 runs >= 0.8: %d/10\n", ly,
    sum(LOO$ARI_P1_hc[LOO$layer == ly] >= 0.8), sum(LOO$ARI_P2_pca_kmeans[LOO$layer == ly] >= 0.8)))
cat("\nReplicate check (PD1-PD2 vs other within-group pairs):\n"); print(REPSUM, row.names = FALSE, digits = 4)
cat("\nOutputs written to", OUT, "and results/metrics.csv\n")