# ======================================================================
# C3 compare arm: does the clustering of the 4,703 DE proteins depend on analysis choices?
# Run from the repo root:  Rscript analysis/C3_coexpression/improved/c3_compare.R
# Quick test (skips slow parts):  C3_QUICK=1 Rscript analysis/C3_coexpression/improved/c3_compare.R
# Pre-registered questions Q-A, Q-B, Q-C are in LOG.md.
# ======================================================================
suppressPackageStartupMessages(library(WGCNA))
options(stringsAsFactors = FALSE)
set.seed(1)
stopifnot(file.exists("CLAUDE.md"))

QUICK       <- Sys.getenv("C3_QUICK") == "1"
RUN_ALLFEAT <- !QUICK
RUN_LOO     <- !QUICK
ALLFEAT_N   <- 8000
POWER       <- 16
MIN_MODULE  <- 100
OUT <- "analysis/C3_coexpression/improved"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
say <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), sprintf(...))
try(enableWGCNAThreads(), silent = TRUE)

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
pct_identical <- function(lab, ref) {
  tab <- table(lab, ref)
  maj <- colnames(tab)[apply(tab, 1, which.max)]
  names(maj) <- rownames(tab)
  mean(maj[as.character(lab)] == as.character(ref))
}
# Authors' pipeline (wgcna_v1.0.R): signed adjacency -> TOM -> average-linkage hclust ->
# cutreeDynamic -> mergeCloseModules. TOM is computed once; the grid only changes deepSplit and merge height.
run_wgcna_grid <- function(expr, power, deepSplits, mergeHeights, minModuleSize) {
  adj  <- adjacency(expr, power = power, type = "signed")
  TOM  <- TOMsimilarity(adj); rm(adj)
  diss <- 1 - TOM; rm(TOM)
  tree <- hclust(as.dist(diss), method = "average")
  out <- list()
  for (ds in deepSplits) {
    mods <- cutreeDynamic(dendro = tree, distM = diss, deepSplit = ds,
                          pamRespectsDendro = FALSE, minClusterSize = minModuleSize, verbose = 0)
    cols <- labels2colors(mods)
    for (mh in mergeHeights) {
      mg <- mergeCloseModules(expr, cols, cutHeight = mh, verbose = 0)
      v <- mg$colors; names(v) <- colnames(expr)
      out[[sprintf("ds%d_mh%.2f", ds, mh)]] <- v
    }
  }
  rm(diss, tree); invisible(gc(verbose = FALSE))
  out
}
# Q-B: per-cluster group means of z-scored intensities, and overlap with the paper's WP-C1
qb_clusters <- function(lab, z, variant, wpc1, de_acc) {
  grp <- sub("[0-9]+$", "", colnames(z))
  cls <- setdiff(unique(lab), c("grey", NA))
  rows <- lapply(cls, function(cl) {
    mem    <- intersect(names(lab)[which(lab == cl)], rownames(z))
    mem_de <- intersect(mem, de_acc)
    m <- colMeans(z[mem, , drop = FALSE])
    g <- tapply(m, grp, mean)
    data.frame(variant = variant, cluster = cl, n_members = length(mem), n_DE = length(mem_de),
               mean_Cortex = unname(g["Cortex"]), mean_PD = unname(g["PD"]), mean_TRK = unname(g["TRK"]),
               rising = isTRUE(g["Cortex"] < g["PD"] && g["PD"] < g["TRK"]),
               jaccard_WPC1 = length(intersect(mem_de, wpc1)) / length(union(mem_de, wpc1)))
  })
  do.call(rbind, rows)
}
save_fig <- function(name, fun, w = 8, h = 6) {
  png(file.path(OUT, paste0(name, ".png")), width = w, height = h, units = "in", res = 150)
  fun(); invisible(dev.off())
  pdf(file.path(OUT, paste0(name, ".pdf")), width = w, height = h)
  fun(); invisible(dev.off())
}
METRICS <- list()
add_metric <- function(variant, metric, value, n) {
  METRICS[[length(METRICS) + 1]] <<- data.frame(claim = "C3", variant = variant, metric = metric,
                                                value = as.numeric(value), n = n)
}

## ---------------- load data ----------------
say("Loading data")
inp <- read.table("JUMPn/JUMPnc/HGG_whl_input.txt", header = TRUE, sep = "\t", quote = "",
                  comment.char = "", check.names = FALSE)
lab <- read.csv("analysis/C3_coexpression/original/labels_with_accession.csv")
stopifnot(nrow(inp) == 4703, identical(lab$Protein_Accession, inp$Protein_Accession))
raw <- read.table("JUMPn/JUMPnk/PhosphoWholeNormalization/hgg_whole_v1_20160420.txt", header = TRUE,
                  sep = "\t", quote = "", comment.char = "", check.names = FALSE)
samples <- c("Cortex1","Cortex2","Cortex3","PD1","PD2","PD3","PD4","TRK1","TRK2","TRK3")
# TMT channel of each sample, from analysis/00_sample_key (verified by data there)
chan <- c(Cortex1="sig130N", Cortex2="sig130C", Cortex3="sig131",
          PD1="sig126", PD2="sig127N", PD3="sig127C", PD4="sig128N",
          TRK1="sig128C", TRK2="sig129N", TRK3="sig129C")
paper   <- setNames(lab$paper,   lab$Protein_Accession)
shipped <- setNames(lab$shipped, lab$Protein_Accession)
wpc1    <- inp$Protein_Accession[inp$Cluster == "WP-C1"]
de_acc  <- inp$Protein_Accession
expr_ctrl <- t(as.matrix(inp[, samples])); colnames(expr_ctrl) <- inp$Protein_Accession
stopifnot(all(expr_ctrl > 0))

# match input proteins to the unscaled matrix, prove the channel mapping, estimate the per-sample rescaling
acc_raw <- raw[[2]]
idx <- match(inp$Protein_Accession, acc_raw); ok <- !is.na(idx)
say("Input proteins not found in hgg_whole_v1 by accession: %d", sum(!ok))
print(inp$Protein_Accession[!ok])
rawm <- as.matrix(raw[idx[ok], unname(chan[samples])]); colnames(rawm) <- samples
storage.mode(rawm) <- "numeric"
rho <- sapply(samples, function(s) cor(inp[ok, s], rawm[, s], method = "spearman"))
sc  <- sapply(samples, function(s) median(inp[ok, s] / rawm[, s]))
cat("\nChannel-mapping proof (Spearman, input vs raw channel) and estimated rescaling factors:\n")
print(round(rbind(spearman = rho, scale_factor = sc), 4))
stopifnot(all(rho > 0.99))
for (s in samples) add_metric("scale", paste0("scale_factor_", s), sc[s], sum(ok))

## ---------------- control + WGCNA grid ----------------
say("Control + grid (TOM computed once)")
grid <- run_wgcna_grid(expr_ctrl, POWER, 1:3, c(0.15, 0.20, 0.25), MIN_MODULE)
names(grid) <- paste0("grid_", names(grid))
control <- grid[["grid_ds2_mh0.20"]]; grid[["grid_ds2_mh0.20"]] <- NULL
ctrl_ari <- ari(control, shipped[names(control)])
cat(sprintf("\nCONTROL CHECK: ARI vs shipped labels = %.4f\n", ctrl_ari))
if (is.na(ctrl_ari) || ctrl_ari < 0.999) stop("Control does not reproduce the shipped labels (ARI < 0.999). Stopping.")

LABS <- c(list(paper = paper, shipped = shipped, control = control), grid)

## ---------------- other variants ----------------
say("Variant: log2 input")
LABS$log2 <- run_wgcna_grid(log2(expr_ctrl), POWER, 2, 0.20, MIN_MODULE)[[1]]

say("Variant: unscaled input (hgg_whole_v1)")
exprU <- t(rawm); colnames(exprU) <- inp$Protein_Accession[ok]
LABS$unscaled <- run_wgcna_grid(exprU, POWER, 2, 0.20, MIN_MODULE)[[1]]

say("Variants: hierarchical (Ward.D2 on 1 - Pearson) and k-means, k = 5")
hc <- hclust(as.dist(1 - cor(expr_ctrl)), method = "ward.D2")
labH <- cutree(hc, k = 5); names(labH) <- colnames(expr_ctrl)
zctrl <- t(scale(expr_ctrl))
set.seed(1); km <- kmeans(zctrl, centers = 5, nstart = 25, iter.max = 100)
labK <- km$cluster; names(labK) <- rownames(zctrl)
LABS$hclust_ward_k5 <- setNames(as.character(labH), names(labH))
LABS$kmeans_k5      <- setNames(as.character(labK), names(labK))
rm(hc)

sil <- NULL
if (requireNamespace("cluster", quietly = TRUE)) {
  say("Silhouette for k-means, k = 2..8")
  dz <- dist(zctrl)
  sil <- data.frame(k = 2:8, mean_silhouette = sapply(2:8, function(k) {
    set.seed(1); kk <- kmeans(zctrl, k, nstart = 10, iter.max = 100)
    mean(cluster::silhouette(kk$cluster, dz)[, "sil_width"])
  }))
  for (i in seq_len(nrow(sil))) add_metric("kmeans", paste0("silhouette_k", sil$k[i]), sil$mean_silhouette[i], 4703)
  write.csv(sil, file.path(OUT, "silhouette.csv"), row.names = FALSE)
  rm(dz)
}

exprA <- NULL; ct_all <- NULL
if (RUN_ALLFEAT) {
  say("Variant: all features (top %d proteins by MAD of log2; per-sample scale factors applied)", ALLFEAT_N)
  keep <- !grepl("^co\\|", acc_raw)                     # drop contaminants
  ra <- as.matrix(raw[keep, unname(chan[samples])]); colnames(ra) <- samples; rownames(ra) <- acc_raw[keep]
  storage.mode(ra) <- "numeric"
  ra <- sweep(ra, 2, sc[samples], "*")
  ra <- ra[complete.cases(ra) & rowSums(ra > 0) == 10, , drop = FALSE]
  mads <- apply(log2(ra), 1, mad)
  top <- names(sort(mads, decreasing = TRUE))[seq_len(min(ALLFEAT_N, nrow(ra)))]
  say("  proteins available %d, used %d, of the 4,703 DE proteins %d are included",
      nrow(ra), length(top), sum(de_acc %in% top))
  exprA <- t(ra[top, ])
  LABS$allfeat <- run_wgcna_grid(exprA, POWER, 2, 0.20, MIN_MODULE)[[1]]
  de_in <- intersect(names(LABS$allfeat), names(paper))
  ct_all <- table(paper = paper[de_in], allfeat_module = LABS$allfeat[de_in])
  write.csv(as.data.frame.matrix(ct_all), file.path(OUT, "allfeat_vs_paper_crosstab.csv"))
  add_metric("allfeat", "n_DE_in_top_proteins", sum(de_acc %in% top), length(top))
  rm(ra)
}

loo <- NULL; n_stable <- NA
if (RUN_LOO) {
  loo <- data.frame(dropped = samples, ari_vs_full = NA_real_)
  for (i in seq_along(samples)) {
    say("Leave-one-out %d/10: dropping %s", i, samples[i])
    li <- run_wgcna_grid(expr_ctrl[-i, ], POWER, 2, 0.20, MIN_MODULE)[[1]]
    loo$ari_vs_full[i] <- ari(li, control[names(li)])
    add_metric(paste0("loo_drop_", samples[i]), "ari_vs_full", loo$ari_vs_full[i], 4703)
  }
  n_stable <- sum(loo$ari_vs_full >= 0.8, na.rm = TRUE)
  add_metric("loo", "n_runs_ari_ge_0.8", n_stable, 10)
  write.csv(loo, file.path(OUT, "loo_stability.csv"), row.names = FALSE)
}

## ---------------- metrics per variant (Q-A, Q-B) ----------------
say("Computing metrics")
VAR <- setdiff(names(LABS), c("paper", "shipped"))
zfun <- function(v) {
  if (v == "log2") return(t(scale(log2(expr_ctrl))))
  if (v == "unscaled") return(t(scale(exprU)))
  if (v == "allfeat") return(t(scale(exprA)))
  zctrl
}
SUMM <- list(); QB <- list(); PROF <- list()
for (v in VAR) {
  l  <- LABS[[v]]
  cm <- intersect(names(l), names(paper))
  cc <- intersect(names(l), names(LABS$control))
  q  <- qb_clusters(l, zfun(v), v, wpc1, de_acc)
  QB[[v]] <- q
  prof <- as.matrix(q[, c("mean_Cortex", "mean_PD", "mean_TRK")]); rownames(prof) <- q$cluster
  colnames(prof) <- c("Cortex", "PD", "TRK"); PROF[[v]] <- prof
  holds  <- any(q$rising & q$jaccard_WPC1 >= 0.8)
  bestJ  <- max(c(0, q$jaccard_WPC1[q$rising]))
  ncl    <- length(setdiff(unique(l), c("grey", NA)))
  nun    <- sum(l == "grey", na.rm = TRUE)
  a_p    <- ari(l[cm], paper[cm]); pid <- pct_identical(l[cm], paper[cm]); a_c <- ari(l[cc], LABS$control[cc])
  SUMM[[v]] <- data.frame(variant = v, n = length(l), clusters = ncl, unassigned = nun,
                          ARI_vs_paper = a_p, pct_identical = pid, ARI_vs_control = a_c,
                          QB_holds = holds, QB_best_rising_Jaccard = bestJ)
  add_metric(v, "n_clusters", ncl, length(l));          add_metric(v, "n_unassigned", nun, length(l))
  add_metric(v, "ari_vs_paper", a_p, length(cm));       add_metric(v, "pct_identical_vs_paper", pid, length(cm))
  add_metric(v, "ari_vs_control", a_c, length(cc))
  add_metric(v, "qB_claim_holds", as.integer(holds), length(l)); add_metric(v, "qB_best_rising_jaccard", bestJ, length(l))
}
summ <- do.call(rbind, SUMM); summ <- summ[order(-summ$ARI_vs_paper), ]
write.csv(summ, file.path(OUT, "summary_variants.csv"), row.names = FALSE)
write.csv(do.call(rbind, QB), file.path(OUT, "qB_clusters.csv"), row.names = FALSE)
wide <- data.frame(Protein_Accession = inp$Protein_Accession, paper = paper[inp$Protein_Accession],
                   shipped = shipped[inp$Protein_Accession])
for (v in VAR) wide[[v]] <- LABS[[v]][wide$Protein_Accession]
write.csv(wide, file.path(OUT, "labels_all_variants.csv"), row.names = FALSE)

M <- do.call(rbind, METRICS); mf <- "results/metrics.csv"
if (file.exists(mf)) {
  old <- read.csv(mf)
  old <- old[!(old$claim == "C3" & old$variant %in% M$variant), , drop = FALSE]
  M <- rbind(old[, names(M)], M)
}
write.csv(M, mf, row.names = FALSE)

## ---------------- figures ----------------
say("Figures")
nm <- c("paper", "shipped", VAR); n <- length(nm)
AM <- matrix(NA_real_, n, n, dimnames = list(nm, nm))
for (i in 1:n) for (j in 1:n) {
  a <- LABS[[nm[i]]]; b <- LABS[[nm[j]]]; cm <- intersect(names(a), names(b))
  AM[i, j] <- ari(a[cm], b[cm])
}
write.csv(round(AM, 4), file.path(OUT, "pairwise_ari.csv"))
save_fig("pairwise_ari_heatmap", function() {
  par(mar = c(9, 9, 3, 1))
  image(1:n, 1:n, AM[, n:1], col = hcl.colors(50, "viridis"), zlim = c(0.7, 1),
        axes = FALSE, xlab = "", ylab = "", main = "C3: pairwise ARI between cluster assignments")
  axis(1, 1:n, nm, las = 2, cex.axis = 0.6); axis(2, 1:n, rev(nm), las = 2, cex.axis = 0.6)
  for (i in 1:n) for (j in 1:n)
    text(i, n - j + 1, sprintf("%.2f", AM[i, j]), cex = 0.45, col = ifelse(AM[i, j] > 0.85, "black", "white"))
}, w = 10, h = 9)

gn <- grep("^grid_", VAR, value = TRUE)
best <- if (length(gn)) gn[which.max(summ$ARI_vs_paper[match(gn, summ$variant)])] else "control"
save_fig("group_profiles", function() {
  show <- intersect(c("control", "grid_ds2_mh0.25", "log2", "unscaled"), names(PROF))
  par(mfrow = c(2, 2), mar = c(3, 4, 3, 1))
  for (v in show) {
    p <- PROF[[v]]; cc <- ifelse(rownames(p) == "yellow", "gold3", rownames(p))
    matplot(1:3, t(p), type = "b", pch = 16, lty = 1, col = cc, xaxt = "n", ylim = c(-1.6, 1.6),
            xlab = "", ylab = "mean z-score", main = v)
    axis(1, 1:3, c("Cortex", "PDGFRA", "NTRK")); abline(h = 0, col = "grey")
    legend("bottomleft", legend = paste0(rownames(p), " (n=", QB[[v]]$n_members, ")"),
           col = cc, lty = 1, pch = 16, cex = 0.6, bty = "n")
  }
}, w = 10, h = 8)
if (!is.null(loo)) save_fig("loo_stability", function() {
  par(mar = c(5, 4, 3, 1))
  barplot(loo$ari_vs_full, names.arg = loo$dropped, las = 2, ylim = c(0, 1),
          ylab = "ARI vs full-data labels", main = "Leave-one-sample-out stability")
  abline(h = 0.8, lty = 2, col = "red")
}, w = 7, h = 5)
if (!is.null(sil)) save_fig("silhouette", function() {
  plot(sil$k, sil$mean_silhouette, type = "b", pch = 16, xlab = "k", ylab = "mean silhouette", main = "k-means silhouette")
}, w = 6, h = 4)

## ---------------- print summary (paste this back) ----------------
cat("\n================ SUMMARY (paste this back) ================\n")
cat(sprintf("Control check ARI vs shipped: %.4f\n", ctrl_ari))
cat("\nVariants, ranked by agreement with the paper's labels (Q-A) and cortex<PD<TRK claim (Q-B):\n")
print(summ, row.names = FALSE, digits = 4)
if (!is.null(loo)) { cat("\nLeave-one-out (Q-C):\n"); print(loo, row.names = FALSE, digits = 4)
                     cat("runs with ARI >= 0.8:", n_stable, "of 10\n") }
if (!is.null(sil)) { cat("\nSilhouette, k-means:\n"); print(sil, row.names = FALSE, digits = 4) }
if (!is.null(ct_all)) { cat("\nAll-features: paper label vs module for the DE proteins:\n"); print(ct_all) }
cat("\nOutputs written to", OUT, "and results/metrics.csv\n")