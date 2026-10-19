# ======================================================================
# C1 part B: is one NTRK sample (TRK3) behind the failures in C1?
# Run from the repo root:  Rscript analysis/C1_exploration/improved/c1b_trk3.R
# Pre-registration (provenance version) is in LOG.md. Companion to c1_exploration.R (same variant grid).
# ======================================================================
options(stringsAsFactors = FALSE)
set.seed(1)
stopifnot(file.exists("CLAUDE.md"))
OUT <- "analysis/C1_exploration/improved"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
say <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), sprintf(...))
ari <- function(a, b) {
  ok <- !is.na(a) & !is.na(b); a <- a[ok]; b <- b[ok]
  tab <- table(a, b); n <- sum(tab); c2 <- function(x) x * (x - 1) / 2
  s <- sum(c2(tab)); sa <- sum(c2(rowSums(tab))); sb <- sum(c2(colSums(tab))); e <- sa * sb / c2(n)
  out <- (s - e) / ((sa + sb) / 2 - e); if (is.nan(out)) NA_real_ else out
}
save_fig <- function(name, fun, w, h) {
  png(file.path(OUT, paste0(name, ".png")), width = w, height = h, units = "in", res = 150); fun(); invisible(dev.off())
  pdf(file.path(OUT, paste0(name, ".pdf")), width = w, height = h); fun(); invisible(dev.off())
}
METRICS <- list()
add_metric <- function(variant, metric, value, n)
  METRICS[[length(METRICS) + 1]] <<- data.frame(claim = "C1", variant = variant, metric = metric, value = as.numeric(value), n = n)

## ---------------- data (same as c1_exploration.R) ----------------
samples <- c("Cortex1","Cortex2","Cortex3","PD1","PD2","PD3","PD4","TRK1","TRK2","TRK3")
group   <- setNames(sub("[0-9]+$", "", samples), samples)
tum     <- samples[group[samples] != "Cortex"]
truth2  <- setNames(ifelse(group == "Cortex", "Cortex", "Tumour"), samples)
map <- c(PDGFRA_1a = "PD1", PDGFRA_1b = "PD2", PDGFRA_2 = "PD3", PDGFRA_3 = "PD4",
         NTRK_1 = "TRK1", NTRK_2 = "TRK2", NTRK_3 = "TRK3", CNTRL_1 = "Cortex1", CNTRL_2 = "Cortex2", CNTRL_3 = "Cortex3")
read_layer <- function(path, drop_contam) {
  d <- read.table(path, header = TRUE, sep = "\t", quote = "", comment.char = "", check.names = FALSE)
  if (drop_contam) d <- d[!grepl("CON_", d[[1]]), ]
  m <- as.matrix(d[, -1]); storage.mode(m) <- "numeric"
  rownames(m) <- make.unique(as.character(d[[1]]))
  stopifnot(!anyNA(map[colnames(m)])); colnames(m) <- map[colnames(m)]
  m[, samples]
}
L <- list(proteome = read_layer("JUMPn/JUMPnk/Integration/whl_hgg_v5_20160513.txt", FALSE),
          phospho  = read_layer("JUMPn/JUMPnk/Integration/pho_hgg_uniq_v5_20160513.txt", TRUE))
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

## ---------------- variant specs (identical grid to c1_exploration.R) ----------------
sp_pca <- expand.grid(layer = names(L), transform = c("log2", "raw"), features = c("all", "top1000"),
                      scaling = c("centred", "scaled"), stringsAsFactors = FALSE)
sp_pca$type <- "pca"; sp_pca$dist <- NA; sp_pca$linkage <- NA
sp_hc <- expand.grid(layer = names(L), transform = c("log2", "raw"), features = c("all", "top1000"),
                     dist = c("euclid_z", "cor"), linkage = c("ward.D", "ward.D2", "average", "complete"), stringsAsFactors = FALSE)
sp_hc$type <- "hc"; sp_hc$scaling <- NA
cn <- c("type", "layer", "transform", "features", "scaling", "dist", "linkage")
V <- rbind(sp_pca[, cn], sp_hc[, cn])
V$name <- ifelse(V$type == "pca", paste("pca", V$layer, V$transform, V$features, V$scaling, sep = "_"),
                 paste("hc", V$layer, V$transform, V$features, V$dist, V$linkage, sep = "_"))
run_spec <- function(s, cols, k) {
  X <- prep(s$layer, s$transform, s$features, cols)
  if (s$type == "pca") {
    p <- prcomp(t(X), center = TRUE, scale. = (s$scaling == "scaled"))
    set.seed(1); km <- kmeans(p$x[, 1:2, drop = FALSE], centers = k, nstart = 25)
    cl <- km$cluster; names(cl) <- rownames(p$x)
  } else { cl <- run_hc(X, s$dist, s$linkage, k); names(cl) <- colnames(X) }
  cl
}
mis <- function(cl, truth) {
  tab <- table(cl, truth[names(cl)])
  maj <- colnames(tab)[apply(tab, 1, which.max)]; names(maj) <- rownames(tab)
  names(cl)[maj[as.character(cl)] != truth[names(cl)]]
}
lab <- function(x) if (length(x)) paste(x, collapse = "+") else "none"

## ---------------- Q-F1: which sample is misgrouped? ----------------
say("F1: misgrouped samples")
F1 <- do.call(rbind, lapply(seq_len(nrow(V)), function(i) {
  s <- V[i, ]; c3 <- run_spec(s, samples, 3); c2 <- run_spec(s, tum, 2)
  data.frame(variant = s$name, layer = s$layer, transform = s$transform,
             ARI_all3 = ari(c3, group[names(c3)]), mis_all3 = lab(mis(c3, group)),
             ARI_tumour2 = ari(c2, group[names(c2)]), mis_tumour2 = lab(mis(c2, group)))
}))
write.csv(F1, file.path(OUT, "c1b_misgrouped.csv"), row.names = FALSE)

## ---------------- Q-F2: tumour vs cortex ----------------
say("F2: tumour vs cortex")
F2 <- do.call(rbind, lapply(seq_len(nrow(V)), function(i) {
  s <- V[i, ]; c2 <- run_spec(s, samples, 2)
  data.frame(variant = s$name, layer = s$layer, transform = s$transform, ARI_tumour_vs_cortex = ari(c2, truth2[names(c2)]))
}))
write.csv(F2, file.path(OUT, "c1b_tumour_vs_cortex.csv"), row.names = FALSE)
F2sum <- aggregate(ARI_tumour_vs_cortex ~ layer + transform, F2, function(x) mean(x > 0.999))

## ---------------- Q-F3: TRK3 profile ----------------
say("F3: TRK3 profile")
inp <- read.table("JUMPn/JUMPnc/HGG_whl_input.txt", header = TRUE, sep = "\t", quote = "", comment.char = "", check.names = FALSE)
cu <- inp[grepl("^cu\\|", inp$Protein_Accession), c("Protein_Accession", samples)]
cu$TRK3_over_mean_TRK1_TRK2 <- cu$TRK3 / rowMeans(cu[, c("TRK1", "TRK2")])
write.csv(cu, file.path(OUT, "c1b_markers.csv"), row.names = FALSE)
POS <- do.call(rbind, lapply(names(L), function(ly) {
  X <- prep(ly, "log2", "all", tum); p <- prcomp(t(X), center = TRUE, scale. = FALSE); s1 <- p$x[, 1]
  pd <- mean(s1[c("PD1", "PD2", "PD3", "PD4")]); idx <- (s1["TRK3"] - pd) / (mean(s1[c("TRK1", "TRK2")]) - pd)
  data.frame(layer = ly, sample = names(s1), PC1 = as.numeric(s1), TRK3_position_index = as.numeric(idx))
}))
write.csv(POS, file.path(OUT, "c1b_tumour_positions.csv"), row.names = FALSE)

## ---------------- Q-F4: drop one sample at a time ----------------
say("F4: drop-one specificity (this is the slow part)")
Vlog <- V[V$transform == "log2", ]; F4 <- NULL
for (d in c("none", samples)) {
  keep <- setdiff(samples, d)
  for (ly in names(L)) {
    sub <- Vlog[Vlog$layer == ly, ]
    a <- sapply(seq_len(nrow(sub)), function(i) { cl <- run_spec(sub[i, ], keep, 3); ari(cl, group[names(cl)]) })
    F4 <- rbind(F4, data.frame(dropped = d, layer = ly, n_variants = length(a), fraction_ARI1 = mean(a > 0.999), mean_ARI = mean(a)))
  }
}
write.csv(F4, file.path(OUT, "c1b_drop_one.csv"), row.names = FALSE)

## ---------------- Q-F5: linkage with and without TRK3 ----------------
say("F5: linkage with and without TRK3")
F5 <- NULL
for (ly in names(L)) for (lk in c("ward.D", "ward.D2", "average", "complete")) for (cond in c("all 10", "without TRK3")) {
  cols <- if (cond == "all 10") samples else setdiff(samples, "TRK3")
  X <- prep(ly, "log2", "top1000", cols); cl <- run_hc(X, "euclid_z", lk, 3)
  F5 <- rbind(F5, data.frame(layer = ly, linkage = lk, samples = cond, ARI = ari(cl, group[colnames(X)])))
}
write.csv(F5, file.path(OUT, "c1b_linkage_with_without_trk3.csv"), row.names = FALSE)

## ---------------- metrics ----------------
for (ly in names(L)) {
  f <- F1[F1$layer == ly & F1$transform == "log2" & F1$ARI_all3 < 0.999, ]
  add_metric(paste0("partB_", ly), "share_failing_log2_variants_only_TRK3_misgrouped", if (nrow(f)) mean(f$mis_all3 == "TRK3") else NA, nrow(f))
  add_metric(paste0("partB_", ly), "TRK3_position_index", POS$TRK3_position_index[POS$layer == ly][1], 7)
}
for (i in seq_len(nrow(F4))) add_metric(paste0("partB_drop_", F4$dropped[i], "_", F4$layer[i]), "fraction_log2_variants_ARI1", F4$fraction_ARI1[i], F4$n_variants[i])
for (i in seq_len(nrow(F2sum))) add_metric(paste0("partB_tumour_vs_cortex_", F2sum$layer[i], "_", F2sum$transform[i]), "fraction_variants_ARI1", F2sum$ARI_tumour_vs_cortex[i], 20)
M <- do.call(rbind, METRICS); mf <- "results/metrics.csv"
if (file.exists(mf)) { old <- read.csv(mf); old <- old[!(old$claim == "C1" & old$variant %in% M$variant), , drop = FALSE]; M <- rbind(old[, names(M)], M) }
write.csv(M, mf, row.names = FALSE)

## ---------------- figures ----------------
save_fig("c1b_drop_one", function() {
  par(mfrow = c(1, 2), mar = c(6, 4.5, 3, 1))
  for (ly in names(L)) {
    d <- F4[F4$layer == ly, ]
    barplot(d$fraction_ARI1, names.arg = d$dropped, las = 2, ylim = c(0, 1),
            col = ifelse(d$dropped == "TRK3", "red3", ifelse(d$dropped == "none", "grey30", "grey75")),
            ylab = "fraction of log2 variants with ARI = 1", main = paste(ly, ": drop one sample"))
  }
}, w = 11, h = 5)
save_fig("c1b_tumour_positions", function() {
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  for (ly in names(L)) {
    d <- POS[POS$layer == ly, ]
    plot(d$PC1, rep(1, nrow(d)), pch = 16, cex = 2, yaxt = "n", ylab = "", main = ly,
         col = ifelse(grepl("^PD", d$sample), "deepskyblue3", "red3"), xlab = "PC1 (tumours only, log2, all features)")
    text(d$PC1, 1, d$sample, pos = 3, cex = 0.8)
  }
}, w = 11, h = 3.5)

## ---------------- summary ----------------
cat("\n================ SUMMARY (paste this back) ================\n")
ff <- F1[F1$ARI_all3 < 0.999, ]
tt <- as.data.frame(table(layer = ff$layer, transform = ff$transform, misgrouped = ff$mis_all3))
cat("\nQ-F1: misgrouped sample(s) in failing variants (all samples, k = 3):\n"); print(tt[tt$Freq > 0, ], row.names = FALSE)
for (ly in names(L)) {
  f <- ff[ff$layer == ly & ff$transform == "log2", ]
  cat(sprintf("  %s, failing log2 variants: %d; only TRK3 misgrouped in %d (%.0f%%)\n", ly, nrow(f), sum(f$mis_all3 == "TRK3"), if (nrow(f)) 100 * mean(f$mis_all3 == "TRK3") else NA))
}
fo <- F1[F1$ARI_tumour2 < 0.999, ]; to <- as.data.frame(table(layer = fo$layer, misgrouped = fo$mis_tumour2))
cat("\nTumour-only (k = 2) misgrouped sample(s):\n"); print(to[to$Freq > 0, ], row.names = FALSE)
cat("\nQ-F2: fraction of variants separating tumour from cortex (ARI = 1):\n"); print(F2sum, row.names = FALSE)
cat("\nQ-F3: transgene/marker rows, and TRK3 relative to TRK1/TRK2:\n")
cup <- cu; cup[, -1] <- round(cup[, -1], 2); print(cup, row.names = FALSE)
cat("\nTRK3 position index (0 = PDGFRA mean, 1 = TRK1/2 mean):\n"); print(unique(POS[, c("layer", "TRK3_position_index")]), row.names = FALSE)
cat("\nQ-F4: fraction of log2 variants with ARI = 1 after dropping each sample ('none' = full data):\n"); print(F4, row.names = FALSE, digits = 3)
cat("\nQ-F5: paper-setup linkage, with and without TRK3:\n"); print(F5, row.names = FALSE, digits = 3)
cat("\nOutputs written to", OUT, "and results/metrics.csv\n")