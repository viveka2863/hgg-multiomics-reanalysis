# ======================================================================
# C1 part C (group-review round 1, Q-Q): sample-level QC of the proteome and phospho layers in data versions A, B, C,
# and marker-panel scores in the proteome and the RNA layer.
# Design and thresholds: LOG.md "PRE-REGISTRATION: group-review amendments, round 1". Run from the repo root:
#   Rscript analysis/C1_exploration/improved/c1c_qc.R
# Needs data/external/csv/{sd1,sd2a,sd3}.csv (c2b_convert_xlsx.R, c11_convert_sd3.R).
# Versions (docs/DATA_VERSIONS.md):  A = repository files hgg_whole_v1 / pho_hgg_uniq_v5 (unscaled);
#   B = published Supp Data 1 / 2a intensities (rescaled by the authors); C = version A log2 with per-sample median centring.
# All outputs are descriptive; no sample is excluded.
# ======================================================================
options(stringsAsFactors = FALSE)
set.seed(1)
source("analysis/C11_rna_layer/improved/c11_lib.R")   # rdx, num, read_sd3, protein_map, load_protein/load_phospho (c2_lib), SAMPLES, GROUP
O <- "analysis/C1_exploration/improved"
msg <- function(...) cat(sprintf(...), "\n")
wr <- function(df, f) write.csv(df, file.path(O, f), row.names = FALSE)

GRP <- c(PD1 = "PDGFRA", PD2 = "PDGFRA", PD3 = "PDGFRA", PD4 = "PDGFRA", TRK1 = "NTRK", TRK2 = "NTRK", TRK3 = "NTRK",
         Cortex1 = "CNTRL", Cortex2 = "CNTRL", Cortex3 = "CNTRL")[SAMPLES]
COL <- c(CNTRL = "#999999", PDGFRA = "#0072B2", NTRK = "#D55E00")

# ---- version C: per-sample median centring (each sample's median is set to the grand median = median of the 10 sample medians) ----
# Input: a log2 matrix (features x samples). A per-sample constant shift only; correlations over features are unchanged by it.
median_centre <- function(L) {
  m <- apply(L, 2, median); g <- median(m)
  sweep(L, 2, m - g)
}

# ---- load the two layers in versions A and B -----------------------------------------------------------------------
SD_ORDER <- c("Cortex1","Cortex2","Cortex3","PD1","PD2","PD3","PD4","TRK1","TRK2","TRK3")   # column order of SD1 / SD2a (verified by values in C2 part B)
A_prot <- log2(load_protein()); A_pho <- log2(load_phospho())
x <- rdx("sd1.csv"); stopifnot(x[1, 1] == "Coexpresion Clusters", x[2, 16] == "NTRK/PDGFRA"); x <- x[-(1:2), ]; x <- x[!is.na(x[, 3]), ]
B_prot <- num(x[, 6:15]); colnames(B_prot) <- SD_ORDER; rownames(B_prot) <- x[, 3]; B_prot <- B_prot[, SAMPLES]
B_prot_gene <- setNames(x[, 2], x[, 3])
y <- rdx("sd2a.csv"); stopifnot(y[1, 1] == "Coexpresion Clusters", y[2, 17] == "NTRK/PDGFRA"); y <- y[-(1:2), ]; y <- y[!is.na(y[, 3]), ]
B_pho <- num(y[, 7:16]); colnames(B_pho) <- SD_ORDER
rownames(B_pho) <- paste(y[, 2], sapply(strsplit(y[, 3], "[|]"), `[`, 2), y[, 4], sep = "_"); B_pho <- B_pho[, SAMPLES]
okrow <- function(M) apply(M, 1, function(v) all(is.finite(v) & v > 0))
msg("rows with non-positive or missing values in B: protein %d of %d, phospho %d of %d",
    sum(!okrow(B_prot)), nrow(B_prot), sum(!okrow(B_pho)), nrow(B_pho))
B_prot <- log2(B_prot[okrow(B_prot), ]); B_pho <- log2(B_pho[okrow(B_pho), ])
L <- list(proteome = list(A = A_prot, B = B_prot, C = median_centre(A_prot)),
          phospho  = list(A = A_pho,  B = B_pho,  C = median_centre(A_pho)))
stopifnot(!anyNA(unlist(lapply(L, function(l) lapply(l, function(m) sum(!is.finite(m)))))) , all(sapply(L, function(l) sapply(l, function(m) sum(!is.finite(m)))) == 0))

# ---- Q1: distribution per sample --------------------------------------------------------------------------------------
Q1 <- NULL
for (ly in names(L)) for (v in names(L[[ly]])) {
  M <- L[[ly]][[v]]; med <- apply(M, 2, median); iqr <- apply(M, 2, IQR); sdv <- apply(M, 2, sd)
  sc <- 1.4826 * median(abs(med - median(med)))                  # 1.4826 x MAD of the 10 sample medians
  z <- if (sc > 1e-9) (med - median(med)) / sc else rep(NA_real_, length(med))   # undefined when all medians are equal (version C by construction)
  Q1 <- rbind(Q1, data.frame(layer = ly, version = v, sample = colnames(M), group = unname(GRP[colnames(M)]), n_features = nrow(M),
                             median_log2 = med, IQR_log2 = iqr, SD_log2 = sdv, z_median = z, distribution_outlier = abs(z) > 3,
                             note = ifelse(sc > 1e-9, "", "MAD of sample medians = 0: z undefined"), row.names = NULL))
}
wr(Q1, "c1c_q1_distribution.csv")

# ---- Q2: gene-centred sample-sample correlations; replicate correlation = mean correlation with the other members of the own group -----
Q2 <- NULL; QR <- NULL; CM <- list()
for (ly in names(L)) for (v in names(L[[ly]])) {
  M <- L[[ly]][[v]]; Mc <- M - rowMeans(M); R <- cor(Mc)           # feature (gene / site) centred across the 10 samples
  CM[[paste(ly, v)]] <- R
  ij <- which(upper.tri(R), arr.ind = TRUE)
  Q2 <- rbind(Q2, data.frame(layer = ly, version = v, sample_i = rownames(R)[ij[, 1]], sample_j = colnames(R)[ij[, 2]], pearson = R[ij]))
  rc <- sapply(colnames(R), function(s) { o <- setdiff(names(GRP)[GRP == GRP[s]], s); mean(R[s, o]) })
  QR <- rbind(QR, data.frame(layer = ly, version = v, sample = names(rc), group = unname(GRP[names(rc)]), replicate_corr = unname(rc),
                             rank_lowest = rank(rc, ties.method = "min"), replicate_outlier = rank(rc, ties.method = "min") == 1, row.names = NULL))
}
wr(Q2, "c1c_q2_sample_correlation.csv"); wr(QR, "c1c_q2_replicate_correlation.csv")

# ---- Q3: marker panels, mean z of log2 intensity per sample --------------------------------------------------------------
PANEL <- list(neuron = c("Rbfox3", "Snap25", "Syt1", "Map2", "Gad1", "Slc17a7"),
              astrocyte = c("Gfap", "Aqp4", "Slc1a3", "Aldh1l1", "S100b"),
              oligodendrocyte = c("Mbp", "Plp1", "Mog", "Mobp"),
              microglia_immune = c("Aif1", "Cx3cr1", "Tmem119", "Ptprc"),
              proliferation = c("Mki67", "Top2a", "Pcna"))
# one row per gene (case-insensitive): the row with the highest mean linear intensity (same rule as C11)
gene_rows <- function(M, genes) {                                   # M: log2 matrix, genes: symbol per row (same order)
  k <- toupper(genes); ok <- !is.na(k) & k != "NA"; idx <- which(ok)
  idx <- idx[order(-rowMeans(2^M[idx, , drop = FALSE]))]; idx <- idx[!duplicated(k[idx])]
  list(idx = idx, key = k[idx])
}
pm <- protein_map()
dsets <- list(
  proteome_A = list(M = L$proteome$A, gene = unname(pm[rownames(L$proteome$A)])),
  proteome_B = list(M = L$proteome$B, gene = unname(B_prot_gene[rownames(L$proteome$B)])),
  proteome_C = list(M = L$proteome$C, gene = unname(pm[rownames(L$proteome$C)])))
sd3 <- read_sd3(); Lr <- log2(sd3$F)
dsets$RNA <- list(M = Lr, gene = sd3$gene)
GRP_RNA <- c(Cortex = "CNTRL", PDGFRA = "PDGFRA", NTRK = "NTRK")[sd3$meta$group]; names(GRP_RNA) <- sd3$meta$sample
Q3 <- NULL; QG <- NULL; Q3found <- NULL
for (d in names(dsets)) {
  M <- dsets[[d]]$M; gr <- gene_rows(M, dsets[[d]]$gene); Z <- t(scale(t(M[gr$idx, , drop = FALSE]))); rownames(Z) <- gr$key
  grp <- if (d == "RNA") GRP_RNA[colnames(M)] else GRP[colnames(M)]
  for (p in names(PANEL)) {
    g <- toupper(PANEL[[p]]); f <- g[g %in% rownames(Z)]
    Q3found <- rbind(Q3found, data.frame(dataset = d, panel = p, gene = PANEL[[p]], found = g %in% rownames(Z),
                                         feature = ifelse(g %in% gr$key, rownames(M)[gr$idx[match(g, gr$key)]], NA)))
    sc <- if (length(f)) colMeans(Z[f, , drop = FALSE]) else rep(NA_real_, ncol(Z))
    Q3 <- rbind(Q3, data.frame(dataset = d, panel = p, sample = colnames(M), group = unname(grp), n_genes = length(f), score = unname(sc), row.names = NULL))
  }
}
QG <- do.call(rbind, lapply(split(Q3, list(Q3$dataset, Q3$panel), drop = TRUE), function(s) {
  gm <- tapply(s$score, s$group, mean)
  data.frame(dataset = s$dataset[1], panel = s$panel[1], group = names(gm), n_samples = as.integer(table(s$group)[names(gm)]),
             mean_score = as.numeric(gm), diff_vs_cortex = as.numeric(gm - gm["CNTRL"]), row.names = NULL)
}))
wr(Q3, "c1c_q3_marker_scores.csv"); wr(QG, "c1c_q3_marker_group_scores.csv"); wr(Q3found, "c1c_q3_marker_genes_found.csv")

# cortex-likeness of TRK3 in neuronal markers (proteome; the RNA samples cannot be linked to TRK3, see sample_key_extended.csv)
CL <- do.call(rbind, lapply(c("proteome_A", "proteome_B", "proteome_C"), function(d) {
  s <- Q3[Q3$dataset == d & Q3$panel == "neuron", ]; v <- setNames(s$score, s$sample)
  pd <- mean(v[c("PD1", "PD2", "PD3", "PD4")]); cx <- mean(v[c("Cortex1", "Cortex2", "Cortex3")])
  data.frame(dataset = d, n_neuron_genes = s$n_genes[1], TRK1 = v[["TRK1"]], TRK2 = v[["TRK2"]], TRK3 = v[["TRK3"]], PDGFRA_mean = pd, Cortex_mean = cx,
             TRK3_above_TRK1_and_TRK2 = v[["TRK3"]] > max(v[["TRK1"]], v[["TRK2"]]), TRK3_above_PDGFRA_mean = v[["TRK3"]] > pd,
             TRK3_cortex_like_in_neuronal_markers = v[["TRK3"]] > max(v[["TRK1"]], v[["TRK2"]]) && v[["TRK3"]] > pd)
}))
wr(CL, "c1c_q3_trk3_cortex_likeness.csv")

# ---- figures (PNG) -------------------------------------------------------------------------------------------------------------
png(file.path(O, "c1c_fig1_boxplots.png"), 15, 8, units = "in", res = 150)
par(mfrow = c(2, 3), mar = c(5, 4.5, 3, 1))
for (ly in names(L)) for (v in names(L[[ly]])) {
  M <- L[[ly]][[v]]; is3 <- colnames(M) == "TRK3"
  boxplot(M, col = COL[GRP[colnames(M)]], border = ifelse(is3, "black", "grey25"), lwd = ifelse(is3, 3, 1), outline = FALSE, las = 2,
          ylab = "log2 intensity", main = sprintf("%s, version %s (%s features)", ly, v, format(nrow(M), big.mark = ",")), xaxt = "n")
  axis(1, at = which(!is3), labels = colnames(M)[!is3], las = 2, col.axis = "grey30")
  axis(1, at = which(is3), labels = "TRK3", las = 2, font = 2)
  abline(h = median(apply(M, 2, median)), lty = 3, col = "grey40")
  if (v == "A") legend("topright", c(names(COL), "TRK3 (heavy outline)"), fill = c(COL, "white"), border = c(rep("grey25", 3), "black"), bty = "n", cex = 0.8)
}
invisible(dev.off())

png(file.path(O, "c1c_fig2_correlation_heatmaps.png"), 15, 9.5, units = "in", res = 150)
par(mfrow = c(2, 3), mar = c(6, 6, 3, 1))
for (ly in names(L)) for (v in names(L[[ly]])) {
  R <- CM[[paste(ly, v)]]; rc <- QR[QR$layer == ly & QR$version == v, ]; rc <- setNames(rc$replicate_corr, rc$sample)[colnames(R)]
  rng <- range(R[upper.tri(R)]); n <- ncol(R)
  image(seq_len(n), seq_len(n), t(R[n:1, ]), col = hcl.colors(30, "Blues 3", rev = TRUE), zlim = rng, axes = FALSE, xlab = "", ylab = "",
        main = sprintf("%s, version %s: gene-centred correlation", ly, v))
  lab <- sprintf("%s (%.2f)", colnames(R), rc)
  axis(1, seq_len(n), lab, las = 2, cex.axis = 0.8, font = 1); axis(2, seq_len(n), rev(lab), las = 1, cex.axis = 0.8)
  for (i in seq_len(n)) for (j in seq_len(n)) text(j, n + 1 - i, sprintf("%.2f", R[i, j]), cex = 0.55, col = ifelse((R[i, j] - rng[1]) / diff(rng) > 0.6, "white", "black"))
  mtext("label = sample (mean correlation with own group)", 1, 5, cex = 0.7)
}
invisible(dev.off())

png(file.path(O, "c1c_fig3_marker_scores.png"), 17, 11, units = "in", res = 150)
par(mfrow = c(4, 5), mar = c(5, 4, 3, 0.5))
for (d in names(dsets)) for (p in names(PANEL)) {
  s <- Q3[Q3$dataset == d & Q3$panel == p, ]; s$is3 <- s$sample == "TRK3"
  bp <- barplot(s$score, col = COL[s$group], border = ifelse(s$is3, "black", NA), names.arg = rep("", nrow(s)), ylim = range(c(0, s$score), na.rm = TRUE) * 1.1,
                main = sprintf("%s: %s (%d genes)", d, p, s$n_genes[1]), ylab = "mean z")
  axis(1, bp, labels = s$sample, las = 2, cex.axis = 0.65, font = 1)
  if (any(s$is3)) points(bp[s$is3], s$score[s$is3] + sign(s$score[s$is3]) * 0.1 * diff(range(c(0, s$score))), pch = 8, cex = 1.4)
  if (d == "proteome_A" && p == "neuron") legend("topleft", c(names(COL), "* TRK3"), fill = c(COL, NA), border = NA, bty = "n", cex = 0.75)
}
invisible(dev.off())

# ---- metrics (append-only; rows with variant prefix amendR1_ are replaced on re-run) -------------------------------------------
mt <- list(); am <- function(v, m, x, n) if (is.finite(x)) mt[[length(mt) + 1]] <<- data.frame(claim = "C1", variant = paste0("amendR1_", v), metric = m, value = x, n = n)
for (ly in names(L)) for (v in names(L[[ly]])) {
  q <- Q1[Q1$layer == ly & Q1$version == v, ]; r <- QR[QR$layer == ly & QR$version == v, ]; tag <- paste0(ly, "_", v)
  am(paste0("Q1_", tag), "max_abs_z_sample_median", if (all(is.na(q$z_median))) NA else max(abs(q$z_median)), 10)
  am(paste0("Q1_", tag), "n_distribution_outliers", if (all(is.na(q$z_median))) NA else sum(q$distribution_outlier), 10)
  am(paste0("Q2_", tag), "min_replicate_corr", min(r$replicate_corr), 10)
  am(paste0("Q2_", tag), "TRK3_replicate_corr", r$replicate_corr[r$sample == "TRK3"], 10)
  am(paste0("Q2_", tag), "TRK3_is_lowest_replicate_corr", as.numeric(r$replicate_outlier[r$sample == "TRK3"]), 10)
}
for (i in seq_len(nrow(CL))) {
  tag <- paste0("Q3_", CL$dataset[i])
  am(tag, "neuron_score_TRK3", CL$TRK3[i], 10); am(tag, "neuron_score_PDGFRA_mean", CL$PDGFRA_mean[i], 10)
  am(tag, "TRK3_cortex_like_neuronal", as.numeric(CL$TRK3_cortex_like_in_neuronal_markers[i]), 10)
}
for (d in names(dsets)) { g <- QG[QG$dataset == d & QG$panel == "neuron" & QG$group == "NTRK", ]
  am(paste0("Q3_", d), "neuron_NTRK_minus_cortex", g$diff_vs_cortex, g$n_samples) }
mt <- do.call(rbind, mt)
mf <- "results/metrics.csv"; ln <- readLines(mf); keep <- !grepl("\"amendR1_Q[1-3]_", ln)   # drop earlier rows of this script (C1 Q-Q only)
if (!all(keep)) writeLines(ln[keep], mf)
stopifnot(substr(readChar(mf, file.info(mf)$size), file.info(mf)$size, file.info(mf)$size) == "\n")
write.table(mt, mf, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE, quote = TRUE)

# ---- printed answers --------------------------------------------------------------------------------------------------------------
cat("\n== Q1 distribution outliers (|z| > 3 on the 10 sample medians) ==\n")
print(aggregate(distribution_outlier ~ layer + version, Q1, sum)); print(Q1[abs(Q1$z_median) > 3 & !is.na(Q1$z_median), c("layer", "version", "sample", "z_median")], row.names = FALSE)
cat("max |z| per layer/version:\n"); print(aggregate(z_median ~ layer + version, Q1, function(z) max(abs(z))), digits = 3)
cat("TRK3 z:\n"); print(Q1[Q1$sample == "TRK3", c("layer", "version", "median_log2", "z_median")], digits = 3, row.names = FALSE)
cat("\n== Q2 replicate outliers (lowest replicate correlation of all 10) ==\n"); print(QR[QR$replicate_outlier, c("layer", "version", "sample", "group", "replicate_corr")], digits = 3, row.names = FALSE)
cat("TRK3 replicate correlation and rank (1 = lowest):\n"); print(QR[QR$sample == "TRK3", c("layer", "version", "replicate_corr", "rank_lowest")], digits = 3, row.names = FALSE)
cat("\n== Q3 TRK3 cortex-likeness in neuronal markers ==\n"); print(CL, digits = 3, row.names = FALSE)
cat("neuron-panel genes found:\n"); print(table(Q3found$dataset[Q3found$panel == "neuron"], Q3found$found[Q3found$panel == "neuron"]))
cat("all panels, genes found / listed:\n"); print(aggregate(found ~ dataset + panel, Q3found, function(x) sprintf("%d/%d", sum(x), length(x))))
writeLines(capture.output(sessionInfo()), "env/sessionInfo_R1.txt")
