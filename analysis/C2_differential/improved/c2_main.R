# C2 Q-B..Q-F + figures + metrics. Run from repo root AFTER original/c2_qa.R.
set.seed(1)
source("analysis/C2_differential/improved/c2_lib.R")

layers <- list(protein = log2(load_protein()), phospho = log2(load_phospho()))
methods <- c("anova_param", "anova_perm", "limma_F"); fdrs <- c(0.05, 0.01)
fcs <- list(none = NA_real_, paper = NULL)   # paper cutoff is layer specific
noT3 <- setdiff(SAMPLES, "TRK3")
metrics <- list(); addm <- function(variant, metric, value, n) metrics[[length(metrics) + 1]] <<-
  data.frame(claim = "C2", variant = variant, metric = metric, value = value, n = n)

lists <- list()   # name -> character vector of DE feature IDs (per layer in name)
flags <- list()   # layer -> data.frame of flags
counts <- NULL
setcount <- function(layer, variant, de, n_feat, settings) {
  lists[[paste(layer, variant, sep = "|")]] <<- de
  counts <<- rbind(counts, data.frame(layer = layer, variant = variant, n_de = length(de), n_features = n_feat, settings = settings))
  addm(paste0(layer, "_", variant), "n_de", length(de), n_feat)
}

# ---- Q-B: method x FC x FDR grid -------------------------------------------
for (ly in names(layers)) {
  x <- layers[[ly]]; g <- GROUP[colnames(x)]
  mfc <- max_pairwise(x, g)
  fl <- data.frame(row.names = rownames(x))
  fcv <- c(none = NA_real_, "1.5" = log2(1.5), "2" = 1)
  for (me in methods) {
    p <- pvals(me, x, g)
    for (fn in names(fcv)) for (fd in fdrs) {
      v <- sprintf("%s_FC%s_FDR%s", me, fn, fd)
      de <- de_flag(p, mfc, fd, fcv[[fn]])
      fl[[v]] <- de
      setcount(ly, v, rownames(x)[de], nrow(x), sprintf("method=%s;FC=%s;FDR=%s;n=10", me, fn, fd))
    }
  }
  # TRK3-dropped limma F at the layer's paper thresholds
  x9 <- x[, noT3]; g9 <- GROUP[noT3]
  de9 <- de_flag(limma_F_p(x9, g9), max_pairwise(x9, g9), 0.05, CUT[[ly]])
  fl[["limma_F_paperFC_FDR0.05_noTRK3"]] <- de9
  setcount(ly, "limma_F_paperFC_FDR0.05_noTRK3", rownames(x)[de9], nrow(x), "method=limma_F;paper thresholds;n=9 (no TRK3)")
  flags[[ly]] <- fl
}
paperkey <- function(ly, me) sprintf("%s|%s_FC%s_FDR0.05", ly, me, if (ly == "protein") "1.5" else "2")

qb <- NULL
for (ly in names(layers)) {
  ref <- lists[[paperkey(ly, "anova_perm")]]
  for (me in c("anova_param", "limma_F")) {
    j <- jaccard(lists[[paperkey(ly, me)]], ref)
    qb <- rbind(qb, data.frame(layer = ly, comparison = paste0(me, "_vs_anova_perm_paper_thresholds"), jaccard = j, agree = j >= 0.8))
    addm(paste0("QB_", ly, "_", me, "_vs_perm"), "jaccard", j, nrow(layers[[ly]]))
  }
  j <- jaccard(lists[[paste0(ly, "|limma_F_paperFC_FDR0.05_noTRK3")]], lists[[paperkey(ly, "limma_F")]])
  qb <- rbind(qb, data.frame(layer = ly, comparison = "limma_F_noTRK3_vs_limma_F_all10", jaccard = j, agree = j >= 0.8))
  addm(paste0("QB_", ly, "_limmaF_noTRK3_vs_all10"), "jaccard", j, nrow(layers[[ly]]))
}
write.csv(qb, file.path(OUT, "c2_qb_method_jaccard.csv"), row.names = FALSE)

# ---- Q-D: genotype contrast ---------------------------------------------------
qd <- NULL; geno <- list()
for (ly in names(layers)) {
  x <- layers[[ly]]
  for (arm in c("all10", "noTRK3")) {
    s <- if (arm == "all10") SAMPLES else noT3
    r <- limma_geno(x[, s], GROUP[s])
    de <- !is.na(r$adjP) & r$adjP < 0.05 & abs(r$logFC) >= CUT[[ly]]
    geno[[paste(ly, arm, sep = "|")]] <- list(r = r, de = de)
    setcount(ly, paste0("geno_NTRKvsPDGFRA_", arm), rownames(x)[de], nrow(x),
             sprintf("limma t NTRK-PDGFRA;adjP<0.05;|log2FC|>=%.3f;%s", CUT[[ly]], arm))
    qd <- rbind(qd, data.frame(layer = ly, arm = arm, n_de = sum(de), share_higher_in_NTRK = mean(r$logFC[de] > 0),
                               n_features = nrow(x)))
  }
  j <- jaccard(lists[[paste0(ly, "|geno_NTRKvsPDGFRA_all10")]], lists[[paste0(ly, "|geno_NTRKvsPDGFRA_noTRK3")]])
  qd <- rbind(qd, data.frame(layer = ly, arm = "jaccard_all10_vs_noTRK3", n_de = NA, share_higher_in_NTRK = j, n_features = nrow(x)))
  addm(paste0("QD_", ly), "jaccard_all10_vs_noTRK3", j, nrow(x))
  addm(paste0("QD_", ly, "_all10"), "share_higher_in_NTRK", qd$share_higher_in_NTRK[qd$layer == ly & qd$arm == "all10"], sum(geno[[paste0(ly, "|all10")]]$de))
}
write.csv(qd, file.path(OUT, "c2_qd_genotype.csv"), row.names = FALSE)
for (ly in names(layers)) {
  a <- geno[[paste0(ly, "|all10")]]; b <- geno[[paste0(ly, "|noTRK3")]]
  write_gz(data.frame(Feature = rownames(a$r), logFC_all10 = round(a$r$logFC, 4), adjP_all10 = signif(a$r$adjP, 4), DE_all10 = a$de,
                      logFC_noTRK3 = round(b$r$logFC, 4), adjP_noTRK3 = signif(b$r$adjP, 4), DE_noTRK3 = b$de),
           file.path(OUT, sprintf("c2_%s_genotype_lists.csv.gz", ly)))
  write_gz(data.frame(Feature = rownames(flags[[ly]]), flags[[ly]]), file.path(OUT, sprintf("c2_%s_DE_flags.csv.gz", ly)))
}

# ---- Q-E: leave-one-sample-out ------------------------------------------------
qe <- NULL
for (ly in names(layers)) {
  x <- layers[[ly]]
  full3 <- lists[[paperkey(ly, "limma_F")]]; fullg <- lists[[paste0(ly, "|geno_NTRKvsPDGFRA_all10")]]
  for (s in SAMPLES) {
    keep <- setdiff(SAMPLES, s); xs <- x[, keep]; gs <- GROUP[keep]
    d3 <- rownames(x)[de_flag(limma_F_p(xs, gs), max_pairwise(xs, gs), 0.05, CUT[[ly]])]
    r <- limma_geno(xs, gs)
    dg <- rownames(x)[!is.na(r$adjP) & r$adjP < 0.05 & abs(r$logFC) >= CUT[[ly]]]
    qe <- rbind(qe, data.frame(layer = ly, left_out = s, list = "limma_F_3group", n_de = length(d3), jaccard = jaccard(d3, full3)),
                    data.frame(layer = ly, left_out = s, list = "NTRK_vs_PDGFRA", n_de = length(dg), jaccard = jaccard(dg, fullg)))
  }
}
write.csv(qe, file.path(OUT, "c2_qe_leave_one_out.csv"), row.names = FALSE)
qe_sum <- aggregate(jaccard ~ layer + list, qe, function(v) c(min = min(v), median = median(v), n_ge_0.8 = sum(v >= 0.8)))
qe_sum <- do.call(data.frame, qe_sum); qe_sum$stable <- qe_sum$jaccard.n_ge_0.8 >= 9
write.csv(qe_sum, file.path(OUT, "c2_qe_summary.csv"), row.names = FALSE)
for (i in seq_len(nrow(qe_sum))) {
  v <- paste0("QE_", qe_sum$layer[i], "_", qe_sum$list[i])
  addm(v, "jaccard_min", qe_sum$jaccard.min[i], 10); addm(v, "jaccard_median", qe_sum$jaccard.median[i], 10)
  addm(v, "runs_jaccard_ge_0.8", qe_sum$jaccard.n_ge_0.8[i], 10)
}

# ---- Q-F: phospho normalised to protein ----------------------------------------
w <- rd("JUMPn/JUMPnk/Integration/whl_hgg_v5_20160513.txt")
wm <- as.matrix(w[, PHOS_COLS]); storage.mode(wm) <- "double"; colnames(wm) <- SAMPLES
gsum <- rowsum(wm, w[[1]]); stopifnot(all(gsum > 0)); lgp <- log2(gsum)
ph <- layers$phospho
tok <- sub("_.*", "", rownames(ph))
ok <- tok != "NA" & tok %in% rownames(lgp)
cat("Q-F: sites", nrow(ph), "; token NA:", sum(tok == "NA"), "; gene unmatched (non-NA token):", sum(tok != "NA" & !tok %in% rownames(lgp)),
    "; matched:", sum(ok), "\n")
xn <- ph[ok, ] - lgp[tok[ok], ]; xu <- ph[ok, ]
g <- GROUP[SAMPLES]
fl <- function(x) de_flag(limma_F_p(x, g), max_pairwise(x, g), 0.05, CUT[["phospho"]])
du <- rownames(xu)[fl(xu)]; dn <- rownames(xn)[fl(xn)]
qf <- data.frame(n_sites_total = nrow(ph), n_token_NA = sum(tok == "NA"), n_unmatched_gene = sum(tok != "NA" & !tok %in% rownames(lgp)),
                 n_matched = sum(ok), n_de_unnormalised = length(du), n_de_normalised = length(dn),
                 n_both = length(intersect(du, dn)), n_only_unnormalised = length(setdiff(du, dn)),
                 n_only_normalised = length(setdiff(dn, du)), jaccard = jaccard(du, dn))
write.csv(qf, file.path(OUT, "c2_qf_normalisation.csv"), row.names = FALSE); print(t(qf))
lists[["phospho|matched_unnormalised"]] <- du; lists[["phospho|matched_normalised"]] <- dn
for (nm in setdiff(names(qf), "n_sites_total")) addm("QF_phospho_norm", nm, qf[[nm]], qf$n_matched)
write_gz(data.frame(Site = rownames(xu), DE_unnormalised = rownames(xu) %in% du, DE_normalised = rownames(xu) %in% dn),
         file.path(OUT, "c2_phospho_normalisation_flags.csv.gz"))

# ---- tidy summaries --------------------------------------------------------------
write.csv(counts, file.path(OUT, "c2_counts_per_variant.csv"), row.names = FALSE)
jm <- list()
for (ly in names(layers)) {
  nm <- grep(paste0("^", ly, "[|]"), names(lists), value = TRUE)
  nm <- nm[!grepl("matched_", nm)]
  J <- outer(nm, nm, Vectorize(function(a, b) jaccard(lists[[a]], lists[[b]])))
  dimnames(J) <- list(sub(".*[|]", "", nm), sub(".*[|]", "", nm)); jm[[ly]] <- J
  write.csv(data.frame(variant = rownames(J), round(J, 4), check.names = FALSE), file.path(OUT, sprintf("c2_%s_jaccard_matrix.csv", ly)), row.names = FALSE)
}
# benchmark comparison for proteins (all grid variants vs the paper's 4,703)
bacc <- rd("JUMPn/JUMPnc/HGG_whl_input.txt")$Protein_Accession
bp <- do.call(rbind, lapply(grep("^protein[|]", names(lists), value = TRUE), function(k) {
  d <- lists[[k]]; tp <- length(intersect(d, bacc))
  data.frame(variant = sub(".*[|]", "", k), n_de = length(d), jaccard = jaccard(d, bacc), precision = if (length(d)) tp / length(d) else NA, recall = tp / length(bacc))
}))
write.csv(bp, file.path(OUT, "c2_protein_vs_paper_list.csv"), row.names = FALSE)

# ---- figures ------------------------------------------------------------------------
pal <- c(anova_param = "#4E79A7", anova_perm = "#F28E2B", limma_F = "#59A14F", other = "#9C755F")
save_fig("c2_fig1_counts", 11, 9, function() {
  par(mfrow = c(2, 1), mar = c(11, 6, 3, 1))
  for (ly in names(layers)) {
    d <- counts[counts$layer == ly, ]
    col <- pal[sub("_FC.*|_paper.*", "", sub("^geno.*", "other", d$variant))]; col[is.na(col)] <- pal["other"]
    bp <- barplot(d$n_de, names.arg = d$variant, las = 2, cex.names = 0.55, col = col, ylab = "DE features",
                  main = sprintf("%s (%d features tested); dashed = paper's %s", ly, d$n_features[1], format(PAPER_N[[ly]], big.mark = ",")))
    abline(h = PAPER_N[[ly]], lty = 2, lwd = 2)
  }
})
save_fig("c2_fig2_jaccard_heatmap", 13, 7, function() {
  par(mfrow = c(1, 2), mar = c(12, 12, 3, 1))
  cols <- colorRampPalette(c("#FFFFFF", "#6BAED6", "#08306B"))(50)
  for (ly in names(layers)) {
    J <- jm[[ly]]; n <- nrow(J)
    image(seq_len(n), seq_len(n), J[, n:1], zlim = c(0, 1), col = cols, axes = FALSE, xlab = "", ylab = "", main = paste("Jaccard,", ly))
    axis(1, seq_len(n), rownames(J), las = 2, cex.axis = 0.45); axis(2, seq_len(n), rev(rownames(J)), las = 2, cex.axis = 0.45)
    box()
  }
})
save_fig("c2_fig3_leave_one_out", 11, 7, function() {
  par(mfrow = c(2, 2), mar = c(6, 4.5, 3, 1))
  for (ly in names(layers)) for (l in c("limma_F_3group", "NTRK_vs_PDGFRA")) {
    d <- qe[qe$layer == ly & qe$list == l, ]
    barplot(d$jaccard, names.arg = d$left_out, las = 2, ylim = c(0, 1), col = ifelse(d$jaccard >= 0.8, "#59A14F", "#E15759"),
            ylab = "Jaccard vs full list", main = paste(ly, l))
    abline(h = 0.8, lty = 2, lwd = 2)
  }
})

mm <- do.call(rbind, metrics)
old <- read.csv("results/metrics.csv", stringsAsFactors = FALSE)
append_metrics(mm, drop_variants = setdiff(old$variant[old$claim == "C2"], c("QA_unscaled", "QA_rescaled", "QC_rescaling")))
writeLines(capture.output(sessionInfo()), "env/sessionInfo_C2.txt")
cat("done\n")
