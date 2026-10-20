# Shared helpers for C2 (pre-registered in LOG.md, 2026-10-18). Run from the repo root.
suppressPackageStartupMessages(library(limma))

OUT <- "analysis/C2_differential"
SAMPLES <- c("PD1","PD2","PD3","PD4","TRK1","TRK2","TRK3","Cortex1","Cortex2","Cortex3")
GROUP <- factor(c(rep("PDGFRA",4), rep("NTRK",3), rep("Cortex",3)), levels = c("Cortex","PDGFRA","NTRK"))
names(GROUP) <- SAMPLES
CHANNEL <- c(PD1="sig126", PD2="sig127N", PD3="sig127C", PD4="sig128N", TRK1="sig128C",
             TRK2="sig129N", TRK3="sig129C", Cortex1="sig130N", Cortex2="sig130C", Cortex3="sig131")
PHOS_COLS <- c(PD1="PDGFRA_1a", PD2="PDGFRA_1b", PD3="PDGFRA_2", PD4="PDGFRA_3",
               TRK1="NTRK_1", TRK2="NTRK_2", TRK3="NTRK_3",
               Cortex1="CNTRL_1", Cortex2="CNTRL_2", Cortex3="CNTRL_3")
CUT <- c(protein = log2(1.5), phospho = log2(2))
PAPER_N <- c(protein = 4703, phospho = 6768)

rd <- function(f) read.delim(f, check.names = FALSE, quote = "", comment.char = "", stringsAsFactors = FALSE,
                             na.strings = character(0))

load_protein <- function() {
  p <- rd("JUMPn/JUMPnk/PhosphoWholeNormalization/hgg_whole_v1_20160420.txt")
  p <- p[!grepl("^co[|]", p[["Protein Accession #"]]), ]
  m <- as.matrix(p[, CHANNEL]); storage.mode(m) <- "double"
  colnames(m) <- SAMPLES; rownames(m) <- p[["Protein Accession #"]]
  stopifnot(!anyNA(m), all(m > 0), !anyDuplicated(rownames(m)))
  m
}
load_phospho <- function() {
  p <- rd("JUMPn/JUMPnk/Integration/pho_hgg_uniq_v5_20160513.txt")
  p <- p[!grepl("CON_", p[[1]]), ]
  m <- as.matrix(p[, PHOS_COLS]); storage.mode(m) <- "double"
  colnames(m) <- SAMPLES; rownames(m) <- make.unique(p[[1]])
  stopifnot(!anyNA(m), all(m > 0))
  m
}

# ---- statistics ------------------------------------------------------------
# one-way ANOVA F for every row; g = factor over columns
anova_F <- function(x, g) {
  g <- droplevels(g); k <- nlevels(g); n <- ncol(x)
  ind <- model.matrix(~ 0 + g); cnt <- colSums(ind)
  gm <- (x %*% ind) %*% diag(1 / cnt)               # group means
  grand <- rowMeans(x)
  ssb <- as.vector((gm - grand)^2 %*% cnt)
  ssw <- rowSums((x - gm %*% t(ind))^2)
  (ssb / (k - 1)) / (ssw / (n - k))
}
anova_param_p <- function(x, g) {
  F <- anova_F(x, g); k <- nlevels(droplevels(g)); n <- ncol(x)
  pf(F, k - 1, n - k, lower.tail = FALSE)
}
# pooled-null permutation p: (1 + #null >= F) / (1 + #null); NaN F treated as NA
anova_perm_p <- function(x, g, nperm = 1000) {
  F <- anova_F(x, g)
  null <- unlist(lapply(seq_len(nperm), function(i) anova_F(x, sample(g))), use.names = FALSE)
  null <- sort(null[is.finite(null) | is.infinite(null) & null > 0])
  N <- length(null)
  cnt_lt <- findInterval(F, null, left.open = TRUE)  # number of null < F
  (1 + N - cnt_lt) / (1 + N)
}
limma_fit <- function(x, g) {
  g <- droplevels(g)
  d <- model.matrix(~ 0 + g); colnames(d) <- levels(g)
  lmFit(x, d)
}
limma_F_p <- function(x, g) {
  fit <- limma_fit(x, g); lv <- levels(droplevels(g))
  cm <- makeContrasts(contrasts = paste0(lv[-1], "-", lv[1]), levels = fit$design)
  eBayes(contrasts.fit(fit, cm))$F.p.value
}
# NTRK vs PDGFRA moderated t (three-group design so variance is pooled over all groups present)
limma_geno <- function(x, g) {
  fit <- limma_fit(x, g)
  fit <- eBayes(contrasts.fit(fit, makeContrasts("NTRK-PDGFRA", levels = fit$design)))
  data.frame(logFC = fit$coefficients[, 1], P = fit$p.value[, 1], adjP = p.adjust(fit$p.value[, 1], "BH"),
             row.names = rownames(x))
}
# max pairwise |log2FC| from group means
max_pairwise <- function(x, g) {
  g <- droplevels(g)
  gm <- sapply(levels(g), function(l) rowMeans(x[, g == l, drop = FALSE]))
  pr <- combn(ncol(gm), 2, function(i) abs(gm[, i[1]] - gm[, i[2]]))
  apply(pr, 1, max)
}
pvals <- function(method, x, g) switch(method, anova_param = anova_param_p(x, g),
                                       anova_perm = anova_perm_p(x, g), limma_F = limma_F_p(x, g))
# DE flag from p-values: adj.P < fdr and (optionally) FC criterion
de_flag <- function(p, mfc, fdr, fc) {
  adj <- p.adjust(p, "BH")
  ok <- !is.na(adj) & adj < fdr
  if (!is.na(fc)) ok <- ok & mfc >= fc
  ok
}

jaccard <- function(a, b) { u <- length(union(a, b)); if (u == 0) NA_real_ else length(intersect(a, b)) / u }

# ---- output helpers --------------------------------------------------------
append_metrics <- function(df, drop_variants = NULL) {
  f <- "results/metrics.csv"
  old <- read.csv(f, stringsAsFactors = FALSE)
  if (!is.null(drop_variants)) old <- old[!(old$claim == "C2" & old$variant %in% drop_variants), ]
  names(df) <- c("claim", "variant", "metric", "value", "n")
  write.csv(rbind(old, df), f, row.names = FALSE)   # quoted, as in the existing file
}
save_fig <- function(name, w, h, fn) {
  png(file.path(OUT, paste0(name, ".png")), w, h, units = "in", res = 150); fn(); dev.off()
  pdf(file.path(OUT, paste0(name, ".pdf")), w, h); fn(); dev.off()
}
write_gz <- function(df, f) { con <- gzfile(f, "w"); write.csv(df, con, row.names = FALSE); close(con) }
