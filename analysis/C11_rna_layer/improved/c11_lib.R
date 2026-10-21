# Shared helpers for C11 (design pre-registered in LOG.md, "PRE-REGISTRATION: C11 RNA layer"). Run from the repo root.
# Needs data/external/csv/{sd1,sd2a,sd3}.csv: Rscript analysis/C2_differential/improved/c2b_convert_xlsx.R and
# analysis/C11_rna_layer/improved/c11_convert_sd3.R create them.
options(stringsAsFactors = FALSE)
stopifnot(file.exists("CLAUDE.md"))
source("analysis/C6_tf_activity/improved/c6_lib.R")   # also sources c4_lib.R and c2_lib.R (limma, ora, summarise_tab, BENCH, build_c6, append_metrics)
C11 <- "analysis/C11_rna_layer"
RC <- c("NTRK_PDGFRA", "NTRK_Cortex", "PDGFRA_Cortex")
RC_PAPER <- c("NTRK/PDGFRA", "NTRK/Cortex", "PDGFRA/Cortex")
SD_ORDER <- c("Cortex1","Cortex2","Cortex3","PD1","PD2","PD3","PD4","TRK1","TRK2","TRK3")   # order of the 10 intensity columns in SD1 / SD2a

rdx <- function(f) read.csv(file.path("data/external/csv", f), header = FALSE, colClasses = "character", na.strings = "")
num <- function(m) { m <- as.matrix(m); storage.mode(m) <- "double"; m }
key <- function(g) toupper(g)                          # case-insensitive gene matching (pre-registered)
save_fig11 <- function(dir, name, w, h, fn) {          # PNG
  png(file.path(dir, paste0(name, ".png")), w, h, units = "in", res = 150); fn(); invisible(dev.off())
}
put_metrics11 <- function(df, prefixes) {              # idempotent: replaces C11 rows whose variant starts with one of `prefixes`
  names(df) <- c("claim", "variant", "metric", "value", "n")
  f <- "results/metrics.csv"; old <- read.csv(f, stringsAsFactors = FALSE)
  drop <- old$claim == "C11" & grepl(paste0("^(", paste(prefixes, collapse = "|"), ")"), old$variant)
  write.csv(old[!drop, ], f, row.names = FALSE); append_metrics(df)
}
m11 <- function(variant, metric, value, n) data.frame(claim = "C11", variant = variant, metric = metric, value = value, n = n)

# ---- ARI and Ward helper (same definitions as analysis/C1_exploration/improved/c1_exploration.R; copied because that file runs on source) ----
ari <- function(a, b) {
  ok <- !is.na(a) & !is.na(b); a <- a[ok]; b <- b[ok]
  tab <- table(a, b); n <- sum(tab); c2 <- function(x) x * (x - 1) / 2
  s <- sum(c2(tab)); sa <- sum(c2(rowSums(tab))); sb <- sum(c2(colSums(tab)))
  e <- sa * sb / c2(n); out <- (s - e) / ((sa + sb) / 2 - e); if (is.nan(out)) NA_real_ else out
}
zdist <- function(X) { Z <- scale(t(X)); Z[is.na(Z)] <- 0; dist(Z) }          # C1 "euclid_z": samples z-scored over features
run_ward <- function(X, k) cutree(hclust(zdist(X), method = "ward.D"), k = k)   # C1 paper setup: ward.D on euclid_z
top_mad <- function(X, n = 1000) { X <- X[apply(X, 1, sd) > 0, , drop = FALSE]; X[order(apply(X, 1, mad), decreasing = TRUE)[seq_len(min(n, nrow(X)))], , drop = FALSE] }

# ---- Supp Data 3 (RNA) ------------------------------------------------------------------------------------------
read_sd3 <- function() {
  x <- rdx("sd3.csv")
  stopifnot(x[1, 1] == "Gene Names", x[1, 2] == "Cortex", x[1, 6] == "PDGFRA Batch1", x[1, 9] == "PDGFRA Batch2",
            x[1, 13] == "NTRK batch 1", x[1, 16] == "NTRK Batch2", x[1, 23] == "B.H. adjusted ANOVA p Values",
            identical(unname(unlist(x[2, 2:19])), as.character(1:18)), identical(unname(unlist(x[2, 20:22])), RC_PAPER))
  x <- x[-(1:2), ]; stopifnot(!anyNA(x[, 1]))
  F <- num(x[, 2:19]); R <- num(x[, 20:22]); colnames(R) <- RC
  meta <- data.frame(sample = paste0("S", 1:18),
                     group = rep(c("Cortex", "PDGFRA", "NTRK"), c(4, 7, 7)),
                     batch = c(rep(NA, 4), rep(1, 3), rep(2, 4), rep(1, 3), rep(2, 4)))
  colnames(F) <- meta$sample
  p <- as.numeric(x[, 23]); gene <- x[, 1]
  stopifnot(!anyNA(F), !anyNA(R), !anyNA(p), min(F) > 0.9)
  # verify the ratio column order by values (each paper column must correlate best with its own recomputation)
  gm <- function(i) rowMeans(F[, i]); L <- log2(cbind(C = gm(1:4), P = gm(5:11), N = gm(12:18)))
  cand <- cbind(NTRK_PDGFRA = L[, "N"] - L[, "P"], NTRK_Cortex = L[, "N"] - L[, "C"], PDGFRA_Cortex = L[, "P"] - L[, "C"])
  idc <- sapply(RC, function(pc) sapply(RC, function(oc) cor(R[, pc], cand[, oc]))); stopifnot(identical(RC[apply(idc, 2, which.max)], RC))
  list(gene = gene, F = F, R = R, p = p, meta = meta, ratio_identity = idc)
}
# duplicated gene names: keep the row with the highest mean value
dedup_idx <- function(gene, score) { o <- order(-score); keep <- o[!duplicated(key(gene)[o])]; sort(keep) }

# ---- protein / phospho (paper ratio columns) and our rescaled all-10 protein ratios ------------------------------
read_sd1 <- function() {
  x <- rdx("sd1.csv"); stopifnot(x[1, 1] == "Coexpresion Clusters", x[2, 16] == "NTRK/PDGFRA", x[2, 17] == "NTRK/Cortex", x[2, 18] == "PDGFRA/Cortex")
  x <- x[-(1:2), ]; x <- x[!is.na(x[, 3]), ]
  R <- num(x[, 16:18]); colnames(R) <- RC; I <- num(x[, 6:15]); colnames(I) <- SD_ORDER
  d <- data.frame(Protein_Accession = x[, 3], gene = x[, 2], mean_int = rowMeans(I), stringsAsFactors = FALSE)
  list(d = d, R = R)
}
read_sd2a <- function() {
  y <- rdx("sd2a.csv"); stopifnot(y[1, 1] == "Coexpresion Clusters", y[2, 17] == "NTRK/PDGFRA", y[2, 18] == "NTRK/Cortex", y[2, 19] == "PDGFRA/Cortex")
  y <- y[-(1:2), ]; y <- y[!is.na(y[, 3]), ]
  R <- num(y[, 17:19]); colnames(R) <- RC
  d <- data.frame(site_id = paste(y[, 2], sapply(strsplit(y[, 3], "[|]"), `[`, 2), y[, 4], sep = "_"), gene = y[, 2], stringsAsFactors = FALSE)
  stopifnot(!anyDuplicated(d$site_id), nrow(d) == 28527)
  list(d = d, R = R)
}
our_protein_ratios <- function() {      # all 10 samples, C2 per-sample rescaling (= the paper's SD1 scale), log2 of ratio of group means
  raw <- load_protein(); fac <- read.csv("analysis/C2_differential/c2_rescale_factors.csv"); fac <- setNames(fac$factor, fac$sample)[SAMPLES]
  m <- sweep(raw, 2, fac, "*"); mm <- function(s) rowMeans(m[, s])
  C <- mm(c("Cortex1","Cortex2","Cortex3")); P <- mm(c("PD1","PD2","PD3","PD4")); N <- mm(c("TRK1","TRK2","TRK3"))
  pm <- protein_map(); a <- intersect(rownames(m), names(pm))
  R <- cbind(NTRK_PDGFRA = log2(N / P), NTRK_Cortex = log2(N / C), PDGFRA_Cortex = log2(P / C))[a, ]
  list(d = data.frame(Protein_Accession = a, gene = unname(pm[a]), mean_int = rowMeans(m[a, ]), stringsAsFactors = FALSE), R = R)
}
