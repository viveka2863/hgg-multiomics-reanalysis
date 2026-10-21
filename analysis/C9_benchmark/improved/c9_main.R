# C9: CRISPR benchmark of the rankings (pre-registered in LOG.md, 2026-10-21). Run from the repo root.
source("analysis/C2_differential/improved/c2_lib.R")
set.seed(1)
OUT9 <- "analysis/C9_benchmark"
dir.create(file.path(OUT9, "original"), showWarnings = FALSE); dir.create(file.path(OUT9, "improved"), showWarnings = FALSE)

# ---- benchmarks ------------------------------------------------------------
B9_hit  <- c("Jun","Myc","Prkaa1","Prkaa2","Eef2k"); B9_non <- c("Brca1","Ep300","Insr","Jund")
B15_hit <- c(B9_hit, "Cd93"); B15_non <- c(B9_non, "Cd74","Epha2","Spry1","Arhgap18","Dab2")
bench <- rbind(data.frame(gene = c(B15_hit, B15_non), hit = rep(c(TRUE, FALSE), c(length(B15_hit), length(B15_non)))))
bench$in_B9 <- bench$gene %in% c(B9_hit, B9_non)
KEY <- toupper(bench$gene)

# ---- helpers ---------------------------------------------------------------
# AUC = P(hit > non-hit), ties 0.5 (Mann-Whitney)
auc_fun <- function(s, hit) {
  ok <- !is.na(s); s <- s[ok]; hit <- hit[ok]
  h <- sum(hit); r <- rank(s)
  (sum(r[hit]) - h * (h + 1) / 2) / (h * sum(!hit))
}
# exact one-sided permutation p over all labelings with the same number of hits
exact_p <- function(s, hit) {
  ok <- !is.na(s); s <- s[ok]; hit <- hit[ok]
  n <- length(s); h <- sum(hit); r <- rank(s)
  obs <- (sum(r[hit]) - h * (h + 1) / 2) / (h * (n - h))
  cmb <- combn(n, h)
  perm <- (colSums(matrix(r[cmb], nrow = h)) - h * (h + 1) / 2) / (h * (n - h))
  c(AUC = obs, p = mean(perm >= obs - 1e-12), min_p = 1 / choose(n, h), n = n, n_hit = h, n_perms = choose(n, h))
}

# best-row-per-gene tables --------------------------------------------------
prot_meta <- rd("JUMPn/JUMPnk/PhosphoWholeNormalization/hgg_whole_v1_20160420.txt")
prot_meta <- prot_meta[!grepl("^co[|]", prot_meta[["Protein Accession #"]]), ]
P <- log2(load_protein()); stopifnot(identical(rownames(P), prot_meta[["Protein Accession #"]]))
prot_gene <- toupper(prot_meta$GN)
Q <- log2(load_phospho())
ph_raw <- rd("JUMPn/JUMPnk/Integration/pho_hgg_uniq_v5_20160513.txt"); ph_raw <- ph_raw[!grepl("CON_", ph_raw[[1]]), ]
stopifnot(nrow(ph_raw) == nrow(Q))
ph_gene <- toupper(sub("_.*", "", ph_raw[[1]]))

best_row <- function(genes, stat, key, decreasing = FALSE) {   # index of best row per key gene, NA if none
  sapply(key, function(k) { i <- which(genes == k); if (!length(i)) return(NA_integer_)
    i[if (decreasing) which.max(stat[i]) else which.min(stat[i])] })
}
# S1-S4 for a sample subset
compute_S14 <- function(keep) {
  g <- droplevels(GROUP[keep])
  stats <- function(M) {
    M <- M[, keep]
    list(Fadj = p.adjust(limma_F_p(M, g), "BH"), con = limma_geno(M, g))
  }
  sp <- stats(P); sq <- stats(Q)
  ip1 <- best_row(prot_gene, sp$Fadj, KEY);        ip2 <- best_row(prot_gene, sp$con$adjP, KEY)
  iq1 <- best_row(ph_gene, sq$Fadj, KEY);          iq2 <- best_row(ph_gene, sq$con$adjP, KEY)
  nl <- function(v, i) ifelse(is.na(i), NA_real_, -log10(v[i]))
  data.frame(gene = bench$gene, S1 = nl(sp$Fadj, ip1), S2 = nl(sp$con$adjP, ip2),
             S3 = nl(sq$Fadj, iq1), S4 = nl(sq$con$adjP, iq2),
             row_S1 = rownames(P)[ip1], row_S2 = rownames(P)[ip2], row_S3 = rownames(Q)[iq1], row_S4 = rownames(Q)[iq2],
             stringsAsFactors = FALSE)
}
all10 <- rep(TRUE, 10); names(all10) <- SAMPLES
no3 <- SAMPLES != "TRK3"
sc <- compute_S14(all10)
ip5 <- best_row(prot_gene, rowMeans(P), KEY, decreasing = TRUE)
sc$S5 <- ifelse(is.na(ip5), NA_real_, rowMeans(P)[ip5]); sc$row_S5 <- rownames(P)[ip5]
sc$S6 <- sapply(KEY, function(k) sum(ph_gene == k))        # 0 = no phosphosite row
sc <- cbind(bench, sc[, -1])
cat("Genes with no protein row:", paste(bench$gene[is.na(ip5)], collapse = ","), "\n")
cat("Genes with no phosphosite row:", paste(bench$gene[sc$S6 == 0], collapse = ","), "\n")
write.csv(sc, file.path(OUT9, "improved", "c9_gene_scores.csv"), row.names = FALSE)

SC <- paste0("S", 1:6)
SC_NAME <- c(S1 = "protein F", S2 = "protein NTRK-PDGFRA", S3 = "phospho F", S4 = "phospho NTRK-PDGFRA",
             S5 = "expression (baseline)", S6 = "# phosphosites (baseline)")

# ---- Q-A / Q-B / Q-D -------------------------------------------------------
auc_table <- function(d, benchmark, scores = SC, arm = "all10") {
  rows <- lapply(scores, function(s) { r <- exact_p(d[[s]], d$hit)
    data.frame(benchmark = benchmark, arm = arm, score = s, score_name = SC_NAME[[s]], AUC = r[["AUC"]], exact_p = r[["p"]],
               min_attainable_p = r[["min_p"]], n_genes = r[["n"]], n_hits = r[["n_hit"]], n_perms = r[["n_perms"]],
               genes_missing = paste(d$gene[is.na(d[[s]])], collapse = ";"), stringsAsFactors = FALSE) })
  do.call(rbind, rows)
}
b9 <- sc[sc$in_B9, ]
qa <- auc_table(b9, "B9"); qb <- auc_table(sc, "B15")
qab <- rbind(qa, qb)
qab$separates <- qab$AUC >= 0.8 & qab$exact_p < 0.05 / 6
base <- setNames(qab$AUC[qab$score == "S5"], qab$benchmark[qab$score == "S5"])
qab$AUC_minus_S5 <- qab$AUC - base[qab$benchmark]
qab$adds_value <- qab$separates & qab$AUC_minus_S5 >= 0.1
qab$exceeds_S5_by_0.1 <- qab$AUC_minus_S5 >= 0.1 & qab$score != "S5"
write.csv(qa, file.path(OUT9, "original", "c9_original_QA_B9_auc.csv"), row.names = FALSE)
write.csv(b9[, c("gene","hit",SC,paste0("row_",SC[1:5]))], file.path(OUT9, "original", "c9_original_B9_scores.csv"), row.names = FALSE)
write.csv(qab, file.path(OUT9, "improved", "c9_QA_QB_QD_auc.csv"), row.names = FALSE)

# per-gene ranks (1 = highest score; ties averaged) with hit labels, per benchmark
rank_tab <- function(d, benchmark) {
  do.call(rbind, lapply(SC, function(s) { ok <- !is.na(d[[s]])
    data.frame(benchmark = benchmark, score = s, gene = d$gene[ok], hit = d$hit[ok], value = d[[s]][ok],
               rank = rank(-d[[s]][ok], ties.method = "average"), n_ranked = sum(ok), stringsAsFactors = FALSE) }))
}
ranks <- rbind(rank_tab(b9, "B9"), rank_tab(sc, "B15"))
write.csv(ranks, file.path(OUT9, "improved", "c9_gene_ranks.csv"), row.names = FALSE)
write.csv(ranks[ranks$benchmark == "B9", ], file.path(OUT9, "original", "c9_original_B9_ranks.csv"), row.names = FALSE)

cat("\n== Minimum attainable p ==\n"); print(qab[, c("benchmark","score","n_genes","n_hits","n_perms","min_attainable_p")], row.names = FALSE)
cat("\n== Q-A/Q-B/Q-D ==\n"); print(qab[, c("benchmark","score","AUC","exact_p","separates","AUC_minus_S5","adds_value","exceeds_S5_by_0.1","genes_missing")], row.names = FALSE, digits = 3)

# ---- Q-C: type-specific rankings (descriptive) ------------------------------
tsp <- list()
add_ts <- function(type, score, genes, vals, hits) {
  d <- data.frame(gene = genes, value = vals, hit = genes %in% hits, stringsAsFactors = FALSE)
  r <- exact_p(d$value, d$hit)
  tsp[[length(tsp) + 1]] <<- data.frame(type = type, score = score, AUC = r[["AUC"]], exact_p = r[["p"]], min_attainable_p = r[["min_p"]],
                                        n_genes = r[["n"]], n_hits = r[["n_hit"]], genes = paste0(d$gene, "=", signif(d$value, 4), collapse = ";"),
                                        stringsAsFactors = FALSE)
}
kin <- c("PRKAA1","PRKAA2","EEF2K","INSR")
for (a in list(c("coded_3_2_3", "analysis/C5_kinase_activity/original/c5_original_limma_coded_3_2_3.csv"),
               c("all10", "analysis/C5_kinase_activity/improved/c5_improved_limma_all10.csv"))) {
  k <- read.csv(a[2], check.names = FALSE); k <- k[match(kin, k$Kinase), ]; stopifnot(!anyNA(k$Kinase))
  hk <- c("PRKAA1","PRKAA2","EEF2K")
  add_ts("K1", paste0("IKAP F -log10 adjP (", a[1], ")"), kin, -log10(k[["adj.P.Val"]]), hk)
  add_ts("K1", paste0("IKAP signed NTRK-PDGFRA effect (", a[1], ")"), kin, k[["group2...group1"]], hk)
}
en <- read.csv(gzfile("analysis/C6_tf_activity/improved/c6_enrichment_tables.csv.gz"), stringsAsFactors = FALSE)
tfs <- c("JUN","MYC","BRCA1","EP300","JUND")
for (u in c("own","U1","U3")) {
  e <- en[en$list == "anova_perm_FC1.5_FDR0.05" & en$method == "ORA" & en$universe == u & en$TF %in% tfs, ]
  e <- e[match(tfs, e$TF), ]; stopifnot(!anyNA(e$TF))
  add_ts("T1", paste0("ORA -log10 p (anova_perm list, ", u, ")"), tfs, -log10(e$p), c("JUN","MYC"))
}
ed <- read.csv("analysis/C8_network/original/c8_accepted_edges_R0.csv", stringsAsFactors = FALSE)
nt <- c("MYC","BRCA1","EP300","JUND")
nacc <- sapply(nt, function(t) sum(ed$TF == t & ed$accept_R0))
stopifnot(identical(unname(nacc), c(6L, 1L, 1L, 2L)))
add_ts("N1", "# R0-accepted edges", nt, as.numeric(nacc), "MYC")
qc <- do.call(rbind, tsp)
write.csv(qc, file.path(OUT9, "improved", "c9_QC_type_specific_auc.csv"), row.names = FALSE)
cat("\n== Q-C ==\n"); print(qc[, 1:7], row.names = FALSE, digits = 3)

# ---- Q-E: TRK3 arm -----------------------------------------------------------
sc3 <- compute_S14(no3)
sc3 <- cbind(bench, sc3[, -1]); sc3$S5 <- sc$S5; sc3$S6 <- sc$S6
write.csv(sc3[, c("gene","hit","in_B9","S1","S2","S3","S4")], file.path(OUT9, "improved", "c9_QE_gene_scores_noTRK3.csv"), row.names = FALSE)
qe <- rbind(auc_table(sc3[sc3$in_B9, ], "B9", SC[1:4], "noTRK3"), auc_table(sc3, "B15", SC[1:4], "noTRK3"),
            qab[qab$score %in% SC[1:4], names(qab)[1:11]])
qe <- qe[order(qe$benchmark, qe$score, qe$arm), ]
write.csv(qe, file.path(OUT9, "improved", "c9_QE_trk3_arm_auc.csv"), row.names = FALSE)
cat("\n== Q-E ==\n"); print(qe[, c("benchmark","score","arm","AUC","exact_p","min_attainable_p","n_genes")], row.names = FALSE, digits = 3)
qe$separates <- qe$AUC >= 0.8 & qe$exact_p < 0.05 / 6
cat("Any S1-S4 arm separating (no TRK3)? ", any(qe$separates[qe$arm == "noTRK3"]), "\n")

# ---- figures ---------------------------------------------------------------
COL <- c(hit = "#C0392B", non = "#7F8C8D")
save_fig9 <- function(name, w, h, fn) {
  f <- file.path(OUT9, "improved", name)
  png(paste0(f, ".png"), w, h, units = "in", res = 150); fn(); dev.off()
}
fig1 <- function() {
  par(mfrow = c(2, 3), mar = c(4, 6, 2.5, 1))
  for (s in SC) {
    r <- ranks[ranks$benchmark == "B15" & ranks$score == s, ]; r <- r[order(r$rank), ]
    plot(r$rank, rev(seq_len(nrow(r))), pch = ifelse(r$hit, 19, 1), col = ifelse(r$hit, COL[["hit"]], COL[["non"]]), cex = 1.3,
         xlim = c(0.5, 15.5), ylim = c(0.5, nrow(r) + 0.5), yaxt = "n", xlab = "rank (1 = highest score)", ylab = "",
         main = paste0(s, ": ", SC_NAME[[s]]), cex.main = 0.9)
    axis(2, at = rev(seq_len(nrow(r))), labels = r$gene, las = 1, cex.axis = 0.65)
    if (s == "S1") legend("bottomleft", c("CRISPR hit", "non-hit"), pch = c(19, 1), col = COL, bty = "n", cex = 0.7)
  }
}
save_fig9("c9_fig1_rank_plot", 11, 7, fig1)
fig2 <- function() {
  par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))
  for (b in c("B9", "B15")) {
    d <- qab[qab$benchmark == b, ]
    x <- barplot(d$AUC, names.arg = d$score, ylim = c(0, 1), col = ifelse(d$score %in% c("S5","S6"), "#BDC3C7", "#2E86C1"),
                 ylab = "AUC (hit vs non-hit)", main = paste0(b, " (n = ", paste(unique(d$n_genes), collapse = "-"), ")"))
    abline(h = 0.8, lty = 2, col = "#C0392B"); abline(h = 0.5, lty = 3)
    abline(h = d$AUC[d$score == "S5"], lty = 1, col = "#7F8C8D")
    text(x, d$AUC + 0.03, sprintf("p=%.2g", d$exact_p), cex = 0.6)
    legend("topright", c("AUC = 0.8", "S5 expression AUC"), lty = c(2, 1), col = c("#C0392B", "#7F8C8D"), bty = "n", cex = 0.7)
  }
}
save_fig9("c9_fig2_auc_bars", 10, 4.5, fig2)
fig3 <- function() {
  par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))
  for (b in c("B9", "B15")) {
    d <- qe[qe$benchmark == b, ]
    m <- rbind(all10 = d$AUC[d$arm == "all10"][order(d$score[d$arm == "all10"])], noTRK3 = d$AUC[d$arm == "noTRK3"][order(d$score[d$arm == "noTRK3"])])
    colnames(m) <- sort(unique(d$score))
    barplot(m, beside = TRUE, ylim = c(0, 1), col = c("#2E86C1", "#F39C12"), ylab = "AUC", main = b)
    abline(h = 0.8, lty = 2, col = "#C0392B")
    legend("topright", c("all 10 samples", "without TRK3"), fill = c("#2E86C1", "#F39C12"), bty = "n", cex = 0.7)
  }
}
save_fig9("c9_fig3_auc_with_without_TRK3", 10, 4.5, fig3)

# ---- metrics.csv -----------------------------------------------------------
m <- function(variant, metric, value, n) data.frame(claim = "C9", variant = variant, metric = metric, value = value, n = n)
mm <- do.call(rbind, c(
  lapply(seq_len(nrow(qab)), function(i) with(qab[i, ], rbind(m(paste0(benchmark, "_", score), "AUC", AUC, n_genes),
                                                              m(paste0(benchmark, "_", score), "exact_p", exact_p, n_genes),
                                                              m(paste0(benchmark, "_", score), "min_attainable_p", min_attainable_p, n_genes)))),
  lapply(seq_len(nrow(qc)), function(i) with(qc[i, ], rbind(m(paste0(type, "_", gsub("[^A-Za-z0-9]+", "_", score)), "AUC", AUC, n_genes),
                                                            m(paste0(type, "_", gsub("[^A-Za-z0-9]+", "_", score)), "min_attainable_p", min_attainable_p, n_genes)))),
  lapply(which(qe$arm == "noTRK3"), function(i) with(qe[i, ], rbind(m(paste0(benchmark, "_", score, "_noTRK3"), "AUC", AUC, n_genes),
                                                                    m(paste0(benchmark, "_", score, "_noTRK3"), "exact_p", exact_p, n_genes))))))
old <- read.csv("results/metrics.csv", stringsAsFactors = FALSE)
write.csv(rbind(old[old$claim != "C9", ], mm), "results/metrics.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "env/sessionInfo_C9.txt")
cat("\nDone.\n")
