# C11 Q-A: concordance of log2 ratios across layers (RNA vs protein, best-site phospho vs protein). Pre-registered in LOG.md.
# Run from the repo root after c2b_convert_xlsx.R and c11_convert_sd3.R: Rscript analysis/C11_rna_layer/original/c11_qa_concordance.R
set.seed(1)
source("analysis/C11_rna_layer/improved/c11_lib.R")
OUT_A <- file.path(C11, "original"); SPCUT <- 0.5

sd3 <- read_sd3(); sd1 <- read_sd1(); sd2 <- read_sd2a(); ours <- our_protein_ratios()
cat("SD3 rows:", length(sd3$gene), " unique (case-insensitive):", length(unique(key(sd3$gene))), "\n")
print(round(sd3$ratio_identity, 3))

# gene-level tables: one row per gene key, contrast columns
gene_level <- function(genes, R, score) {            # keep the highest-score row per gene
  i <- dedup_idx(genes, score); list(key = key(genes[i]), R = R[i, , drop = FALSE], idx = i)
}
rna <- gene_level(sd3$gene, sd3$R, rowMeans(sd3$F))
pro <- gene_level(sd1$d$gene, sd1$R, sd1$d$mean_int)
pro$acc <- sd1$d$Protein_Accession[pro$idx]
prv <- gene_level(ours$d$gene, ours$R, ours$d$mean_int)
prv$acc <- ours$d$Protein_Accession[prv$idx]
# phospho: per gene and contrast the site with the largest |ratio|
ph <- list(key = unique(key(sd2$d$gene))); ph$R <- matrix(NA_real_, length(ph$key), 3, dimnames = list(ph$key, RC)); ph$site <- ph$R; mode(ph$site) <- "character"
kk <- key(sd2$d$gene)
for (ct in RC) {
  o <- order(-abs(sd2$R[, ct])); b <- o[!duplicated(kk[o])]            # first = largest |ratio| per gene
  ph$R[kk[b], ct] <- sd2$R[b, ct]; ph$site[kk[b], ct] <- sd2$d$site_id[b]
}
cat(sprintf("genes: RNA %d, protein(paper) %d, protein(ours) %d, phospho %d\n", length(rna$key), length(pro$key), length(prv$key), length(ph$key)))

pair_stats <- function(a, b, label_a, label_b, source_b) {
  do.call(rbind, lapply(RC, function(ct) {
    g <- intersect(rownames(a), rownames(b)); x <- a[g, ct]; y <- b[g, ct]; ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
    big <- abs(x) >= 1 & abs(y) >= 1
    data.frame(layer_a = label_a, layer_b = label_b, protein_source = source_b, contrast = ct, n_genes = length(x),
               spearman = cor(x, y, method = "spearman"), n_both_abs_ge1 = sum(big), sign_agree_abs_ge1 = if (any(big)) mean(sign(x[big]) == sign(y[big])) else NA_real_)
  }))
}
M <- function(o) { m <- o$R; rownames(m) <- o$key; m }
RNA <- M(rna); PRO <- M(pro); PRV <- M(prv); PH <- ph$R
S <- rbind(pair_stats(RNA, PRO, "RNA", "protein", "paper_SD1"), pair_stats(PH, PRO, "phospho_best_site", "protein", "paper_SD1"),
           pair_stats(RNA, PRV, "RNA", "protein", "ours_rescaled_all10"), pair_stats(PH, PRV, "phospho_best_site", "protein", "ours_rescaled_all10"),
           pair_stats(RNA, PH, "RNA", "phospho_best_site", "paper_SD2a"))
S$concordance_rule_applies <- S$contrast %in% c("NTRK_Cortex", "PDGFRA_Cortex") & S$layer_b == "protein"
S$meets_0.5 <- S$spearman >= SPCUT
rownames(S) <- NULL
write.csv(S, file.path(OUT_A, "c11_qa_concordance.csv"), row.names = FALSE); print(S, digits = 3)
V <- do.call(rbind, lapply(unique(S$protein_source[S$layer_b == "protein"]), function(ps) {
  d <- S[S$layer_b == "protein" & S$protein_source == ps & S$concordance_rule_applies, ]
  data.frame(protein_source = ps, n_tests = nrow(d), n_meeting_0.5 = sum(d$meets_0.5), concordant = all(d$meets_0.5)) }))
write.csv(V, file.path(OUT_A, "c11_qa_decision.csv"), row.names = FALSE); print(V)

# tidy long table (genes present in >= 2 layers), with accessions / site IDs
allk <- Reduce(union, list(rownames(RNA), rownames(PRO), rownames(PH)))
lk <- function(m, ct) m[match(allk, rownames(m)), ct]      # NA for genes absent from a layer
long <- do.call(rbind, lapply(RC, function(ct) {
  data.frame(gene_key = allk, contrast = ct, rna_log2ratio = lk(RNA, ct), protein_log2ratio_paper = lk(PRO, ct),
             Protein_Accession_paper = pro$acc[match(allk, pro$key)], protein_log2ratio_ours = lk(PRV, ct),
             Protein_Accession_ours = prv$acc[match(allk, prv$key)], phospho_log2ratio = lk(PH, ct),
             phospho_site_id = lk(ph$site, ct), row.names = NULL)
}))
n_layers <- rowSums(!is.na(long[, c("rna_log2ratio", "protein_log2ratio_paper", "phospho_log2ratio")]))
long <- long[n_layers >= 2, ]
for (cn in grep("log2ratio", names(long), value = TRUE)) long[[cn]] <- round(long[[cn]], 4)
write_gz(long, file.path(OUT_A, "c11_qa_gene_ratios.csv.gz"))

# ---- figure 1 ------------------------------------------------------------------------------------------------------
save_fig11(OUT_A, "c11_fig1_ratio_scatter", 11, 7, function() {
  par(mfrow = c(2, 3), mar = c(4.2, 4.2, 3, 1))
  pl <- list(list(RNA, "RNA", "rna", "RNA"), list(PH, "phospho (best site)", "phospho", "phospho"))
  for (a in pl) for (ct in RC) {
    g <- intersect(rownames(a[[1]]), rownames(PRO)); x <- PRO[g, ct]; y <- a[[1]][g, ct]; ok <- is.finite(x) & is.finite(y)
    r <- cor(x[ok], y[ok], method = "spearman")
    lim <- c(-6, 6); xy <- unique(round(cbind(pmax(pmin(x[ok], 6), -6), pmax(pmin(y[ok], 6), -6)), 2))   # display only: rounded, de-duplicated (keeps the PDF small)
    plot(xy[, 1], xy[, 2], pch = 16, cex = 0.3, col = rgb(0, 0, 0, 0.25), xlim = lim, ylim = lim,
         xlab = "protein log2 ratio (SD1); display clipped at +/-6", ylab = paste(a[[2]], "log2 ratio"), main = gsub("_", "/", ct), cex.main = 1)
    abline(h = 0, v = 0, col = "grey70", lty = 3); abline(0, 1, col = "#D55E00", lty = 2)
    legend("topleft", sprintf("Spearman %.2f\nn = %d genes", r, sum(ok)), bty = "n", cex = 0.9)
  }
})
# ---- metrics ------------------------------------------------------------------------------------------------------
mm <- list()
for (i in seq_len(nrow(S))) { v <- paste("QA", S$layer_a[i], "vs", S$layer_b[i], S$protein_source[i], S$contrast[i], sep = "_")
  mm[[length(mm) + 1]] <- m11(v, "spearman", S$spearman[i], S$n_genes[i])
  mm[[length(mm) + 1]] <- m11(v, "sign_agree_abs_ge1", S$sign_agree_abs_ge1[i], S$n_both_abs_ge1[i]) }
for (i in seq_len(nrow(V))) mm[[length(mm) + 1]] <- m11(paste0("QA_decision_", V$protein_source[i]), "n_meeting_spearman_0.5_of_4", V$n_meeting_0.5[i], 4)
put_metrics11(do.call(rbind, mm), "QA")
dir.create("env", showWarnings = FALSE); writeLines(capture.output(sessionInfo()), "env/sessionInfo_C11.txt")
