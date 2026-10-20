# C2 Q-A: paper-style pipeline (permutation ANOVA + BH 0.05 + FC>=1.5 in any pairwise) on proteins,
# unscaled and per-sample rescaled, versus the paper's 4,703 list. Run from repo root.
set.seed(1)
source("analysis/C2_differential/improved/c2_lib.R")

raw <- load_protein()
bench <- rd("JUMPn/JUMPnc/HGG_whl_input.txt")
bacc <- bench$Protein_Accession
stopifnot(length(bacc) == 4703, !anyDuplicated(bacc))
unmatched <- setdiff(bacc, rownames(raw))
cat("Benchmark accessions not in raw file:", length(unmatched), "->", unmatched, "\n")

# per-sample factor: median over matched accessions of HGG_whl_input / raw
mt <- intersect(bacc, rownames(raw))
bm <- as.matrix(bench[match(mt, bacc), SAMPLES]); storage.mode(bm) <- "double"
fac <- sapply(SAMPLES, function(s) median(bm[, s] / raw[mt, s]))
cat("Rescaling factors (n matched =", length(mt), "):\n"); print(round(fac, 4))
resc <- sweep(raw, 2, fac, "*")

paper_style <- function(m) {
  x <- log2(m)
  p <- anova_perm_p(x, GROUP[colnames(x)], 1000)
  mfc <- max_pairwise(x, GROUP[colnames(x)])
  list(p = p, adjP = p.adjust(p, "BH"), mfc = mfc, fc_pass = mfc >= CUT["protein"],
       de = de_flag(p, mfc, 0.05, CUT["protein"]))
}
res <- list(unscaled = paper_style(raw), rescaled = paper_style(resc))

rows <- list(); sm <- NULL
for (v in names(res)) {
  de <- rownames(raw)[res[[v]]$de]
  tp <- length(intersect(de, bacc))
  sm <- rbind(sm, data.frame(variant = v, n_de = length(de), n_paper = length(bacc), tp = tp,
                             jaccard = jaccard(de, bacc), precision = tp / length(de), recall = tp / length(bacc),
                             n_fc_pass = sum(res[[v]]$fc_pass), n_adjP_lt_0.05 = sum(res[[v]]$adjP < 0.05),
                             n_proteins_tested = nrow(raw)))
}
write.csv(sm, file.path(OUT, "c2_qa_summary.csv"), row.names = FALSE)
print(sm)

# Q-C pieces (computed here because both pipelines are in memory)
ru <- res$unscaled; rr <- res$rescaled
de_u <- rownames(raw)[ru$de]; de_r <- rownames(raw)[rr$de]
qc <- data.frame(
  jaccard_lists_resc_vs_unscaled = jaccard(de_u, de_r),
  n_fc_pass_only_rescaled = sum(rr$fc_pass & !ru$fc_pass),
  n_fc_pass_only_unscaled = sum(!rr$fc_pass & ru$fc_pass),
  n_de_only_rescaled = length(setdiff(de_r, de_u)), n_de_only_unscaled = length(setdiff(de_u, de_r)),
  n_de_only_rescaled_due_to_fc = sum(rr$de & !ru$de & rr$adjP < 0.05 & ru$adjP < 0.05 & rr$fc_pass & !ru$fc_pass),
  n_de_only_unscaled_due_to_fc = sum(ru$de & !rr$de & ru$adjP < 0.05 & rr$adjP < 0.05 & ru$fc_pass & !rr$fc_pass))
write.csv(qc, file.path(OUT, "c2_qc_rescaling.csv"), row.names = FALSE)
print(t(qc))

write.csv(data.frame(sample = SAMPLES, factor = as.numeric(fac)), file.path(OUT, "c2_rescale_factors.csv"), row.names = FALSE)
write_gz(data.frame(Protein_Accession = rownames(raw), in_paper_list = rownames(raw) %in% bacc,
                    paper_style_unscaled = ru$de, paper_style_rescaled = rr$de,
                    adjP_unscaled = signif(ru$adjP, 4), adjP_rescaled = signif(rr$adjP, 4),
                    maxLog2FC_unscaled = round(ru$mfc, 4), maxLog2FC_rescaled = round(rr$mfc, 4)),
         file.path(OUT, "c2_qa_protein_lists.csv.gz"))

m <- rbind(
  data.frame(claim = "C2", variant = paste0("QA_", sm$variant), metric = "n_de", value = sm$n_de, n = sm$n_proteins_tested),
  data.frame(claim = "C2", variant = paste0("QA_", sm$variant), metric = "jaccard_vs_paper", value = sm$jaccard, n = sm$n_paper),
  data.frame(claim = "C2", variant = paste0("QA_", sm$variant), metric = "precision", value = sm$precision, n = sm$n_de),
  data.frame(claim = "C2", variant = paste0("QA_", sm$variant), metric = "recall", value = sm$recall, n = sm$n_paper),
  data.frame(claim = "C2", variant = "QC_rescaling", metric = "jaccard_resc_vs_unscaled", value = qc$jaccard_lists_resc_vs_unscaled, n = nrow(raw)),
  data.frame(claim = "C2", variant = "QC_rescaling", metric = "n_fc_pass_only_rescaled", value = qc$n_fc_pass_only_rescaled, n = nrow(raw)),
  data.frame(claim = "C2", variant = "QC_rescaling", metric = "n_fc_pass_only_unscaled", value = qc$n_fc_pass_only_unscaled, n = nrow(raw)),
  data.frame(claim = "C2", variant = "QC_rescaling", metric = "n_de_only_rescaled_due_to_fc", value = qc$n_de_only_rescaled_due_to_fc, n = nrow(raw)),
  data.frame(claim = "C2", variant = "QC_rescaling", metric = "n_de_only_unscaled_due_to_fc", value = qc$n_de_only_unscaled_due_to_fc, n = nrow(raw)))
append_metrics(m, drop_variants = unique(m$variant))
