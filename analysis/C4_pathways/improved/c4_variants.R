# C4 improved arm: Q-B variants V1-V7, Q-C (normalisation), Q-D (method agreement), figures 1-2. Run after original/c4_main.R.
set.seed(1)
source("analysis/C4_pathways/improved/c4_lib.R")
ctx <- build_ctx(); sets <- ctx$sets; ortho <- ctx$ortho; bg <- ctx$bg; stat <- ctx$stat
HAVE_GSEA <- requireNamespace("fgsea", quietly = TRUE)

# sanity: hypergeometric tail == one-sided Fisher exact test
{ tt <- matrix(c(7, 3, 5, 85), 2)
  stopifnot(isTRUE(all.equal(fisher.test(tt, alternative = "greater")$p.value, phyper(7 - 1, 12, 88, 10, lower.tail = FALSE)))) }

pf <- read.csv(gzfile("analysis/C2_differential/c2_protein_DE_flags.csv.gz"), check.names = FALSE, stringsAsFactors = FALSE)
hf <- read.csv(gzfile("analysis/C2_differential/c2_phospho_DE_flags.csv.gz"), check.names = FALSE, stringsAsFactors = FALSE)
nf <- read.csv(gzfile("analysis/C2_differential/c2_phospho_normalisation_flags.csv.gz"), stringsAsFactors = FALSE)
deP <- genes_of(pf$Feature[pf[["anova_perm_FC1.5_FDR0.05"]]], ctx$maps$protein)
deH <- genes_of(hf$Feature[hf[["anova_perm_FC2_FDR0.05"]]], ctx$maps$phospho)
sbP <- sets_in_bg(sets, bg$protein, "table", ortho); sbH <- sets_in_bg(sets, bg$phospho, "table", ortho)

tabs <- list(); cls <- list(); summ <- NULL
run <- function(v, tP, tH, thr = 0.05, desc) {
  tabs[[v]] <<- rbind(tag(tP, layer = "proteome"), tag(tH, layer = "phospho"))
  cl <- classify_wide(tP, tH, thr); cls[[v]] <<- cl
  summ <<- rbind(summ, data.frame(variant = v, description = desc, threshold = thr,
                                  n_enriched_all_proteome = sum(!is.na(tP$BH) & tP$BH < thr),
                                  n_enriched_all_phospho = sum(!is.na(tH$BH) & tH$BH < thr), claim_summary(cl)))
}
run("main", ora(deP, bg$protein, sbP), ora(deH, bg$phospho, sbH), desc = "ORA, own backgrounds, BH<0.05, anova_perm lists")

# V1 common background
common <- intersect(bg$protein, bg$phospho)
sbC <- sets_in_bg(sets, common, "table", ortho)
run("V1_common_bg", ora(intersect(deP, common), common, sbC), ora(intersect(deH, common), common, sbC),
    desc = "background = genes measured in both layers; lists restricted to it")
# V2 size-matched top-N
N <- min(length(deP), length(deH)); topN <- function(s) names(sort(s, decreasing = TRUE))[seq_len(N)]
run("V2_size_matched", ora(topN(stat$protein), bg$protein, sbP), ora(topN(stat$phospho), bg$phospho, sbH),
    desc = sprintf("top-%d genes by -log10(limma F adj.P) in each layer", N))
# V3 normalised phospho list (background unchanged = all phospho genes, as pre-registered)
mapN <- phospho_map(nf$Site)
deHn <- genes_of(nf$Site[nf$DE_normalised], ctx$maps$phospho)
run("V3_phospho_normalised", ora(deP, bg$protein, sbP), ora(deHn, bg$phospho, sbH),
    desc = "phospho list = DE after normalisation to protein (limma F, FC2, adj.P<0.05)")
# V4 FCS
fP <- fcs(stat$protein, sbP); fH <- fcs(stat$phospho, sbH)
run("V4_FCS", fP, fH, desc = "Wilcoxon rank-sum FCS on -log10(limma F adj.P)")
# V5 GSEA
if (HAVE_GSEA) { gP <- gsea(stat$protein, sbP); gH <- gsea(stat$phospho, sbH)
  run("V5_GSEA", gP, gH, desc = "fgsea preranked (scoreType pos), same statistic") }
# V6 BH < 0.2
run("V6_BH0.2", ora(deP, bg$protein, sbP), ora(deH, bg$phospho, sbH), thr = 0.2, desc = "main-setting ORA, BH<0.2")
# V7 upper-case symbols
sbPu <- sets_in_bg(sets, bg$protein, "upper"); sbHu <- sets_in_bg(sets, bg$phospho, "upper")
run("V7_uppercase", ora(deP, bg$protein, sbPu), ora(deH, bg$phospho, sbHu), desc = "mouse symbols upper-cased, matched to human symbols")

all_tabs <- do.call(rbind, Map(function(t, v) tag(t, variant = v), tabs, names(tabs)))
write_gz(all_tabs[, c("variant", "layer", "set", "size", "overlap", "n_list", "n_bg", "p", "BH")], file.path(OUT, "c4_variant_enrichment_tables.csv.gz"))
class_long <- do.call(rbind, Map(function(c, v) tag(c, variant = v), cls, names(cls)))
write.csv(class_long, file.path(OUT, "c4_S_wide_classification_per_variant.csv"), row.names = FALSE)
write.csv(summ, file.path(OUT, "c4_qb_variant_summary.csv"), row.names = FALSE)
print(summ[, -2]); 
nvar <- nrow(summ); nhold <- sum(summ$claim_holds)
cat(sprintf("\nQ-B: claim holds in %d of %d variants (%.2f); holds in main: %s; GSEA run: %s\n", nhold, nvar, nhold / nvar, summ$claim_holds[1], HAVE_GSEA))
# descriptive: S_core + cell cycle per variant
d <- all_tabs[all_tabs$set %in% c(S_CORE, CTRL), c("variant", "layer", "set", "size", "overlap", "p", "BH")]
write.csv(d, file.path(OUT, "c4_S_core_and_cell_cycle_per_variant.csv"), row.names = FALSE)
cc <- d[d$set == CTRL & d$variant %in% c("main", "V4_FCS"), ]; print(cc)

# ---- Q-C: normalisation at pathway level (same matched-site basis for both lists) ----
duG <- genes_of(nf$Site[nf$DE_unnormalised], ctx$maps$phospho)
tU <- ora(duG, bg$phospho, sbH); tN <- ora(deHn, bg$phospho, sbH); stopifnot(identical(tU$set, tN$set))
eU <- enriched(tU); eN <- enriched(tN)
qc <- data.frame(n_genes_unnormalised = length(duG), n_genes_normalised = length(deHn), n_sets_tested = nrow(tU),
                 n_enriched_unnormalised = length(eU), n_enriched_normalised = length(eN), n_both = length(intersect(eU, eN)),
                 jaccard = jaccard(eU, eN), spearman_neglog10p = cor(-log10(tU$p), -log10(tN$p), method = "spearman"))
qc$agree <- !is.na(qc$jaccard) && qc$jaccard >= 0.8
write.csv(qc, file.path(OUT, "c4_qc_normalisation.csv"), row.names = FALSE)
write.csv(data.frame(set = tU$set, size = tU$size, p_unnorm = tU$p, BH_unnorm = tU$BH, p_norm = tN$p, BH_norm = tN$BH),
          file.path(OUT, "c4_qc_per_set.csv"), row.names = FALSE)
cat("\nQ-C:\n"); print(t(qc))

# ---- Q-D: method agreement (ORA vs FCS, ORA vs GSEA), all sets ----------------
qd <- NULL
for (ly in c("proteome", "phospho")) {
  base <- tabs$main[tabs$main$layer == ly, ]; eo <- enriched(base)
  for (m in c("V4_FCS", if (HAVE_GSEA) "V5_GSEA")) {
    o <- tabs[[m]][tabs[[m]]$layer == ly, ]; e <- enriched(o)
    qd <- rbind(qd, data.frame(layer = ly, comparison = paste0("ORA_vs_", sub("^V[0-9]_", "", m)), n_sets_ora = nrow(base), n_sets_other = nrow(o),
      n_enriched_ora = length(eo), n_enriched_other = length(e), n_both = length(intersect(eo, e)), jaccard = jaccard(eo, e),
      spearman_neglog10p = cor(-log10(base$p[match(o$set, base$set)]), -log10(o$p), method = "spearman")))
  }
  if (HAVE_GSEA) { e1 <- enriched(tabs$V4_FCS[tabs$V4_FCS$layer == ly, ]); e2 <- enriched(tabs$V5_GSEA[tabs$V5_GSEA$layer == ly, ])
    qd <- rbind(qd, data.frame(layer = ly, comparison = "FCS_vs_GSEA", n_sets_ora = NA, n_sets_other = NA, n_enriched_ora = length(e1),
      n_enriched_other = length(e2), n_both = length(intersect(e1, e2)), jaccard = jaccard(e1, e2), spearman_neglog10p = NA)) }
}
qd$agree <- !is.na(qd$jaccard) & qd$jaccard >= 0.8
write.csv(qd, file.path(OUT, "c4_qd_method_agreement.csv"), row.names = FALSE); cat("\nQ-D:\n"); print(qd)

# ---- figures ----------------------------------------------------------------------
cols <- c("phospho-only" = "#D55E00", "proteome-only" = "#0072B2", both = "#009E73", neither = "#E8E8E8")
lab <- c("phospho-only" = "P", "proteome-only" = "R", both = "B", neither = "")
vs <- summ$variant
save_fig("c4_fig1_S_wide_classification", 9, 5.5, function() {
  par(mar = c(11, 17, 3, 1))
  plot(NA, xlim = c(0.5, length(vs) + 0.5), ylim = c(0.5, 9.5), axes = FALSE, xlab = "", ylab = "", yaxs = "i", xaxs = "i")
  for (i in seq_along(vs)) for (j in seq_along(S_WIDE)) {
    k <- cls[[vs[i]]]$class[j]; y <- 10 - j
    rect(i - .5, y - .5, i + .5, y + .5, col = cols[k], border = "white", lwd = 2)
    text(i, y, lab[k], cex = 0.9, font = 2, col = "white")
  }
  axis(1, at = seq_along(vs), labels = vs, las = 2, tick = FALSE, cex.axis = 0.8, line = -0.5)
  axis(2, at = 9:1, labels = sub("KEGG_", "", S_WIDE), las = 1, tick = FALSE, cex.axis = 0.8)
  title("S_wide pathway classification per variant (BH cut-off as registered)", cex.main = 0.9)
  legend(x = length(vs) / 2 + 0.5, y = -4.5, xjust = 0.5, xpd = TRUE, horiz = TRUE, bty = "n", cex = 0.8, fill = cols,
         legend = c("phospho-only (P)", "proteome-only (R)", "both (B)", "neither"))
})
bars <- rbind(data.frame(lab = "Q-C unnormalised\nvs normalised\n(phospho)", j = qc$jaccard, grp = "Q-C"),
              data.frame(lab = paste0("Q-D ", qd$layer, "\n", sub("_", " ", qd$comparison)), j = qd$jaccard, grp = "Q-D"))
save_fig("c4_fig2_jaccard_QC_QD", 10, 5, function() {
  par(mar = c(7, 4.5, 3, 1)); v <- ifelse(is.na(bars$j), 0, bars$j)
  b <- barplot(v, col = ifelse(bars$grp == "Q-C", "#CC79A7", "#0072B2"), ylim = c(0, 1), ylab = "Jaccard of enriched sets", names.arg = NA)
  text(b, -0.03, bars$lab, srt = 45, adj = 1, xpd = TRUE, cex = 0.7)
  text(b, v + 0.03, ifelse(is.na(bars$j), "undefined\n(0 vs 0 sets)", sprintf("%.2f", bars$j)), cex = 0.7)
  abline(h = 0.8, lty = 2, lwd = 1.5); text(max(b), 0.83, "agree >= 0.8", adj = 1, cex = 0.8)
  title("Pathway-level agreement across normalisation (Q-C) and methods (Q-D)", cex.main = 0.9)
})

# ---- metrics ------------------------------------------------------------------------
m <- NULL
for (i in seq_len(nrow(summ))) { v <- paste0("QB_", summ$variant[i])
  m <- rbind(m, mrow(v, "S_wide_enriched_either", summ$n_enriched_either[i], 9), mrow(v, "S_wide_phospho_only_share", summ$phospho_only_share[i], summ$n_enriched_either[i]),
             mrow(v, "claim_holds", as.integer(summ$claim_holds[i]), 9)) }
m <- rbind(m, mrow("QB_all", "n_variants_claim_holds", nhold, nvar),
           mrow("QC_normalisation", "jaccard_enriched_sets", qc$jaccard, qc$n_sets_tested),
           mrow("QC_normalisation", "spearman_neglog10p", qc$spearman_neglog10p, qc$n_sets_tested))
for (i in seq_len(nrow(qd))) m <- rbind(m, mrow(paste0("QD_", qd$layer[i], "_", qd$comparison[i]), "jaccard_enriched_sets", qd$jaccard[i], qd$n_sets_other[i]))
put_metrics(m)
writeLines(capture.output(sessionInfo()), "env/sessionInfo_C4.txt")
