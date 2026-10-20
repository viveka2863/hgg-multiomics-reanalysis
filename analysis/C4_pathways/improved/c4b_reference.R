# C4 part B: Q-G1 (detection control), Q-G2 (reference variants V8/V9), Q-G3 (normalisation under V8).
# Pre-registered in LOG.md ("PRE-REGISTRATION: C4 part B"). Run from the repo root. Exploratory.
set.seed(1)
source("analysis/C4_pathways/improved/c4_lib.R")
ctx <- build_ctx(); sets <- ctx$sets; ortho <- ctx$ortho; bg <- ctx$bg

U1 <- sort(unique(ortho$V2))
sbU1 <- sets_in_bg(sets, U1, "table", ortho)
U2 <- sort(unique(unlist(sbU1, use.names = FALSE)))
cat(sprintf("U1 = %d mouse symbols; U2 = %d mouse genes in %d KEGG sets\n", length(U1), length(U2), length(sets)))
REF <- list(U1 = U1, U2 = U2)
SB  <- list(U1 = sbU1, U2 = sets_in_bg(sets, U2, "table", ortho))

# detected gene lists
detP <- bg$protein; detH <- bg$phospho
cat(sprintf("detected: proteome %d, phospho %d (pre-registered: 11,717 and 5,421)\n", length(detP), length(detH)))

# DE lists exactly as in the C4 main setting
pf <- read.csv(gzfile("analysis/C2_differential/c2_protein_DE_flags.csv.gz"), check.names = FALSE, stringsAsFactors = FALSE)
hf <- read.csv(gzfile("analysis/C2_differential/c2_phospho_DE_flags.csv.gz"), check.names = FALSE, stringsAsFactors = FALSE)
nf <- read.csv(gzfile("analysis/C2_differential/c2_phospho_normalisation_flags.csv.gz"), stringsAsFactors = FALSE)
deP <- genes_of(pf$Feature[pf[["anova_perm_FC1.5_FDR0.05"]]], ctx$maps$protein)
deH <- genes_of(hf$Feature[hf[["anova_perm_FC2_FDR0.05"]]], ctx$maps$phospho)
deHn <- genes_of(nf$Site[nf$DE_normalised], ctx$maps$phospho)
deHu <- genes_of(nf$Site[nf$DE_unnormalised], ctx$maps$phospho)

# reproduce the committed main classification before relying on the shared code path
sbP0 <- sets_in_bg(sets, bg$protein, "table", ortho); sbH0 <- sets_in_bg(sets, bg$phospho, "table", ortho)
cl_main <- classify_wide(ora(deP, bg$protein, sbP0), ora(deH, bg$phospho, sbH0))
old <- read.csv("analysis/C4_pathways/improved/c4_S_wide_classification_per_variant.csv", stringsAsFactors = FALSE)
old <- old[old$variant == "main", ]
stopifnot(identical(cl_main$class, old$class[match(cl_main$set, old$set)]))

# ---- list coverage in the reference ------------------------------------------
cover <- do.call(rbind, lapply(names(REF), function(r) data.frame(reference = r, n_ref = length(REF[[r]]),
  list = c("detected_proteome", "detected_phospho", "DE_proteome", "DE_phospho", "DE_phospho_unnorm", "DE_phospho_norm"),
  n_list = c(length(detP), length(detH), length(deP), length(deH), length(deHu), length(deHn)),
  n_in_ref = sapply(list(detP, detH, deP, deH, deHu, deHn), function(l) length(intersect(l, REF[[r]]))))))
write.csv(cover, "analysis/C4_pathways/improved/c4b_list_coverage.csv", row.names = FALSE); print(cover)

# ---- Q-G1: detection control -------------------------------------------------
tabs <- list(); g1 <- NULL
for (r in names(REF)) for (l in c("phospho", "proteome")) {
  det <- if (l == "phospho") detH else detP
  t <- ora(det, REF[[r]], SB[[r]]); tabs[[paste(r, l)]] <- tag(t, reference = r, list = paste0("detected_", l))
  g1 <- rbind(g1, data.frame(reference = r, list = paste0("detected_", l), n_list_in_ref = t$n_list[1], n_ref = length(REF[[r]]),
              n_sets_tested = nrow(t), n_S_wide_tested = sum(S_WIDE %in% t$set), n_S_wide_enriched = sum(S_WIDE %in% enriched(t)),
              n_all_enriched = length(enriched(t)), S_wide_enriched = paste(sub("KEGG_", "", intersect(S_WIDE, enriched(t))), collapse = ";")))
}
write.csv(do.call(rbind, tabs), "analysis/C4_pathways/improved/c4b_qg1_detected_enrichment_tables.csv", row.names = FALSE)
write.csv(g1, "analysis/C4_pathways/improved/c4b_qg1_counts.csv", row.names = FALSE); print(g1[, -9])
nH <- g1$n_S_wide_enriched[g1$reference == "U1" & g1$list == "detected_phospho"]
nP <- g1$n_S_wide_enriched[g1$reference == "U1" & g1$list == "detected_proteome"]
qg1 <- nH >= 4 && nP < nH / 2
cat(sprintf("\nQ-G1 (U1): phospho %d/9, proteome %d/9 -> detection bias supported: %s\n", nH, nP, qg1))

# ---- Q-G2: classification under main, V8 (U1), V9 (U2) -----------------------
cls <- list(main = cl_main); enr <- list()
for (v in c("V8_U1", "V9_U2")) {
  r <- sub(".*_", "", v)
  tP <- ora(deP, REF[[r]], SB[[r]]); tH <- ora(deH, REF[[r]], SB[[r]])
  enr[[v]] <- rbind(tag(tP, layer = "proteome", reference = r), tag(tH, layer = "phospho", reference = r))
  cls[[v]] <- classify_wide(tP, tH)
}
write.csv(do.call(rbind, enr), "analysis/C4_pathways/improved/c4b_qg2_enrichment_tables.csv", row.names = FALSE)
cl_long <- do.call(rbind, Map(function(c, v) tag(c, variant = v), cls, names(cls)))
# S_wide sets with < 5 genes in the reference are dropped by ora(); record them as NA
miss <- setdiff(S_WIDE, Reduce(intersect, lapply(enr, function(e) e$set[e$layer == "phospho"])))
cat("S_wide sets dropped (<5 genes in a reference):", if (length(miss)) miss else "none", "\n")
write.csv(cl_long, "analysis/C4_pathways/improved/c4b_qg2_S_wide_classification.csv", row.names = FALSE)
qg2 <- do.call(rbind, Map(function(c, v) tag(claim_summary(c), variant = v), cls, names(cls)))
write.csv(qg2, "analysis/C4_pathways/improved/c4b_qg2_summary.csv", row.names = FALSE); print(qg2)
# all-set counts per layer, for context
allc <- do.call(rbind, lapply(names(enr), function(v) { e <- enr[[v]]
  do.call(rbind, lapply(c("proteome", "phospho"), function(l) { x <- e[e$layer == l, ]
    data.frame(variant = v, layer = l, n_list = x$n_list[1], n_sets = nrow(x), n_enriched = length(enriched(x))) })) }))
write.csv(allc, "analysis/C4_pathways/improved/c4b_qg2_all_set_counts.csv", row.names = FALSE); print(allc)

# ---- Q-G3: normalisation at pathway level under V8 ---------------------------
tU <- ora(deHu, U1, SB$U1); tN <- ora(deHn, U1, SB$U1); stopifnot(identical(tU$set, tN$set))
eU <- enriched(tU); eN <- enriched(tN)
qg3 <- data.frame(reference = "U1", n_genes_unnormalised = tU$n_list[1], n_genes_normalised = tN$n_list[1], n_sets_tested = nrow(tU),
                  n_enriched_unnormalised = length(eU), n_enriched_normalised = length(eN), n_both = length(intersect(eU, eN)),
                  jaccard = jaccard(eU, eN), spearman_neglog10p = cor(-log10(tU$p), -log10(tN$p), method = "spearman"))
qg3$agree <- !is.na(qg3$jaccard) && qg3$jaccard >= 0.8
write.csv(qg3, "analysis/C4_pathways/improved/c4b_qg3_normalisation.csv", row.names = FALSE)
write.csv(data.frame(set = tU$set, size = tU$size, p_unnorm = tU$p, BH_unnorm = tU$BH, p_norm = tN$p, BH_norm = tN$BH),
          "analysis/C4_pathways/improved/c4b_qg3_per_set.csv", row.names = FALSE)
cat("\nQ-G3:\n"); print(t(qg3))

# ---- figures -----------------------------------------------------------------
cols <- c("phospho-only" = "#D55E00", "proteome-only" = "#0072B2", both = "#009E73", neither = "#E8E8E8")
lab <- c("phospho-only" = "P", "proteome-only" = "R", both = "B", neither = "")
vs <- names(cls)
save_fig("c4b_fig1_S_wide_classification", 6.5, 5.5, function() {
  par(mar = c(11, 17, 3, 1))
  plot(NA, xlim = c(0.5, length(vs) + 0.5), ylim = c(0.5, 9.5), axes = FALSE, xlab = "", ylab = "", yaxs = "i", xaxs = "i")
  for (i in seq_along(vs)) for (j in seq_along(S_WIDE)) {
    k <- cls[[vs[i]]]$class[cls[[vs[i]]]$set == S_WIDE[j]]; y <- 10 - j
    rect(i - .5, y - .5, i + .5, y + .5, col = cols[k], border = "white", lwd = 2)
    text(i, y, lab[k], cex = 0.9, font = 2, col = "white")
  }
  axis(1, at = seq_along(vs), labels = vs, las = 2, tick = FALSE, cex.axis = 0.8, line = -0.5)
  axis(2, at = 9:1, labels = sub("KEGG_", "", S_WIDE), las = 1, tick = FALSE, cex.axis = 0.8)
  title("S_wide classification: main vs V8 (U1) vs V9 (U2)", cex.main = 0.9)
  legend(x = length(vs) / 2 + 0.5, y = -4.5, xjust = 0.5, xpd = TRUE, ncol = 2, bty = "n", cex = 0.8, fill = cols,
         legend = c("phospho-only (P)", "proteome-only (R)", "both (B)", "neither"))
})
save_fig("c4b_fig2_detected_S_wide", 6.5, 5, function() {
  par(mar = c(5, 4.5, 3, 1))
  v <- g1$n_S_wide_enriched; names(v) <- paste0(sub("detected_", "", g1$list), "\n", g1$reference)
  b <- barplot(v, ylim = c(0, 9.8), col = ifelse(g1$list == "detected_phospho", "#D55E00", "#0072B2"),
               ylab = "S_wide sets enriched (BH < 0.05) of 9")
  text(b, v + 0.4, v); abline(h = 4, lty = 2)
  title("All detected genes vs reference: S_wide enrichment", cex.main = 0.9)
})

# ---- metrics -----------------------------------------------------------------
m <- NULL
for (i in seq_len(nrow(g1))) m <- rbind(m, mrow(paste0("partB_QG1_", g1$list[i], "_", g1$reference[i]), "S_wide_enriched", g1$n_S_wide_enriched[i], g1$n_S_wide_tested[i]))
m <- rbind(m, mrow("partB_QG1", "detection_bias_supported", as.integer(qg1), 9))
for (i in seq_len(nrow(qg2))) { v <- paste0("partB_QG2_", qg2$variant[i])
  m <- rbind(m, mrow(v, "S_wide_enriched_either", qg2$n_enriched_either[i], 9), mrow(v, "S_wide_phospho_only_share", qg2$phospho_only_share[i], qg2$n_enriched_either[i]),
             mrow(v, "claim_holds", as.integer(qg2$claim_holds[i]), 9)) }
m <- rbind(m, mrow("partB_QG3_V8_normalisation", "jaccard_enriched_sets", qg3$jaccard, qg3$n_sets_tested),
           mrow("partB_QG3_V8_normalisation", "spearman_neglog10p", qg3$spearman_neglog10p, qg3$n_sets_tested))
put_metrics(m)
writeLines(capture.output(sessionInfo()), "env/sessionInfo_C4b.txt")
