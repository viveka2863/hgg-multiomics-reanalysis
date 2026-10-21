# ======================================================================
# C5 part B (supplement round 1, Q-B): reconcile the repo's kinase tables with the paper's Supplementary Data 4e and Source Data Fig. 4a.
# Design and thresholds: LOG.md "PRE-REGISTRATION: supplement round 1". Run from the repo root AFTER
# analysis/C2_differential/improved/c2b_convert_xlsx.R:   Rscript analysis/C5_kinase_activity/c5b_reconcile_paper_tables.R
# ======================================================================
options(stringsAsFactors = FALSE)
set.seed(1)
stopifnot(file.exists("CLAUDE.md"))
OI <- "analysis/C5_kinase_activity/improved"
D  <- "JUMPn/JUMPnk"
rdx <- function(f) read.csv(file.path("data/external/csv", f), header = FALSE, colClasses = "character", na.strings = "")
wr <- function(df, f) write.csv(df, file.path(OI, f), row.names = FALSE)
msg <- function(...) cat(sprintf(...), "\n")
jac <- function(a, b) if (length(union(a, b)) == 0) NA_real_ else length(intersect(a, b)) / length(union(a, b))

# ---- paper tables ----------------------------------------------------------------------------------------
e <- rdx("t4e.csv"); stopifnot(e[1, 1] == "Kinases", e[1, 5] == "ANOVA P values", e[1, 6] == "BH FDR")
e <- e[-1, ]; e <- e[!is.na(e[, 1]), ]
T4e <- data.frame(Kinase = e[, 1], P_paper = as.numeric(e[, 5]), FDR_paper = as.numeric(e[, 6]))
stopifnot(nrow(T4e) == 90, !anyDuplicated(T4e$Kinase))
s4 <- rdx("src4a.csv"); stopifnot(s4[1, 1] == "Kinase")
paper41 <- s4[-1, 1]; stopifnot(length(paper41) == 41, !anyDuplicated(paper41))
msg("Sheet 4e: %d kinases; BH FDR < 0.05: %d; (BH recomputed from its P values: %d). Source 4a: %d kinases.", nrow(T4e), sum(T4e$FDR_paper < 0.05),
    sum(p.adjust(T4e$P_paper, "BH") < 0.05), length(paper41))

# ---- repo tables --------------------------------------------------------------------------------------------
pub <- read.delim(file.path(D, "KinaseActivityHeatmap/ALL_LIMMA_ANOVA_HGG_Result_20160708.txt"), check.names = FALSE)
all10 <- read.csv(file.path(OI, "c5_improved_limma_all10.csv"), check.names = FALSE)
ik <- read.delim(file.path(D, "KinaseActivityHeatmap/IKAP_result_v2_20160708.txt"), check.names = FALSE, row.names = 1)
ph <- read.delim(file.path(D, "PhosphoWholeNormalization/hgg_phosphosite_matrix_normalized_whole_v1_20160420.txt"), check.names = FALSE)
stopifnot(nrow(ph) == 984)

# ---- B1: sheet 4e vs shipped coded-list table and vs our all-10 table ---------------------------------------------
nm <- data.frame(quantity = c("n_4e", "n_shipped", "n_all10", "n_4e_in_shipped", "n_4e_in_all10", "n_4e_in_IKAP_rownames", "n_paper41_in_4e", "n_paper41_in_IKAP_rownames"),
                 value = c(nrow(T4e), nrow(pub), nrow(all10), sum(T4e$Kinase %in% rownames(pub)), sum(T4e$Kinase %in% all10$Kinase), sum(T4e$Kinase %in% rownames(ik)),
                           sum(paper41 %in% T4e$Kinase), sum(paper41 %in% rownames(ik))))
wr(nm, "c5b_b1_name_matching.csv"); print(nm, row.names = FALSE)
if (any(!T4e$Kinase %in% rownames(pub))) msg("4e kinases not in the shipped table: %s", paste(setdiff(T4e$Kinase, rownames(pub)), collapse = ", "))
b1 <- NULL
for (nmv in c("shipped_coded", "all10")) {
  R <- if (nmv == "shipped_coded") data.frame(Kinase = rownames(pub), P = pub$P.Value, adjP = pub$adj.P.Val) else data.frame(Kinase = all10$Kinase, P = all10$P.Value, adjP = all10$adj.P.Val)
  m <- merge(T4e, R, by = "Kinase")
  b1 <- rbind(b1, data.frame(table = nmv, n_common = nrow(m),
    spearman_P = cor(m$P_paper, m$P, method = "spearman"), spearman_FDR = cor(m$FDR_paper, m$adjP, method = "spearman"),
    median_abs_log10P_diff = median(abs(log10(m$P_paper) - log10(m$P))), max_abs_log10P_diff = max(abs(log10(m$P_paper) - log10(m$P))),
    median_abs_log10FDR_diff = median(abs(log10(m$FDR_paper) - log10(m$adjP))),
    n_BH_lt_0.05_paper_4e = sum(m$FDR_paper < 0.05), n_BH_lt_0.05_repo_table = sum(m$adjP < 0.05),
    n_both = sum(m$FDR_paper < 0.05 & m$adjP < 0.05), jaccard_BH_lt_0.05 = jac(m$Kinase[m$FDR_paper < 0.05], m$Kinase[m$adjP < 0.05])))
  wr(data.frame(m, log10P_diff = log10(m$P_paper) - log10(m$P)), sprintf("c5b_b1_pvalues_4e_vs_%s.csv", nmv))
}
wr(b1, "c5b_b1_4e_vs_repo_tables.csv"); print(b1, row.names = FALSE, digits = 4)

# ---- B2: substrate-count filter ---------------------------------------------------------------------------------------
sub <- lapply(seq_len(nrow(ph)), function(i) trimws(strsplit(ph$KINASE[i], "[,/]")[[1]]))
pairs <- unique(data.frame(site = rep(ph$Phosphosite, lengths(sub)), kinase = unlist(sub)))
pairs <- pairs[pairs$kinase != "" & !is.na(pairs$kinase), ]
nsub <- table(pairs$kinase)
T4e$n_substrates <- as.integer(nsub[T4e$Kinase]); T4e$n_substrates[is.na(T4e$n_substrates)] <- 0L
T4e$BH_lt_0.05 <- T4e$FDR_paper < 0.05
T4e$pass_ge3 <- T4e$BH_lt_0.05 & T4e$n_substrates >= 3
T4e$in_source_4a_41 <- T4e$Kinase %in% paper41
msg("normalised matrix: %d site-kinase pairs, %d distinct kinases; 4e kinases with no substrate row: %d", nrow(pairs), length(nsub), sum(T4e$n_substrates == 0))
wr(T4e, "c5b_b2_4e_kinases_with_substrate_counts.csv")
k_chain <- function(k) sum(T4e$BH_lt_0.05 & T4e$n_substrates >= k)
chain <- data.frame(stage = c("kinases in sheet 4e", "BH FDR < 0.05", "BH < 0.05 and >= 3 quantified substrates", "Source Data 4a (paper's Fig. 4a)", "paper's text: significant", "paper's text: >= 3 substrates"),
                    n = c(nrow(T4e), sum(T4e$BH_lt_0.05), sum(T4e$pass_ge3), length(paper41), 47, 41),
                    source = c("data", "data", "data", "data", "paper (as stated in the pre-registration)", "paper (as stated in the pre-registration)"))
wr(chain, "c5b_b2_filter_chain.csv"); print(chain, row.names = FALSE)
sweep_k <- data.frame(min_substrates = 0:10, n_BH_lt_0.05_and_ge_k = sapply(0:10, k_chain))
wr(sweep_k, "c5b_b2_substrate_threshold_sweep.csv"); print(t(sweep_k), quote = FALSE)

# ---- B3: filtered set vs the 41 kinases of Source Data 4a ---------------------------------------------------------------------
filt <- T4e$Kinase[T4e$pass_ge3]; bh <- T4e$Kinase[T4e$BH_lt_0.05]
mk <- function(lbl, ours) { tp <- length(intersect(ours, paper41))
  data.frame(set = lbl, n_ours = length(ours), n_paper = 41, tp = tp, precision = tp / length(ours), recall = tp / 41, jaccard = jac(ours, paper41), reconciled_jaccard_ge_0.8 = jac(ours, paper41) >= 0.8) }
b3 <- rbind(mk("4e BH<0.05 and >=3 substrates", filt), mk("4e BH<0.05 (no substrate filter)", bh),
            mk("shipped coded table adj.P<0.05", rownames(pub)[pub$adj.P.Val < 0.05]), mk("our all-10 table adj.P<0.05", all10$Kinase[all10$adj.P.Val < 0.05]))
wr(b3, "c5b_b3_vs_source_4a.csv"); print(b3, row.names = FALSE, digits = 3)
d <- data.frame(Kinase = union(filt, paper41), in_filtered_4e = union(filt, paper41) %in% filt, in_source_4a = union(filt, paper41) %in% paper41)
d <- d[!(d$in_filtered_4e & d$in_source_4a), ]
d$n_substrates <- as.integer(nsub[d$Kinase]); d$n_substrates[is.na(d$n_substrates)] <- 0L
d$FDR_4e <- T4e$FDR_paper[match(d$Kinase, T4e$Kinase)]
wr(d, "c5b_b3_disagreeing_kinases.csv"); msg("kinases in only one of {filtered 4e set, Source 4a}: %d", nrow(d)); print(d, row.names = FALSE)

# ---- figure -------------------------------------------------------------------------------------------------------------------------
save_fig2 <- function(name, w, h, fn) {
  png(file.path(OI, paste0(name, ".png")), w, h, units = "in", res = 150); fn(); invisible(dev.off())
}
save_fig2("c5b_fig1_filter_chain", 11, 5, function() {
  par(mfrow = c(1, 2), mar = c(9, 4.5, 3, 1))
  v <- chain$n[1:4]; names(v) <- c("kinases in 4e", "BH < 0.05", "+ >= 3 substrates", "Source Data 4a")
  bp <- barplot(v, col = c("grey80", "#4E79A7", "#4E79A7", "#F28E2B"), las = 2, ylab = "kinases", ylim = c(0, 100), main = "Filter chain (data)")
  text(bp, v, v, pos = 3)
  abline(h = 47, lty = 2, col = "red", lwd = 2); abline(h = 41, lty = 3, col = "red", lwd = 2)
  legend("topright", bty = "n", lty = c(2, 3), col = "red", lwd = 2, legend = c("paper: 47", "paper: 41"), cex = 0.9)
  plot(sweep_k$min_substrates, sweep_k$n_BH_lt_0.05_and_ge_k, type = "b", pch = 16, xlab = "minimum quantified substrates", ylab = "kinases with BH < 0.05 passing the filter",
       main = "Substrate-threshold sweep", ylim = c(0, max(sweep_k$n_BH_lt_0.05_and_ge_k, 50)))
  abline(h = 47, lty = 2, col = "red", lwd = 2); abline(h = 41, lty = 3, col = "red", lwd = 2); abline(v = 3, lty = 2, col = "grey40")
  points(3, k_chain(3), pch = 21, bg = "#F28E2B", cex = 2)
})

# ---- metrics ----------------------------------------------------------------------------------------------------------------------------
M <- NULL
add <- function(variant, metric, value, n) M <<- rbind(M, data.frame(claim = "C5", variant = variant, metric = metric, value = value, n = n))
for (i in seq_len(nrow(b1))) { v <- paste0("suppB_B1_4e_vs_", b1$table[i])
  add(v, "spearman_P", b1$spearman_P[i], b1$n_common[i]); add(v, "median_abs_log10P_diff", b1$median_abs_log10P_diff[i], b1$n_common[i])
  add(v, "n_BH_lt_0.05_repo", b1$n_BH_lt_0.05_repo_table[i], b1$n_common[i]); add(v, "jaccard_BH_lt_0.05", b1$jaccard_BH_lt_0.05[i], b1$n_common[i]) }
add("suppB_B2_4e", "n_BH_lt_0.05", sum(T4e$BH_lt_0.05), 90); add("suppB_B2_4e", "n_BH_lt_0.05_and_ge3_substrates", sum(T4e$pass_ge3), 90)
for (i in seq_len(nrow(b3))) { v <- paste0("suppB_B3_", c("filtered", "BH_only", "shipped", "all10")[i])
  for (mm in c("jaccard", "precision", "recall")) add(v, mm, b3[[mm]][i], b3$n_ours[i]) }
mf <- "results/metrics.csv"; old <- read.csv(mf, stringsAsFactors = FALSE)
old <- old[!(old$claim == "C5" & grepl("^suppB_", old$variant)), ]
write.csv(rbind(old, M[, names(old)]), mf, row.names = FALSE)
cat("done\n")
