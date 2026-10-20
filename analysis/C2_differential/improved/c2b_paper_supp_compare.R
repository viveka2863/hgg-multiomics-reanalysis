# ======================================================================
# C2 part B (supplement round 1, Q-A): our C2 lists / statistics vs the paper's Supplementary Data 1 (proteins) and 2a (phosphosites).
# Design and thresholds: LOG.md "PRE-REGISTRATION: supplement round 1". Run from the repo root AFTER c2b_convert_xlsx.R:
#   Rscript analysis/C2_differential/improved/c2b_paper_supp_compare.R
# ======================================================================
options(stringsAsFactors = FALSE)
set.seed(1)
source("analysis/C2_differential/improved/c2_lib.R")   # sets OUT, SAMPLES, GROUP, CUT, loaders, p-value helpers, jaccard
O <- "analysis/C2_differential/improved"
dir.create("env", showWarnings = FALSE)
for (f in c("sd1.csv", "sd2a.csv")) stopifnot(file.exists(file.path("data/external/csv", f)))
SD_ORDER <- c("Cortex1","Cortex2","Cortex3","PD1","PD2","PD3","PD4","TRK1","TRK2","TRK3")   # sheet order: Cortex 1-3, PDGFRA 4-7, NTRK 8-10 (mapping to repo names is verified below by values)
rdx <- function(f) read.csv(file.path("data/external/csv", f), header = FALSE, colClasses = "character", na.strings = "")
num <- function(m) { m <- as.matrix(m); storage.mode(m) <- "double"; m }
wr <- function(df, f) write.csv(df, file.path(O, f), row.names = FALSE)
msg <- function(...) cat(sprintf(...), "\n")

# ---- Supp Data 1 (proteins) ---------------------------------------------------------------------
x <- rdx("sd1.csv"); stopifnot(x[1, 1] == "Coexpresion Clusters", x[1, 19] == "ANOVA p Values", x[2, 16] == "NTRK/PDGFRA")
x <- x[-(1:2), ]; x <- x[!is.na(x[, 3]), ]
P1 <- data.frame(Protein_Accession = x[, 3], cluster = x[, 1], p_paper = as.numeric(x[, 19]),
                 r_NTRK_PDGFRA = as.numeric(x[, 16]), r_NTRK_Cortex = as.numeric(x[, 17]), r_PDGFRA_Cortex = as.numeric(x[, 18]))
I1 <- num(x[, 6:15]); rownames(I1) <- x[, 3]          # sample numbers 1-10: Cortex 1-3, PDGFRA 4-7, NTRK 8-10 (sheet order)
stopifnot(!anyDuplicated(P1$Protein_Accession))
msg("SD1: %d proteins, %d labelled (clusters: %s)", nrow(P1), sum(grepl("^WP-C", P1$cluster)), paste(names(table(P1$cluster)), table(P1$cluster), collapse = ", "))
paper_prot <- P1$Protein_Accession[grepl("^WP-C", P1$cluster)]   # unlabelled rows carry the text "NA" in the cluster column
stopifnot(length(paper_prot) == 4703)

raw <- load_protein(); fac <- read.csv(file.path(OUT, "c2_rescale_factors.csv")); fac <- setNames(fac$factor, fac$sample)[SAMPLES]
resc <- sweep(raw, 2, fac, "*")
bench <- rd("JUMPn/JUMPnc/HGG_whl_input.txt")
msg("repo benchmark: %d proteins; raw matrix: %d proteins", nrow(bench), nrow(raw))

# ---- A0: benchmark vs Supp Data 1 list ------------------------------------------------------------
a0 <- data.frame(quantity = c("n_repo_benchmark", "n_SD1_labelled", "n_intersect", "jaccard", "n_SD1_all", "n_repo_raw_matrix", "n_SD1_all_in_repo_raw"),
                 value = c(nrow(bench), length(paper_prot), length(intersect(bench$Protein_Accession, paper_prot)),
                           jaccard(bench$Protein_Accession, paper_prot), nrow(P1), nrow(raw), sum(P1$Protein_Accession %in% rownames(raw))))
cl <- merge(bench[, c("Protein_Accession", "Cluster")], P1[, c("Protein_Accession", "cluster")], by = "Protein_Accession")
a0 <- rbind(a0, data.frame(quantity = "cluster_label_agreement_on_intersect", value = mean(cl$Cluster == cl$cluster)))
wr(a0, "c2b_a0_benchmark_vs_SD1.csv"); print(a0, row.names = FALSE)

# ---- sample mapping and scaling checks (protein) ---------------------------------------------------
lb <- log2(as.matrix(bench[match(paper_prot, bench$Protein_Accession), SAMPLES]))
lb <- matrix(as.numeric(lb), nrow(lb), dimnames = list(paper_prot, SAMPLES))
# SD1 sample number j (1-10) vs repo columns: median |log2 ratio| on the labelled proteins
map_tab <- function(Ipaper, Ours, label) {
  M <- sapply(colnames(Ours), function(s) apply(Ipaper, 2, function(v) median(abs(log2(v) - Ours[, s]), na.rm = TRUE)))
  best <- colnames(M)[apply(M, 1, which.min)]
  data.frame(matrix = label, paper_sample_no = 1:10, best_repo_column = best, median_abs_log2_diff_best = apply(M, 1, min),
             expected_repo_column = SD_ORDER, as_expected = best == SD_ORDER)
}
mp1 <- map_tab(I1[paper_prot, ], lb, "SD1_vs_HGG_whl_input")
mp2 <- map_tab(I1[intersect(rownames(I1), rownames(raw)), ], log2(resc)[intersect(rownames(I1), rownames(raw)), ], "SD1_vs_rescaled_raw")
wr(rbind(mp1, mp2), "c2b_a_sample_mapping_protein.csv"); print(rbind(mp1, mp2), row.names = FALSE)
stopifnot(all(mp1$as_expected), all(mp2$as_expected))
# is SD1 raw or rescaled? per-sample median of SD1 / raw, labelled vs unlabelled, next to the C2 rescale factors
cm <- intersect(rownames(I1), rownames(raw)); lab <- cm %in% paper_prot
Iv <- I1[cm, ]; colnames(Iv) <- SD_ORDER; Iv <- Iv[, SAMPLES]   # reorder to repo order only after the mapping check passed
sc <- do.call(rbind, lapply(SAMPLES, function(s) data.frame(sample = s, c2_rescale_factor = unname(fac[s]),
        median_SD1_over_raw_labelled = median(Iv[lab, s] / raw[cm[lab], s]), median_SD1_over_raw_unlabelled = median(Iv[!lab, s] / raw[cm[!lab], s]),
        sd_log2_SD1_over_raw_unlabelled = sd(log2(Iv[!lab, s] / raw[cm[!lab], s])))))
wr(sc, "c2b_a_sd1_scaling_check.csv"); print(sc, row.names = FALSE, digits = 4)

# ---- statistics: recompute p-values on unscaled and rescaled matrices -------------------------------------
cat("computing p-values (permutation null: 1000 permutations, seed 1) ...\n")
pstats <- function(m) { x <- log2(m); g <- GROUP[colnames(x)]
  list(p = data.frame(anova_param = anova_param_p(x, g), anova_perm = anova_perm_p(x, g, 1000), limma_F = limma_F_p(x, g), row.names = rownames(x)),
       mfc = max_pairwise(x, g), mean_log2 = sapply(levels(g), function(l) rowMeans(x[, g == l, drop = FALSE]))) }
S <- list(unscaled = pstats(raw), rescaled = pstats(resc))
grid_lists <- function(st, rn) {     # same 18-variant grid as c2_main.R
  fcv <- c(none = NA_real_, "1.5" = log2(1.5), "2" = 1); out <- list()
  for (me in colnames(st$p)) for (fn in names(fcv)) for (fd in c(0.05, 0.01))
    out[[sprintf("%s_FC%s_FDR%s", me, fn, fd)]] <- rn[de_flag(st$p[[me]], st$mfc, fd, fcv[[fn]])]
  out
}
L_prot <- list(unscaled = grid_lists(S$unscaled, rownames(raw)), rescaled = grid_lists(S$rescaled, rownames(raw)))
# consistency with the committed C2 flags (unscaled)
fl <- read.csv(file.path(OUT, "c2_protein_DE_flags.csv.gz"), check.names = FALSE)
chk <- do.call(rbind, lapply(names(L_prot$unscaled), function(v) data.frame(variant = v, n_committed = sum(fl[[v]]), n_recomputed = length(L_prot$unscaled[[v]]),
          jaccard = jaccard(fl$Feature[fl[[v]]], L_prot$unscaled[[v]]))))
wr(chk, "c2b_check_recomputed_vs_committed_flags.csv"); msg("recomputed vs committed unscaled protein flags: min Jaccard %.4f", min(chk$jaccard, na.rm = TRUE))
qa <- read.csv(file.path(OUT, "c2_qa_summary.csv"))
msg("rescaled paper-style list: recomputed %d (committed c2_qa_summary: %d)", length(L_prot$rescaled$anova_perm_FC1.5_FDR0.05), qa$n_de[qa$variant == "rescaled"])

# ---- A1: protein lists vs SD1 list ------------------------------------------------------------------------------------
cmp_list <- function(ours, paper) { tp <- length(intersect(ours, paper))
  data.frame(n_ours = length(ours), n_paper = length(paper), tp = tp, precision = if (length(ours)) tp / length(ours) else NA,
             recall = tp / length(paper), jaccard = jaccard(ours, paper)) }
a1 <- do.call(rbind, lapply(names(L_prot), function(mx) do.call(rbind, lapply(names(L_prot[[mx]]), function(v)
  cbind(layer = "protein", matrix = mx, variant = v, cmp_list(L_prot[[mx]][[v]], paper_prot))))))
a1$reproduced_jaccard_ge_0.8 <- a1$jaccard >= 0.8
a1$paper_settings <- a1$variant == "anova_perm_FC1.5_FDR0.05"
wr(a1, "c2b_a1_protein_lists_vs_SD1.csv")
print(a1[a1$paper_settings | a1$variant %in% c("anova_param_FC1.5_FDR0.05", "limma_F_FC1.5_FDR0.05"), ], row.names = FALSE, digits = 3)
msg("A1 protein: best Jaccard %.3f (%s, %s); variants with Jaccard >= 0.8: %d of %d", max(a1$jaccard), a1$matrix[which.max(a1$jaccard)],
    a1$variant[which.max(a1$jaccard)], sum(a1$reproduced_jaccard_ge_0.8), nrow(a1))

# ---- A3 (protein): p-values and log2 ratios ------------------------------------------------------------------
rownames(P1) <- P1$Protein_Accession
sp_tab <- NULL; rat_tab <- NULL
rat_cols <- c(NTRK_PDGFRA = "r_NTRK_PDGFRA", NTRK_Cortex = "r_NTRK_Cortex", PDGFRA_Cortex = "r_PDGFRA_Cortex")
ours_ratio <- function(st) cbind(NTRK_PDGFRA = st$mean_log2[, "NTRK"] - st$mean_log2[, "PDGFRA"], NTRK_Cortex = st$mean_log2[, "NTRK"] - st$mean_log2[, "Cortex"],
                                 PDGFRA_Cortex = st$mean_log2[, "PDGFRA"] - st$mean_log2[, "Cortex"])
# which of our three contrasts does each paper column correspond to? (correlation table, unscaled)
R_u <- ours_ratio(S$unscaled)
idcheck <- do.call(rbind, lapply(names(rat_cols), function(pc) data.frame(paper_column = pc, t(sapply(colnames(R_u),
             function(oc) cor(P1[cm, rat_cols[[pc]]], R_u[cm, oc], use = "complete.obs"))), check.names = FALSE)))
idcheck$best_match <- colnames(R_u)[apply(idcheck[, -1], 1, which.max)]
wr(idcheck, "c2b_a3_ratio_column_identity_protein.csv"); print(idcheck, row.names = FALSE, digits = 3)
stopifnot(identical(idcheck$best_match, idcheck$paper_column))
# floor: ratios recomputed from the paper's own intensities (mean of log2) vs the paper's rounded ratios
grp_ratio <- function(L) { gm <- sapply(levels(GROUP), function(l) rowMeans(L[, names(GROUP)[GROUP == l], drop = FALSE]))
  cbind(NTRK_PDGFRA = gm[, "NTRK"] - gm[, "PDGFRA"], NTRK_Cortex = gm[, "NTRK"] - gm[, "Cortex"], PDGFRA_Cortex = gm[, "PDGFRA"] - gm[, "Cortex"]) }
Rp <- grp_ratio(log2(Iv))
grp_ratio_lm <- function(I) { gm <- sapply(levels(GROUP), function(l) log2(rowMeans(I[, names(GROUP)[GROUP == l], drop = FALSE])))
  cbind(NTRK_PDGFRA = gm[, "NTRK"] - gm[, "PDGFRA"], NTRK_Cortex = gm[, "NTRK"] - gm[, "Cortex"], PDGFRA_Cortex = gm[, "PDGFRA"] - gm[, "Cortex"]) }
Rpl <- grp_ratio_lm(Iv)
for (mx in c("unscaled", "rescaled", "paper_intensities", "paper_intensities_log2_of_mean")) for (sub in c("all", "labelled")) for (ct in names(rat_cols)) {
  R <- switch(mx, unscaled = ours_ratio(S$unscaled), rescaled = ours_ratio(S$rescaled), paper_intensities = Rp, paper_intensities_log2_of_mean = Rpl)
  ids <- if (sub == "all") cm else intersect(cm, paper_prot)
  d <- R[ids, ct] - P1[ids, rat_cols[[ct]]]
  rat_tab <- rbind(rat_tab, data.frame(layer = "protein", matrix = mx, subset = sub, contrast = ct, n = length(ids), median_abs_diff = median(abs(d), na.rm = TRUE),
                   mean_abs_diff = mean(abs(d), na.rm = TRUE), median_signed_diff = median(d, na.rm = TRUE),
                   spearman = cor(R[ids, ct], P1[ids, rat_cols[[ct]]], method = "spearman", use = "complete.obs")))
}
for (mx in c("unscaled", "rescaled")) for (me in colnames(S[[mx]]$p)) {
  sp_tab <- rbind(sp_tab, data.frame(layer = "protein", matrix = mx, method = me, n = length(cm),
                   spearman_p = cor(S[[mx]]$p[cm, me], P1[cm, "p_paper"], method = "spearman", use = "complete.obs"),
                   share_paper_p_le_0.005 = mean(P1[cm, "p_paper"] <= 0.005, na.rm = TRUE)))
}

# ---- Supp Data 2a (phosphosites) --------------------------------------------------------------------------------------
y <- rdx("sd2a.csv"); stopifnot(y[1, 1] == "Coexpresion Clusters", y[1, 20] == "ANOVA p Values", y[2, 17] == "NTRK/PDGFRA")
y <- y[-(1:2), ]; y <- y[!is.na(y[, 3]), ]
acc2 <- sapply(strsplit(y[, 3], "[|]"), `[`, 2)
P2 <- data.frame(key = paste(y[, 2], acc2, y[, 4], sep = "_"), cluster = y[, 1], p_paper = as.numeric(y[, 20]),
                 r_NTRK_PDGFRA = as.numeric(y[, 17]), r_NTRK_Cortex = as.numeric(y[, 18]), r_PDGFRA_Cortex = as.numeric(y[, 19]))
I2 <- num(y[, 7:16]); rownames(I2) <- P2$key
stopifnot(!anyDuplicated(P2$key), nrow(P2) == 28527)
paper_pho <- P2$key[grepl("^PP-C", P2$cluster)]; stopifnot(length(paper_pho) == 6768)
rownames(P2) <- P2$key
ph <- load_phospho(); pid <- rownames(ph); pkey <- sub("\\.[0-9]+$", "", pid)   # strip make.unique suffix added in load_phospho
first <- !duplicated(pkey)    # first occurrence in file order represents a key
msg("SD2a: %d sites, %d labelled; ours: %d rows, %d unique keys (%d rows share a key with an earlier row)", nrow(P2), length(paper_pho), nrow(ph), length(unique(pkey)), sum(!first))
mt <- data.frame(quantity = c("n_SD2a_sites", "n_SD2a_labelled", "n_ours_rows", "n_ours_unique_keys", "n_ours_rows_duplicate_key",
                              "SD2a_all_matched_in_ours", "SD2a_all_match_rate", "SD2a_labelled_matched_in_ours", "SD2a_labelled_match_rate",
                              "ours_keys_matched_in_SD2a", "ours_keys_match_rate", "ours_keys_unmatched", "SD2a_all_unmatched", "SD2a_labelled_unmatched"),
                 value = c(nrow(P2), length(paper_pho), nrow(ph), length(unique(pkey)), sum(!first),
                           sum(P2$key %in% pkey), mean(P2$key %in% pkey), sum(paper_pho %in% pkey), mean(paper_pho %in% pkey),
                           sum(unique(pkey) %in% P2$key), mean(unique(pkey) %in% P2$key), sum(!unique(pkey) %in% P2$key),
                           sum(!P2$key %in% pkey), sum(!paper_pho %in% pkey)))
wr(mt, "c2b_a2_phospho_match_rates.csv"); print(mt, row.names = FALSE)
um <- setdiff(P2$key, pkey)    # diagnostic only: do unmatched SD2a sites exist under another accession of the same gene?
gs <- function(k) sub("^([^_]+)_[^_]+_(.*)$", "\\1_\\2", k)
msg("diagnostic: of %d unmatched SD2a sites, %d have the same gene+site under another accession in our file", length(um), sum(gs(um) %in% gs(pkey)))
wr(data.frame(key_unmatched_SD2a = um, labelled = um %in% paper_pho), "c2b_a2_unmatched_SD2a_keys.csv")

# sample mapping and scaling (phospho), on matched keys (first occurrence)
cm2 <- P2$key[P2$key %in% pkey]; ridx <- match(cm2, pkey)
Rm <- log2(ph[ridx, ]); rownames(Rm) <- cm2
mp3 <- map_tab(I2[cm2, ], Rm, "SD2a_vs_pho_hgg_uniq_v5"); wr(mp3, "c2b_a_sample_mapping_phospho.csv"); print(mp3, row.names = FALSE)
stopifnot(all(mp3$as_expected))
Iv2 <- I2[cm2, ]; colnames(Iv2) <- SD_ORDER; Iv2 <- Iv2[, SAMPLES]
wr(data.frame(sample = SAMPLES, median_SD2a_over_raw = sapply(SAMPLES, function(s) median(Iv2[, s] / ph[ridx, s])),
              sd_log2_SD2a_over_raw = sapply(SAMPLES, function(s) sd(log2(Iv2[, s] / ph[ridx, s])))), "c2b_a_sd2a_scaling_check.csv")

# do the matched rows carry the same intensities in both files? (Cortex 1-3 are never rescaled, so they must be equal if it is the same row)
same <- rowSums(abs(log2(Iv2[, c("Cortex1", "Cortex2", "Cortex3")]) - log2(as.matrix(ph[ridx, c("Cortex1", "Cortex2", "Cortex3")]))) < 1e-6) == 3
dupk <- cm2 %in% pkey[duplicated(pkey)]
wr(data.frame(quantity = c("n_matched_keys", "n_identical_cortex_intensities", "n_different", "n_matched_keys_with_duplicate_rows_in_ours",
                           "n_different_among_duplicate_keys", "n_different_among_unique_keys"),
              value = c(length(cm2), sum(same), sum(!same), sum(dupk), sum(!same & dupk), sum(!same & !dupk))), "c2b_a2_phospho_intensity_identity.csv")
msg("phospho matched keys with identical Cortex intensities in both files: %d of %d (different: %d, of which %d are keys duplicated in our file)", sum(same), length(cm2), sum(!same), sum(!same & dupk))

# ---- A2: phospho lists vs SD2a list --------------------------------------------------------------------------------------
fp <- read.csv(file.path(OUT, "c2_phospho_DE_flags.csv.gz"), check.names = FALSE, colClasses = c(Feature = "character"))
stopifnot(identical(fp$Feature, pid))
a2 <- NULL
for (v in setdiff(names(fp), "Feature")) {
  ok <- unique(pkey[fp[[v]]])
  r <- cmp_list(ok, paper_pho)
  ur <- intersect(unique(pkey), P2$key)           # secondary: restrict both lists to sites present in both files
  r2 <- cmp_list(intersect(ok, ur), intersect(paper_pho, ur))
  a2 <- rbind(a2, cbind(layer = "phospho", matrix = "unscaled", variant = v, r, n_ours_matched_key = length(intersect(ok, P2$key)),
                        jaccard_restricted_to_common_universe = r2$jaccard, precision_restricted = r2$precision, recall_restricted = r2$recall,
                        n_paper_in_common_universe = r2$n_paper))
}
a2$reproduced_jaccard_ge_0.8 <- a2$jaccard >= 0.8; a2$paper_settings <- a2$variant == "anova_perm_FC2_FDR0.05"
wr(a2, "c2b_a2_phospho_lists_vs_SD2a.csv")
print(a2[a2$paper_settings | a2$variant %in% c("anova_param_FC2_FDR0.05", "limma_F_FC2_FDR0.05"), c("variant", "n_ours", "tp", "precision", "recall", "jaccard", "jaccard_restricted_to_common_universe")], row.names = FALSE, digits = 3)

# ---- A3 (phospho) ------------------------------------------------------------------------------------------------------------
cat("computing phospho p-values ...\n")
S$phospho <- pstats(ph)
R2 <- ours_ratio(S$phospho)
idcheck2 <- do.call(rbind, lapply(names(rat_cols), function(pc) data.frame(paper_column = pc, t(sapply(colnames(R2),
              function(oc) cor(P2[cm2, rat_cols[[pc]]], R2[ridx, oc], use = "complete.obs"))), check.names = FALSE)))
idcheck2$best_match <- colnames(R2)[apply(idcheck2[, -1], 1, which.max)]
wr(idcheck2, "c2b_a3_ratio_column_identity_phospho.csv"); print(idcheck2, row.names = FALSE, digits = 3)
stopifnot(identical(idcheck2$best_match, idcheck2$paper_column))
Rp2 <- grp_ratio(log2(Iv2)); Rpl2 <- grp_ratio_lm(Iv2)
for (mx in c("unscaled", "paper_intensities", "paper_intensities_log2_of_mean")) for (sub in c("all", "labelled", "identical_intensity")) for (ct in names(rat_cols)) {
  ids <- switch(sub, all = cm2, labelled = intersect(cm2, paper_pho), identical_intensity = cm2[same])
  ours <- switch(mx, unscaled = R2[match(ids, pkey), ct], paper_intensities = Rp2[ids, ct], paper_intensities_log2_of_mean = Rpl2[ids, ct])
  d <- ours - P2[ids, rat_cols[[ct]]]
  rat_tab <- rbind(rat_tab, data.frame(layer = "phospho", matrix = mx, subset = sub, contrast = ct, n = length(ids), median_abs_diff = median(abs(d), na.rm = TRUE),
                   mean_abs_diff = mean(abs(d), na.rm = TRUE), median_signed_diff = median(d, na.rm = TRUE),
                   spearman = cor(ours, P2[ids, rat_cols[[ct]]], method = "spearman", use = "complete.obs")))
}
for (me in colnames(S$phospho$p)) for (sub in c("all", "identical_intensity")) { sel <- if (sub == "all") rep(TRUE, length(cm2)) else same
  sp_tab <- rbind(sp_tab, data.frame(layer = "phospho", matrix = "unscaled", method = me, n = sum(sel),
    spearman_p = cor(S$phospho$p[ridx[sel], me], P2[cm2[sel], "p_paper"], method = "spearman", use = "complete.obs"), share_paper_p_le_0.005 = mean(P2[cm2[sel], "p_paper"] <= 0.005, na.rm = TRUE))) }
sp_tab$subset <- ifelse(sp_tab$layer == "phospho" & duplicated(sp_tab[, c("layer", "matrix", "method")], fromLast = TRUE), "all", ifelse(sp_tab$layer == "phospho", "identical_intensity", "all"))
wr(sp_tab, "c2b_a3_pvalue_spearman.csv"); wr(rat_tab, "c2b_a3_log2ratio_differences.csv")
print(sp_tab, row.names = FALSE, digits = 3); print(rat_tab[rat_tab$subset == "all", ], row.names = FALSE, digits = 3)

# decision for "paper used rescaled data" (protein, all matched proteins, three contrasts pooled as the mean of the per-contrast medians)
pool <- function(mx, sub = "all") mean(rat_tab$median_abs_diff[rat_tab$layer == "protein" & rat_tab$matrix == mx & rat_tab$subset == sub])
dec <- data.frame(subset = c("all", "labelled"), pooled_median_abs_diff_unscaled = sapply(c("all", "labelled"), function(s) pool("unscaled", s)),
                  pooled_median_abs_diff_rescaled = sapply(c("all", "labelled"), function(s) pool("rescaled", s)))
dec$ratio_unscaled_over_rescaled <- dec$pooled_median_abs_diff_unscaled / dec$pooled_median_abs_diff_rescaled
dec$rescaled_supported <- dec$ratio_unscaled_over_rescaled >= 2 | dec$pooled_median_abs_diff_rescaled < 0.02
wr(dec, "c2b_a3_rescaled_decision.csv"); print(dec, row.names = FALSE, digits = 4)

# ---- figures ---------------------------------------------------------------------------------------------------------------------------
save_fig2 <- function(name, w, h, fn) {
  png(file.path(O, paste0(name, ".png")), w, h, units = "in", res = 150); fn(); invisible(dev.off())
  pdf(file.path(O, paste0(name, ".pdf")), w, h); fn(); invisible(dev.off())
}
save_fig2("c2b_fig1_log2ratios_vs_paper", 11, 11, function() {
  par(mfrow = c(3, 3), mar = c(4, 4, 3, 1))
  panels <- list(list("protein", "unscaled", function(ct) cbind(R_u[cm, ct], P1[cm, rat_cols[[ct]]])),
                 list("protein", "rescaled", function(ct) cbind(ours_ratio(S$rescaled)[cm, ct], P1[cm, rat_cols[[ct]]])),
                 list("phospho", "unscaled", function(ct) cbind(R2[ridx, ct], P2[cm2, rat_cols[[ct]]])))
  for (pn in panels) for (ct in names(rat_cols)) {
    xy <- pn[[3]](ct); lim <- range(xy, na.rm = TRUE)
    plot(xy[, 2], xy[, 1], pch = 16, cex = 0.25, col = rgb(0.2, 0.3, 0.6, 0.25), xlim = lim, ylim = lim, xlab = "paper log2 ratio (Supp Data)", ylab = "our log2 ratio",
         main = sprintf("%s, %s: %s", pn[[1]], pn[[2]], gsub("_", "/", ct)), cex.main = 0.9)
    abline(0, 1, col = "red")
    rr <- rat_tab[rat_tab$layer == pn[[1]] & rat_tab$matrix == pn[[2]] & rat_tab$subset == "all" & rat_tab$contrast == ct, ]
    legend("topleft", bty = "n", cex = 0.8, legend = sprintf("median |diff| %.3f\nn = %d", rr$median_abs_diff, rr$n))
  }
})
save_fig2("c2b_fig2_jaccard_vs_paper", 13, 11, function() {
  par(mfrow = c(3, 1), mar = c(11, 4.5, 3, 1))
  pcol <- c(anova_param = "#4E79A7", anova_perm = "#F28E2B", limma_F = "#59A14F")
  for (pn in list(list(a1[a1$matrix == "unscaled", ], "Proteins, unscaled (vs 4,703 in Supp Data 1)"),
                  list(a1[a1$matrix == "rescaled", ], "Proteins, rescaled with the C2 factors (vs 4,703)"),
                  list(a2, "Phosphosites (vs 6,768 in Supp Data 2a; full lists)"))) {
    d <- pn[[1]]; me <- sub("_FC.*", "", d$variant)
    bp <- barplot(d$jaccard, names.arg = d$variant, las = 2, cex.names = 0.7, ylim = c(0, 1), col = pcol[me], border = ifelse(d$paper_settings, "black", NA),
                  ylab = "Jaccard vs paper list", main = pn[[2]])
    abline(h = 0.8, lty = 2, lwd = 2)
    legend("topright", bty = "n", fill = pcol, legend = names(pcol), cex = 0.8, horiz = TRUE)
  }
})

# ---- metrics.csv + sessionInfo -----------------------------------------------------------------------------------------------------------
pv <- a1[a1$paper_settings, ]
M <- data.frame(claim = "C2", variant = "suppB_A0_benchmark_vs_SD1", metric = "jaccard", value = a0$value[a0$quantity == "jaccard"], n = length(paper_prot))
add <- function(variant, metric, value, n) M <<- rbind(M, data.frame(claim = "C2", variant = variant, metric = metric, value = value, n = n))
for (i in seq_len(nrow(pv))) for (mm in c("jaccard", "precision", "recall")) add(paste0("suppB_A1_protein_", pv$matrix[i], "_paper_settings"), mm, pv[[mm]][i], pv$n_ours[i])
add("suppB_A1_protein_best_variant", "jaccard", max(a1$jaccard), length(paper_prot))
add("suppB_A1_protein_n_variants_jaccard_ge_0.8", "count", sum(a1$jaccard >= 0.8), nrow(a1))
pp <- a2[a2$paper_settings, ]
for (mm in c("jaccard", "precision", "recall", "jaccard_restricted_to_common_universe")) add("suppB_A2_phospho_paper_settings", mm, pp[[mm]], pp$n_ours)
add("suppB_A2_phospho_best_variant", "jaccard", max(a2$jaccard), length(paper_pho))
add("suppB_A2_phospho_n_variants_jaccard_ge_0.8", "count", sum(a2$jaccard >= 0.8), nrow(a2))
add("suppB_A2_phospho_match", "SD2a_labelled_match_rate", mt$value[mt$quantity == "SD2a_labelled_match_rate"], length(paper_pho))
add("suppB_A2_phospho_match", "SD2a_all_match_rate", mt$value[mt$quantity == "SD2a_all_match_rate"], nrow(P2))
for (i in seq_len(nrow(sp_tab))) add(sprintf("suppB_A3_%s_%s_%s%s", sp_tab$layer[i], sp_tab$matrix[i], sp_tab$method[i], ifelse(sp_tab$subset[i] == "all", "", "_identical_intensity")), "spearman_p_vs_paper", sp_tab$spearman_p[i], sp_tab$n[i])
ra_all <- rat_tab[rat_tab$subset == "all" & !grepl("^paper_intensities", rat_tab$matrix), ]
for (i in seq_len(nrow(ra_all))) add(sprintf("suppB_A3_%s_%s_%s", ra_all$layer[i], ra_all$matrix[i], ra_all$contrast[i]), "median_abs_log2ratio_diff", ra_all$median_abs_diff[i], ra_all$n[i])
for (i in seq_len(nrow(dec))) { add(sprintf("suppB_A3_rescaled_decision_%s", dec$subset[i]), "ratio_unscaled_over_rescaled", dec$ratio_unscaled_over_rescaled[i], nrow(raw)) }
mf <- "results/metrics.csv"; old <- read.csv(mf, stringsAsFactors = FALSE)
old <- old[!(old$claim == "C2" & grepl("^suppB_", old$variant)), ]
write.csv(rbind(old, M[, names(old)]), mf, row.names = FALSE)
writeLines(capture.output(sessionInfo()), "env/sessionInfo_suppB.txt")
cat("done: c2b outputs in", O, "\n")
