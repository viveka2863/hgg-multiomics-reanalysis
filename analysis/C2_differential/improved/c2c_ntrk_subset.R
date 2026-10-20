# ======================================================================
# C2 part C (exploratory): which NTRK samples did the paper's supplement use?
# Design and thresholds: LOG.md "PRE-REGISTRATION: C2 part C" (Q-H1..Q-H4). Run from the repo root AFTER c2b_convert_xlsx.R:
#   Rscript analysis/C2_differential/improved/c2c_ntrk_subset.R
# ======================================================================
options(stringsAsFactors = FALSE)
set.seed(1)
source("analysis/C2_differential/improved/c2_lib.R")   # OUT, SAMPLES, GROUP, CUT, loaders, p-value helpers, jaccard
O <- "analysis/C2_differential/improved"
dir.create("env", showWarnings = FALSE)
for (f in c("sd1.csv", "sd2a.csv", "t4h.csv")) stopifnot(file.exists(file.path("data/external/csv", f)))
SD_ORDER <- c("Cortex1","Cortex2","Cortex3","PD1","PD2","PD3","PD4","TRK1","TRK2","TRK3")   # sheet order of the 10 intensity columns
rdx <- function(f) read.csv(file.path("data/external/csv", f), header = FALSE, colClasses = "character", na.strings = "")
num <- function(m) { m <- as.matrix(m); storage.mode(m) <- "double"; m }
wr <- function(df, f) write.csv(df, file.path(O, f), row.names = FALSE)
msg <- function(...) cat(sprintf(...), "\n")
EXACT <- 0.01; JCUT <- 0.8; SPCUT <- 0.95

# sample order was verified by values in c2b (committed mapping tables): re-assert, do not recompute
for (f in c("c2b_a_sample_mapping_protein.csv", "c2b_a_sample_mapping_phospho.csv")) {
  mp <- read.csv(file.path(O, f)); stopifnot(all(mp$as_expected), all(mp$expected_repo_column == rep(SD_ORDER, nrow(mp) / 10))) }   # protein file stacks two matrices

# ---- Supp Data 1 (proteins) and 2a (phosphosites) as in c2b ------------------------------------------------
x <- rdx("sd1.csv"); stopifnot(x[1, 1] == "Coexpresion Clusters", x[1, 19] == "ANOVA p Values", x[2, 16] == "NTRK/PDGFRA", x[2, 17] == "NTRK/Cortex", x[2, 18] == "PDGFRA/Cortex")
x <- x[-(1:2), ]; x <- x[!is.na(x[, 3]), ]
RC <- c("NTRK_PDGFRA", "NTRK_Cortex", "PDGFRA_Cortex")
P1 <- data.frame(id = x[, 3], cluster = x[, 1], p_paper = as.numeric(x[, 19])); R1 <- num(x[, 16:18]); colnames(R1) <- RC
I1 <- num(x[, 6:15]); colnames(I1) <- SD_ORDER; rownames(I1) <- rownames(R1) <- P1$id; stopifnot(!anyDuplicated(P1$id))
lab1 <- P1$id[grepl("^WP-C", P1$cluster)]; stopifnot(length(lab1) == 4703)
y <- rdx("sd2a.csv"); stopifnot(y[1, 1] == "Coexpresion Clusters", y[1, 20] == "ANOVA p Values", y[2, 17] == "NTRK/PDGFRA", y[2, 18] == "NTRK/Cortex", y[2, 19] == "PDGFRA/Cortex")
y <- y[-(1:2), ]; y <- y[!is.na(y[, 3]), ]
P2 <- data.frame(id = paste(y[, 2], sapply(strsplit(y[, 3], "[|]"), `[`, 2), y[, 4], sep = "_"), gene = y[, 2], site = y[, 4], cluster = y[, 1], p_paper = as.numeric(y[, 20]))
R2 <- num(y[, 17:19]); colnames(R2) <- RC
I2 <- num(y[, 7:16]); colnames(I2) <- SD_ORDER; rownames(I2) <- rownames(R2) <- P2$id; stopifnot(!anyDuplicated(P2$id), nrow(P2) == 28527)
lab2 <- P2$id[grepl("^PP-C", P2$cluster)]; stopifnot(length(lab2) == 6768)
L <- list(protein = list(I = I1, R = R1, P = P1, lab = lab1, fc = log2(1.5)), phospho = list(I = I2, R = R2, P = P2, lab = lab2, fc = 1))
for (ly in names(L)) msg("%s: %d rows, %d labelled; non-positive intensities: %d, NA: %d", ly, nrow(L[[ly]]$I), length(L[[ly]]$lab), sum(L[[ly]]$I <= 0, na.rm = TRUE), sum(is.na(L[[ly]]$I)))

# ---- NTRK subsets (Q-H1) ----------------------------------------------------------------------------------------
SUBS <- list("all3" = c("TRK1", "TRK2", "TRK3"), "TRK1+TRK2" = c("TRK1", "TRK2"), "TRK1+TRK3" = c("TRK1", "TRK3"), "TRK2+TRK3" = c("TRK2", "TRK3"),
             "TRK1" = "TRK1", "TRK2" = "TRK2", "TRK3" = "TRK3")
CTX <- c("Cortex1", "Cortex2", "Cortex3"); PDG <- c("PD1", "PD2", "PD3", "PD4")
ratios <- function(I, trk) {    # log2 of ratio of group means on the raw intensity scale
  m <- function(cols) { v <- rowMeans(I[, cols, drop = FALSE]); v[!is.finite(v) | v <= 0] <- NA; v }
  n <- m(trk); p <- m(PDG); c0 <- m(CTX)
  cbind(NTRK_PDGFRA = log2(n / p), NTRK_Cortex = log2(n / c0), PDGFRA_Cortex = log2(p / c0))
}

# identity of the ratio columns by values (against the all-3 recomputation): each paper column must correlate best with its own name
for (ly in names(L)) { Rc <- ratios(L[[ly]]$I, SUBS$all3)
  idc <- sapply(RC, function(pc) sapply(RC, function(oc) cor(L[[ly]]$R[, pc], Rc[, oc], use = "complete.obs"))); print(round(idc, 3))
  stopifnot(identical(RC[apply(idc, 2, which.max)], RC)) }

h1 <- NULL
for (ly in names(L)) for (sb in names(SUBS)) { Rc <- ratios(L[[ly]]$I, SUBS[[sb]])
  for (sset in c("all", "labelled")) { ids <- if (sset == "all") rownames(Rc) else L[[ly]]$lab
    for (ct in RC) { d <- Rc[ids, ct] - L[[ly]]$R[ids, ct]
      h1 <- rbind(h1, data.frame(layer = ly, subset_rows = sset, ratio_column = ct, ntrk_subset = sb, n = sum(!is.na(d)), median_abs_diff = median(abs(d), na.rm = TRUE),
                                 max_abs_diff = max(abs(d), na.rm = TRUE), spearman = cor(Rc[ids, ct], L[[ly]]$R[ids, ct], method = "spearman", use = "complete.obs"))) } } }
h1$exact_lt_0.01 <- h1$median_abs_diff < EXACT
best <- do.call(rbind, lapply(split(h1, list(h1$layer, h1$subset_rows, h1$ratio_column), drop = TRUE), function(d) { i <- which.min(d$median_abs_diff)
  data.frame(layer = d$layer[1], subset_rows = d$subset_rows[1], ratio_column = d$ratio_column[1], best_ntrk_subset = d$ntrk_subset[i], min_median_abs_diff = d$median_abs_diff[i],
             second_best = d$ntrk_subset[order(d$median_abs_diff)[2]], second_median_abs_diff = sort(d$median_abs_diff)[2], exact = d$median_abs_diff[i] < EXACT) }))
rownames(best) <- NULL
wr(h1, "c2c_h1_ratio_difference_by_subset.csv"); wr(best, "c2c_h1_best_subset.csv"); print(best, row.names = FALSE, digits = 3)

# ---- Q-H2 p-values ----------------------------------------------------------------------------------------------------
sp <- NULL
for (ly in names(L)) { lg <- log2(L[[ly]]$I); lg[!is.finite(lg)] <- NA; pp <- L[[ly]]$P$p_paper
  runs <- c(list(none_dropped = SD_ORDER), setNames(lapply(SD_ORDER, function(s) setdiff(SD_ORDER, s)), paste0("drop_", SD_ORDER)), list("NTRK=TRK1+TRK2" = setdiff(SD_ORDER, "TRK3")))
  for (nm in names(runs)) { cols <- runs[[nm]]; g <- GROUP[ifelse(cols == "PD1" | cols == "PD2" | cols == "PD3" | cols == "PD4", cols, cols)]
    g <- factor(c(Cortex = "Cortex", PD = "PDGFRA", TRK = "NTRK")[sub("[0-9]+$", "", cols)], levels = c("Cortex", "PDGFRA", "NTRK"))
    p <- anova_param_p(lg[, cols], g)
    sp <- rbind(sp, data.frame(layer = ly, run = nm, n_samples = length(cols), n = sum(is.finite(p) & !is.na(pp)), spearman_p = cor(p, pp, method = "spearman", use = "complete.obs"))) } }
sp$consistent_ge_0.95 <- NA
for (ly in names(L)) { i <- sp$layer == ly; sp$consistent_ge_0.95[i] <- sp$spearman_p[i] == max(sp$spearman_p[i]) & max(sp$spearman_p[i]) >= SPCUT }
wr(sp, "c2c_h2_pvalue_spearman.csv"); print(sp, row.names = FALSE, digits = 4)

# ---- Q-H3 DE lists ----------------------------------------------------------------------------------------------------
cmp_list <- function(ours, paper) { tp <- length(intersect(ours, paper)); data.frame(n_ours = length(ours), n_paper = length(paper), tp = tp,
  precision = if (length(ours)) tp / length(ours) else NA, recall = tp / length(paper), jaccard = jaccard(ours, paper)) }
maxpw <- function(gm) apply(gm, 1, function(v) max(v) - min(v))          # max pairwise |difference| of group-level log2 values
de_from <- function(I, cols, fc, fcmode) {
  lg <- log2(I[, cols]); lg[!is.finite(lg)] <- NA
  g <- factor(c(Cortex = "Cortex", PD = "PDGFRA", TRK = "NTRK")[sub("[0-9]+$", "", cols)], levels = c("Cortex", "PDGFRA", "NTRK"))
  p <- anova_param_p(lg, g)
  gm <- sapply(levels(g), function(l) { cc <- cols[g == l]; if (fcmode == "log2_of_mean") log2(rowMeans(I[, cc, drop = FALSE])) else rowMeans(lg[, cc, drop = FALSE]) })
  de_flag(p, maxpw(gm), 0.05, fc) & !is.na(p)
}
h3 <- NULL; H3L <- list()
for (ly in names(L)) for (fm in c("log2_of_mean", "mean_of_log2")) for (sb in c("all3", "TRK1+TRK2")) {
  cols <- setdiff(SD_ORDER, if (sb == "all3") character(0) else "TRK3")
  fl <- de_from(L[[ly]]$I, cols, L[[ly]]$fc, fm); ids <- rownames(L[[ly]]$I)[fl]; H3L[[paste(ly, fm, sb)]] <- ids
  h3 <- rbind(h3, cbind(layer = ly, arm = "sheet_intensities", fc_mode = fm, ntrk = sb, cmp_list(ids, L[[ly]]$lab))) }

# repo arm (proteins): committed C2 pipeline on the rescaled repo matrix, TRK3 excluded (and all-3 for reference)
raw <- load_protein(); fac <- read.csv(file.path(OUT, "c2_rescale_factors.csv")); fac <- setNames(fac$factor, fac$sample)[SAMPLES]
resc <- sweep(raw, 2, fac, "*")
cat("repo arm: p-values (permutation null: 1000 permutations, seed 1) ...\n")
for (sb in c("all3", "TRK1+TRK2")) { cols <- setdiff(SAMPLES, if (sb == "all3") character(0) else "TRK3"); x <- log2(resc[, cols]); g <- GROUP[cols]
  mfc <- max_pairwise(x, g)
  for (me in c("anova_perm", "limma_F")) { set.seed(1); p <- pvals(me, x, g); ids <- rownames(x)[de_flag(p, mfc, 0.05, CUT["protein"])]
    h3 <- rbind(h3, cbind(layer = "protein", arm = paste0("repo_rescaled_", me), fc_mode = "mean_of_log2", ntrk = sb, cmp_list(ids, lab1))) } }
h3$reproduced_jaccard_ge_0.8 <- h3$jaccard >= JCUT
wr(h3, "c2c_h3_de_lists_vs_labelled.csv"); print(h3, row.names = FALSE, digits = 3)

# ---- Q-H4 sheet 4h ----------------------------------------------------------------------------------------------------
t4 <- rdx("t4h.csv"); stopifnot(t4[1, 1] == "Phosphosites", grepl("^Fi PDGFRA", t4[1, 3]), grepl("^Fi NTRK", t4[1, 4])); t4 <- t4[-1, ]
Fi <- data.frame(entry = t4[, 1], Fi_PDGFRA = as.numeric(t4[, 3]), Fi_NTRK = as.numeric(t4[, 4]))
Fi$gene <- sub("[.].*", "", Fi$entry); Fi$sites <- strsplit(sub("^[^.]*[.]", "", Fi$entry), "/")
RA <- ratios(I2, SUBS$all3); RB <- ratios(I2, SUBS[["TRK1+TRK2"]])
h4 <- do.call(rbind, lapply(seq_len(nrow(Fi)), function(i) {
  rows <- sapply(Fi$sites[[i]], function(s) { m <- which(P2$gene == Fi$gene[i] & P2$site == s); if (length(m) == 1) m else NA_integer_ })
  rows <- unlist(rows); ok <- !is.na(rows); r <- rows[ok]
  data.frame(entry = Fi$entry[i], n_sites = length(rows), n_matched = sum(ok), Fi_PDGFRA = Fi$Fi_PDGFRA[i], Fi_NTRK = Fi$Fi_NTRK[i],
             sheet_NTRK_Cortex = mean(R2[r, "NTRK_Cortex"]), sheet_PDGFRA_Cortex = mean(R2[r, "PDGFRA_Cortex"]),
             recomp_all3_NTRK_Cortex = mean(RA[r, "NTRK_Cortex"]), recomp_TRK12_NTRK_Cortex = mean(RB[r, "NTRK_Cortex"]),
             recomp_PDGFRA_Cortex = mean(RA[r, "PDGFRA_Cortex"]))
}))
stopifnot(all(h4$n_matched == h4$n_sites))   # every constituent site must be found by gene + site (checked below; a stop here means the matching rule needs a decision)
h4$d_Fi_vs_sheet <- abs(h4$Fi_NTRK - h4$sheet_NTRK_Cortex); h4$d_Fi_vs_all3 <- abs(h4$Fi_NTRK - h4$recomp_all3_NTRK_Cortex); h4$d_Fi_vs_TRK12 <- abs(h4$Fi_NTRK - h4$recomp_TRK12_NTRK_Cortex)
h4$d_FiPD_vs_sheet <- abs(h4$Fi_PDGFRA - h4$sheet_PDGFRA_Cortex)
wr(h4, "c2c_h4_sheet4h_vs_SD2a.csv")
h4s <- do.call(rbind, lapply(c("all_entries", "single_site_entries"), function(sel) { d <- if (sel == "all_entries") h4 else h4[h4$n_sites == 1, ]
  data.frame(entries = sel, n = nrow(d), median_abs_diff_Fi_NTRK_vs_sheet_column = median(d$d_Fi_vs_sheet), median_abs_diff_vs_recomputed_all3 = median(d$d_Fi_vs_all3),
             median_abs_diff_vs_recomputed_TRK1TRK2 = median(d$d_Fi_vs_TRK12), median_abs_diff_Fi_PDGFRA_vs_sheet_column = median(d$d_FiPD_vs_sheet),
             spearman_Fi_vs_sheet = cor(d$Fi_NTRK, d$sheet_NTRK_Cortex, method = "spearman"), spearman_Fi_vs_all3 = cor(d$Fi_NTRK, d$recomp_all3_NTRK_Cortex, method = "spearman"),
             spearman_Fi_vs_TRK12 = cor(d$Fi_NTRK, d$recomp_TRK12_NTRK_Cortex, method = "spearman")) }))
wr(h4s, "c2c_h4_summary.csv"); print(t(h4s), digits = 4)

# ---- figures --------------------------------------------------------------------------------------------------------------
save_fig2 <- function(name, w, h, fn) { png(file.path(O, paste0(name, ".png")), w, h, units = "in", res = 150); fn(); invisible(dev.off())
  pdf(file.path(O, paste0(name, ".pdf")), w, h); fn(); invisible(dev.off()) }
FLOOR <- 1e-6   # exact zeros cannot be drawn on a log axis
cc <- c(NTRK_PDGFRA = "#4E79A7", NTRK_Cortex = "#F28E2B", PDGFRA_Cortex = "#59A14F")
save_fig2("c2c_fig1_ratio_difference_by_ntrk_subset", 11, 5, function() {
  par(mfrow = c(1, 2), mar = c(7, 4.5, 3, 1))
  for (ly in names(L)) { d <- h1[h1$layer == ly & h1$subset_rows == "all", ]
    plot(NA, xlim = c(0.5, length(SUBS) + 0.5), ylim = range(c(pmax(d$median_abs_diff, FLOOR), EXACT)) * c(0.7, 1.4), log = "y", xaxt = "n", xlab = "", ylab = "median |recomputed - sheet| log2 ratio (floored at 1e-6)",
         main = sprintf("%s (all %d rows)", ly, nrow(L[[ly]]$I)))
    axis(1, at = seq_along(SUBS), labels = paste("NTRK =", names(SUBS)), las = 2, cex.axis = 0.8)
    for (k in seq_along(RC)) { dd <- d[d$ratio_column == RC[k], ]; points(match(dd$ntrk_subset, names(SUBS)) + (k - 2) * 0.12, pmax(dd$median_abs_diff, FLOOR), pch = 16, col = cc[RC[k]], cex = 1.3) }
    abline(h = EXACT, lty = 2, lwd = 2); legend("topright", bty = "n", pch = 16, col = cc, legend = gsub("_", "/", RC), cex = 0.8) } })
j <- h3[h3$fc_mode == "log2_of_mean" | h3$arm != "sheet_intensities", ]
j$label <- ifelse(j$arm == "sheet_intensities", paste0("sheet intensities, NTRK=", j$ntrk), paste0(sub("repo_rescaled_", "repo ", j$arm), ", NTRK=", j$ntrk))
j$lay <- j$layer
save_fig2("c2c_fig2_jaccard_recomputed_vs_labelled", 10, 6, function() {
  par(mfrow = c(1, 2), mar = c(11, 4.5, 3, 1))
  for (ly in names(L)) { d <- j[j$lay == ly, ]
    barplot(d$jaccard, names.arg = d$label, las = 2, cex.names = 0.75, ylim = c(0, 1), col = ifelse(grepl("TRK1\\+TRK2", d$label), "#F28E2B", "#4E79A7"), ylab = "Jaccard vs labelled (DE) list",
            main = sprintf("%s (labelled: %d)", ly, length(L[[ly]]$lab)))
    abline(h = JCUT, lty = 2, lwd = 2) } })

# ---- verdict rule (pre-registered) ---------------------------------------------------------------------------------------
vr <- sapply(names(L), function(ly) {
  b <- best[best$layer == ly & best$subset_rows == "all" & best$ratio_column %in% c("NTRK_PDGFRA", "NTRK_Cortex"), ]
  h1ok <- all(b$best_ntrk_subset == "TRK1+TRK2" & b$exact)
  s <- h3[h3$layer == ly & h3$arm == "sheet_intensities" & h3$fc_mode == "log2_of_mean", ]
  h3ok <- s$jaccard[s$ntrk == "TRK1+TRK2"] >= JCUT && s$jaccard[s$ntrk == "all3"] < JCUT
  c(H1_TRK12_exact_for_NTRK_columns = h1ok, H3_J_TRK12_ge_0.8_and_J_all3_lt_0.8 = h3ok) })
print(vr); msg("VERDICT RULE (supplement computed without TRK3): %s", all(vr))

# ---- metrics.csv + sessionInfo --------------------------------------------------------------------------------------------
M <- NULL; add <- function(variant, metric, value, n) M <<- rbind(M, data.frame(claim = "C2", variant = variant, metric = metric, value = value, n = n))
for (i in which(h1$subset_rows == "all")) add(sprintf("suppC_H1_%s_%s_NTRK_%s", h1$layer[i], h1$ratio_column[i], h1$ntrk_subset[i]), "median_abs_log2ratio_diff", h1$median_abs_diff[i], h1$n[i])
for (i in seq_len(nrow(sp))) add(sprintf("suppC_H2_%s_%s", sp$layer[i], sp$run[i]), "spearman_p_vs_paper", sp$spearman_p[i], sp$n[i])
for (i in seq_len(nrow(h3))) for (mm in c("jaccard", "precision", "recall")) add(sprintf("suppC_H3_%s_%s_%s_NTRK_%s", h3$layer[i], h3$arm[i], h3$fc_mode[i], h3$ntrk[i]), mm, h3[[mm]][i], h3$n_ours[i])
for (i in seq_len(nrow(h4s))) for (mm in setdiff(names(h4s), c("entries", "n"))) add(sprintf("suppC_H4_%s", h4s$entries[i]), mm, h4s[[mm]][i], h4s$n[i])
mf <- "results/metrics.csv"; old <- read.csv(mf, stringsAsFactors = FALSE)
old <- old[!(old$claim == "C2" & grepl("^suppC_", old$variant)), ]
write.csv(rbind(old, M[, names(old)]), mf, row.names = FALSE)
writeLines(capture.output(sessionInfo()), "env/sessionInfo_C2c.txt")
cat("done: c2c outputs in", O, "\n")
