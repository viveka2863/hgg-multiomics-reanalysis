# Build and PROVE the sample key across all input files. Run from the repo root:
#   Rscript analysis/00_sample_key/make_sample_key.R
# Outputs (in this folder): sample_key.csv, mapping_evidence.csv, transgene_check.csv, input_vs_channel_spearman.csv
# Nothing under JUMPn/ or JUMPp/ is written. No RNG is used.

out <- "analysis/00_sample_key/"
J   <- "JUMPn/"
PW  <- paste0(J, "JUMPnk/PhosphoWholeNormalization/")
INT <- paste0(J, "JUMPnk/Integration/")
rd  <- function(f) read.delim(f, check.names = FALSE, quote = "", stringsAsFactors = FALSE)

inp  <- rd(paste0(J, "JUMPnc/HGG_whl_input.txt"))
v1   <- rd(paste0(PW, "hgg_whole_v1_20160420.txt"))
v5   <- rd(paste0(INT, "whl_hgg_v5_20160513.txt"))
p5   <- rd(paste0(INT, "pho_hgg_uniq_v5_20160513.txt"))
k5   <- rd(paste0(INT, "ikap_hgg_v5_20160513.txt"))
kres <- rd(paste0(J, "JUMPnk/KinaseActivityHeatmap/IKAP_result_v2_20160708.txt"))
raw  <- rd(paste0(PW, "id_all_site_quan_header_removed.txt"))
jraw <- rd(paste0(PW, "jpeng_hgg_phosphosite_matrix_v1_20160420.txt"))
jnrm <- rd(paste0(PW, "jpeng_phosphosite_matrix_normalized_v1_20160420.txt"))
nwh  <- rd(paste0(PW, "hgg_phosphosite_matrix_normalized_whole_v1_20160420.txt"))
wnrm <- rd(paste0(PW, "hgg_whole_normalized_v1_20160420.txt.txt"))

ch <- c("sig126","sig127N","sig127C","sig128N","sig128C","sig129N","sig129C","sig130N","sig130C","sig131")
stopifnot(identical(names(raw)[24:33], ch), identical(names(v1)[3:12], ch))

ev <- data.frame()  # tidy evidence log
add <- function(check, a, b, stat, value, note = "")
  ev <<- rbind(ev, data.frame(check = check, a = a, b = b, statistic = stat, value = value, note = note))

# ---- 1. whl_hgg_v5 columns vs hgg_whole_v1 channels: exact value match --------------------------
v5m <- as.matrix(v5[, 2:11]); v1m <- as.matrix(v1[, 3:12])
v5_channel <- character(10)
for (j in 1:10) {
  frac <- sapply(1:10, function(i) mean(v5m[, j] %in% v1m[, i]))
  v5_channel[j] <- ch[which.max(frac)]
  add("exact_v5_vs_whole_v1", names(v5)[1 + j], ch[which.max(frac)], "fraction_of_values_found", max(frac),
      paste("next best:", signif(sort(frac, TRUE)[2], 3)))
}
stopifnot(identical(v5_channel, ch))  # v5 column j == channel j, in order

# ---- 2. pho_hgg_uniq_v5 columns vs raw site file: exact value match -----------------------------
p5m <- as.matrix(p5[, 2:11]); rawm <- as.matrix(raw[, ch])
p5_channel <- character(10)
for (j in 1:10) {
  frac <- sapply(1:10, function(i) mean(p5m[, j] %in% rawm[, i]))
  p5_channel[j] <- ch[which.max(frac)]
  add("exact_pho_v5_vs_raw_sites", names(p5)[1 + j], ch[which.max(frac)], "fraction_of_values_found", max(frac),
      paste("next best:", signif(sort(frac, TRUE)[2], 3)))
}
stopifnot(identical(p5_channel, ch))
stopifnot(identical(names(p5)[2:11], names(v5)[2:11]))

# ---- 3. jpeng raw matrix Data1..10 == raw site sig126..131 (exact, row by row) ------------------
stopifnot(nrow(jraw) == nrow(raw))
for (j in 1:10)
  add("exact_jpeng_raw_matrix_vs_raw_sites", paste0("Data", j), ch[j], "all_rows_identical",
      as.numeric(isTRUE(all.equal(jraw[[paste0("Data", j)]], raw[[ch[j]]], tolerance = 0))))

# ---- 4. HGG_whl_input vs channels --------------------------------------------------------------
im <- as.matrix(inp[, 5:14]); colnames(im) <- names(inp)[5:14]
# 4a. exact attempt, row-aware (same accession): fraction of proteins whose value is identical
m <- match(inp$Protein_Accession, v1[["Protein Accession #"]]); ok <- !is.na(m)
exact <- t(sapply(1:10, function(i) sapply(1:10, function(j) mean(round(im[ok, i], 2) == round(v1m[m[ok], j], 2)))))
dimnames(exact) <- list(colnames(im), ch)
# Tumour columns have no identical values: test for a constant per-column scale factor instead.
# prop = fraction of proteins whose input/channel ratio is within 1e-4 (relative) of that pair's median ratio.
prop <- t(sapply(1:10, function(i) sapply(1:10, function(j) {
  r <- im[ok, i] / v1m[m[ok], j]; mean(abs(r / median(r, na.rm = TRUE) - 1) < 1e-4, na.rm = TRUE) })))
dimnames(prop) <- list(colnames(im), ch)
pr_best <- apply(prop, 1, which.max)
for (i in 1:10) {
  e <- sort(prop[i, ], TRUE); rr <- im[ok, i] / v1m[m[ok], pr_best[i]]
  add("proportional_input_vs_channel", colnames(im)[i], ch[pr_best[i]], "fraction_proteins_constant_ratio", e[1],
      sprintf("runner-up %s %.4f; median ratio %.4f; exact-identical fraction %.4f", names(e)[2], e[2], median(rr, na.rm = TRUE), exact[i, pr_best[i]]))
}
# 4b. Spearman fallback over proteins shared by accession (v1 carries accession, v5 does not;
#     v5 == v1 proven exact above, so v1 channel == v5 column)
add("input_proteins_matched_by_accession_to_v1", "HGG_whl_input", "hgg_whole_v1", "n_matched", sum(!is.na(m)),
    paste("unmatched:", paste(inp$Protein_Accession[is.na(m)], collapse = ";")))
S <- t(sapply(1:10, function(i) sapply(1:10, function(j) cor(im[ok, i], v1m[m[ok], j], method = "spearman"))))  # rows = input samples, cols = channels
dimnames(S) <- list(colnames(im), ch)
write.csv(round(S, 4), paste0(out, "input_vs_channel_spearman.csv"))
best <- apply(S, 1, which.max)
for (i in 1:10) {
  s <- sort(S[i, ], TRUE)
  add("spearman_input_vs_channel", colnames(im)[i], ch[best[i]], "spearman_rho_best", s[1],
      sprintf("runner-up %s rho=%.4f; margin=%.4f", names(s)[2], s[2], s[1] - s[2]))
}
print(round(prop, 3))
stopifnot(!any(duplicated(best)), identical(unname(best), unname(pr_best)))  # one-to-one; constant-ratio and Spearman agree

# 4c. transgene/vector proteins inside the input file itself (columns Cortex|PD|TRK)
for (a in grep("^cu\\|", inp$Protein_Accession, value = TRUE)) {
  x <- im[inp$Protein_Accession == a, ]; g <- sub("[0-9]+$", "", names(x))
  add("transgene_protein_input", a, "Cortex | PD | TRK mean", "mean_intensity", NA,
      paste(sprintf("%.0f", tapply(x, g, mean)[c("Cortex", "PD", "TRK")]), collapse = " | "))
}

# ---- 5. Transgene check on raw site file ---------------------------------------------------------
grp_data  <- c(rep("PDGFRA", 4), rep("NTRK", 3), rep("CNTRL", 3))        # layout implied by file headers
grp_paper <- c(rep("Cortex", 3), rep("PDGFRA", 4), rep("NTRK", 3))       # layout stated in the task (Fig 1c/d), channels 1-10
tg <- raw[grepl("cu\\|(PDGFRA_D842V|TMP3_NTRK1)\\|", raw[[2]]), ]
tgm <- as.matrix(tg[, ch])
gm <- function(x, g) tapply(x, g, mean)
rows <- lapply(seq_len(nrow(tg)), function(r) {
  a <- gm(tgm[r, ], grp_data); b <- gm(tgm[r, ], grp_paper)
  data.frame(site = tg[[2]][r], transgene = sub("cu\\|([^|]+)\\|.*", "\\1", tg[[2]][r]),
             top_channels = paste(ch[order(tgm[r, ], decreasing = TRUE)[1:3]], collapse = ","),
             PDGFRA_over_CNTRL_headerlayout = a["PDGFRA"] / a["CNTRL"],
             NTRK_over_CNTRL_headerlayout = a["NTRK"] / a["CNTRL"],
             PDGFRA_over_rest_headerlayout = a["PDGFRA"] / mean(tgm[r, 5:10]),
             NTRK_over_rest_headerlayout = a["NTRK"] / mean(tgm[r, c(1:4, 8:10)]),
             PDGFRA_over_rest_paperlayout = b["PDGFRA"] / mean(tgm[r, c(1:3, 8:10)]),
             NTRK_over_rest_paperlayout = b["NTRK"] / mean(tgm[r, 1:7]),
             row.names = NULL)
})
tgt <- do.call(rbind, rows)
tgt[, -(1:3)] <- round(tgt[, -(1:3)], 2)
write.csv(tgt, paste0(out, "transgene_check.csv"), row.names = FALSE)
# marker proteins in the whole-proteome file (GFP / mCherry vectors)
mk <- v1[grepl("^cu\\|(GFP|mCherry|PDGFRA_D842V|TMP3_NTRK1)\\|", v1[["Protein Accession #"]]), ]
for (r in seq_len(nrow(mk))) {
  x <- as.numeric(mk[r, ch]); g <- gm(x, grp_data)
  add("transgene_protein_whole_v1", mk[["Protein Accession #"]][r], "ch1-4 | ch5-7 | ch8-10 mean",
      "mean_intensity", NA, sprintf("%.0f | %.0f | %.0f", g["PDGFRA"], g["NTRK"], g["CNTRL"]))
}
# pooled call from the two kinase-domain sites that discriminate cleanly
pd <- tgm[grep("PDGFRA_D842V:Y849", tg[[2]]), ]; nt <- tgm[grep("TMP3_NTRK1:S501", tg[[2]]), ]
add("transgene_PDGFRA_Y849_site", "mean ch1-4", "mean ch5-10", "ratio", mean(pd[1:4]) / mean(pd[5:10]))
add("transgene_NTRK1_S501_site", "mean ch5-7", "mean other 7", "ratio", mean(nt[5:7]) / mean(nt[c(1:4, 8:10)]))

# ---- 6. Kinase activity files (ikap) -------------------------------------------------------------
# ikap_hgg_v5 has header PDGFRA1,PDGFRA2,PDGFRA2,PDGFRA4 -> map by exact values to IKAP_result
stopifnot(setequal(k5$Kinase, kres$Kinase))
kk <- kres[match(k5$Kinase, kres$Kinase), ]
for (j in 2:11) {
  hit <- sapply(2:11, function(i) isTRUE(all.equal(as.numeric(k5[[j]]), as.numeric(kk[[i]]), tolerance = 0)))
  add("exact_ikap_v5_vs_IKAP_result", sprintf("ikap_hgg_v5 pos %d ('%s')", j - 1, names(k5)[j]),
      names(kres)[1 + which(hit)], "all_90_kinases_identical", as.numeric(any(hit)),
      if (names(k5)[j] == "PDGFRA2" && j == 4) "duplicated header: this 3rd column is IKAP_result PDGFRA3" else "")
}
# Which Data column does each IKAP column track? Spearman across kinases of IKAP vs mean of
# normalised substrate-site log-ratios per kinase (supporting evidence only, not a proof).
nm <- jnrm[jnrm$KINASE != "NA" & !is.na(jnrm$KINASE), ]
dcols <- paste0("Data", 1:10)
km <- aggregate(nm[, dcols], by = list(Kinase = nm$KINASE), FUN = mean)
common <- intersect(km$Kinase, kres$Kinase)
add("ikap_support_kinases_with_substrate_data", "jpeng_normalized", "IKAP_result", "n_kinases", length(common))
for (i in 2:11) {
  r <- sapply(dcols, function(d) cor(kres[match(common, kres$Kinase), i], km[match(common, km$Kinase), d], method = "spearman"))
  add("ikap_support_spearman", names(kres)[i], dcols[which.max(r)], "spearman_rho_best", max(r),
      paste("runner-up", names(sort(r, TRUE))[2], signif(sort(r, TRUE)[2], 3)))
}

# ---- 7. Normalised matrices: Data1..10 order -----------------------------------------------------
# step2 normalises against columns '7,8,9' (0-based -> Data8-10). Controls should centre on ~0.
cm <- colMeans(as.matrix(jnrm[, dcols]), na.rm = TRUE)
for (j in 1:10) add("normalized_phospho_column_mean_log2", dcols[j], "", "mean", cm[j],
                    "controls (Data8-10) expected ~0 if step2 control columns are Data8-10")
cw <- colMeans(as.matrix(wnrm[, dcols]), na.rm = TRUE)
for (j in 1:10) add("normalized_whole_column_mean_log2", dcols[j], "", "mean", cw[j])
# phospho-normalised-by-whole vs phospho-normalised: same column should correlate best
cc <- match(nwh$Phosphosite, jnrm$Phosphosite); okc <- !is.na(cc)
for (j in 1:10) {
  r <- sapply(dcols, function(d) cor(nwh[[dcols[j]]][okc], jnrm[[d]][cc[okc]], use = "complete.obs", method = "spearman"))
  add("spearman_whole_normalized_vs_phospho_normalized", dcols[j], dcols[which.max(r)], "spearman_rho_best", max(r))
}
# which raw channel is each Data column? (jpeng raw matrix == raw site file, section 3)

write.csv(ev, paste0(out, "mapping_evidence.csv"), row.names = FALSE)

# ---- 8. Assemble the key ---------------------------------------------------------------------------
# Channel for each canonical sample comes from the DATA (section 4), group from the v5 header (section 1).
canon <- colnames(im)                         # Cortex1-3, PD1-4, TRK1-3
chan  <- best[canon]
idx   <- unname(as.integer(chan))
grp   <- c(Cortex = "CNTRL", PD = "PDGFRA", TRK = "NTRK")[sub("[0-9]+$", "", canon)]
stopifnot(identical(unname(grp), unname(grp_data[idx])))  # input-derived group agrees with v5-header group
ord <- order(idx)
key <- data.frame(
  canonical_id = canon, group = unname(grp), tmt_channel = ch[idx], channel_index = idx,
  HGG_whl_input = canon,
  hgg_whole_v1 = ch[idx],
  whl_hgg_v5 = names(v5)[1 + idx],
  pho_hgg_uniq_v5 = names(p5)[1 + idx],
  ikap_hgg_v5 = names(k5)[1 + idx], ikap_hgg_v5_position = idx,
  IKAP_result = NA_character_,
  normalized_matrix = paste0("Data", idx),
  raw_site_file = ch[idx], stringsAsFactors = FALSE)
# IKAP_result columns: CNTRL1-3, PDGFRA1-4, NTRK1-3 (group-blocked, different order). Map by exact values.
kcol <- sapply(seq_len(nrow(key)), function(r) {
  pos <- key$ikap_hgg_v5_position[r]
  hit <- sapply(2:11, function(i) isTRUE(all.equal(as.numeric(k5[[1 + pos]]), as.numeric(kk[[i]]), tolerance = 0)))
  names(kres)[1 + which(hit)]
})
key$IKAP_result <- kcol
key <- key[order(key$channel_index), ]
write.csv(key, paste0(out, "sample_key.csv"), row.names = FALSE)
print(key); print(round(S, 3)); print(tgt)
print(ev[ev$check %in% c("exact_input_vs_whole_v1_and_v5", "input_proteins_matched_by_accession_to_v1"), ])
dir.create("env", showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), "env/sessionInfo_00_sample_key.txt")
