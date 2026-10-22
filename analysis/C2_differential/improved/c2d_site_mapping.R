# ======================================================================
# C2 part D (group-review round 1, Q-M): can repository phospho rows be carried across to the 28,527 Supp Data 2a sites?
# Design and thresholds: LOG.md "PRE-REGISTRATION: group-review amendments, round 1" (M1-M3). Run from the repo root AFTER c2b_convert_xlsx.R:
#   Rscript analysis/C2_differential/improved/c2d_site_mapping.R
# Keys tried (SD2a site -> repository rows, single-site repository rows unless stated):
#   exact = Gene_ACCESSION_SITE as written; acc_site = ACCESSION+SITE; accnoiso_site = ACCESSION without isoform suffix (-n) + SITE;
#   gene_site = GENE+SITE; sorted_multisite = GENE_ACCESSION_ + sorted site list (matters only for multi-site repository rows).
# Difference (M3): max over the 10 samples of |log2(SD2a) - log2(repo x per-sample phospho rescale factor)|; factors = median SD2a/raw per sample
#   from c2b_a_sd2a_scaling_check.csv (the protein factors of c2_rescale_factors.csv are the wrong scale for phospho; shown as a sensitivity column).
#   Where several repository rows share a key, "safe" uses the largest difference over those rows (all candidates must agree).
# ======================================================================
options(stringsAsFactors = FALSE)
set.seed(1)
source("analysis/C2_differential/improved/c2_lib.R")   # rd, load_phospho, SAMPLES, PHOS_COLS
O <- "analysis/C2_differential/improved"
SD_ORDER <- c("Cortex1","Cortex2","Cortex3","PD1","PD2","PD3","PD4","TRK1","TRK2","TRK3")
msg <- function(...) cat(sprintf(...), "\n")
rdx <- function(f) read.csv(file.path("data/external/csv", f), header = FALSE, colClasses = "character", na.strings = "")
num <- function(m) { m <- as.matrix(m); storage.mode(m) <- "double"; m }
wr <- function(df, f) write.csv(df, file.path(O, f), row.names = FALSE)
TOL <- 0.01

# ---- repository phospho rows (contaminants removed, as load_phospho) -------------------------------------------------------
ph <- rd("JUMPn/JUMPnk/Integration/pho_hgg_uniq_v5_20160513.txt"); ph <- ph[!grepl("CON_", ph[[1]]), ]
M  <- as.matrix(ph[, PHOS_COLS]); storage.mode(M) <- "double"; colnames(M) <- SAMPLES
stopifnot(all(M > 0), identical(unname(M), unname(load_phospho())))
rid <- ph[[1]]; tk <- strsplit(rid, "_", fixed = TRUE)
r <- data.frame(id = rid, gene = sapply(tk, `[`, 1), sites = sapply(tk, function(v) v[length(v)]),
                acc = sapply(tk, function(v) paste(v[-c(1, length(v))], collapse = "_")), stringsAsFactors = FALSE)
r$n_sites <- lengths(strsplit(r$sites, ",", fixed = TRUE)); r$multi <- r$n_sites > 1
r$accnoiso <- sub("-[0-9]+$", "", r$acc)
srt <- function(s) vapply(strsplit(s, ",", fixed = TRUE), function(v) paste(v[order(as.integer(sub("^[A-Z]", "", v)), v)], collapse = ","), "")
r$sorted_sites <- srt(r$sites)
msg("repository rows %d; multi-site rows %d; ids with 3 / 4+ underscore parts: %d / %d", nrow(r), sum(r$multi), sum(lengths(tk) == 3), sum(lengths(tk) > 3))

# ---- SD2a ------------------------------------------------------------------------------------------------------------------------
y <- rdx("sd2a.csv"); stopifnot(y[1, 1] == "Coexpresion Clusters" || y[1, 1] == "Coexpresion Clusters ", y[2, 17] == "NTRK/PDGFRA")
y <- y[-(1:2), ]; y <- y[!is.na(y[, 3]), ]; stopifnot(nrow(y) == 28527)
s <- data.frame(gene = y[, 2], acc = sapply(strsplit(y[, 3], "[|]"), `[`, 2), site = y[, 4], labelled = grepl("^PP-C", y[, 1]), stringsAsFactors = FALSE)
s$accnoiso <- sub("-[0-9]+$", "", s$acc); s$site_id <- paste(s$gene, s$acc, s$site, sep = "_")
stopifnot(!anyDuplicated(s$site_id), sum(s$labelled) == 6768, !any(grepl(",", s$site)))
I2 <- num(y[, 7:16]); colnames(I2) <- SD_ORDER; I2 <- I2[, SAMPLES]; stopifnot(all(I2 > 0))

# ---- candidate rows under each key (single-site repository rows for the site keys) ------------------------------------------------
sing <- which(!r$multi)
cand <- function(rkey, skey, rows) { sp <- split(rows, rkey[rows]); lapply(skey, function(k) if (k %in% names(sp)) sp[[k]] else integer(0)) }
C <- list(
  exact            = cand(r$id, s$site_id, seq_len(nrow(r))),
  acc_site         = cand(paste(r$acc, r$sites), paste(s$acc, s$site), sing),
  accnoiso_site    = cand(paste(r$accnoiso, r$sites), paste(s$accnoiso, s$site), sing),
  gene_site        = cand(paste(r$gene, r$sites), paste(s$gene, s$site), sing),
  sorted_multisite = cand(paste(r$gene, r$acc, r$sorted_sites, sep = "_"), s$site_id, seq_len(nrow(r))))
nC <- sapply(C, lengths)
# diagnostic only (not a recoverable key): SD2a site contained in a multi-site repository row of the same gene and accession
mr <- which(r$multi); msp <- split(mr, paste(r$gene[mr], r$acc[mr], sep = "_"))
contained <- sapply(seq_len(nrow(s)), function(i) { v <- msp[[paste(s$gene[i], s$acc[i], sep = "_")]]; !is.null(v) && any(sapply(strsplit(r$sites[v], ",", fixed = TRUE), function(z) s$site[i] %in% z)) })

# ---- M1: status and non-match classes ----------------------------------------------------------------------------------------------
KEYS <- colnames(nC)
first_key <- apply(nC > 0, 1, function(b) if (any(b)) KEYS[which(b)[1]] else NA_character_)
status <- ifelse(is.na(first_key), NA, ifelse(first_key == "exact", "exact", "recoverable"))
acc_in_repo <- s$acc %in% r$acc | s$accnoiso %in% r$accnoiso
nm <- ifelse(!is.na(first_key) & first_key == "exact", "exact",
      ifelse(nC[, "gene_site"] > 0 & nC[, "acc_site"] == 0 & nC[, "accnoiso_site"] == 0, "gene_site_matches_accession_differs",
      ifelse(nC[, "acc_site"] > 0 | nC[, "accnoiso_site"] > 0, "accession_site_matches_gene_label_differs",
      ifelse(acc_in_repo, "accession_present_site_differs", "absent"))))
status[is.na(status)] <- ifelse(acc_in_repo[is.na(status)], "not_recoverable_accession_present", "not_recoverable_absent")
msg("exact %d; recoverable %d; not recoverable %d (of 28,527)", sum(status == "exact"), sum(status == "recoverable"), sum(grepl("^not", status)))

# ---- M2: duplicated repository ids ---------------------------------------------------------------------------------------------------
dup_ids <- unique(r$id[duplicated(r$id)])
dg <- lapply(dup_ids, function(d) M[r$id == d, , drop = FALSE])
ident <- sapply(dg, function(m) all(apply(m, 1, function(v) all(v == m[1, ]))))
in_sd2a <- dup_ids %in% s$site_id
DUP <- data.frame(quantity = c("n_repo_rows", "n_distinct_ids", "n_rows_sharing_an_id", "n_duplicated_ids", "n_dup_ids_identical_values", "n_dup_ids_distinct_values",
                               "n_dup_ids_in_SD2a", "n_dup_ids_in_SD2a_identical", "n_dup_ids_in_SD2a_distinct", "n_SD2a_sites_with_exact_dup_candidates",
                               "n_multisite_repo_rows", "n_SD2a_sites_contained_in_multisite_repo_row_diagnostic"),
                  value = c(nrow(r), length(unique(r$id)), sum(r$id %in% dup_ids), length(dup_ids), sum(ident), sum(!ident), sum(in_sd2a), sum(in_sd2a & ident), sum(in_sd2a & !ident),
                            sum(nC[, "exact"] > 1), sum(r$multi), sum(contained)))
wr(DUP, "c2d_duplicate_summary.csv"); print(DUP, row.names = FALSE)

# ---- M3: differences after rescaling -------------------------------------------------------------------------------------------------------
fz <- read.csv(file.path(O, "c2b_a_sd2a_scaling_check.csv")); fp <- setNames(fz$median_SD2a_over_raw, fz$sample)[SAMPLES]
fq <- read.csv("analysis/C2_differential/c2_rescale_factors.csv"); fq <- setNames(fq$factor, fq$sample)[SAMPLES]
Lp <- log2(I2)
d_rows <- function(rows, i, f) sapply(rows, function(k) max(abs(Lp[i, ] - (log2(M[k, ]) + log2(f)))))
sel <- first_key; mind <- maxd <- maxd_unsc <- maxd_prot <- rep(NA_real_, nrow(s)); ncand <- integer(nrow(s)); rep_ids <- rep(NA_character_, nrow(s))
for (i in which(!is.na(sel))) {
  rows <- C[[sel[i]]][[i]]; ncand[i] <- length(rows)
  d <- d_rows(rows, i, fp); mind[i] <- min(d); maxd[i] <- max(d)
  maxd_unsc[i] <- max(d_rows(rows, i, rep(1, 10))); maxd_prot[i] <- max(d_rows(rows, i, fq))
  rep_ids[i] <- paste(head(r$id[rows], 3), collapse = ";")
}
safe <- !is.na(maxd) & maxd < TOL
ver <- data.frame(site_id = s$site_id, gene = s$gene, accession = s$acc, site = s$site, labelled_DE = s$labelled, status = status, first_key = first_key,
                  nonmatch_class = nm, n_candidate_rows = ncand, n_candidate_accessions = sapply(seq_len(nrow(s)), function(i) if (ncand[i]) length(unique(r$acc[C[[first_key[i]]][[i]]])) else 0L),
                  repo_ids_first3 = rep_ids, key_exact = nC[, "exact"] > 0, key_acc_site = nC[, "acc_site"] > 0, key_accnoiso_site = nC[, "accnoiso_site"] > 0,
                  key_gene_site = nC[, "gene_site"] > 0, key_sorted_multisite = nC[, "sorted_multisite"] > 0, contained_in_multisite_row = contained,
                  max_abs_log2_diff_min = mind, max_abs_log2_diff_max = maxd, max_abs_log2_diff_unscaled = maxd_unsc, max_abs_log2_diff_protein_factors = maxd_prot,
                  safe_to_carry = safe)
write_gz(ver, file.path(O, "c2d_site_mapping_verified.csv.gz"))

# ---- summaries ---------------------------------------------------------------------------------------------------------------------------------
both <- function(f) rbind(data.frame(subset = "all_28527", f(rep(TRUE, nrow(s)))), data.frame(subset = "labelled_DE_6768", f(s$labelled)))
REC <- do.call(rbind, lapply(KEYS, function(k) both(function(m) data.frame(key = k, n_sd2a_sites_with_key_match = sum(nC[m, k] > 0),
          n_additional_to_exact = sum(nC[m, k] > 0 & nC[m, "exact"] == 0), n_with_more_than_one_candidate_row = sum(nC[m, k] > 1)))))
REC <- rbind(REC, both(function(m) data.frame(key = "any_key", n_sd2a_sites_with_key_match = sum(!is.na(first_key[m])), n_additional_to_exact = sum(!is.na(first_key[m]) & nC[m, "exact"] == 0), n_with_more_than_one_candidate_row = NA)))
wr(REC, "c2d_recoverable_counts.csv")
NMC <- do.call(rbind, lapply(sort(unique(nm)), function(k) data.frame(nonmatch_class = k, n_all = sum(nm == k), n_labelled = sum(nm == k & s$labelled))))
NMC$recoverable_by <- c(accession_present_site_differs = "none", absent = "none", exact = "(already exact)",
                        accession_site_matches_gene_label_differs = "acc_site and/or accnoiso_site", gene_site_matches_accession_differs = "gene_site")[NMC$nonmatch_class]
wr(NMC, "c2d_nonmatch_classes.csv")
DIF <- do.call(rbind, lapply(c("exact", "recoverable"), function(st) both(function(m) { k <- m & status == st & is.finite(maxd)
  data.frame(status = st, n_matched = sum(k), n_safe = sum(safe & k), share_safe = if (sum(k)) mean(safe[k]) else NA,
             median_max_diff = median(maxd[k]), q90_max_diff = unname(quantile(maxd[k], 0.9)), n_zero_diff_unscaled = sum(maxd_unsc[k] < 1e-6),
             n_safe_with_protein_factors = sum(maxd_prot[k] < TOL), n_safe_best_candidate = sum(mind[k] < TOL), n_single_candidate_safe = sum(safe & ncand == 1 & k))
})))
wr(DIF, "c2d_difference_summary.csv"); print(DIF, row.names = FALSE, digits = 3); print(REC, row.names = FALSE); print(NMC, row.names = FALSE)

# ---- figure ------------------------------------------------------------------------------------------------------------------------------------
png(file.path(O, "c2d_fig1_site_mapping_summary.png"), 12, 5, units = "in", res = 150)
par(mfrow = c(1, 2), mar = c(7, 4.5, 3, 1))
lv <- c("exact", "recoverable", "not_recoverable_accession_present", "not_recoverable_absent")
tb <- rbind(all = table(factor(status, lv)), labelled_DE = table(factor(status[s$labelled], lv)))
tb <- tb / rowSums(tb) * 100
barplot(tb, beside = TRUE, col = c("#999999", "#0072B2"), names.arg = c("exact", "recoverable\n(other key)", "not recov.:\naccession present", "not recov.:\nabsent"), las = 1, cex.names = 0.75,
        ylab = "% of Supp Data 2a sites", main = "Status of the 28,527 sites in the repository file", legend.text = c("all 28,527", "6,768 DE"), args.legend = list(bty = "n", x = "topright"))
k <- is.finite(maxd); h <- pmax(maxd[k], 1e-4)
hist(log10(h), breaks = 40, col = "#999999", border = "white", main = "Max |log2 diff| over 10 samples (matched sites)", xlab = "log10 max |log2 difference| (0.01 = -2)", ylab = "SD2a sites")
abline(v = log10(TOL), col = "#D55E00", lwd = 2); legend("topright", sprintf("safe (< 0.01): %d of %d matched", sum(safe), sum(k)), bty = "n", text.col = "#D55E00")
invisible(dev.off())

# ---- metrics (append-only; rows with variant prefix amendR1_QM_ are replaced on re-run) ---------------------------------------------------------------
mt <- rbind(
  data.frame(variant = "amendR1_QM_status", metric = c("n_exact", "n_recoverable_other_key", "n_not_recoverable", "n_exact_labelled_DE"), value = c(sum(status == "exact"), sum(status == "recoverable"), sum(grepl("^not", status)), sum(status == "exact" & s$labelled))),
  data.frame(variant = "amendR1_QM_duplicates", metric = c("n_duplicated_repo_ids", "n_dup_ids_distinct_values"), value = c(length(dup_ids), sum(!ident))),
  data.frame(variant = "amendR1_QM_difference", metric = c("n_safe_to_carry", "share_safe_of_matched", "median_max_abs_log2_diff"), value = c(sum(safe), mean(safe[k]), median(maxd[k]))))
mt <- data.frame(claim = "C2", mt, n = nrow(s))
mf <- "results/metrics.csv"; ln <- readLines(mf); keep <- !grepl("\"amendR1_QM_", ln); if (!all(keep)) writeLines(ln[keep], mf)
write.table(mt, mf, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE, quote = TRUE)
writeLines(capture.output(sessionInfo()), "env/sessionInfo_R1.txt")   # same session as c1c_qc.R (both run with the same R and packages)
