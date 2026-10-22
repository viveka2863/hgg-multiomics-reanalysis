# Group-review round 1, Q-X: extended sample table (MS samples and RNA samples, with batch and pairing status).
# Design: LOG.md "PRE-REGISTRATION: group-review amendments, round 1". Run from the repo root:
#   Rscript analysis/00_sample_key/make_sample_key_extended.R
# Reads sample_key.csv (unchanged) and Supp Data 3 (data/external/csv/sd3.csv; made by analysis/C11_rna_layer/improved/c11_convert_sd3.R).
# Output: analysis/00_sample_key/sample_key_extended.csv. No RNG is used. Nothing under JUMPn/ or JUMPp/ is written.
options(stringsAsFactors = FALSE)
stopifnot(file.exists("CLAUDE.md"))
source("analysis/C11_rna_layer/improved/c11_lib.R")   # read_sd3() verifies the Supp Data 3 header and returns meta (sample, group, batch)

sk  <- read.csv("analysis/00_sample_key/sample_key.csv", check.names = FALSE)
sd3 <- read_sd3()$meta
stopifnot(nrow(sk) == 10, nrow(sd3) == 18, identical(sk$channel_index, 1:10))

# ---- MS rows (one per TMT channel) --------------------------------------------------------------
sd_order <- c("Cortex1","Cortex2","Cortex3","PD1","PD2","PD3","PD4","TRK1","TRK2","TRK3")   # column order of SD1 / SD2a (verified by values in C2 part B)
grp_ms <- sk$group
rep_ms <- ave(seq_len(10), grp_ms, FUN = seq_along)
ms <- data.frame(
  sample_id = sk$canonical_id, modality = "MS", group = grp_ms, replicate_index = rep_ms, tmt_channel = sk$tmt_channel,
  col_HGG_whl_input = sk$HGG_whl_input, col_hgg_whole_v1 = sk$hgg_whole_v1, col_whl_hgg_v5 = sk$whl_hgg_v5,
  col_pho_hgg_uniq_v5 = sk$pho_hgg_uniq_v5, col_ikap_hgg_v5 = sk$ikap_hgg_v5, col_IKAP_result = sk$IKAP_result,
  supplement_col_SD1_SD2a = match(sk$canonical_id, sd_order), sd3_col = NA_integer_, rna_batch = NA_character_,
  stringsAsFactors = FALSE)

# ---- RNA rows (Supp Data 3 columns 1-18) ----------------------------------------------------------
grp_rna <- c(Cortex = "CNTRL", PDGFRA = "PDGFRA", NTRK = "NTRK")[sd3$group]
rep_rna <- ave(seq_len(18), grp_rna, FUN = seq_along)
rna <- data.frame(
  sample_id = paste0("RNA_", sd3$sample), modality = "RNA", group = unname(grp_rna), replicate_index = rep_rna, tmt_channel = NA_character_,
  col_HGG_whl_input = NA, col_hgg_whole_v1 = NA, col_whl_hgg_v5 = NA, col_pho_hgg_uniq_v5 = NA, col_ikap_hgg_v5 = NA, col_IKAP_result = NA,
  supplement_col_SD1_SD2a = NA_integer_, sd3_col = seq_len(18),
  rna_batch = ifelse(is.na(sd3$batch), "none_stated", paste0("batch", sd3$batch)), stringsAsFactors = FALSE)
stopifnot(identical(rna$rna_batch[rna$group == "PDGFRA"], c(rep("batch1", 3), rep("batch2", 4))),
          identical(rna$rna_batch[rna$group == "NTRK"],   c(rep("batch1", 3), rep("batch2", 4))))

# ---- naive pairing candidate: same group, same replicate index (an assumption, NOT evidence) ----------------
key_ms  <- paste(ms$group,  ms$replicate_index);  key_rna <- paste(rna$group, rna$replicate_index)
ms$pairing_candidate  <- ifelse(key_ms %in% key_rna, rna$sample_id[match(key_ms, key_rna)], NA)
rna$pairing_candidate <- ifelse(key_rna %in% key_ms, ms$sample_id[match(key_rna, key_ms)], NA)
ms$pairing_status  <- ifelse(is.na(ms$pairing_candidate),  "no_rna_candidate", "unverified")
rna$pairing_status <- ifelse(is.na(rna$pairing_candidate), "no_ms_candidate",  "unverified")
x <- rbind(ms, rna)
write.csv(x, "analysis/00_sample_key/sample_key_extended.csv", row.names = FALSE, na = "")
cat(sprintf("rows: %d (MS %d, RNA %d)\n", nrow(x), sum(x$modality == "MS"), sum(x$modality == "RNA")))
print(table(modality = x$modality, pairing_status = x$pairing_status))
print(table(x$modality, x$group))
