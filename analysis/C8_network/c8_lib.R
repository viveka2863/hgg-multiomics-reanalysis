# C8 shared functions: candidate set, three layers, rule R0. Pre-registered in LOG.md (2026-10-20).
set.seed(1)
S <- c(paste0("PD",1:4), paste0("TRK",1:3), paste0("Cortex",1:3))
grp_of <- function(s) ifelse(grepl("^PD", s), "P", ifelse(grepl("^TRK", s), "N", "C"))
ordered_NPC <- function(N, P, C) !is.na(N) & !is.na(P) & !is.na(C) & N > P & P > C
gmeans <- function(m, samples) {   # m: matrix rows x S (colnames S)
  g <- grp_of(samples)
  sapply(c(N = "N", P = "P", C = "C"), function(k) rowMeans(m[, samples[g == k], drop = FALSE]))
}

load_all <- function() {
  d <- "JUMPn/"
  # kinase-substrate (header has 10 names, rows 11 fields)
  ks <- read.delim(paste0(d, "JUMPnk/PhosphoWholeNormalization/Mouse_Kinase_Substrate.txt"), header = FALSE,
                   skip = 1, quote = "", stringsAsFactors = FALSE, comment.char = "")
  names(ks) <- c("KinaseGeneName","KinaseUniprot","KinaseOrganism","HomologUniprot","HomologLocation",
                 "OrigSubstrateUniprot","OrigSubstrateGeneName","OrigSubstrateOrganism","OrigSubstrateLocation",
                 "OrigSubstrateSeq","Type")
  stopifnot(ks[1, 1] == "Pak2", ks[1, 2] == "Q64303", ks[1, 3] == "rat", ks[1, 4] == "P31938", ks[1, 5] == "S298")
  ks$K <- toupper(ks$KinaseGeneName)
  # phospho
  ph <- read.delim(paste0(d, "JUMPnk/Integration/pho_hgg_uniq_v5_20160513.txt"), quote = "", stringsAsFactors = FALSE, check.names = FALSE)
  ph <- ph[!grepl("CON_", ph$Gene), ]
  m <- regmatches(ph$Gene, regexec("^(.*)_([A-Z0-9]+)(-[0-9]+)?_([STY][0-9]+(,[STY][0-9]+)*)$", ph$Gene))
  ok <- lengths(m) > 0
  stopifnot(mean(ok) > 0.99)
  ph <- ph[ok, ]; m <- m[ok]
  pm <- as.matrix(ph[, -1]); mode(pm) <- "numeric"
  pm[!is.na(pm) & pm <= 0] <- NA
  pm <- log2(pm)
  colnames(pm) <- c(paste0("PD", 1:4), paste0("TRK", 1:3), paste0("Cortex", 1:3))
  info <- data.frame(id = ph$Gene, gene = sapply(m, `[`, 2), acc = sapply(m, `[`, 3),
                     sites = sapply(m, `[`, 5), stringsAsFactors = FALSE)
  # one row per (id, site)
  sl <- strsplit(info$sites, ",")
  site_long <- data.frame(row = rep(seq_len(nrow(info)), lengths(sl)), site = unlist(sl), stringsAsFactors = FALSE)
  site_long$gene <- info$gene[site_long$row]; site_long$acc <- info$acc[site_long$row]
  site_long$key <- paste(site_long$acc, site_long$site, sep = "_")
  # IKAP
  ik <- read.delim(paste0(d, "JUMPnk/KinaseActivityHeatmap/IKAP_result_v2_20160708.txt"), stringsAsFactors = FALSE, check.names = FALSE)
  km <- as.matrix(ik[, -1]); rownames(km) <- toupper(ik$Kinase)
  km <- km[, c(paste0("PDGFRA", 1:4), paste0("NTRK", 1:3), paste0("CNTRL", 1:3))]
  colnames(km) <- S
  stopifnot(!anyDuplicated(rownames(km)))
  # proteome
  pr <- read.delim(paste0(d, "JUMPnk/PhosphoWholeNormalization/hgg_whole_v1_20160420.txt"), quote = "", stringsAsFactors = FALSE, check.names = FALSE)
  pr <- pr[!grepl("^co\\|", pr[["Protein Accession #"]]), ]
  ch <- c("sig126","sig127N","sig127C","sig128N","sig128C","sig129N","sig129C","sig130N","sig130C","sig131")
  prm <- as.matrix(pr[, ch]); mode(prm) <- "numeric"; colnames(prm) <- S
  prm[!is.na(prm) & prm <= 0] <- NA
  prm <- log2(prm); rownames(prm) <- NULL
  gn <- pr$GN
  keep <- !is.na(gn) & gn != "NA" & rowSums(is.na(prm)) == 0
  prm <- prm[keep, ]; gn <- gn[keep]
  # duplicate gene symbols: keep the row with the highest mean log2 intensity (documented in LOG)
  o <- order(-rowMeans(prm)); prm <- prm[o, ]; gn <- gn[o]
  dup <- duplicated(toupper(gn)); n_dup <- sum(dup)
  prm <- prm[!dup, ]; rownames(prm) <- toupper(gn[!dup])
  # TF targets
  hm <- read.delim(paste0(d, "JUMPnp/example_data/hs_mm_homo_r66.txt"), header = FALSE, stringsAsFactors = FALSE, quote = "")
  names(hm) <- c("hs", "mm"); hm$hs <- toupper(hm$hs); hm$mm <- toupper(hm$mm)
  fs <- list.files(paste0(d, "JUMPnt/TF_DB"), pattern = "^tf_tg_proximal_filtered_", full.names = TRUE)
  tfs <- sub("\\.txt$", "", sub("^tf_tg_proximal_filtered_", "", basename(fs)))
  tg <- lapply(fs, function(f) { x <- readLines(f)[-1]; unique(toupper(trimws(x[nzchar(trimws(x))]))) })
  names(tg) <- toupper(tfs)
  tgm <- lapply(tg, function(x) { mm <- unique(hm$mm[hm$hs %in% x]); mm[mm %in% rownames(prm)] })
  list(ks = ks, pm = pm, info = info, site_long = site_long, km = km, prm = prm, hm = hm, tg = tg, tgm = tgm,
       n_dup = n_dup, n_tf = length(tg))
}

# candidate edges: kinase (in IKAP) -> TF (one of 115) via a measured site
candidates <- function(D) {
  tfmm <- lapply(names(D$tg), function(t) unique(c(t, D$hm$mm[D$hm$hs == t])))
  names(tfmm) <- names(D$tg)
  sl <- D$site_long
  sl$TF <- NA_character_
  for (t in names(tfmm)) sl$TF[toupper(sl$gene) %in% tfmm[[t]] & is.na(sl$TF)] <- t
  sl <- sl[!is.na(sl$TF), ]
  ks <- D$ks[D$ks$K %in% rownames(D$km), ]
  ks$key <- paste(ks$HomologUniprot, ks$HomologLocation, sep = "_")
  e <- merge(ks[, c("K", "key")], unique(sl[, c("key", "TF", "gene", "acc", "site")]), by = "key")
  e <- unique(e)
  e <- e[order(e$K, e$TF, e$site), ]; rownames(e) <- NULL
  names(e)[names(e) == "K"] <- "Kinase"; names(e)[names(e) == "gene"] <- "TF_phospho_gene"
  e
}

# layers for a given sample subset; returns data.frames of group means + flags
layers <- function(D, samples, E = NULL) {
  L1 <- gmeans(D$km, samples)
  # site-level: mean over IDs containing the site
  sl <- D$site_long
  v <- D$pm[sl$row, samples, drop = FALSE]
  agg <- rowsum(v, sl$key, na.rm = FALSE) / as.vector(table(sl$key)[rownames(rowsum(v, sl$key))])
  L2 <- gmeans(agg, samples)
  # TF target layer: per-protein z across the included samples
  z <- t(scale(t(D$prm[, samples, drop = FALSE])))
  z <- z[!is.na(z[, 1]), , drop = FALSE]
  g <- grp_of(samples)
  tgs <- names(D$tgm)
  n_det <- sapply(D$tgm, length)
  L3 <- t(sapply(tgs, function(t) {
    if (n_det[[t]] < 10) return(c(N = NA, P = NA, C = NA))
    s <- colMeans(z[D$tgm[[t]], , drop = FALSE])
    c(N = mean(s[g == "N"]), P = mean(s[g == "P"]), C = mean(s[g == "C"]))
  }))
  list(L1 = L1, L2 = L2, L3 = L3, n_det = n_det)
}

edge_table <- function(D, E, lay) {
  f <- function(M, k) { r <- M[k, , drop = FALSE]; r }
  l1 <- lay$L1[E$Kinase, , drop = FALSE]
  l2 <- lay$L2[E$key, , drop = FALSE]
  l3 <- lay$L3[E$TF, , drop = FALSE]
  out <- cbind(E[, c("Kinase", "TF", "TF_phospho_gene", "acc", "site")],
               L1_N = l1[, "N"], L1_P = l1[, "P"], L1_C = l1[, "C"],
               L2_N = l2[, "N"], L2_P = l2[, "P"], L2_C = l2[, "C"],
               L3_N = l3[, "N"], L3_P = l3[, "P"], L3_C = l3[, "C"])
  out$ord_L1 <- ordered_NPC(out$L1_N, out$L1_P, out$L1_C)
  out$ord_L2 <- ordered_NPC(out$L2_N, out$L2_P, out$L2_C)
  out$ord_L3 <- ordered_NPC(out$L3_N, out$L3_P, out$L3_C)
  out$n_ordered <- out$ord_L1 + out$ord_L2 + out$ord_L3
  out$accept_R0 <- out$n_ordered == 3
  out$accept_V1 <- out$n_ordered >= 2
  rownames(out) <- NULL
  out
}
nodes_of <- function(et, flag) {
  a <- et[et[[flag]], ]
  unique(c(paste0("K:", a$Kinase), paste0("T:", a$TF)))
}
jacc <- function(a, b) if (length(union(a, b)) == 0) NA_real_ else length(intersect(a, b)) / length(union(a, b))
append_metrics <- function(rows) {
  f <- "results/metrics.csv"
  write.table(rows, f, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE, qmethod = "double",
              quote = c(1, 2, 3))
}
