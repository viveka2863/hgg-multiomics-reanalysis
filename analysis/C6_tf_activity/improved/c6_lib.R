# Shared helpers for C6 TF activity (design pre-registered in LOG.md, "PRE-REGISTRATION: C6 TF activity"). Run from the repo root.
source("analysis/C4_pathways/improved/c4_lib.R")    # also sources c2_lib.R (limma, de_flag, ora, fcs, ortholog + protein map helpers)
C6 <- "analysis/C6_tf_activity"
TF_DIR <- "JUMPn/JUMPnt/TF_DB"
BENCH <- c("MYC","JUND","JUN","EP300","BRCA1","GTF2B","GTF2F1","TBP","TAF1","TAF7","NFKB1","BCL3","E2F4","STAT2",
           "SIN3A","RAD21","CEBPB","IRF3","ZBTB33","BCLAF1","KAT2A","SMARCB1","REST")
stopifnot(length(BENCH) == 23, !anyDuplicated(BENCH))
CRISPR <- c(JUN = "hit", MYC = "hit", JUND = "non-hit", EP300 = "non-hit", BRCA1 = "non-hit")
UNIVERSES <- c("own", "U1", "U3")

read_tf_sets <- function() {
  fs <- list.files(TF_DIR, "^tf_tg_proximal_filtered_.*\\.txt$", full.names = TRUE)
  stopifnot(length(fs) == 115)
  s <- lapply(fs, function(f) { l <- readLines(f, warn = FALSE); v <- trimws(l[-1]); unique(v[nzchar(v)]) })
  names(s) <- sub("^tf_tg_proximal_filtered_(.*)\\.txt$", "\\1", basename(fs))
  stopifnot(all(BENCH %in% names(s))); s
}
# human symbols -> mouse symbols (all orthologs in the table), independent of any universe
map_tf_sets <- function(sets, ortho) {
  hm <- split(ortho$V2, ortho$V1)
  lapply(sets, function(s) unique(unlist(hm[intersect(s, names(hm))], use.names = FALSE)))
}
restrict <- function(mapped, u) lapply(mapped, function(s) intersect(s, u))

# Mann-Whitney probability that a benchmark TF scores better than a non-benchmark TF (ties 0.5)
auc_rank <- function(score, is_b) {            # larger score = better
  r <- rank(score); nb <- sum(is_b); no <- sum(!is_b)
  (sum(r[is_b]) - nb * (nb + 1) / 2) / (nb * no)
}
auc_p <- function(p, is_b) auc_rank(-p, is_b)
summarise_tab <- function(tab) {               # tab: ORA/FCS table with benchmark flag
  b <- tab$benchmark; enr <- !is.na(tab$BH) & tab$BH < 0.05
  nT <- nrow(tab); nB <- sum(b); nE <- sum(enr); x <- sum(b & enr)
  data.frame(n_tested = nT, n_bench_tested = nB, n_enriched = nE, n_bench_enriched = x,
             AUC = auc_p(tab$p, b), AUC_size_only = auc_rank(tab$size, b),
             p_overlap_bench_enriched = if (nE > 0) phyper(x - 1, nB, nT - nB, nE, lower.tail = FALSE) else NA_real_)
}
run_ora <- function(de, u, mapped) {
  t <- ora(de, u, restrict(mapped, u)); t$benchmark <- t$set %in% BENCH; names(t)[1] <- "TF"; t
}
c6_put_metrics <- function(df) {               # idempotent for C6
  names(df) <- c("claim", "variant", "metric", "value", "n")
  f <- "results/metrics.csv"; old <- read.csv(f, stringsAsFactors = FALSE)
  write.csv(old[old$claim != "C6", ], f, row.names = FALSE)
  append_metrics(df)
}
m6 <- function(variant, metric, value, n) data.frame(claim = "C6", variant = variant, metric = metric, value = value, n = n)

# context: gene lists, universes, mapped TF sets
build_c6 <- function() {
  ortho <- read_ortho(); sets <- read_tf_sets(); mapped <- map_tf_sets(sets, ortho)
  pm <- protein_map()
  u <- list(own = sort(unique(unname(pm))), U1 = sort(unique(ortho$V2)), U3 = sort(unique(unlist(mapped, use.names = FALSE))))
  fl <- read.csv(file.path("analysis/C2_differential", "c2_protein_DE_flags.csv.gz"), check.names = FALSE, stringsAsFactors = FALSE)
  gl <- read.csv(file.path("analysis/C2_differential", "c2_protein_genotype_lists.csv.gz"), check.names = FALSE, stringsAsFactors = FALSE)
  lists <- list(anova_perm_FC1.5_FDR0.05 = fl$Feature[fl[["anova_perm_FC1.5_FDR0.05"]]],
                DE_all10 = gl$Feature[gl$DE_all10], DE_noTRK3 = gl$Feature[gl$DE_noTRK3])
  genes <- lapply(lists, genes_of, map = pm)
  list(ortho = ortho, sets = sets, mapped = mapped, pm = pm, u = u, lists = lists, genes = genes, fl = fl)
}
