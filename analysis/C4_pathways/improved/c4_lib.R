# Shared helpers for C4 (design pre-registered in LOG.md, 2026-10-19). Run from the repo root.
source("analysis/C2_differential/improved/c2_lib.R")
C4 <- "analysis/C4_pathways"
OUT <- file.path(C4, "improved")        # used by save_fig() from c2_lib.R

S_WIDE <- c("KEGG_GLIOMA", "KEGG_NEUROTROPHIN_SIGNALING_PATHWAY", "KEGG_MTOR_SIGNALING_PATHWAY",
            "KEGG_MAPK_SIGNALING_PATHWAY", "KEGG_ERBB_SIGNALING_PATHWAY", "KEGG_INSULIN_SIGNALING_PATHWAY",
            "KEGG_PATHWAYS_IN_CANCER", "KEGG_P53_SIGNALING_PATHWAY", "KEGG_VEGF_SIGNALING_PATHWAY")
S_CORE <- S_WIDE[1:4]
CTRL <- "KEGG_CELL_CYCLE"
WPC5 <- c("KEGG_CITRATE_CYCLE_TCA_CYCLE", "KEGG_PYRUVATE_METABOLISM", "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
          "KEGG_OXIDATIVE_PHOSPHORYLATION", "KEGG_AXON_GUIDANCE", "KEGG_CALCIUM_SIGNALING_PATHWAY")
MIN_SET <- 5

# ---- gene sets and orthologs ------------------------------------------------
read_sets <- function() {
  fs <- list.files("JUMPn/JUMPnp/kegg_DB", "_gmt\\.txt$", full.names = TRUE)
  ln <- lapply(fs, readLines, warn = FALSE)
  stopifnot(length(fs) == 186, all(vapply(ln, function(l) grepl("^http", l[1]), TRUE)))
  s <- lapply(ln, function(l) { v <- trimws(l[-1]); unique(v[nzchar(v)]) })
  names(s) <- sub("_gmt\\.txt$", "", basename(fs)); s
}
read_ortho <- function() {
  h <- read.delim("JUMPn/JUMPnp/example_data/hs_mm_homo_r66.txt", header = FALSE, quote = "",
                  stringsAsFactors = FALSE, na.strings = character(0))
  stopifnot(ncol(h) == 2); h
}
# gene set -> members of the background (mouse symbols), via the ortholog table or upper-casing
sets_in_bg <- function(sets, bg, mode = c("table", "upper"), ortho = NULL) {
  mode <- match.arg(mode)
  if (mode == "table") {
    hm <- split(ortho$V2, ortho$V1)
    lapply(sets, function(s) intersect(unique(unlist(hm[intersect(s, names(hm))], use.names = FALSE)), bg))
  } else {
    bgU <- toupper(bg); lapply(sets, function(s) bg[bgU %in% s])
  }
}

# ---- feature -> gene maps, backgrounds --------------------------------------
protein_map <- function() {
  p <- rd("JUMPn/JUMPnk/PhosphoWholeNormalization/hgg_whole_v1_20160420.txt")
  p <- p[!grepl("^co[|]", p[["Protein Accession #"]]) & p$GN != "NA", ]
  setNames(p$GN, p[["Protein Accession #"]])
}
phospho_map <- function(sites) { g <- sub("_.*", "", sites); names(g) <- sites; g[g != "NA"] }
genes_of <- function(features, map) unique(unname(map[intersect(features, names(map))]))

# gene-level statistic: -log10(limma F adj.P), best feature per gene
gene_stat <- function(x, map) {
  p <- limma_F_p(x, GROUP[colnames(x)]); adj <- p.adjust(p, "BH"); names(adj) <- rownames(x)
  ok <- intersect(names(map), names(adj))
  tapply(-log10(adj[ok]), map[ok], max)
}

# ---- enrichment methods -----------------------------------------------------
# one-sided Fisher exact test == hypergeometric upper tail
ora <- function(de, bg, sb) {
  de <- intersect(de, bg); N <- length(bg); n <- length(de)
  sb <- sb[lengths(sb) >= MIN_SET]
  K <- lengths(sb); a <- vapply(sb, function(s) sum(s %in% de), 1L)
  p <- phyper(a - 1, K, N - K, n, lower.tail = FALSE)
  data.frame(set = names(sb), size = K, overlap = a, n_list = n, n_bg = N, p = p, BH = p.adjust(p, "BH"),
             row.names = NULL, stringsAsFactors = FALSE)
}
fcs <- function(stat, sb) {
  sb <- sb[lengths(sb) >= MIN_SET]
  p <- vapply(sb, function(s) { i <- names(stat) %in% s
    suppressWarnings(wilcox.test(stat[i], stat[!i], alternative = "greater", exact = FALSE)$p.value) }, 1)
  data.frame(set = names(sb), size = lengths(sb), overlap = NA_integer_, n_list = NA_integer_, n_bg = length(stat),
             p = p, BH = p.adjust(p, "BH"), row.names = NULL, stringsAsFactors = FALSE)
}
gsea <- function(stat, sb) {
  sb <- sb[lengths(sb) >= MIN_SET]
  set.seed(1)
  r <- suppressWarnings(fgsea::fgsea(pathways = sb, stats = stat, minSize = MIN_SET, maxSize = Inf, scoreType = "pos"))
  r <- as.data.frame(r)[, c("pathway", "size", "pval")]
  r <- r[match(names(sb), r$pathway), ]
  data.frame(set = names(sb), size = lengths(sb), overlap = NA_integer_, n_list = NA_integer_, n_bg = length(stat),
             p = r$pval, BH = p.adjust(r$pval, "BH"), row.names = NULL, stringsAsFactors = FALSE)
}
enriched <- function(tab, thr = 0.05) tab$set[!is.na(tab$BH) & tab$BH < thr]

# ---- Q-B classification -----------------------------------------------------
classify_wide <- function(tabP, tabH, thr = 0.05, sets = S_WIDE) {
  eP <- sets %in% enriched(tabP, thr); eH <- sets %in% enriched(tabH, thr)
  data.frame(set = sets, proteome = eP, phospho = eH,
             class = ifelse(eP & eH, "both", ifelse(eH, "phospho-only", ifelse(eP, "proteome-only", "neither"))),
             stringsAsFactors = FALSE)
}
claim_summary <- function(cl) {
  either <- sum(cl$class != "neither"); po <- sum(cl$class == "phospho-only")
  share <- if (either > 0) po / either else NA_real_
  data.frame(n_enriched_either = either, n_phospho_only = po, n_proteome_only = sum(cl$class == "proteome-only"),
             n_both = sum(cl$class == "both"), phospho_only_share = share,
             claim_holds = either >= 3 && !is.na(share) && share >= 0.5)
}

# ---- output helpers ---------------------------------------------------------
put_metrics <- function(df) {      # idempotent: replaces C4 rows with the same variant names, then uses append_metrics()
  names(df) <- c("claim", "variant", "metric", "value", "n")
  f <- "results/metrics.csv"; old <- read.csv(f, stringsAsFactors = FALSE)
  write.csv(old[!(old$claim == "C4" & old$variant %in% df$variant), ], f, row.names = FALSE)
  append_metrics(df)
}
mrow <- function(variant, metric, value, n) data.frame(claim = "C4", variant = variant, metric = metric, value = value, n = n)
tag <- function(tab, ...) { a <- list(...); for (k in names(a)) tab[[k]] <- a[[k]]; tab }

# context: log2 matrices, maps, backgrounds, gene-level statistics (full 10 samples)
build_ctx <- function() {
  L <- list(protein = log2(load_protein()), phospho = log2(load_phospho()))
  maps <- list(protein = protein_map(), phospho = phospho_map(rownames(L$phospho)))
  bg <- lapply(maps, function(m) sort(unique(unname(m))))
  stat <- list(protein = gene_stat(L$protein, maps$protein), phospho = gene_stat(L$phospho, maps$phospho))
  stopifnot(setequal(names(stat$protein), bg$protein), setequal(names(stat$phospho), bg$phospho))
  list(L = L, maps = maps, bg = bg, stat = stat, sets = read_sets(), ortho = read_ortho())
}
