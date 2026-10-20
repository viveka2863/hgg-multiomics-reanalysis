# Shared helpers for C10 (design pre-registered in LOG.md, 2026-10-21). Run from the repo root.
source("analysis/C4_pathways/improved/c4_lib.R")      # also sources c2_lib.R; defines rd, SAMPLES, GROUP, append_metrics
suppressPackageStartupMessages(library(Rcpp))
C10 <- "analysis/C10_function_prediction"
NNULL <- 1000

sourceCpp(code = '
#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;
// leave-one-gene-out neighbour voting; mem = 0-based member indices, ism = membership flag per gene.
// returns per-member AUC against non-members (ties 0.5)
// [[Rcpp::export]]
NumericVector gba_member_auc(NumericMatrix W, IntegerVector mem, LogicalVector ism) {
  int N = W.nrow(), k = mem.size();
  std::vector<double> total(N, 0.0);
  for (int m = 0; m < k; m++) for (int j = 0; j < N; j++) total[j] += W(j, mem[m]);
  int nn = 0; for (int j = 0; j < N; j++) if (!ism[j]) nn++;
  NumericVector out(k); std::vector<double> v(nn);
  for (int m = 0; m < k; m++) {
    int g = mem[m], c = 0;
    for (int j = 0; j < N; j++) if (!ism[j]) v[c++] = total[j] - W(j, g);
    std::sort(v.begin(), v.end());
    double sg = total[g] - W(g, g);
    double lo = std::lower_bound(v.begin(), v.end(), sg) - v.begin();
    double hi = std::upper_bound(v.begin(), v.end(), sg) - v.begin();
    out[m] = (lo + 0.5 * (hi - lo)) / nn;
  }
  return out;
}
// [[Rcpp::export]]
NumericVector score_member_auc(NumericVector s, IntegerVector mem, LogicalVector ism) {
  int N = s.size(), k = mem.size(), nn = 0;
  std::vector<double> v; for (int j = 0; j < N; j++) if (!ism[j]) v.push_back(s[j]);
  nn = v.size(); std::sort(v.begin(), v.end());
  NumericVector out(k);
  for (int m = 0; m < k; m++) {
    double lo = std::lower_bound(v.begin(), v.end(), s[mem[m]]) - v.begin();
    double hi = std::upper_bound(v.begin(), v.end(), s[mem[m]]) - v.begin();
    out[m] = (lo + 0.5 * (hi - lo)) / nn;
  }
  return out;
}')

# ---- universe, expression ---------------------------------------------------
load_universe <- function() {
  lab <- read.csv("analysis/C3_coexpression/original/labels_with_accession.csv", stringsAsFactors = FALSE)
  raw <- load_protein(); fac <- read.csv("analysis/C2_differential/c2_rescale_factors.csv")
  fac <- setNames(fac$factor, fac$sample)[SAMPLES]
  resc <- sweep(raw, 2, fac, "*")
  stopifnot(!anyDuplicated(lab$Protein_Accession))
  nomatrix <- setdiff(lab$Protein_Accession, rownames(resc))        # the cu| fusion entry is not in the whole-proteome file
  stopifnot(identical(nomatrix, "cu|TPM3_NTRK1|TMP3_NTRK1"))
  map <- protein_map()
  acc <- intersect(lab$Protein_Accession, names(map))               # drops GN == "NA" and the fusion entry
  d <- data.frame(acc = acc, gene = unname(map[acc]), ab = rowMeans(resc[acc, ]), stringsAsFactors = FALSE)
  d <- d[order(-d$ab), ]; d <- d[!duplicated(d$gene), ]            # one row per gene, highest mean abundance
  d <- d[order(d$gene), ]
  X <- log2(resc[d$acc, ]); rownames(X) <- d$gene
  wp <- setNames(lab$paper, lab$Protein_Accession)[d$acc]; names(wp) <- d$gene
  list(X = X, wp = wp, n_labelled_input = nrow(lab), n_dropped_GN_NA = nrow(lab) - length(acc),
       n_dropped_dup_gene = length(acc) - nrow(d), acc = setNames(d$acc, d$gene))
}
get_sets <- function(genes, min_size = 10) {
  sb <- sets_in_bg(read_sets(), genes, "table", read_ortho())
  list(all = sb, main = sb[lengths(sb) >= min_size])
}

# ---- networks ---------------------------------------------------------------
net <- function(X, method = "pearson") {
  W <- cor(t(X), method = method); W[W < 0] <- 0; diag(W) <- 0; W
}
within_group_X <- function(X) {
  g <- GROUP[colnames(X)]
  for (l in levels(g)) { i <- g == l; X[, i] <- X[, i] - rowMeans(X[, i, drop = FALSE]) }
  X
}

# ---- scorers: return per-member AUC for a set (indices into universe) -------
mk_ism <- function(N, idx) { z <- logical(N); z[idx] <- TRUE; z }
auc_gba <- function(W, idx) gba_member_auc(W, as.integer(idx - 1L), mk_ism(nrow(W), idx))
auc_deg <- function(deg, idx) score_member_auc(deg, as.integer(idx - 1L), mk_ism(length(deg), idx))
auc_clu <- function(clus, idx) {                       # share of the other set members in the gene's WP cluster
  k <- length(idx); cl_tot <- tabulate(clus, 5); cnt <- tabulate(clus[idx], 5); nm <- cl_tot - cnt
  out <- numeric(5)
  for (c in which(cnt > 0)) {
    other <- cnt - (seq_len(5) == c); sc <- other / (k - 1); sg <- sc[c]
    out[c] <- (sum(nm[sc < sg]) + 0.5 * sum(nm[sc == sg])) / sum(nm)
  }
  rep(out, cnt)    # one entry per member, grouped by cluster (mean is what is used)
}
set_auc <- function(scorer, idx) mean(scorer(idx))

# null sets, shared by all networks/baselines: 1000 random sets per distinct size, seed 1
make_null_sets <- function(sizes, N) {
  set.seed(1)
  setNames(lapply(sort(unique(sizes)), function(k) replicate(NNULL, sort(sample.int(N, k)), simplify = FALSE)),
           sort(unique(sizes)))
}
# observed set AUCs + null summary + p for a scorer
run_scorer <- function(scorer, idxs, nulls, cores = 10) {
  obs <- vapply(idxs, function(i) set_auc(scorer, i), 1)
  sz <- lengths(idxs)
  nullv <- parallel::mclapply(names(nulls), function(k) vapply(nulls[[k]], function(i) set_auc(scorer, i), 1), mc.cores = cores)
  names(nullv) <- names(nulls)
  nm <- nullv[as.character(sz)]
  p <- (1 + mapply(function(o, n) sum(n >= o), obs, nm)) / (1 + NNULL)
  data.frame(set = names(idxs), size = sz, auc = obs, null_mean = vapply(nm, mean, 1), null_sd = vapply(nm, sd, 1),
             null_q95 = vapply(nm, quantile, 1, probs = .95), p = p, BH = p.adjust(p, "BH"), row.names = NULL, stringsAsFactors = FALSE)
}

mrow <- function(variant, metric, value, n) data.frame(claim = "C10", variant = variant, metric = metric, value = value, n = n)
put_metrics <- function(df, variants) {
  names(df) <- c("claim", "variant", "metric", "value", "n")
  f <- "results/metrics.csv"; old <- read.csv(f, stringsAsFactors = FALSE)
  write.csv(rbind(old[!(old$claim == "C10" & old$variant %in% variants), ], df), f, row.names = FALSE)
}
fig <- function(dir, name, w, h, fn) {
  png(file.path(dir, paste0(name, ".png")), w, h, units = "in", res = 150); fn(); dev.off()
  pdf(file.path(dir, paste0(name, ".pdf")), w, h); fn(); dev.off()
}
