# C4 original arm: mapping coverage, Q-A (proteome clusters) and the main-setting ORA (Q-B main). Run from repo root.
set.seed(1)
source("analysis/C4_pathways/improved/c4_lib.R")
O <- file.path(C4, "original")
ctx <- build_ctx(); sets <- ctx$sets; ortho <- ctx$ortho
cat("proteome genes:", length(ctx$bg$protein), " phospho genes:", length(ctx$bg$phospho),
    " in both:", length(intersect(ctx$bg$protein, ctx$bg$phospho)), "\n")

# ---- mapping coverage --------------------------------------------------------
hm <- split(ortho$V2, ortho$V1)
cov <- do.call(rbind, lapply(names(sets), function(s) {
  g <- sets[[s]]; mapped <- g %in% names(hm)
  m <- unique(unlist(hm[g[mapped]], use.names = FALSE))
  upp <- lapply(ctx$bg, function(b) sum(g %in% toupper(b)))
  data.frame(set = s, n_human = length(g), n_human_mapped = sum(mapped), n_mouse = length(m),
             n_in_proteome_bg = length(intersect(m, ctx$bg$protein)), n_in_phospho_bg = length(intersect(m, ctx$bg$phospho)),
             n_upper_proteome_bg = upp$protein, n_upper_phospho_bg = upp$phospho)
}))
cov$in_S_wide <- cov$set %in% S_WIDE
write.csv(cov, file.path(O, "c4_mapping_coverage_per_set.csv"), row.names = FALSE)
allm <- unique(unlist(hm[unique(unlist(sets))], use.names = FALSE))
bgcov <- do.call(rbind, lapply(names(ctx$bg), function(l) data.frame(layer = l, n_bg = length(ctx$bg[[l]]),
  bg_in_ortholog_table = mean(ctx$bg[[l]] %in% ortho$V2), bg_in_any_mapped_set = mean(ctx$bg[[l]] %in% allm),
  bg_upper_in_human_symbols = mean(toupper(ctx$bg[[l]]) %in% unlist(sets)))))
cs <- function(d) c(share_human_mapped = sum(d$n_human_mapped) / sum(d$n_human),
                    share_mapped_in_proteome_bg = sum(d$n_in_proteome_bg) / sum(d$n_mouse),
                    share_mapped_in_phospho_bg = sum(d$n_in_phospho_bg) / sum(d$n_mouse),
                    share_human_matched_upper_proteome = sum(d$n_upper_proteome_bg) / sum(d$n_human),
                    share_human_matched_upper_phospho = sum(d$n_upper_phospho_bg) / sum(d$n_human))
setcov <- rbind(data.frame(group = "all_186_sets", t(cs(cov))), data.frame(group = "S_wide", t(cs(cov[cov$in_S_wide, ]))))
write.csv(bgcov, file.path(O, "c4_mapping_coverage_background.csv"), row.names = FALSE)
write.csv(setcov, file.path(O, "c4_mapping_coverage_summary.csv"), row.names = FALSE)
print(setcov); print(bgcov)

# ---- Q-A: ORA of WP-C1 / WP-C5 (background = proteome genes) -----------------
lab <- read.csv("analysis/C3_coexpression/original/labels_with_accession.csv", stringsAsFactors = FALSE)
sb_p <- sets_in_bg(sets, ctx$bg$protein, "table", ortho)
qa <- NULL; qa_n <- NULL
for (cl in c("WP-C1", "WP-C5")) {
  acc <- lab$Protein_Accession[lab$paper == cl]
  g <- genes_of(acc, ctx$maps$protein)
  qa_n <- rbind(qa_n, data.frame(cluster = cl, n_accessions = length(acc), n_acc_with_gene = sum(acc %in% names(ctx$maps$protein)),
                                 n_genes = length(g)))
  qa <- rbind(qa, tag(ora(g, ctx$bg$protein, sb_p), cluster = cl))
}
write.csv(qa, file.path(O, "c4_qa_ora_tables.csv"), row.names = FALSE)
write.csv(qa_n, file.path(O, "c4_qa_cluster_sizes.csv"), row.names = FALSE)
get <- function(cl, s) qa[qa$cluster == cl & qa$set == s, ]
cc <- get("WP-C1", CTRL); th <- do.call(rbind, lapply(WPC5, function(s) get("WP-C5", s)))
qa_sum <- data.frame(item = c("cell_cycle_in_WP-C1", WPC5),
                     cluster = c("WP-C1", rep("WP-C5", 6)), size = c(cc$size, th$size), overlap = c(cc$overlap, th$overlap),
                     p = c(cc$p, th$p), BH = c(cc$BH, th$BH), enriched = c(cc$BH, th$BH) < 0.05)
n_th <- sum(qa_sum$enriched[-1])
cat("\nQ-A:\n"); print(qa_sum); cat("cell cycle enriched in WP-C1:", qa_sum$enriched[1], "; WP-C5 themes enriched:", n_th, "/ 6\n")
write.csv(qa_sum, file.path(O, "c4_qa_summary.csv"), row.names = FALSE)
cat("Q-A consistent with Fig. 3 (pre-registered rule):", qa_sum$enriched[1] && n_th >= 3, "\n")

# ---- main setting: ORA on anova_perm gene-level DE lists, own backgrounds ----
pf <- read.csv(gzfile("analysis/C2_differential/c2_protein_DE_flags.csv.gz"), check.names = FALSE, stringsAsFactors = FALSE)
hf <- read.csv(gzfile("analysis/C2_differential/c2_phospho_DE_flags.csv.gz"), check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(setequal(pf$Feature, rownames(ctx$L$protein)), setequal(hf$Feature, rownames(ctx$L$phospho)))
deP <- genes_of(pf$Feature[pf[["anova_perm_FC1.5_FDR0.05"]]], ctx$maps$protein)
deH <- genes_of(hf$Feature[hf[["anova_perm_FC2_FDR0.05"]]], ctx$maps$phospho)
cat("\nmain-setting DE features: protein", sum(pf[["anova_perm_FC1.5_FDR0.05"]]), "-> genes", length(deP),
    "; phospho", sum(hf[["anova_perm_FC2_FDR0.05"]]), "-> genes", length(deH), "\n")
sbH <- sets_in_bg(sets, ctx$bg$phospho, "table", ortho)
tP <- ora(deP, ctx$bg$protein, sb_p); tH <- ora(deH, ctx$bg$phospho, sbH)
write.csv(rbind(tag(tP, layer = "proteome", method = "ORA"), tag(tH, layer = "phospho", method = "ORA")),
          file.path(O, "c4_main_ora_tables.csv"), row.names = FALSE)
cl <- classify_wide(tP, tH); cs0 <- claim_summary(cl)
write.csv(cl, file.path(O, "c4_main_S_wide_classification.csv"), row.names = FALSE)
print(cl); print(cs0)
desc <- function(sets_) { x <- rbind(tag(tP[tP$set %in% sets_, ], layer = "proteome"), tag(tH[tH$set %in% sets_, ], layer = "phospho"))
  x[, c("layer", "set", "size", "overlap", "n_list", "p", "BH")] }
write.csv(desc(c(S_CORE, CTRL)), file.path(O, "c4_main_S_core_and_cell_cycle.csv"), row.names = FALSE)
print(desc(c(S_CORE, CTRL)))
cat("sets tested (>=5 in bg): proteome", nrow(tP), " phospho", nrow(tH), "; enriched proteome", length(enriched(tP)), " phospho", length(enriched(tH)), "\n")

m <- rbind(mrow("QA", "cell_cycle_enriched_WP-C1", as.integer(qa_sum$enriched[1]), qa_n$n_genes[1]),
           mrow("QA", "WP-C5_themes_enriched_of6", n_th, qa_n$n_genes[2]),
           mrow("main_ORA", "n_DE_genes_proteome", length(deP), length(ctx$bg$protein)),
           mrow("main_ORA", "n_DE_genes_phospho", length(deH), length(ctx$bg$phospho)),
           mrow("main_ORA", "n_enriched_sets_proteome", length(enriched(tP)), nrow(tP)),
           mrow("main_ORA", "n_enriched_sets_phospho", length(enriched(tH)), nrow(tH)),
           mrow("main_ORA", "S_wide_enriched_either", cs0$n_enriched_either, 9),
           mrow("main_ORA", "S_wide_phospho_only_share", cs0$phospho_only_share, cs0$n_enriched_either),
           mrow("main_ORA", "claim_holds", as.integer(cs0$claim_holds), 9),
           mrow("mapping", "share_S_wide_human_mapped", setcov$share_human_mapped[2], sum(cov$n_human[cov$in_S_wide])))
put_metrics(m)
saveRDS(list(deP = deP, deH = deH), file.path("scratch", "c4_main_lists.rds"))
writeLines(capture.output(sessionInfo()), "env/sessionInfo_C4.txt")
