# C5 part B, post hoc (not pre-registered): compare the filtered kinase set (BH < 0.05 and >= 3 quantified substrates)
# with the kinase labels printed in the paper's Fig. 4a (typed from the PDF) and with Source Data sheet 4a.
# Run from the repo root: Rscript analysis/C5_kinase_activity/improved/c5b_posthoc_printed_fig4a.R
options(stringsAsFactors = FALSE)
stopifnot(file.exists("CLAUDE.md"))
D <- "analysis/C5_kinase_activity/improved"
lab <- toupper(read.csv(file.path(D, "c5b_fig4a_labels_from_pdf.csv"))$kinase)
t <- read.csv(file.path(D, "c5b_b2_4e_kinases_with_substrate_counts.csv"), check.names = FALSE)
t$k <- toupper(t$Kinase)
ours <- t$k[t$FDR_paper < 0.05 & t$n_substrates >= 3]
src  <- t$k[as.logical(t$in_source_4a_41)]
stopifnot(length(lab) == 41, length(unique(lab)) == 41)
jac <- function(a, b) length(intersect(a, b)) / length(union(a, b))
res <- data.frame(
  reference = c("Source Data 4a (pre-registered benchmark)", "Fig. 4a labels printed in the PDF"),
  n_reference = c(length(src), length(lab)), n_ours = length(ours),
  shared = c(length(intersect(ours, src)), length(intersect(ours, lab))),
  jaccard = c(jac(ours, src), jac(ours, lab)))
res$precision <- res$shared / res$n_ours
res$recall <- res$shared / res$n_reference
write.csv(res, file.path(D, "c5b_posthoc_fig4a_printed_vs_source.csv"), row.names = FALSE)
d <- data.frame(kinase = sort(unique(c(ours, lab, src))))
d$ours <- d$kinase %in% ours; d$printed_fig4a <- d$kinase %in% lab; d$source_data_4a <- d$kinase %in% src
d <- d[!(d$ours & d$printed_fig4a & d$source_data_4a), ]
write.csv(d, file.path(D, "c5b_posthoc_fig4a_differences.csv"), row.names = FALSE)
print(res, digits = 3, row.names = FALSE)
cat("\nKinases that differ between the printed figure and Source Data 4a:",
    paste(setdiff(union(lab, src), intersect(lab, src)), collapse = ", "), "\n")
M <- data.frame(claim = "C5", variant = c("suppB_posthoc_source4a", "suppB_posthoc_printed_fig4a"),
                metric = "jaccard", value = res$jaccard, n = res$n_reference)
mf <- "results/metrics.csv"; old <- read.csv(mf)
old <- old[!(old$claim == "C5" & old$variant %in% M$variant), , drop = FALSE]
write.csv(rbind(old[, names(M)], M), mf, row.names = FALSE)
cat("Outputs written to", D, "and results/metrics.csv\n")
