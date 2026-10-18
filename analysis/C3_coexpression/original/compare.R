# Base R only. Run from repo root after run_control.sh. Rows align by order (4,703).
ari <- function(a, b) {
  t <- table(a, b); c2 <- function(x) x * (x - 1) / 2
  s <- sum(c2(t)); sa <- sum(c2(rowSums(t))); sb <- sum(c2(colSums(t))); n <- c2(sum(t))
  (s - sa * sb / n) / (0.5 * (sa + sb) - sa * sb / n)
}
inp <- read.table("JUMPn/JUMPnc/HGG_whl_input.txt", header = TRUE, sep = "\t", quote = "")
shp <- read.table("scratch/JUMPnc/shipped_hgg_whl_cluster_output.txt", header = TRUE, sep = "\t")
rer <- read.table("scratch/JUMPnc/hgg_whl_cluster_output.txt", header = TRUE, sep = "\t")
cat("rows: input", nrow(inp), "shipped", nrow(shp), "rerun", nrow(rer), "\n")
stopifnot(nrow(inp) == nrow(shp), nrow(shp) == nrow(rer))
cat("intensities identical input vs shipped:", isTRUE(all.equal(unname(as.matrix(inp[, 5:14])), unname(as.matrix(shp[, 1:10])), check.attributes = FALSE)), "\n")
paper <- inp$Cluster; shipped <- shp$group; rerun <- rer$group
out <- "analysis/C3_coexpression/original/"
xt <- list(paper_vs_shipped = table(paper, shipped), shipped_vs_rerun = table(shipped, rerun), paper_vs_rerun = table(paper, rerun))
for (n in names(xt)) { print(xt[[n]]); write.csv(as.data.frame.matrix(xt[[n]]), paste0(out, "crosstab_", n, ".csv")) }
cat("sizes shipped:\n"); print(table(shipped)); cat("sizes rerun:\n"); print(table(rerun))
res <- data.frame(comparison = c("shipped_vs_rerun", "paper_vs_rerun", "paper_vs_shipped"),
  ARI = c(ari(shipped, rerun), ari(paper, rerun), ari(paper, shipped)))
print(res); write.csv(res, paste0(out, "ari.csv"), row.names = FALSE)
