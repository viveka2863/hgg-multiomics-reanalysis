# C8 original arm: candidate set, rule R0, Q-A recovery. Run from repo root.
set.seed(1)
source("analysis/C8_network/c8_lib.R")
out <- "analysis/C8_network/original/"
D <- load_all()
E <- candidates(D)
cat("duplicate proteome gene symbols dropped:", D$n_dup, "\n")
cat("candidate edges:", nrow(E), " kinases:", length(unique(E$Kinase)), " TFs:", length(unique(E$TF)), "\n")
lay <- layers(D, S)
et <- edge_table(D, E, lay)
write.csv(et, paste0(out, "c8_candidate_edges.csv"), row.names = FALSE)
acc <- et[et$accept_R0, ]
write.csv(acc, paste0(out, "c8_accepted_edges_R0.csv"), row.names = FALSE)
cat("accepted R0 edges:", nrow(acc), "\n"); print(acc[, c("Kinase", "TF", "site")])
nodes <- nodes_of(et, "accept_R0")
paper <- c("T:MYC", "T:JUND", "T:EP300", "T:BRCA1", "K:AKT1", "K:MAPK3", "K:CDK5")
rec <- data.frame(node = paper, in_candidate = paper %in% c(paste0("K:", et$Kinase), paste0("T:", et$TF)),
                  accepted_R0 = paper %in% nodes)
akt <- et[et$Kinase == "AKT1" & et$TF == "BRCA1", ]
rec_edge <- nrow(akt[akt$site == "S686" & akt$accept_R0, ]) > 0
print(rec); print(akt[, c("Kinase","TF","site","ord_L1","ord_L2","ord_L3")])
write.csv(rec, paste0(out, "c8_qa_recovery.csv"), row.names = FALSE)
n_rec <- sum(rec$accepted_R0)
qa <- if (n_rec == 7) "recovered" else if (n_rec >= 4) "partly" else "not recovered"
cat("Q-A:", n_rec, "/7 ->", qa, "; AKT1->BRCA1 S686 accepted:", rec_edge, "\n")
append_metrics(data.frame("C8", c("candidate","candidate","candidate","R0","R0","QA","QA"),
  c("n_edges","n_kinases","n_TFs","n_accepted_edges","n_accepted_nodes","paper_nodes_recovered","AKT1_BRCA1_S686_accepted"),
  c(nrow(E), length(unique(E$Kinase)), length(unique(E$TF)), nrow(acc), length(nodes), n_rec, as.numeric(rec_edge)),
  c(nrow(E), nrow(E), nrow(E), nrow(E), nrow(E), 7, 1)))
library(igraph)
g <- graph_from_data_frame(data.frame(from = paste0("K:", et$Kinase), to = paste0("T:", et$TF), acc = et$accept_R0), directed = FALSE)
vn <- V(g)$name; isK <- startsWith(vn, "K:")
draw <- function() {
  par(mar = c(1, 1, 2, 1))
  lay2 <- layout_as_bipartite(g, types = isK)
  plot(g, layout = lay2[, 2:1], vertex.label = sub("^[KT]:", "", vn), vertex.size = 14,
       vertex.shape = ifelse(isK, "circle", "square"), vertex.color = ifelse(vn %in% nodes, "#d95f02", "grey85"),
       vertex.label.cex = 0.6, edge.color = ifelse(E(g)$acc, "#d95f02", "grey70"), edge.width = ifelse(E(g)$acc, 3, 1),
       main = "C8 candidate kinase-TF network (circles kinases, squares TFs; orange = accepted under R0)")
}
png(paste0(out, "c8_fig1_candidate_network.png"), 1400, 800, res = 130); draw(); dev.off()
dir.create("env", showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), "env/sessionInfo_C8.txt")
