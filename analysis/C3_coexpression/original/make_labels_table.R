inp  <- read.table("JUMPn/JUMPnc/HGG_whl_input.txt", header=TRUE, sep="\t", quote="")
ship <- read.table("scratch/JUMPnc/shipped_hgg_whl_cluster_output.txt", header=TRUE, sep="\t")
new  <- read.table("scratch/JUMPnc/hgg_whl_cluster_output.txt", header=TRUE, sep="\t")
num  <- c("Cortex1","Cortex2","Cortex3","PD1","PD2","PD3","PD4","TRK1","TRK2","TRK3")
stopifnot(nrow(inp)==nrow(ship), nrow(ship)==nrow(new),
          isTRUE(all.equal(as.matrix(inp[,num]), as.matrix(ship[,num]), check.attributes=FALSE)),
          isTRUE(all.equal(as.matrix(inp[,num]), as.matrix(new[,num]),  check.attributes=FALSE)))
labels <- data.frame(Protein_Accession=inp$Protein_Accession, paper=inp$Cluster,
                     shipped=ship$group, rerun=new$group)
write.csv(labels, "analysis/C3_coexpression/original/labels_with_accession.csv", row.names=FALSE)
cat("wrote", nrow(labels), "rows\n")
