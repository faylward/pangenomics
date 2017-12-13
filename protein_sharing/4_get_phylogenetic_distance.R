library(ape)
library(reshape2)
tree <- read.tree(file="mlst_tree.nw")
m <- cophenetic.phylo(tree)
melted <- melt(m)


new_names = list()
for(i in 1:length(melted$Var1)) {
	name1 <- melted$Var1[i]
	name2 <- melted$Var2[i]
	name3 <- paste(name1, name2, sep="||")
	print(name3)
	new_names[i] <- name3
}
new_names <- as.character(new_names)



row.names(melted) <- new_names
colnames(melted) <- c("genome1", "genome2", "pd")

melted$pd <- melted$pd / max(melted$pd)
melted$pd <- 1- melted$pd


#final <- data.frame(cbind(new_names, melted))
#final$new_names <- as.character(final$new_names)

g <- read.table("genome_distances.txt", row.names=1, sep="\t", header=T)
#final$new_names %in% row.names(g)

merged <- cbind(g, melted[row.names(g),])


pdf(file="genomic_distance_vs_phylogenetic_distance.pdf", width=8, height=8)
plot(merged$pd, merged$av_prop, type="p", col="blue", xlim=c(0, 1), ylim=c(0, 1))
abline(0, 1)
dev.off()

