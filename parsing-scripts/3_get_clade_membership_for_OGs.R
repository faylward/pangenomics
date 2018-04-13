# open the protienotho output and place it into a data.frame called "x". The first line becomes the column names. 
x <- read.table(file="final_set_dir.proteinortho.cluster.subset", sep="\t", header=T, row.names=1)
subset <- x[,4:221]
pa <- subset != "*"
counts <- rowSums(pa)

genomes <- gsub("_", ".", colnames(pa))
colnames(pa) <- genomes

# get EPI and DEEP clades
fbio <- read.table(file="log10_tpm_all1.final.txt", header=T, row.names=1, sep="\t")
names <- gsub("-", ".", row.names(fbio))
row.names(fbio) <- names
freq <- rowSums(fbio > 0)
mfreq <- names(freq)[which(freq >=5)]
bio <- fbio[mfreq,]
biocor <- cor(t(bio), method="pearson")
bn <- as.dist(1-biocor)
bclust <- hclust(bn, method="average")
t <- cutree(bclust, h=1.1)

fpa <- pa[names(counts[counts > 10]),]
EPI <-  fpa[,names(t[t==2])]
DEEP <- fpa[,names(t[t==1])]

fepi <- rowSums(EPI)
fdeep <- rowSums(DEEP)

epi_prop <- fepi / dim(EPI)[2]
deep_prop <- fdeep / dim(DEEP)[2]

epi_sum <- as.numeric(sum(fepi))
deep_sum <- as.numeric(sum(fdeep))

pvals <- list()
for(i in 1:length(fepi)) {
	one <- as.numeric(fepi[i])
	two <- as.numeric(fdeep[i])
	mat <- matrix(c(one, two, epi_sum, deep_sum), nrow=2)
	test <- fisher.test(mat, alternative="two.sided")
	pvals[i] <- test$p.value
}

fdr <- p.adjust(as.numeric(pvals), method="fdr")
names(fdr) <- names(fepi)

annot <- read.table(file="full_annotation.tsv", header=T, row.names=1, sep="\t", quote="")
subannot <- cbind(annot[names(fdr),], epi_prop, deep_prop, as.numeric(pvals), fdr)
colnames(subannot) <- c("representative", "length", "COG", "COG_Score", "COG_Annotation", "Pfam", "Pfam_score", "Pfam_Annot", "TIGRFam", "TIGRfam_score", "TIGR_Annote", "epi_prop", "deep_prop", "pvals", "fdr")

sig <- subannot[which(subannot$fdr < 0.01),]
#write.table(subannot, file="full_annotation.tsv", quote=F, sep="\t")
write.table(sig, file="enriched_clusters_epi_vs_deep.tsv", quote=F, sep="\t", col.names=NA)




# density plot of OG membership
#plot(sort(counts, decreasing=T), log="y", type="h", lwd=2, lend=2, col="dodgerblue")

# now get clade information
c <- read.table(file="clade_designations.tbl", header=T, row.names=1, sep="\t")
shared <- intersect(colnames(subset), row.names(c))
clades <- c[shared,]


library(naturalsort)
clade_names <- naturalsort(names(table(clades$newclade)[1:10]))

clade_table <- data.frame(matrix(ncol=10, nrow=length(counts)))
for(i in 1:10) {
clade_genomes <- row.names(clades)[which(clades$newclade == clade_names[i])]
clade_subset <- pa[,clade_genomes]
clade_counts <- rowSums(clade_subset)
clade_table[,i] <- clade_counts
}
colnames(clade_table) <- clade_names
row.names(clade_table) <- row.names(pa)

clade_counts <- as.numeric(table(clades$newclade)[colnames(clade_table)])
clade_prop_table <- scale(clade_table, scale=clade_counts, center=F)


# we can also cluster the protein families to find ones that have similar distributions across genomes. 
a <- pa[counts>50,]
d <- dist(a)
clust <- hclust(d, method="average")

# and plot it
plot(clust)
pdf("protein_family_cluster.pdf", width=10, height=5)
plot(clust)
dev.off()



########################################################################
######## marker gene clade-by-clade deep vs epi comparison #############
########################################################################
x <- read.table(file="final_set_dir.proteinortho.cluster.subset", sep="\t", header=T, row.names=1, quote="", check.names=F)
subset <- x[,4:221]
pa <- subset != "*"
cnames <- gsub("_", ".", colnames(pa))
cnames2 <- gsub("-", ".", cnames)
colnames(pa) <- cnames2

cyt <- c("c_00251", "c_02339", "c_05556", "c_03893", "c_04345", "c_04811", "c_06063", "c_03197", "c_05598", "c_03308", "c_09525", "c_07488", "c_10324", "c_07623", "c_09526", "c_09582", "c_10347", "c_06129", "c_09579")
cyt_names <- c("Cytochrome B subunit", "cytochrome c", "Cytochrome c biogenesis protein", "cytochrome C family protein", "Cytochrome oxidase subunit IV", "cytochrome C family protein", "cytochrome-b5 reductase activity, acting on NAD(P)H", "Planctomycete cytochrome C", "cytochrome C family protein", "Cytochrome C oxidase subunit I", "Cytochrome C Oxidase, Subunit III", "cytochrome C peroxidase", "Ni Fe-hydrogenase, b-type cytochrome subunit", "succinate dehydrogenase, cytochrome", "caa(3)-type oxidase subunit IV", "cytochrome C, class I", "Cytochrome c", "cytC biogenesis protein", "cytochrome C, class I")

markers <- c('c_04683', 'c_04584', 'c_04699', 'c_05548', 'c_04685',   'c_10254', 'c_09534', 'c_06818', cyt)
names1  <- c('c_04683', 'c_04584', 'c_04699', 'c_05548', 'c_04685',   'c_10254', 'c_09534', 'c_06818', cyt, 'vNUO', 'NQR')
labels <-  c('PR',      'PL1',     'PL2',     'PL3',     'PL4',       'NarH',    'NarG',    'PFOR', cyt_names,  'vNUO', 'NQR')

# get EPI and DEEP clades
fbio <- read.table(file="log10_tpm_all1.final.txt", header=T, row.names=1, sep="\t")
names <- gsub("-", ".", row.names(fbio))
row.names(fbio) <- names
freq <- rowSums(fbio > 0)
mfreq <- names(freq)[which(freq >=5)]
bio <- fbio[mfreq,]
biocor <- cor(t(bio), method="pearson")
bn <- as.dist(1-biocor)
bclust <- hclust(bn, method="average")
t <- cutree(bclust, h=1.1)

c <- read.table(file="clade_designations.tbl", header=T, row.names=1, sep="\t")
clades <- c("Clade_1", "Clade_2", "Clade_5", "Clade_6", "Clade_7")

gnames <- gsub("_", ".", row.names(c))
row.names(c) <- gnames

epi <-  pa[,names(t[t==2])]
deep <- pa[,names(t[t==1])]

cluster_list <- list()
clade_list <- list()
diff_list <- list()
index <- 0
for(i in 1:length(markers)) {
	cluster <- markers[[i]]
	epi_occur <- pa[cluster, names(t[t==2])]
	dep_occur <- pa[cluster, names(t[t==1])]

	for(j in 1:5) {
		clade <- clades[j]
		genomes <- row.names(c)[which(c$newclade == clade)]

		epi_genomes <- intersect(genomes, names(epi_occur))
		dep_genomes <- intersect(genomes, names(dep_occur))

		epi_frac <- sum(pa[cluster, epi_genomes]) / length(epi_genomes)
		dep_frac <- sum(pa[cluster, dep_genomes]) / length(dep_genomes)

		diff <- dep_frac - epi_frac
		
		print(paste(cluster, clade, diff), sep="\t")
		index <- index + 1
		cluster_list[index] <- cluster
		clade_list[index] <- clade
		diff_list[index] <- diff
	}
}
final1 <- data.frame(cbind(as.character(cluster_list), as.character(clade_list), as.numeric(diff_list)))
colnames(final1) <- c("marker_gene", "clade", "difference")


####### and now for the respiratory complexes

nqr <- c('c_00081', 'c_01163', 'c_00176', 'c_00175', 'c_00174')
cnuo <- c('c_00320', 'c_00950', 'c_00951', 'c_00952',                                  'c_00743', 'c_00741', 'c_00740', 'c_00739', 'c_00742') #, 'c_00738', 'c_01162')
vnuo <- c('c_09360', 'c_01747', 'c_05728', 'c_03199', 'c_07489', 'c_07436', 'c_07415', 'c_04089', 'c_04093', 'c_09361', 'c_04088', 'c_03381') #,            'c_02919')

vnuo_matrix <- pa[vnuo,]
vnuo_sum <- colSums(vnuo_matrix)
vNUO <- as.numeric(vnuo_sum >= 6)
names(vNUO) <- names(vnuo_sum)

cnuo_matrix <- pa[cnuo,]
cnuo_sum <- colSums(cnuo_matrix)
cNUO <- as.numeric(cnuo_sum >= 5)
names(cNUO) <- names(cnuo_sum)

nqr_matrix <- pa[nqr,]
nqr_sum <- colSums(nqr_matrix)
NQR <- as.numeric(nqr_sum > 3)
names(NQR) <- names(nqr_sum)

resp <- rbind(vNUO, cNUO, NQR)
sets <- c('vNUO', 'NQR')
cluster_list <- list()
clade_list <- list()
diff_list <- list()
index <- 0
for(i in 1:2) {
	cluster <- sets[i]
	epi_occur <- resp[cluster, names(t[t==2])]
	dep_occur <- resp[cluster, names(t[t==1])]

	for(j in 1:5) {
		clade <- clades[j]
		genomes <- row.names(c)[which(c$newclade == clade)]

		epi_genomes <- intersect(genomes, names(epi_occur))
		dep_genomes <- intersect(genomes, names(dep_occur))

		epi_frac <- sum(resp[cluster, epi_genomes]) / length(epi_genomes)
		dep_frac <- sum(resp[cluster, dep_genomes]) / length(dep_genomes)

		#print(paste(clade, epi_frac, dep_frac, sep=" ")

		diff <- dep_frac - epi_frac
		
		print(paste(cluster, clade, epi_frac, dep_frac, diff), sep="\t")
		index <- index + 1
		cluster_list[index] <- cluster
		clade_list[index] <- clade
		diff_list[index] <- diff
	}
}
final2 <- data.frame(cbind(as.character(cluster_list), as.character(clade_list), as.numeric(diff_list)))
colnames(final2) <- c("marker_gene", "clade", "difference")

final <- rbind(final1, final2)
final$difference <- as.numeric(as.character(final$difference))

library(ggplot2)
pdf("marker_genes_by_clade.pdf", height=16, width=8)
ggplot(final, aes(x=marker_gene, y=difference, fill=clade, group=clade)) + geom_bar(alpha=0.75, stat="identity", position = position_dodge(width=0.7), width=0.7) + coord_flip() + scale_y_continuous(limits=c(-1, 1)) + theme_bw() + scale_fill_manual(values=c("#0000FF", "#CC0000", "#CC6600", "#336600", "#6600CC")) + scale_x_discrete(limit=names1, labels=labels)
dev.off()































