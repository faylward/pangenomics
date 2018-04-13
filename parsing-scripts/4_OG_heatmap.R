
# labels2colors
#
#---------------------------------------------------------------------------------------------------------
# This function converts integer numerical labels labels into color names in the order either given by
# colorSeq,
# or (if colorSeq==NULL) by standardColors(). If GreyIsZero == TRUE, labels 0 will be assigned
# the color grey; otherwise presence of labels below 1 will trigger an error.
# dimensions of labels (if present) are preserved.

labels2colors = function(labels, zeroIsGrey = TRUE, colorSeq = NULL, naColor = "grey",
                         commonColorCode = TRUE)
{
  if (is.null(colorSeq)) colorSeq = standardColors();

  if (is.numeric(labels))
  {
    if (zeroIsGrey) minLabel = 0 else minLabel = 1
    if (any(labels<0, na.rm = TRUE)) minLabel = min(c(labels), na.rm = TRUE)
    nLabels = labels;
  } else {
    
    if (commonColorCode)
    {
      factors = factor(c(as.matrix(as.data.frame(labels))))
      nLabels = as.numeric(factors)
      dim(nLabels)= dim(labels);
    } else {
      labels = as.matrix(as.data.frame(labels));
      factors = list();
      for (c in 1:ncol(labels))
        factors[[c]] = factor(labels[, c]);
      nLabels = sapply(factors, as.numeric)
    }
  }
      
  if (max(nLabels, na.rm = TRUE) > length(colorSeq))
  {
     nRepeats = as.integer((max(labels)-1)/length(colorSeq)) + 1;
     warning(paste("labels2colors: Number of labels exceeds number of avilable colors.", 
                   "Some colors will be repeated", nRepeats, "times."))
     extColorSeq = colorSeq;
     for (rep in 1:nRepeats) 
       extColorSeq = c(extColorSeq, paste(colorSeq, ".", rep, sep=""));
  } else {
     nRepeats = 1;
     extColorSeq = colorSeq;
  }
  colors = rep("grey", length(nLabels));
  fin = !is.na(nLabels);
  colors[!fin] = naColor;
  finLabels = nLabels[fin];
  colors[fin][finLabels!=0] = extColorSeq[finLabels[finLabels!=0]];
  if (!is.null(dim(labels)))
    dim(colors) = dim(labels);
  
  colors;
}

##########################################################
################# Get Phylogeny ##########################
##########################################################
library(dendextend)
library(ape)

tree <- read.tree(file="proteins_for_phylogeny.faa.final_tree.nw")
new_tree <- drop.tip(tree, c("GCA.000754365.1.ASM75436v1", "GCA.000146505.1.ASM14650v1", "GCA.000009925.1.ASM992v1"))

c <- chronos(new_tree)
tree_dendrogram <- as.dendrogram(c)


# open the protienotho output and place it into a data.frame called "x". The first line becomes the column names. 
x <- read.table(file="final_set_dir.proteinortho.cluster", sep="\t", header=T, row.names=1)
subset <- x[,4:227]
pa <- subset != "*"
counts <- rowSums(pa)
counts1 <- colSums(pa)

to_include <- names(counts[counts > 10])
final <- t(pa[to_include,])

#corr <- cor(final, method="pearson")
#n <- as.dist(1-corr)
n <- dist(t(final))
clust <- hclust(n, method="complete")

# plot heatmap
library(RColorBrewer)
hmcol = colorRampPalette(brewer.pal(9, "Blues"))(100)
library(gplots)
pdf(file="OG_clustering2.pdf", width=16, height=10)
heatmap.2(as.matrix(final*1), Rowv=tree_dendrogram, Colv=as.dendrogram(clust), trace="none", dendrogram="both", scale="none", col=hmcol, key=F)
dev.off()


# get core and variable genomes
pan <- cutree(clust, h=14)
pangenome <- gsub("1", "core_pangenome", pan)
pangenome <- gsub("2", "variable_pangenome", pangenome)
names(pangenome) <- names(pan)


annot <- read.table(file="full_annotation.tsv", sep="\t", header=T, row.names=1, quote="") 
pan_stats <- cbind(pangenome, annot[names(pangenome),])
write.table(pan_stats, file="pangenome_stats_and_annotation.tsv", quote=F, sep="\t")


#####################################################################################################################################
############### tanglegram comparing functional clustering of genomes to biogeographic clustering of genomes ########################
#####################################################################################################################################

# now cluster genomes by their distribution
fbio <- read.table(file="log10_tpm_all1.txt", header=T, row.names=1, sep="\t")
names <- gsub("-", ".", row.names(fbio))
row.names(fbio) <- names
freq <- rowSums(fbio > 0)
mfreq <- names(freq)[which(freq >=5)]

cgen <- names(counts1[counts1>00])
most_freq <- intersect(mfreq, cgen)
bio <- fbio[most_freq,]

biocor <- cor(t(bio), method="pearson")
bn <- as.dist(1-biocor)
bclust <- hclust(bn, method="average")
t <- cutree(bclust, h=1.1)

# cluster genomes by functional profiles
to_include <- names(counts[counts > 10])
names <- gsub("_", ".", colnames(pa))
colnames(pa) <- names
final <- pa[to_include, most_freq]

#library(philentropy)
#n <- distance(t(final), method="jaccard")
corr <- cor(final, method="pearson")
n <- as.dist(1-corr)
fclust <- hclust(n, method="average")
#fclust <- hclust(as.dist(n), method="average")


# plot heatmap
library(RColorBrewer)
hmcol = colorRampPalette(brewer.pal(9, "Blues"))(100)
library(gplots)
pdf(file="functional_clustering_by_biogeography.pdf", width=6, height=7)
heatmap.2(as.matrix(bio), Rowv=as.dendrogram(fclust), Colv=T, trace="none", dendrogram="both", scale="none", col=hmcol, RowSideColors=labels2colors(t, colorSeq=c("#5e3c99", "#e66101")), labRow=c(rep("", dim(bio)[1])), labCol=c(rep("", dim(bio)[2])), key=F)
dev.off()









# plot heatmap
library(RColorBrewer)
hmcol = colorRampPalette(brewer.pal(9, "Blues"))(100)
library(gplots)
pdf(file="genome_clustering.pdf", width=16, height=10)
heatmap.2(as.matrix(final*1), Rowv=T, Colv=T, trace="none", dendrogram="both", scale="none", col=hmcol, key=F)
dev.off()






tanglegram(bclust, fclust)










