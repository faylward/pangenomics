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
cyt_names <- c("Cytochrome B subunit", "cytochrome c", "Cytochrome c biogenesis protein", "Cytochrome C family protein", "Cytochrome oxidase subunit IV", "cytochrome C family protein", "Cytochrome-b5 reductase", "Cytochrome C", "Cytochrome C family protein", "Cytochrome C oxidase subunit I", "Cytochrome C Oxidase, Subunit III", "Cytochrome C peroxidase", "B-type cytochrome subunit", "Succinate dehydrogenase, cytochrome subunit", "caa(3)-type oxidase subunit IV", "cytochrome C, class I", "Cytochrome c", "cytC biogenesis protein", "cytochrome C, class I")

nqr_names <- c('NqrA', 'NqrB', 'NqrC', 'NqrD', 'NqrE')
vnuo_names <- c('NuoA', 'NuoB', 'NuoC', 'NuoD', 'NuoE', 'NuoF', 'NuoG', 'NuoH', 'NuoI','NuoJ','NuoK', 'NuoL')
cnuo_names <- c('NuoA', 'NuoB', 'NuoC',                         'NuoG', 'NuoH', 'NuoI','NuoJ','NuoK', 'NuoL')
pfor_names <- c('PFOR', 'PFOR', 'PFOR', 'PFOR')

markers <- c('c_04683', 'c_04584', 'c_04699', 'c_05548', 'c_04685',   'c_10254', 'c_09534', 'c_06818', cyt)
nqr <- c('c_00081', 'c_01163', 'c_00176', 'c_00175', 'c_00174')
vnuo <- c('c_09360', 'c_01747', 'c_05728', 'c_03199', 'c_07489', 'c_07436', 'c_07415', 'c_04089', 'c_04093', 'c_09361', 'c_04088')
cnuo <- c('c_00320', 'c_00950', 'c_00951', 'c_00952',                                  'c_00743', 'c_00741', 'c_00740', 'c_00739', 'c_00742')
pfor <- c('c_06818', 'c_01788', 'c_09510', 'c_04962')
narg <- c("c_09534", "c_14563")
narh <- c("c_10254")

labels <-  c(cnuo_names, nqr_names, vnuo_names, 'yaa',     'PR',      'PL1',     'PL2',     'PL3',     'PL4',     'Pstarv',  cyt_names, pfor_names, 'NarH',    'NarG')
names1  <- c(cnuo,       nqr,       vnuo,       'c_01489', 'c_04683', 'c_04584', 'c_04699', 'c_05548', 'c_04685', 'c_08238', cyt,       pfor,       narg, narh)
#labels <-  c('PR',      'PL1',     'PL2',     'PL3',     'PL4',       'NarH',    'NarG',    'PFOR', cyt_names,  'vNUO', 'NQR')

pfa <- 1 * as.matrix(pa[names1,bclust$labels])

# now get clade information
c <- read.table(file="clade_designations.tbl", header=T, row.names=1, sep="\t")
row.names(c) <- gsub("_", ".", row.names(c))
shared <- intersect(bclust$labels, row.names(c))
clades <- c[shared,]

clade_colors <- gsub('[*]', '#', clades[names(t),3])

# plot heatmap
library(RColorBrewer)
hmcol = colorRampPalette(brewer.pal(9, "Blues"))(100)
library(gplots)
pdf(file="functional_gene_heatmap.pdf", width=8, height=10)
heatmap.2(as.matrix(pfa), Rowv=F, Colv=as.dendrogram(bclust), trace="none", dendrogram="column", scale="none", col=hmcol, ColSideColors=labels2colors(t, colorSeq=c("#5e3c99", "#e66101")), key=F, margins=c(12,6)) #, labCol=c(rep("", dim(bio)[2])), key=F)
dev.off()


# plot heatmap
library(RColorBrewer)
hmcol = colorRampPalette(brewer.pal(9, "Blues"))(100)
library(gplots)
pdf(file="functional_gene_heatmap_cladecolors.pdf", width=8, height=12)
heatmap.2(as.matrix(pfa), Rowv=F, Colv=as.dendrogram(bclust), trace="none", dendrogram="column", scale="none", col=hmcol, ColSideColors=clade_colors, key=F, margins=c(12,6)) #, labCol=c(rep("", dim(bio)[2])), key=F)
dev.off()

































