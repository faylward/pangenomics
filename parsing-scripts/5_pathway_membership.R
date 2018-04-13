

###########################################################
################# Nuo membership ##########################
###########################################################

# variable Nuo cluster
# Nuo genes     A         B          C         D          E          F           G          H         I          J          K         L          M          N
vnuo <- c('COG0838', 'COG0377', 'COG0852', 'COG0649', 'COG1905', 'COG1894', 'COG1034', 'COG1005', 'COG1143','COG0839',  'COG0713', 'COG1009', 'COG1008', 'COG1007')
vnuo <- c('c_09360', 'c_01747', 'c_05728', 'c_03199', 'c_07489', 'c_07436', 'c_07415', 'c_04089', 'c_04093', 'c_09361', 'c_04088', 'c_03381') #,            'c_02919')


# core Nuo cluster
# Nuo genes        A             B          C         D          E          F           G          H         I          J          K         L          M          N
cnuo <-     c('COG0838', 'COG0377', 'COG0852', 'COG0649', 'COG1905', 'COG1894', 'COG1034', 'COG1005', 'COG1143', 'COG0839', 'COG0713', 'COG1009', 'COG1008', 'COG1007')
cnuo <-     c('c_00320', 'c_00950', 'c_00951', 'c_00952',                                  'c_00743', 'c_00741', 'c_00740', 'c_00739', 'c_00742') #, 'c_00738', 'c_01162')



# in TOBG_MED-645 NuoK appears to be encoded by 36907_27 in cluster_02940


# Nqr cluster
#           A         B         C          D          E
#           COG1726   COG1805   COG2869    COG1347    COG2209
nqr <- c('c_00081', 'c_01163', 'c_00176', 'c_00175', 'c_00174')


pr <- c('c_04683')
nosz <- c('c_07394')
plya <- c('c_04584')
cytb <- c('c_00251')
cytc <- c('c_09582')

narg <- c('c_09534')
narh <- c('c_10254', 'c_10390')
narj <- c('c_10346')
nari <- c('c_10345')

pfor1 <- c('c_04962') # COG1013
pfor2 <- c('c_06818') # COG1014
pfor3 <- c('c_09510') # COG1014

# open the protienotho output and place it into a data.frame called "x". The first line becomes the column names. 
x <- read.table(file="final_set_dir.proteinortho.cluster2", sep="\t", header=T, row.names=1, quote="", check.names=F)
subset <- x[,4:227]
pa <- subset != "*"
cnames <- gsub("_", ".", colnames(pa))
colnames(pa) <- cnames

vnuo_matrix <- pa[vnuo[1:12],]
vnuo_sum <- colSums(vnuo_matrix)
vnuo_binary <- as.numeric(vnuo_sum >= 6)
names(vnuo_binary) <- names(vnuo_sum)

cnuo_matrix <- pa[cnuo,]
cnuo_sum <- colSums(cnuo_matrix)
cnuo_binary <- as.numeric(cnuo_sum >= 5)
names(cnuo_binary) <- names(cnuo_sum)

nqr_matrix <- pa[nqr,]
nqr_sum <- colSums(nqr_matrix)
nqr_binary <- as.numeric(nqr_sum > 3)
names(nqr_binary) <- names(nqr_sum)

#write.table(vnuo_binary, file="datasets/vnuo.txt", sep="\t", quote=F)
#write.table(cnuo_binary, file="datasets/cnuo.txt", sep="\t", quote=F)
#write.table(nqr_binary, file="datasets/nqr.txt", sep="\t", quote=F)



pr_val <- as.numeric(pa[pr,])
names(pr_val) <- colnames(pa)
#write.table(pr_val, file="datasets/pr.txt", sep="\t", quote=F)

nosz_val <- as.numeric(pa[nosz,])
names(nosz_val) <- colnames(pa)
#write.table(nosz_val, file="datasets/nosz.txt", sep="\t", quote=F)

plya_val <- as.numeric(pa[plya,])
names(plya_val) <- colnames(pa)
#write.table(plya_val, file="datasets/plya.txt", sep="\t", quote=F)

pfor_val <- as.numeric(pa[pfor,])
names(pfor_val) <- colnames(pa)
#write.table(pfor_val, file="datasets/pfor.txt", sep="\t", quote=F)

cytb_val <- as.numeric(pa[cytb,])
names(cytb_val) <- colnames(pa)
#write.table(cytb_val, file="datasets/cytb.txt", sep="\t", quote=F)

cytc_val <- as.numeric(pa[cytc,])
names(cytc_val) <- colnames(pa)
#write.table(cytc_val, file="datasets/cytc.txt", sep="\t", quote=F)


pfor1_val <- as.numeric(pa[pfor1,])
names(pfor1_val) <- colnames(pa)
write.table(pfor1_val, file="datasets/pfor1.txt", sep="\t", quote=F)

pfor2_val <- as.numeric(pa[pfor2,])
names(pfor2_val) <- colnames(pa)
write.table(pfor2_val, file="datasets/pfora2.txt", sep="\t", quote=F)

pfor3_val <- as.numeric(pa[pfor3,])
names(pfor3_val) <- colnames(pa)
write.table(pfor3_val, file="datasets/pfor3.txt", sep="\t", quote=F)



narg_val <- as.numeric(pa[narg,])
names(narg_val) <- colnames(pa)
#write.table(narg_val, file="datasets/narg.txt", sep="\t", quote=F)

narh_val <- as.numeric(colSums((pa[narh,])))
narh_val2 <- as.numeric(narh_val > 0)
names(narh_val2) <- colnames(pa)
#write.table(narh_val2, file="datasets/narh.txt", sep="\t", quote=F)

narj_val <- as.numeric(pa[narj,])
names(narj_val) <- colnames(pa)
#write.table(narj_val, file="datasets/narj.txt", sep="\t", quote=F)

nari_val <- as.numeric(pa[nari,])
names(nari_val) <- colnames(pa)
#write.table(nari_val, file="datasets/nari.txt", sep="\t", quote=F)











