setwd('/zfs/tillers/panicle/ssapkot/HDP')

library(data.table)
library(sommer)

M <- fread('data/Genomic/SAP_WGS_filtered.012')
ind <- read.table('data/Genomic/SAP_WGS_filtered.012.indv',header=F)

M <- as.matrix(M[,-1])-1
rownames(M) <- ind$V1

### get additive matrix for inbreds

#G <- A.mat(M)

D <-D.mat(M)

saveRDS(D, file='data/unphased_snps/D_matrix_WGS_inbred.rds')
#rm(G)
### get hybrid genotypes from the inbred parents

#males <- read.table('unphased_snps/male_geno.txt',header=F)
#males <- males$V1
#females <- read.table('unphased_snps/female_geno.txt',header=F)
#females <- females$V1
#
#M1 <- M[females,]
#M2 <- M[males,]
#rm(M)
#gc()
#
#HMM <- build.HMM(M1,M2)
#summary(HMM)
#summary(HMM$HMM.add)
#H <- HMM$HMM.add
#H[1:5,1:5]
#
#A <-A.mat(H)
#
##saveRDS(A, file='/zfs/tillers/panicle/ssapkot/git_repo/GP_Sorghumhybrids/data/A-matrix_hybrids.rds')
#
#
#H <- HMM$HMM.dom
#rm(HMM)
#D <- D.mat(H)
#
#save(G,A,D, file='/zfs/tillers/panicle/ssapkot/git_repo/GP_Sorghumhybrids/data/GRM_unphased_HDP.RData')

