
library(sommer)
library(tidyverse)

######## Creating hybrid genotypic matrix and relationship matrices

## using a thinned matrix to ease computation load
M <- fread('/scratch1/ssapkot/HDP/data/SAP_1Kb_thinned_filtered.012', header=F)
indv <- read.table('/scratch1/ssapkot/HDP/data/SAP_1Kb_thinned_filtered.012.indv', header=F)
##adjust the first column due to fread and convert into -1,0,1 from 0,1,2
M <- M[,-1]-1
M <- cbind(indv$V1,M) ##use first row to subset later by genotype names because data.table doesn't support rownames

## Make M1 genotypic matrix for female genotype and M2 for male genotype

males <- read.table('/scratch1/ssapkot/HDP/data/hom_snps/male_geno.txt',header=F)
males <- males$V1

females <- read.table('/scratch1/ssapkot/HDP/data/hom_snps/female_geno.txt',header=F)
females <- females$V1

M1 <- subset(M, V1 %in% females)[,-1]
M2 <- subset(M, V1 %in% males)[,-1]

HMM <- build.HMM(M1,M2) ## sommer function to build hybrid genotypes

hybrids <- HMM$data.used
hybrids <- hybrids[,c(3,1,2)]
colnames(hybrids) <- c('CUSo_Geno','Female','Male')
hybrids[,1] <- gsub(':','.',hybrids[,1])

######### Phenotypic analysis and BLUPs calculation
hdp <- read.csv('/zfs/tillers/panicle/ssapkot/HDP/data/Phenotypes/HDP_compiled_Phenotypes_spread.csv', header=T)

traits <- c('DTA','Height','Anthracnose','Lodging','Deoxynivalenol','Fumonisinm.','IVSD','Prolamin','TGW','Yield','Seed_No','Starch','Amylose','Protein','KCal_lb','Fat','Nitrogen')

fac.var <- colnames(hdp[,c(5,9,11,15:18,20:21)])
hdp[,fac.var] <- lapply(hdp[,fac.var],factor)

#df <- hdp[,c(fac.var,traits)]
Year <- c('2019','2020')

VComp <- c()
GenoBLUPs <- c()
for (i in Year){
    varcomp <- c()
    BLUPs <- unique(hdp[,c(5,9,11)])
    
    for (j in traits){
        blups <- c()
        df <- hdp[Year==i,c(fac.var,j)]
        colnames(df) <- c(fac.var,'trait')

        fit1 <- mmer(trait~Rep,
            random= ~ CUSo_Geno + BlockA,
            rcov=~units,
            data=df)
        blups <- as.data.frame(fit1$U$CUSo_Geno)
        blups <- cbind(rownames(blups),blups)
        colnames(blups) <- c('CUSo_Geno',j)
        blups[,1] <- gsub('CUSo_Geno','',blups[,1])

        BLUPs <- BLUPs %>% left_join(blups, by='CUSo_Geno')
        varcomp <- rbind(varcomp, summary(fit1)$varcomp)
        rownames(varcomp) <- gsub('trait',j,rownames(varcomp))
    }
    #write.table(varcomp, file=paste0('/scratch1/ssapkot/HDP/Results/HDP_adjust_varcomp-',i,'.txt'), row.names=T)
    VComp <- rbind(VComp,cbind(varcomp,i))
    GenoBLUPs <- rbind(GenoBLUPs, cbind(BLUPs,i))
}

