
#setwd('/zfs/tillers/panicle/ssapkot/git_repo/GP_SorghumHybrids')

library(tidyverse)
library(reshape2)
library(sommer)

hdp <- read.csv('results/output/Heterosis_BLUPs_byYear.csv', header=T)


fac.var <- colnames(hdp[1:5])
#num.var <- colnames(hdp[22:41])

hdp[,fac.var] <- lapply(hdp[,fac.var],factor)
#hdp[num.var] <- lapply(hdp[num.var],as.numeric)

#traits <- c('DTA','Height','Yield','TGW','Seed_No','Amylose','Starch','Protein','Fat','KCal_lb')
traits <- levels(hdp$Trait)

hyb <- read.table('/scratch1/ssapkot/HDP/data/unphased_snps/hybrid_geno.txt', header=T)

A <- readRDS('/scratch1/ssapkot/HDP/data/unphased_snps/G-matrix_WGS.rds')

hyb_df <- hdp[which(hdp$CUSo_Geno %in% hyb$Female.Male),]
hyb_df$Male_CUSo <- as.factor(as.character(hyb_df$Male_CUSo)) ## to make sure the variables has the right levels otherwise there will be error in running the model with G matrix

VarComp <- c()
GCA_main <- GCA_19 <- GCA_20 <- data.frame(unique(hyb_df$Male_CUSo))
colnames(GCA_main) <- colnames(GCA_19) <- colnames(GCA_20) <- 'CUSo_Geno'

for (i in traits){

    df <- hyb_df[hyb_df$Trait==i,]
    
    E <- diag(length(unique(df$Year))) 
    rownames(E) <- colnames(E) <- unique(df$Year) 
    EA <- kronecker(E,A, make.dimnames = TRUE)

    fit1 <- mmer(MPH~Female_CUSo,
                random=~vs(Male_CUSo) + vs(ds(Year),Male_CUSo),
                rcov=~vs(units),
                data=df)

    varcomp <- summary(fit1)$varcomp
    rownames(varcomp) <- gsub('MPH',i,rownames(varcomp))
    VarComp <- rbind(VarComp, varcomp)

    gca_m <- data.frame(fit1$U$`u:Male_CUSo`)
    gca_m <- cbind(rownames(gca_m),gca_m)
    colnames(gca_m) <- c('CUSo_Geno',i)
    
    gca_19 <- data.frame(fit1$U$`2019:Male_CUSo`)
    gca_19 <- cbind(rownames(gca_19),gca_19)
    colnames(gca_19) <- c('CUSo_Geno',i)

    gca_20 <- data.frame(fit1$U$`2020:Male_CUSo`)
    gca_20 <- cbind(rownames(gca_20),gca_20)
    colnames(gca_20) <- c('CUSo_Geno',i)
    
    GCA_main <- GCA_main %>% left_join(gca_m, by='CUSo_Geno')
    GCA_19 <- GCA_19 %>% left_join(gca_19, by='CUSo_Geno')
    GCA_20 <- GCA_20 %>% left_join(gca_20, by='CUSo_Geno')

}

GCA_males <- list(GCA_main,GCA_19,GCA_20)
save(VarComp, GCA_males, file='results/output/VarComp_GCA_MPH.RData')
