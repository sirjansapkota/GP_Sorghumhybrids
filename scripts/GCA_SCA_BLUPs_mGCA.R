library(reshape2)
library(tidyverse)
library(sommer)

### read in phenotypic data
hdp <- read.csv('/zfs/tillers/panicle/ssapkot/git_repo/GP_Sorghumhybrids/data/Phenotypes/HDP_Heterosis_ByBlock.csv', header=T)
fac.var <- colnames(hdp[1:12])
hdp[,fac.var] <- lapply(hdp[,fac.var],factor)
hdp <- hdp[,c(colnames(hdp)[1:12],'Hyb_v')] ## select only the column with values to be modeled (13 = hybrid values)

traits <- c('DTA','PH','TGW','GY','GNP')

### For GBLUP hybrids
load('/zfs/tillers/panicle/ssapkot/git_repo/GP_Sorghumhybrids/data/GRM_phased_thinned.RData')
#rownames(A) <- colnames(A) <- gsub(":",".",rownames(A)) ##replace : with . to keep consistent with phenotype file
#rownames(D) <- colnames(D) <- gsub(":",".",rownames(D)) ##replace : with . to keep consistent with phenotype file
geno_id <- rownames(A)

### get the list of male and female lines
males <- read.table('/zfs/tillers/panicle/ssapkot/HDP/data/Genomic/male_geno.txt',header=F)
males <- males$V1
females <- read.table('/zfs/tillers/panicle/ssapkot/HDP/data/Genomic/female_geno.txt',header=F)
females <- females$V1

Gm <- G[order(males),order(males)]
rownames(Gm) <- colnames(Gm) <- males
Gf <- G[order(females),order(females)]
rownames(Gf) <- colnames(Gf) <- females
Gfm <- kronecker(Gf,Gm, make.dimnames=TRUE)
#### make kronecker product for geno x year covariance matrix
E <- diag(length(unique(hdp$Year)))
rownames(E) <- colnames(E) <- unique(hdp$Year)
EGm <- kronecker(E,Gm, make.dimnames = TRUE)
EGf <- kronecker(E,Gf, make.dimnames = TRUE)
EGfm <- kronecker(E,Gfm, make.dimnames=TRUE)

## subset phenotypic data to include only the lines with genomic data (for GBLUP)
pheno <- hdp[hdp$CUSo_Geno %in% geno_id,]
pheno$CUSo_Geno <- as.factor(as.character(pheno$CUSo_Geno)) ##to change values in levels(pheno$CUSo_Geno)
pheno$Male_CUSo <- as.factor(as.character(pheno$Male_CUSo))
pheno$Female_CUSo <- as.factor(as.character(pheno$Female_CUSo))

varcomp_gca <- c()
varcomp_sca <- c()
GBLUPs_GCA <- c()
GBLUPs_SCA <- c()
model_gca <- c()
model_sca <- c()

###### For combining abilities without genomic data gca , gca + sca model
for (j in traits){
	blups <- c()

	## for gca only model
	df_a <- pheno[pheno$Trait==j,c(1:13)]
	
	fit_a <- mmer(Hyb_v~Year + Female_CUSo,
#                random= ~ vs(Male_CUSo) + vs(Year:Male_CUSo) + vs(Year:Rep:BlockA),
		random= ~ vs(Male_CUSo,Gu=Gm)  + vs(Year:Male_CUSo,Gu=EGm) + vs(Year:Rep:BlockA),
                rcov=~units,
                data=df_a)
	
	#fit_a <- mkConv(fit_a)        ##for asreml models to converge if not converged already

	blups <- predict(fit_a, classify='Male_CUSo')$pval
	blups$trait <- gsub('Hyb_v',j,blups$trait)
	blups$pev <- blups[,'standard.error']^2
	Vgca <- fit_a$sigma$`u:Male_CUSo`[1]
	blups$rel <- 1-(blups$pev/Vgca)
	GBLUPs_GCA <- bind_rows(GBLUPs_GCA,blups)

	#### variance component analysis
	suma <- summary(fit_a)$varcomp
	
	varcomp_gca <- bind_rows(varcomp_gca, suma)
	rownames(varcomp_gca) <- gsub('Hyb_v',j,rownames(varcomp_gca))
	
	par <- cbind(j,fit_a$AIC,fit_a$BIC)
	model_gca <- rbind(model_gca,par)	

	## for gca + sca models
	blups <- c()
	fit_b <- mmer(Hyb_v~ Year + Female_CUSo,
#        random= ~ vs(Male_CUSo)  + vs(Year:Male_CUSo) + vs(Female_CUSo:Male_CUSo) + vs(Year:Female_CUSo:Male_CUSo) + vs(Year:Rep:BlockA),
	random= ~ vs(Male_CUSo,Gu=Gm) + vs(Year:Male_CUSo,Gu=EGm) + vs(Female_CUSo:Male_CUSo,Gu=Gfm) + vs(Year:Female_CUSo:Male_CUSo,Gu=EGfm) + vs(Year:Rep:BlockA),        
	rcov=~units,
        data=df_a)
	
	blups <- predict(fit_b, classify='Male_CUSo')$pval
        blups$trait <- gsub('Hyb_v',j,blups$trait)
        blups$pev <- blups[,'standard.error']^2
        Vgca <- fit_b$sigma$`u:Male_CUSo`[1]
        blups$rel <- 1-(blups$pev/Vgca)
        GBLUPs_SCA<- bind_rows(GBLUPs_SCA,blups)

        #### variance component analysis
        suma <- summary(fit_b)$varcomp

        varcomp_sca <- bind_rows(varcomp_sca, suma)
        rownames(varcomp_sca) <- gsub('Hyb_v',j,rownames(varcomp_sca))

	par <- cbind(j,fit_b$AIC,fit_b$BIC)
        model_sca <- rbind(model_sca,par)
}

colnames(model_gca) <- colnames(model_sca) <- c('Trait','AIC','BIC')
## making a list of stuff to save
keep_list <- c('varcomp_sca','varcomp_gca','GBLUPs_GCA','GBLUPs_SCA','model_gca','model_sca')
del_list <- setdiff(ls(),keep_list)
rm(list=del_list)

save(list=ls(), file='analysis2/output/GCA_SCA_mGCA_LxT_Hyb_v_GBLUPs.RData')
