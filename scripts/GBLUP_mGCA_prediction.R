library(reshape2)
library(tidyverse)
library(sommer)

### read in phenotypic blups and genomic matrices
hdp <- read.csv('/zfs/tillers/panicle/ssapkot/git_repo/GP_Sorghumhybrids/data/Phenotypes/HDP_Heterosis_ByBlock.csv', header=T)
fac.var <- colnames(hdp[1:12])
hdp[,fac.var] <- lapply(hdp[,fac.var],factor)
hdp <- hdp[,c(colnames(hdp)[1:12],'Hyb_v')]


traits <- c('DTA','PH','TGW','GY','GNP')	

load('/zfs/tillers/panicle/ssapkot/git_repo/GP_Sorghumhybrids/data/GRM_phased_thinned.RData')

hyb_geno <- rownames(A)
pheno <- hdp[hdp$CUSo_Geno %in% hyb_geno,]
pheno$CUSo_Geno <- as.factor(as.character(gsub('\\.',':',pheno$CUSo_Geno))) ##to change values in levels(pheno$CUSo_Geno)
pheno$Male_CUSo <- as.factor(as.character(pheno$Male_CUSo))
pheno$Female_CUSo <- as.factor(as.character(pheno$Female_CUSo))
hyb_geno <- gsub('\\.',':',hyb_geno)


### /BLUPs_all_hybrids_hybridvalues.RData')get the list of male and female lines
males <- read.table('/zfs/tillers/panicle/ssapkot/HDP/data/Genomic/male_geno.txt',header=F)
males <- males$V1
females <- read.table('/zfs/tillers/panicle/ssapkot/HDP/data/Genomic/female_geno.txt',header=F)
females <- females$V1

###  make genomic matrics
Gm <- G[sort(males),sort(males)]
rownames(Gm) <- colnames(Gm) <- males
Gf <- G[sort(females),sort(females)]
rownames(Gf) <- colnames(Gf) <- females
Gfm <- kronecker(Gf,Gm, make.dimnames=TRUE)
#### make kronecker product for geno x year covariance matrix
E <- diag(length(unique(hdp$Year)))
rownames(E) <- colnames(E) <- unique(hdp$Year)
EGm <- kronecker(E,Gm, make.dimnames = TRUE)
EGf <- kronecker(E,Gf, make.dimnames = TRUE)
EGfm <- kronecker(E,Gfm, make.dimnames=TRUE)


## load BLUPs of genotypes and reliability measures to calculate accuracy
load('analysis2/output/GCA_SCA_mGCA_LxT_Hyb_v_BLUPs.RData')
Reliability <- read.csv('analysis2/output/Reliability_mGCA_YieldComp.csv', header=T)

Pred_value <- c()
for (j in 1:length(traits)){
        trait <- traits[j]
        df <- pheno[pheno$Trait==trait,c(colnames(pheno)[c(1:12)],'Hyb_v')]
	rel_gca <- Reliability$Rel_r[Reliability$Trait==trait&Reliability$Model=='GCA'&Reliability$Value=='hyb']
        rel_sca <- Reliability$Rel_r[Reliability$Trait==trait&Reliability$Model=='SCA'&Reliability$Value=='hyb']
        blup_gca <- BLUPs_GCA[BLUPs_GCA$trait==trait,]
        blup_sca <- BLUPs_SCA[BLUPs_SCA$trait==trait,]

	accuracy <- c()	
	for (i in 123:222){
		set.seed(i)
		test <- sample(males,60,replace=FALSE)
	
		test_geno <- c(paste0(females[1],':',test),paste0(females[2],':',test))
	
		pred <- c()
		ability <- c()
		acc <- c()
	
	        #Make training (TRN) and testing (TST) dfs
	  	y.trn <- df
		y.trn[which(y.trn$Male_CUSo %in% test), "Hyb_v"] <- NA
	        
		## gblup model for gca only model
	  	ans <- mmer(Hyb_v~Year + Female_CUSo,
			random= ~ vs(Male_CUSo,Gu=Gm)  + vs(Year:Male_CUSo,Gu=EGm) + vs(Year:Rep:BlockA),
			rcov=~units,
			data=y.trn)
		
		gblups <- predict(ans, classify='Male_CUSo')$pval
		
		ability <- cor(gblups$predicted.value[which(gblups$Male_CUSo %in% test)],blup_gca$predicted.value[which(blup_gca$Male_CUSo %in% test)], use='complete')	
		acc <- ability/rel_gca
		pred <- cbind(i,trait,'GCA',ability,acc)
		accuracy  <- rbind(accuracy,pred)

		## gblup model for gca+sca model
		pred <- c()
		ans <- mmer(Hyb_v~Year + Female_CUSo,
                	random= ~ vs(Male_CUSo, Gu=Gm) + vs(Year:Male_CUSo,Gu=EGm) + vs(Female_CUSo:Male_CUSo,Gu=Gfm) + vs(Year:Female_CUSo:Male_CUSo,Gu=EGfm) + vs(Year:Rep:BlockA),
                	rcov=~units,
                data=y.trn)

		gblups <- predict(ans, classify='Male_CUSo')$pval

                ability <- cor(gblups$predicted.value[which(gblups$Male_CUSo %in% test)],blup_sca$predicted.value[which(blup_sca$Male_CUSo %in% test)], use='complete')
                acc <- ability/rel_sca
                pred <- cbind(i,trait,'GCA_SCA',ability,acc)
                accuracy  <- rbind(accuracy,pred)
		
		}
	    Pred_value <- rbind(Pred_value,accuracy)
	}

Pred_value <- data.frame(Pred_value)
colnames(Pred_value) <- c('Seed','Trait','Model','Predictive_ability','Prediction_accuracy')
Pred_value$Predictive_ability <- as.numeric(Pred_value$Predictive_ability)
Pred_value$Prediction_accuracy <- as.numeric(Pred_value$Prediction_accuracy)
write.csv(Pred_value, file='analysis2/results/GCA_SCA_maleGCA_Prediction_Hyb_v.csv', row.names=F,quote=F)
