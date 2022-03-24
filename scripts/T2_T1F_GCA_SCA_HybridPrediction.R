library(reshape2)
library(tidyverse)
library(sommer)

### read in phenotypic data
hdp <- read.csv('/zfs/tillers/panicle/ssapkot/git_repo/GP_Sorghumhybrids/data/Phenotypes/HDP_Heterosis_ByBlock.csv', header=T)
fac.var <- colnames(hdp[1:12])
hdp[,fac.var] <- lapply(hdp[,fac.var],factor)
hdp <- hdp[,c(colnames(hdp)[1:12],'Hyb_v')] ## select only the column with values to be modeled (13 = hybrid values)

traits <- c('DTA','PH','TGW','GY','GNP')

### get the list of male and female lines

males <- read.table('/zfs/tillers/panicle/ssapkot/HDP/data/Genomic/male_geno.txt',header=F)
males <- males$V1
females <- read.table('/zfs/tillers/panicle/ssapkot/HDP/data/Genomic/female_geno.txt',header=F)
females <- females$V1

## subset genomic matrices for available males and female lines (for GBLUP)
load('/zfs/tillers/panicle/ssapkot/git_repo/GP_Sorghumhybrids/data/GRM_phased_thinned.RData')
geno_id <- rownames(A)

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

####### subset phenotypic data for lines with genomic data
pheno <- hdp[hdp$CUSo_Geno %in% geno_id,]
pheno$CUSo_Geno <- as.factor(as.character(gsub('\\.',':',pheno$CUSo_Geno))) ##to change values in levels(pheno$CUSo_Geno)
pheno$Male_CUSo <- as.factor(as.character(pheno$Male_CUSo))
pheno$Female_CUSo <- as.factor(as.character(pheno$Female_CUSo))

geno_id <- gsub('\\.',':',geno_id)
## load BLUPs and reliability measures
Reliability <- read.csv('analysis2/Reliability_GCA_SCA_LxT_YieldComp.csv', header=T)
load('analysis2/GCA_SCA_models_LxT_hybridValue_BLUPs.RData')

Pred_value <- c()

for (i in 123:222){
	set.seed(i)
	##### Sampling Plan for T2 and T1F #############

	T1F_Tst_males <- sample(males,60,replace=FALSE)
	T1F_Tst_geno <- c(paste0(females[1],':',T1F_Tst_males),paste0(females[2],':',T1F_Tst_males))
	
#	T2_Tst_males <- sample(males,120,replace=FALSE)
#	T2_Tst1_males <- sample(T2_Tst_males,60,replace=FALSE)
#	T2_Tst2_males <- setdiff(T2_Tst_males,T2_Tst1_males)
	
#	T2_Tst_geno <- c(paste0(females[1],':',T2_Tst1_males),paste0(females[2],':',T2_Tst2_males))
	###############################################
	
#	test <- T2_Tst_geno
	test <- T1F_Tst_geno
	train <- setdiff(geno_id,test)
	for (j in 1:length(traits)){
		trait <- traits[j]
		df <- pheno[pheno$Trait==trait,c(colnames(pheno)[c(1:12)],'Hyb_v')]
	        ## make cross validation folds
		rel_gca <- Reliability$Rel_r[Reliability$Trait==trait&Reliability$Model=='GCA'&Reliability$Value=='hyb']
                rel_sca <- Reliability$Rel_r[Reliability$Trait==trait&Reliability$Model=='SCA'&Reliability$Value=='hyb']
		blup_gca <- BLUPs_GCA[BLUPs_GCA$Trait==trait,]
		blup_sca <- BLUPs_SCA[BLUPs_SCA$Trait==trait,] 	
		pred <- c()
		ability <- c()
		acc <- c()
	
	        #Make training (TRN) and testing (TST) dfs
		
		df[which(df$CUSo_Geno %in% test), "Hyb_v"] <- NA
	        
		## gblup model for gca only model
	  	tryCatch({
		ans <- mmer(Hyb_v~Year + Female_CUSo,
		random= ~ vs(Male_CUSo,Gu=Gm)  + vs(Year:Male_CUSo,Gu=EGm) + vs(Year:Rep:BlockA),
		rcov=~units,
		data=df)
		
		u1 <- ans$U$`u:Male_CUSo`$Hyb_v
       		u2 <- c(ans$Beta[3,'Estimate'],0)
		uf <- gsub('Female_CUSo','',ans$Beta[3,'Effect'])
		names(u2) <- c(uf,setdiff(females,uf))
       		beta <- ans$Beta[1,"Estimate"]
       	 	blups <- data.frame(expand.grid(names(u2),names(u1)),expand.grid(u2,u1))
       	 	colnames(blups)[1:4] <- c('Female','Male','Uf','Um')
       	 	blups$Geno <- paste0(blups$Female,':',blups$Male)
       	 	blups <- blups[,c(ncol(blups),1:4)]
       	 	blups$Beta <- beta
       	 	blups$hyb <- blups$Um + blups$Uf + blups$Beta
		ability <- cor(blups$hyb[which(blups$Geno %in% test)],blup_gca$hyb[which(blup_gca$Geno %in% test)], use='complete')
		acc <- ability/rel_gca		
		pred <- cbind(i,trait,'T1F','GCA',ability,acc)
		Pred_value  <- rbind(Pred_value,pred)

		## gblup model for gca+sca model
		ans <- mmer(Hyb_v~Year + Female_CUSo,
                random= ~ vs(Male_CUSo,Gu=Gm)   + vs(Year:Male_CUSo,Gu=EGm) + vs(Female_CUSo:Male_CUSo,Gu=Gfm) + vs(Year:Female_CUSo:Male_CUSo,Gu=EGfm) + vs(Year:Rep:BlockA),
                rcov=~units,
                data=df)

                u1 <- ans$U$`u:Male_CUSo`$Hyb_v
        	u2 <- c(ans$Beta[3,'Estimate'],0)
                uf <- gsub('Female_CUSo','',ans$Beta[3,'Effect'])
                names(u2) <- c(uf,setdiff(females,uf))
        	sca <- data.frame(ans$U$`u:Female_CUSo:Male_CUSo`$Hyb_v)
        	sca$Geno <- rownames(sca)
        	colnames(sca)[1] <- 'Usca'
        	beta <- ans$Beta[1,"Estimate"]
        	blups <- data.frame(expand.grid(names(u2),names(u1)),expand.grid(u2,u1))
        	colnames(blups)[1:4] <- c('Female','Male','Uf','Um')
        	blups$Geno <- paste0(blups$Female,':',blups$Male)
        	blups <- blups[,c(ncol(blups),1:4)]
        	blups <- blups %>% left_join(sca, by='Geno')
        	blups$Beta <- beta
        	blups$hyb <- blups$Um + blups$Uf + blups$Usca + blups$Beta
                ability <- cor(blups$hyb[which(blups$Geno %in% test)],blup_sca$hyb[which(blup_sca$Geno %in% test)], use='complete')
		acc <- ability/rel_sca
		pred <- cbind(i,trait,'T1F','GCA+SCA',ability,acc)
                Pred_value  <- rbind(Pred_value,pred)
		}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})	
	
	}
}

Pred_value <- data.frame(Pred_value)
colnames(Pred_value) <- c('Seed','Trait','Method','Model','Predictive_ability','Prediction_accuracy')
Pred_value$Predictive_ability <- as.numeric(Pred_value$Predictive_ability)
Pred_value$Prediction_accuracy <- as.numeric(Pred_value$Prediction_accuracy)

write.csv(Pred_value, file='analysis2/GCA_SCA_LxT_hybridValue_T1F.csv', row.names=F,quote=F)
