
library(tidyverse)

## to get summary of reliability and accuracy

GCA_Rel <- BLUPs_GCA %>% group_by(trait) %>% summarize(Reliability=round(mean(rel),3), Pheno_Accu=round(sqrt(mean(rel)),3))

SCA_Rel <- BLUPs_SCA %>% group_by(trait) %>% summarize(Reliability=round(mean(rel),3), Pheno_Accu=round(sqrt(mean(rel)),3))


### calculation of reliability (correlation between GBLUP and BLUPs of GCA when full ranked)

## for gca only model

blups <- BLUPs_GCA[,c(colnames(BLUPs_GCA)[1:3],'Trait','hyb')]
blups <- blups %>% spread(.,Trait,hyb)

gblups <- GBLUPs_GCA[,c(colnames(GBLUPs_GCA)[1:3],'Trait','hyb')]
gblups <- gblups %>% spread(.,Trait,hyb)

## for male GCA
blups <- BLUPs[,1:3] %>% spread(.,trait,predicted.value)

gblups <- GBLUPs[,1:3] %>% spread(.,trait,predicted.value)


df <- c()

for (i in 2:6){
	corr <- round(cor(blups[,i],gblups[,i]),3)
	df <- rbind(df,corr)
}

df <- data.frame(cbind(colnames(gblups)[-1],df))
colnames(df) <- c('Trait','Rel_r')
df$Model <- 'GCA'
df$Value <- 'het'
#Rel <- c()
Rel <- rbind(Rel,df)

write.csv(Rel, 'analysis2/output/Reliability_mGCA_YieldComp.csv', row.names=F, quote=F)




