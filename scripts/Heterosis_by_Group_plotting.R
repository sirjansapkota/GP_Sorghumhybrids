
library(tidyverse)

setwd('zfs/tillers/panicle/ssapkot/git_repo/GP_Sorghumhybrids')
load('data/HDP_data_HetSummary_GRM.RData')

race_df <- data.frame(rbind(c(1,'Durra'),c(2,'Kafir'),c(3,'Caudatum'),c(4,'Guinea'),c(5,'Milo'),c(6,'Mixed')))
colnames(race_df) <- c('K.Cluster','SubPop')
race_df$K.Cluster <- as.integer(race_df$K.Cluster)


hdp <- hdp %>% left_join(race_df, by='K.Cluster')
hdp$SubPop <- as.factor(hdp$SubPop)
hdp <- hdp %>% filter(!SubPop=='NA')

## function to remove outliers
calc_boxplot_stat <- function(x) {
  coef <- 1.5
  n <- sum(!is.na(x))
  # calculate quantiles
  stats <- quantile(x, probs = c(0.0, 0.25, 0.5, 0.75, 1.0))
  names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
  iqr <- diff(stats[c(2, 4)])
  # set whiskers
  outliers <- x < (stats[2] - coef * iqr) | x > (stats[4] + coef * iqr)
  if (any(outliers)) {
    stats[c(1, 5)] <- range(c(stats[2:4], x[!outliers]), na.rm = TRUE)
  }
  return(stats)
}

##### plotting boxplot for heterotic values by subpopulation
gg <- hdp %>% 
    group_by(Trait, SubPop) %>% 
    ggplot(.,aes(Female,MPH,fill=SubPop)) + stat_summary(fun.data=calc_boxplot_stat, geom='boxplot', position=position_dodge()) +
  facet_wrap(~Trait, scale='free') + theme_bw() + scale_fill_brewer(palette = 'Dark2')

gg <- gg + theme(
  axis.title.x=element_blank(),
  axis.text.x=element_text(size=15),
  axis.ticks.x=element_blank(),
  axis.title.y=element_text(size=20),
  axis.text.y=element_text(size=15),
  strip.text.x=element_text(size=20),
  panel.background=element_blank(),
  legend.title=element_blank(),
  legend.text=element_text(size=15),
  legend.position = c(0.8,0.3),
  panel.spacing.x=unit(2,"line"))
gg <- gg + labs(x='',y=' % mid-parent heterosis')

ggsave(gg, file='analysis2/results/Boxplots_MPH_byRace.png', units='in', width=8, height=5, dpi=600)

###### plotting barplot for mean and standard deviation of heterotic values
hdp %>% 
  filter(Female=='ATx2928') %>%
  group_by(Trait, SubPop) %>%
  summarize(Mean=mean(MPH, na.rm=T), SD= sd(MPH, na.rm=T)) %>%
  ggplot(., aes(SubPop,MPH,fill=SubPop)) + 
  geom_bar(stat='identity', position=position_dodge(0.9)) +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), colour = "black", width=.07, position = position_dodge(0.9)) +
  facet_wrap(~Trait, scale='free') + 
  theme_bw() + scale_fill_brewer(palette = 'Dark2')

########################################################################

## Calculate genetic distance between female parent and male parents and run correlation with heterotic values of hybrid combinations

##read in genotypic data
library(data.table)
library(sommer)
library(ggpubr)

M <- fread('/zfs/tillers/panicle/ssapkot/HDP/data/Genomic/SAP_WGS_filtered.012')
ind <- read.table('/zfs/tillers/panicle/ssapkot/HDP/data/Genomic/SAP_WGS_filtered.012.indv',header=F)

M <- as.matrix(M[,-1])-1
rownames(M) <- ind$V1

D.euc <- as.matrix(dist(M))
saveRDS(D.euc, file='data/EuclideanDist_parents.rds')
D.euc <- readRDS('data/EuclideanDist_parents.rds')
D.euc <- scale(D.euc)

dist_df <- c()

### repeat for each female parent
atx <- data.frame(D.euc[,colnames(D.euc)=='CUSo09108'])
df <- data.frame(cbind(paste0('CUSo09108.',rownames(atx)),atx))
rownames(df) <- NULL
colnames(df) <- c('CUSo_Geno','Dist')
dist_df <- rbind(dist_df,df)

## merge the genetic distance metric for each hybrid to hdp data frame
hdp <- hdp %>% left_join(dist_df, by='CUSo_Geno')

### calculate correlation between parental genetic distance and mid parent heterosis by Female, Trait, Year, and Rep

Corr_het_GD <- hdp %>% group_by(Female,Year,Rep,Trait) %>% summarize(Corr=cor(MPH,Dist,use='complete'))
Corr_het_GD %>% ggplot(.,aes(Female,Corr,color=Trait)) + geom_point() + theme_classic()
hdp %>% group_by(Female,Trait) %>% summarize(Corr = cor(Dist, MPH, use='complete', method='pearson'))

## group wise correlation test
Corr_het_GD <- hdp %>% group_by(Female,Year,Rep,Trait) %>% summarize(p.value = cor.test(Dist,MPH)$p.value,
                                             estimate = cor.test(Dist,MPH)$estimate)

### plot
gg <- hdp %>% group_by(Female,Trait) %>% ggplot(.,aes(Dist,MPH,color=Female)) 
gg <- gg + geom_point(size=1.2) + geom_smooth(method='lm') + facet_wrap(~Trait, scale='free') 
gg <- gg + theme_classic() + scale_color_manual(values=c("#E69F00", "#56B4E9"))
gg <- gg + theme(
  axis.title.x=element_text(size=20),
  axis.text.x=element_text(size=15),
  axis.title.y=element_text(size=20),
  axis.text.y=element_text(size=15),
  strip.text.x=element_text(size=20),
  panel.background=element_blank(),
  legend.title=element_blank(),
  legend.text=element_text(size=15),
  legend.position = c(0.8,0.3),
  panel.spacing.x=unit(2,"line"))
gg <- gg + labs(x='genetic distance between parents',y=' % mid-parent heterosis')

ggsave(gg, file='analysis2/results/Correlation_GenDist_MPH_ByFemale.png',units='in',height=5,width=8,dpi=600)

### save all
save.image(file='data/HDP_data_HetSummary_GRM.RData')
