
install.packages("ggcorrplot")

library(tidyverse)
library(ggcorrplot)

#### plotting boxplots by category
# geno is subsetted to include all cuso hybrids (SAP with genomic data) and the check and the parental lines of the hybrids
pheno <- read.csv('data/Phenotypes/HDP_compiled_Phenotypes_spread.csv', header=T)

pheno <- pheno %>% filter(CUSo_Geno %in% geno)
pheno <- pheno[,c(5,13,22:24,29:30)]
colnames(pheno)[4:7] <- c('GNP','GY','DTA','PH')
df.m <- melt(pheno)

df.m$Category <- gsub('Hybrid','F1',df.m$Category)
df.m$Category <- gsub('Female','P1',df.m$Category)
df.m$Category <- gsub('Male','P2',df.m$Category)


library(ggbeeswarm)
library(ggpubr)

gg <- ggplot(df.m, aes(Category, value, color=Category)) + geom_boxplot(position=position_dodge(), width=0.9) + facet_wrap(~variable, scale='free') + theme_classic() + scale_color_manual(values=c("#E69F00", "#56B4E9","#000000","#800080"))
gg <- gg + labs(x='',y='') + geom_quasirandom(size=0.3, shape='.')

## stat comparison
my_comparisons <- list(c('Check','F1'),c('F1','P1'), c('F1','P2'),c('P1','P2'))

gg <- gg + stat_compare_means(comparisons=my_comparisons,label="p.signif", vjust=0.3)
####

gg <- gg + theme(
                axis.title.x=element_blank(),
                axis.text.x=element_text(size=15),
#                axis.ticks.x=element_blank(),
                axis.title.y=element_text(size=20),
                axis.text.y=element_text(size=15),
                strip.text.x=element_text(size=20),
#                strip.text.y=element_text(size=20),
                panel.background=element_blank(),
                legend.title=element_blank(),
                legend.text=element_text(size=15),
		legend.position = "none",
                panel.spacing.x=unit(2,"line"))

ggsave(gg, file='analysis2/results/Boxplots_BLUPs_parents_hybrid.png', units='in', width=12, height=6, dpi=600)














### processing BLUPs for correlation plot

load('analysis/BLUPs_all_hybrids_geneticValue.RData')
head(BLUPs)

BLUPs <- BLUPs %>% filter(CUSo_Geno %in% geno)
blups <- BLUPs[,1:3] %>% spread(.,trait,predicted.value)
blups <- blups[,c(1,3,5,8,10,11)]
colnames(blups) <- c('CUSo_Geno','DTA','PH','GNP','TGW','GY')
traits <- colnames(blups)[-1]

M <- cor(blups[,-1])


gg <- ggcorrplot(M, hc.order = FALSE, type = "lower", lab = TRUE, lab_size=6)
gg <- gg   + theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.title.x=element_blank(),
  axis.text.x=element_text(size=17),
#  axis.ticks.x=element_blank(),
  axis.title.y=element_blank(),
  axis.text.y=element_text(size=17),
  legend.text=element_text(size=17),
  legend.title=element_blank(),
  legend.position='none'
)
ggsave(gg, file='analysis2/results/Correlation_blups_HybValue_onlyCUSo_traits.png',units='in',height=5,width=5,dpi=600)


###### heterosis plot

hdp <- read.csv('data/Phenotypes/HDP_Heterosis_ByBlock.csv', header=T)
hdp <- hdp %>% filter(Trait %in% traits)
df <- hdp[,c(1,2,4,7,8,12,21)]

gg <- ggplot(df, aes(x=Year,y=MPH, fill=Year)) + geom_violin(trim=TRUE) 
gg <- gg + facet_grid(vars(Female), vars(Trait), scale="free")
gg <- gg + geom_boxplot(width=0.1) + scale_fill_manual(values=c("#E69F00", "#56B4E9"))
gg <- gg + theme(
		axis.title.x=element_blank(),
                axis.text.x=element_text(size=15),
                axis.ticks.x=element_blank(),
                axis.title.y=element_text(size=20),
                axis.text.y=element_text(size=15),
                strip.text.x=element_text(size=20),
                strip.text.y=element_text(size=20),
                panel.background=element_blank(),
                legend.title=element_blank(),
                legend.position = "none",
                panel.spacing.x=unit(0,"line"))
gg <- gg + ylab("Mid-parent heterosis (%)")
ggsave(gg, file='analysis2/results/Heterosis_byYear_byFemale.png', units='in',width=10,height=6,dpi=300)


### prediction results
t2 <- read.csv('analysis2/output/GCA_SCA_LxT_heterosis_T2.csv', header=T)
t1f <- read.csv('analysis2/output/GCA_SCA_LxT_heterosis_T1F.csv', header=T)
data <- rbind(t2,t1f)
data$Method2 <- paste0(data$Method,'_',data$Model)

gg <- ggplot(data, aes(Trait,Prediction_accuracy, fill=Model)) + geom_violin(position = position_dodge(0.9), trim=FALSE)
gg <- gg + facet_wrap(~Type)
gg <- gg + geom_boxplot(width=0.2,position = position_dodge(0.9)) + scale_fill_manual(values=c("#E69F00", "#56B4E9"))
gg <- gg + theme_light() + theme(
                axis.title.x=element_blank(),
                axis.text.x=element_text(size=15,angle=45,hjust=0.5,vjust=0.7),
                axis.ticks.x=element_blank(),
                axis.title.y=element_text(size=15),
                axis.text.y=element_text(size=13),
		strip.text.x=element_text(size=15),
                legend.title=element_blank(),
		legend.text=element_text(size=15),
                legend.position = "top",
                panel.spacing.x=unit(2,"line"))
gg <- gg + ylab("Prediction accuracy") 
ggsave(gg, file='analysis2/results/Prediction_mGCA_all.png', units='in',width=6,height=4,dpi=600)



## split violin for Reps (see functions in "make_split_violin.R"
gg = ggplot(df, aes(Year,MPH,fill=Rep)) + geom_split_violin() + facet_grid(vars(Female), vars(Trait), scale="free") +  stat_summary(fun.data=data_summary)
 gg = gg + scale_fill_brewer(palette="Dark2") + theme_minimal() + theme(axis.text.x = element_text(angle = 30),legend.position=c(0.87,0.2)) +labs(x='',y='')
 gg = gg +  theme(strip.background = element_blank(), strip.placement = "outside",panel.spacing=unit(20,"pt"))
 gg

