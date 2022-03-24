library(tidyverse)


load('analysis2/output/GCA_SCA_mGCA_LxT_Hyb_v_GBLUPs.RData')

blups <- GBLUPs_GCA %>% filter(trait=='GY')

blups <- blups[order(blups$predicted.value, decreasing=T),]

cuso2pi <- read.csv('/zfs/tillers/panicle/HDP_data/CuSo_to_PI.csv', header=F)
colnames(cuso2pi) <- c('Male_CUSo','PI')
blups <- blups %>% left_join(cuso2pi, by='Male_CUSo')
#blups$predicted.value <- blups$predicted.value*0.0673 ## to converts bu/ac to t/ha
#blups$standard.error <- blups$standard.error*0.0673

df <- blups[1:10,]
df$PI <- factor(df$PI, levels=df$PI)

gg <- ggplot(df, aes(PI,predicted.value)) + geom_bar(position=position_dodge(), stat='identity',fill='tan3',width=0.6)
gg <- gg + geom_errorbar(aes(ymin=predicted.value-standard.error, ymax=predicted.value+standard.error), color='black', width=.07, position=position_dodge(0.9))
gg <- gg + theme_minimal() + labs(y='% heterosis')
gg <- gg + theme_light() + theme(
                axis.title.x=element_blank(),
                axis.text.x=element_text(size=15,angle=45,hjust=0.5,vjust=0.5),
                axis.ticks.x=element_blank(),
                axis.title.y=element_text(size=15),
                axis.text.y=element_text(size=13))


ggsave(gg, file='analysis2/results/Barplot_MaleGCA_GBLUPs_heterosis.png', units='in',width=6,height=4,dpi=600)
