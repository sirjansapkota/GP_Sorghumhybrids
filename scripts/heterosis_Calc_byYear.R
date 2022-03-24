library(tidyverse)
library(ggplot2)
library(reshape2)

hdp <- read.csv('/zfs/tillers/panicle/ssapkot/git_repo/GP_Sorghumhybrids/results/output/BLUPs_FieldEffectsAdj_all.csv', header=T)

traits <- traits <- c('DTA','Height','Yield','TGW','Seed_No','Amylose','Starch','Protein','Fat','KCal_lb')

##subset by year and rep and then run loop by traits to calculate heterosis
Female <- c('CUSo09108','CUSo07750')
year <- c('2019','2020')

het_all <- c()

for (a in year){
    hdp_x <- hdp %>% filter(Year==a)
    df.all <- c()
    for (i in 1:length(traits)){
      df.het <- c()
      trait <- traits[i]
       print(i)
      for (j in 1:2){
        female <- Female[j]
        df <- hdp_x[,c('CUSo_Geno','Female_CUSo','Male_CUSo','Year',trait)]
        colnames(df) <- c('CUSo_Geno','Female_CUSo','Male_CUSo','Year','Hyb_v')
        df.hy <- df[nchar(df$CUSo_Geno)>15,] %>% filter(Female_CUSo==female)
        nonmales <- c('CUSo09108','CUSo07750','FILL','83P17')
        df.males <- df[!df$CUSo_Geno %in% nonmales,] ## remove females, checks and fills
        df.males <- df.males[nchar(df.males$CUSo_Geno)<15,] ## remove all the hybrids from df

        df.female <- df[df$CUSo_Geno==female,]
        df.check <- df[df$CUSo_Geno=='83P17',]
        df.hy$Trait <- trait
        df.hy$Male_v <- df.males$Hyb_v[match(df.hy$Male_CUSo,df.males$CUSo_Geno)]
        df.hy$Female_v <- df.female$Hyb_v
        df.hy$Check_v <- df.check$Hyb_v
        df.hy$MidPar_v <- (df.hy$Male_v + df.hy$Female_v)/2
        df.hy$HighPar_v <- apply(df.hy[,c('Male_v','Female_v')], 1, max)
        df.hy$FemalePH <- 100 * (df.hy$Hyb_v - df.hy$Female_v)/df.hy$Female_v
        df.hy$MalePH <- 100 * (df.hy$Hyb_v - df.hy$Male_v)/df.hy$Male_v
        df.hy$MPH <- 100 *(df.hy$Hyb_v - df.hy$MidPar_v)/df.hy$MidPar_v
        df.hy$BPH <- 100 * (df.hy$Hyb_v - df.hy$HighPar_v)/df.hy$HighPar_v
        df.hy$OverCheck <- 100 * (df.hy$Hyb_v - df.hy$Check_v)/df.hy$Check_v

        df.het <- rbind(df.het,df.hy)
        }

     df.all <- rbind(df.all,df.het)

    }
    het_all <- rbind(het_all, df.all)
}

het_all <- het_all[,c(1:4,6,5,7:ncol(het_all))]

het_all <- het_all %>% mutate_if(is.numeric, round, digits=3)

write.csv(het_all, file='results/output/Heterosis_BLUPs_byYear.csv', row.names=FALSE, quote=FALSE)
