library(tidyverse)

df_cv1 <- read.csv('manuscript_GP1/output/Accuracy_n300_CV1_GBLUP_MPH_hybrids.csv', header=T)
df_642 <- read.csv('manuscript_GP1/output/Accuracy_n300_CV2_CV3_CUSo07750_MPH.csv', header=T)
df_2928 <- read.csv('manuscript_GP1/output/Accuracy_n300_CV2_CV3_CUSo09108_MPH.csv', header=T)

df_2928$Method <- paste0(df_2928$Method,'_Tx2928')
df_642$Method <- paste0(df_642$Method,'_Tx642')

df <- df_cv1 ## first make output dataframe for cv1 alone and then combine CV2-CV3 later
## make a combined dataframe for the two prediction results
df_642$Fparent <- 'ATx642'
df_2928$Fparent <- 'ATx2928'
df <- rbind(df_642,df_2928)

df[,'Trait'] <- gsub('Height','PH',df[,'Trait'])
df[,'Trait'] <- gsub('Seed_No','GNP',df[,'Trait'])
df[,'Trait'] <- gsub('Yield','GY',df[,'Trait'])

## aggregate stats and summary into tables
output <- df %>% group_by(Trait) %>% summarize( CV1=paste0(round(mean(Prediction_accuracy),3),' (',round(sd(Prediction_accuracy),3),')'))

output2 <- df %>% group_by(Trait,Method) %>% summarize( value=paste0(round(mean(Prediction_accuracy),3),' (',round(sd(Prediction_accuracy),3),')'))
output2 <- spread(output2,Method,value)

sum_df <- output %>% left_join(output2, by="Trait")

write.csv(sum_df, 'manuscript_GP1/output/Mean_SD_CV1-CV3_MPH.csv', quote=F, row.names=F)
## mean comparison for phenotypic data

traits <- unique(df$Trait)

output <- c()

for (i in traits){
	data = df[df$Trait==i,]
	test = t.test(MPH~Female, data=data)
        test_r = cbind(i,round(test$estimate[1],3),round(test$estimate[2],3),test$p.value)
        output <- rbind(output,test_r)
}
output <- as.data.frame(output)
colnames(output) <- c('Trait','ATx2928_mean','ATx642_mean','p_value')
rownames(output) <- NULL
output



## mean comparison between results from prediction of two females and two methods in factorial combination
traits <- unique(df$Trait)
methods <- unique(df$Method)
parents <- unique(df$Fparent)

output_parents <- c()
for (i in traits){
	data = df[df$Trait==i,]
	
	for (j in methods){
		data2 = data[data$Method==j,]
		test = t.test(Prediction_accuracy~Fparent, data=data2)
		test_r = cbind(i,j,round(test$estimate[1],3),round(test$estimate[2],3),test$p.value)
		output_parents <- rbind(output_parents,test_r)
	}
}

output_parents <- as.data.frame(output_parents)
colnames(output_parents) <- c('Trait','Method','ATx2928_mean','ATx642_mean','p_value')
rownames(output_parents) <- NULL

traits <- unique(df$Trait)
models <- unique(df$Model)

output_method <- c()
for (i in traits){
        data = df[df$Trait==i,]
	for (j in models){
		data2 <- data[data$Model==j,]
        	test = t.test(Prediction_accuracy~Method, data=data2)
        	test_r = cbind(i,j,round(test$estimate[1],3),round(test$estimate[2],3),test$p.value)
        	output_method <- rbind(output_method,test_r)
        }
}

output_method <- as.data.frame(output_method)
colnames(output_method) <- c('Trait','Model','T1F','T2','p_value')
rownames(output_method) <- NULL
output_method

method <- unique(df$Method)
output_model <- c()
for (i in traits){
        data = df[df$Trait==i,]
        for (j in method){
                data2 <- data[data$Method==j,]
                test = t.test(Prediction_accuracy~Model, data=data2)
                test_r = cbind(i,j,round(test$estimate[1],3),round(test$estimate[2],3),test$p.value)
                output_model <- rbind(output_model,test_r)
        }
}

output_model <- as.data.frame(output_model)
colnames(output_model) <- c('Trait','Method','GCA','GCA+SCA','p_value')
rownames(output_model) <- NULL
output_model


