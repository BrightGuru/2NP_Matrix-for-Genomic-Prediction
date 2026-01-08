#Packages
#This is a list of all packages used in the analysis.
library(data.table)
library(magrittr)
library(BGLR)
library(lme4)
library(foreach)
library(doParallel)
library(doFuture)
## Loading required package: future
library(doMC)
library(naniar)
library(ggpubr)
## Loading required package: ggplot2
library(lme4)
library(plyr)

library(dplyr)
library(tidyverse)

BLUES_Year=list()
Years = c("2018","2019","2020","2021","2022","2023")
trait <- commandArgs(trailingOnly = T)[1]

##BLues for Hybrid for each year
Pheno_data_without_Outlier<- Pheno_data_without_Outlier %>%
  mutate(Hybrid_Field = paste(Hybrid, "_",Field.Location)) 

unique_hybrid_Field <- as.character(unique(Pheno_data_without_Outlier$Hybrid_Field))

BLUES_All <- BLUES_All %>%
  mutate(Hybrid_Field = paste(Hybrid, "_",Field.Location))%>% 
  dplyr::select(Hybrid_Field,Hybrid,Field.Location, Year,Env,traits)%>% 
  filter(Hybrid_Field %in% unique_hybrid_Field)

for (k in 1:length(Years)){
  
  print(paste0(' Year:', Years[k]))
  data = filter(BLUES_All, Year %in% Years[k]) 
  
  data$Hybrid<- as.factor(data$Hybrid)
  
  data$Field.Location <- as.factor(data$Field.Location)
  
  data[[traits]] <- as.numeric(data[[traits]])
  
  BLUES_FUN1 <- function(traits, data= ".") {
    Model <- as.data.frame(fixef (lmer(paste0(traits, "~0 + Hybrid + (1|Field.Location)")
                                       , data=data)))
  }
  
  BLUES_F <- lapply(traits, BLUES_FUN1, data=data)
  BLUES_All_F <- map2(BLUES_F, traits, ~set_names(..1,..2) %>% rownames_to_column(var = "Hybrid"))%>% reduce(full_join)
  
  mean_df = data.frame(BLUES_All_F, Years[k])
  
  colnames(mean_df) =c("Hybrid",traits, "Year" )
  
  BLUES_Year[[k]] = cbind(mean_df)
}

file_name2 <- paste0(trait,"_BLUES_Year", ".csv")

BLUES_Year=ldply(BLUES_Year, rbind)

write.table(BLUES_Year,file = file_name2,col.names = T,row.names = F,sep = '\t',quote = F)

unique(Yield_Mg_ha_BLUES_Year$Hybrid)
summary(as.factor(Yield_Mg_ha_BLUES_Year$Year))
