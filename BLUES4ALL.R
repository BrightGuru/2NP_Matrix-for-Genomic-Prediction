### Packages 
library(tidyverse)
library(data.table)
library(naniar)
library(ggpubr)
library(lme4)
library(plyr)
library(dplyr)

#Load in the phenotype data

Pheno_data_2014_23 <- fread(file = "/usr/users/osatohanmwen/G2F_competition_2024/Training_data/1_Training_Trait_Data_2014_2023.csv")

Pheno_data_2014_23 <- as_tibble(Pheno_data_2014_23)


Pheno_data_2018_2023 <- filter(Pheno_data_2014_23, Year %in% c("2018", "2019", "2020", "2021","2022","2023"))


## Model Hybrid and replicate for outlier identification and remover using the residuals and

### replace _ with . ####
Pheno_data_2018_2023 <- as_tibble(Pheno_data_2018_2023)
COL_names <- gsub("_", ".", colnames(Pheno_data_2018_2023))

colnames(Pheno_data_2018_2023) <- COL_names

Pheno_data_2018_2023 <- Pheno_data_2018_2023 %>%
  mutate(Hybrid_Field = paste(Hybrid, "_",Field.Location)) %>%
  mutate(Hybrid_Env = paste(Hybrid, "_",Env))  %>%
  dplyr::select("Hybrid_Env","Hybrid","Env","Year","Field.Location","Experiment","Replicate",
                "Block","Yield.Mg.ha","Pollen.DAP.days","Silk.DAP.days",
                "Plant.Height.cm","Ear.Height.cm","Yield.Mg.ha","Grain.Moisture","Twt.kg.m3")


#### Preparing for modeling ####

Pheno_data_2018_2023$Year <- as.factor(Pheno_data_2018_2023$Year)
Pheno_data_2018_2023$Hybrid<- as.factor(Pheno_data_2018_2023$Hybrid)
Pheno_data_2018_2023$Replicate <- as.factor(Pheno_data_2018_2023$Replicate)
Pheno_data_2018_2023$Block <- as.factor(Pheno_data_2018_2023$Block)
Pheno_data_2018_2023$Field.Location <- as.factor(Pheno_data_2018_2023$Field.Location)
Pheno_data_2018_2023$Hybrid_Env <- as.factor(Pheno_data_2018_2023$Hybrid_Env)



###Yield without outliers

outlier.resid <- function(data, trait){

  data<- data[complete.cases(data[,paste(trait)]),]

  #Fit a linear model for the residuals
  Model  <- lm(paste0(trait, "~  Hybrid_Env + Replicate"), data=data)

  ##Residual
  data$RESID <- resid(Model)

  #finding out what two standard deviations is and save it

  SD2 <- 2 * sd(resid(Model))

  # Flag the outliers as 1
  out_flag<-ifelse(data$RESID > SD2, 1, 0)

  # Combine the orginal data
  out_data<-cbind(data, out_flag) # Combine the orginal data

  file_name <- paste0(trait,"_outlier", ".csv")

  write.table(out_data,file = file_name,col.names = T,row.names = F,sep = '\t',quote = F)

  #Make a new data.frame with no outliers
  out_data<- out_data[!out_data$out_flag,]

  file_name1 <- paste0(trait,"_without_outlier", ".csv")

  write.table(out_data,file = file_name1,col.names = T,row.names = F,sep = '\t',quote = F)

  return(out_data)
}


trait <- commandArgs(trailingOnly = T)[1]

print(trait)

Pheno_data_2018_2023[[trait]] <- as.numeric(Pheno_data_2018_2023[[trait]])

Pheno_data_without_outlier<- outlier.resid(Pheno_data_2018_2023, trait = trait)

### Summary statistics table of the yield data across years and location####
Pheno_data_without_outlier_summary <- Pheno_data_without_outlier %>%
  dplyr::select(trait,"RESID")%>%
  summarise_all(list(min=min,median=median,max=max,mean=mean,sd=sd), na.rm=TRUE) %>%
  gather(stat,val)%>%
  separate(stat, into = c("phenotype","stat"),sep = "_") %>%
  spread(stat,val)

Pheno_data_without_outlier_summary$CV <- Pheno_data_without_outlier_summary[,6] / Pheno_data_without_outlier_summary[,3] *100
Pheno_data_without_outlier_summary$St.err <-  Pheno_data_without_outlier_summary[,6] / sqrt(length(Pheno_data_without_outlier$Yield.Mg.ha))


trait <- commandArgs(trailingOnly = T)[1]

print(paste("Trait:", trait))

file <- paste0("/home/uni08/osatohanmwen/2NP_final/", trait, "_without_outlier",".csv")

Pheno_data_without_outlier <- data.frame(fread(file))%>%
    dplyr::select("Hybrid_Env","Hybrid","Env","Year","Field.Location","Experiment","Replicate",
                  "Block",trait)

Pheno_data_without_outlier[[trait]] <- as.numeric(Pheno_data_without_outlier[[trait]])

#### Preparing for modeling ####

Pheno_data_without_outlier$Year <- as.factor(Pheno_data_without_outlier$Year)
Pheno_data_without_outlier$Hybrid<- as.factor(Pheno_data_without_outlier$Hybrid)
Pheno_data_without_outlier$Replicate <- as.factor(Pheno_data_without_outlier$Replicate)
Pheno_data_without_outlier$Block <- as.factor(Pheno_data_without_outlier$Block)
Pheno_data_without_outlier$Field.Location <- as.factor(Pheno_data_without_outlier$Field.Location)
Pheno_data_without_outlier$Hybrid_Env <- as.factor(Pheno_data_without_outlier$Hybrid_Env)


BLUES_All_Year=list()

Years = c("2018","2019","2020","2021","2022","2023")

j = 1

for (k in 1:length(Years)){
  
  print(paste0(' Year:', Years[k]))
  data = filter(Pheno_data_without_outlier, Year %in% Years[k]) 
  
  for (f in unique(data$Field.Location)) {
    
    if (length(unique(data[data$Field.Location == f &
                           !is.na(data$Field.Location), 'Replicate'])) > 1){
      print(paste0(f,' : at least 2 reps'))
      
      mod1 = lm(paste0(trait, " ~ 0 +  Hybrid + Replicate "), 
                data=data[data$Field.Location == f & 
                            !is.na(data$Field.Location),])
      
      mean_df = data.frame(Years[k], f,head(gsub("Hybrid", "", row.names(as.data.frame(coef(mod1)))), -1),head(as.data.frame(coef(mod1))[,1], -1),paste0(f, "_",Years[k]))
      
      colnames(mean_df) =c("Year", "Field.Location", "Hybrid",trait, "Env" )
      
      
      BLUES_All_Year[[j]] = cbind(mean_df)
      j = 1 + j
    }
    
    else{
      print(paste0(f,' : only 1 rep'))
      
      BLUES_All_Year[[j]] = data[data$Field.Location == f &!is.na(data[[trait]]),
                                 c("Year", "Field.Location", "Hybrid",trait,"Env")]
    }
    j = 1 + j
  }
}

file_name2 <- paste0(trait,"_BLUES_All", ".csv")

BLUES_All=ldply(BLUES_All_Year, rbind)

write.table(BLUES_All,file = file_name2,col.names = T,row.names = F,sep = '\t',quote = F)

# data frame  parameters
unique_hybrid_env <- as.character(unique(Pheno_data_without_outlier$Hybrid_Env))

##BLUES WITH Env, Hybrid As Random EFFECTS
BLUES_All <- BLUES_All %>%
  mutate(Hybrid_Env = paste(Hybrid, "_",Env))%>%
  dplyr::select(Hybrid_Env,Hybrid,Field.Location, Year,Env,trait)%>%
  filter(Hybrid_Env %in% unique_hybrid_env)

BLUES_All$Hybrid<- as.factor(BLUES_All$Hybrid)
BLUES_All$Env <- as.factor(BLUES_All$Env)

BLUES_All[[trait]] <- as.numeric(BLUES_All[[trait]])


BLUES_FUN1 <- function(trait, data= ".") {
  Model <- as.data.frame(fixef (lmer(paste0(trait, "~0 + Hybrid + (1|Env)") , data=data)))
}

BLUES_1 <- lapply(trait, BLUES_FUN1, data=BLUES_All)

BLUES_All_1 <- map2(BLUES_1, trait, ~set_names(..1,..2) %>%
                      rownames_to_column(var = "Hybrid"))%>%reduce(full_join)

file_name3 <- paste0(trait,"_BLUES_All_final", ".csv")

write.table(BLUES_All_1,file = file_name3,col.names = T,row.names = F,sep = '\t',quote = F)

##BLues for Hybrid for each year
BLUES_Year=list()

for (k in 1:length(Years)){
  
  print(paste0(' Year:', Years[k]))
  
  data = filter(BLUES_All, Year %in% Years[k]) 
  
  data$Hybrid<- as.factor(data$Hybrid)
  
  data$Field.Location <- as.factor(data$Field.Location)
  
  data[[trait]] <- as.numeric(data[[trait]])
  
  BLUES_FUN1 <- function(trait, data= ".") {
    Model <- as.data.frame(fixef (lmer(paste0(trait, "~0 + Hybrid + (1|Field.Location)")
                                       , data=data)))
  }
  
  BLUES_F <- lapply(trait, BLUES_FUN1, data=data)
  BLUES_All_F <- map2(BLUES_F, trait, ~set_names(..1,..2) %>% rownames_to_column(var = "Hybrid"))%>% reduce(full_join)
  
  mean_df = data.frame(BLUES_All_F, Years[k])
  
  colnames(mean_df) =c("Hybrid",trait, "Year" )
  
  BLUES_Year[[k]] = cbind(mean_df)
}

file_name4 <- paste0(trait,"_BLUES_Year", ".csv")

BLUES_Year=ldply(BLUES_Year, rbind)

write.table(BLUES_Year,file = file_name4,col.names = T,row.names = F,sep = '\t',quote = F)

#### Heritability and Variance ####

#calculating variance and heritability on line mean basis

####Create empty data frame

Reshape_da <- data.frame()
DataVar <- data.frame()
DataVarOut <- data.frame()
Heritability_cullis <- data.frame(Trait=NULL,H2=NULL)
Heritability_F_M    <- data.frame(Trait=NULL,H2=NULL)

##what to drop from the variance components fine
drops <- c("sdcor", "var1", "var2")

#### Heritability for all phenotype ####

summary(unique(Pheno_data_without_outlier$Env))
summary(unique(Pheno_data_without_outlier$Replicate))

Model<- lmer(paste0(trait, "~ (1|Hybrid) + Replicate + (1|Env) +  (1|Hybrid:Env)"), data=Pheno_data_without_outlier)

##Broad sense heritability assuming a balanced data (all genotype are equally observed in every environment and year)
##Falconer and Mackay, 1996

varComp <- as.data.frame(VarCorr(Model))
varComp <-varComp[,!(names(varComp) %in% drops)]
varComp$Traits <- trait
DataVar <- rbind (DataVar, varComp)

#Reshaping the variance component to easily run the heritability script
Reshape_da <- reshape(varComp, idvar ="Traits", timevar ="grp", direction = "wide" )

##Broad sense heritability assumming a balanced data (all genotype are equally observed in every environment and year)
##Falconer and Mackay, 1996
summary(Model)

H2 <- ((Reshape_da[,3]) / ((Reshape_da[,3]) + (Reshape_da[,5]/436)))

H2<- data.frame(Trait= paste(trait), H2 = H2)

Heritability_F_M <- rbind(Heritability_F_M, H2)

##Broad sense heritability assuming an unbalanced data using BLUPs (all genotype are not equally observed in every environment and year)
## Cullis et all,2006

library(arm)
ses<- se.ranef(Model)$'Hybrid'
v_BLUP<- ses^2
sigma2_g=VarCorr(Model, comp="Variance")$'Hybrid'[1]

Reliability <- 1- (v_BLUP/ (2 *sigma2_g))

H2<- data.frame(Trait= paste(trait), H2 = round(mean(Reliability),3))
Heritability_cullis <- rbind(Heritability_cullis,H2)

Heritability <- data.frame(Traits=Heritability_F_M$Trait, Heritability_cullis=Heritability_cullis$H2,Heritability_F_M=Heritability_F_M$H2 )

#Reshaping the variance component 
DataVarOut <- reshape(DataVar, idvar ="Traits", timevar ="grp", direction = "wide" )

## adding both heritability estimate to the data frame of variance  
DataVarOut <- merge(DataVarOut,Heritability)

###merge with the the other summary table
Data_her_sum <- merge(DataVarOut, Pheno_data_without_outlier_summary, by.x ="Traits", by.y = "phenotype")
## add co-efficient of variation
#Data_her_sum$CV <- Data_her_sum[,13] / Data_her_sum[,10]
#Data_her_sum$CV1 <- (sqrt(Data_her_sum[,4])) / (Data_her_sum[,10])


#write out a data frame of the data summary , variance and heritability estimates
file_name4 <- paste0(trait,"_h2_sum_var", ".csv")

write.table(Data_her_sum,file = file_name4,col.names = T,row.names = F,sep = '\t',quote = F)








