library(tidyverse)
library(data.table)
library(foreach)
library(doParallel)
library(doFuture)
library(doMC)
library(plyr)
library(dplyr)
library(AGHmatrix)
library(EMMREML)

trait <- commandArgs(trailingOnly = T)[1]

print(paste("Trait:", trait))

file <- paste0("/home/uni08/osatohanmwen/2NP_final/",trait ,"_BLUES_Year",".csv")

Pheno_Blues_2019_2023 <- fread(file)

Pheno_Blues_2019_2023 <- Pheno_Blues_2019_2023 %>%
  mutate(Hybrid=gsub("Hybrid", "",Hybrid))

print(summary(unique(Pheno_Blues_2019_2023$Hybrid)))

##load in filtered and recoded genotype data and sort with the hybrid from phenotype data

SNP_data <- data.frame(fread(file = "/home/uni08/osatohanmwen/G2F_maize/2NP_matrix_4_GP/LGBM_GBLUP_2NP/Genotype_True.csv"))%>%
  filter(Hybrid %in% Pheno_Blues_2019_2023$Hybrid)

Pheno_Blues_2019_2023 <- Pheno_Blues_2019_2023 %>%
  filter(Hybrid %in% SNP_data$Hybrid)

###Join year and Genotype
Pheno_data <- Pheno_Blues_2019_2023 %>%
  mutate(Hybrid_Year= paste(Hybrid, "_",Year)) %>%
  dplyr::select(c("Hybrid_Year", "Hybrid","Year",trait))

Geno_data <- Pheno_data %>% select(c("Hybrid_Year", "Hybrid"))
Geno_data <- left_join(Geno_data,SNP_data)

SNP_data_G <- Geno_data %>% select(-Hybrid_Year, -Hybrid)%>%
  mutate_if(is.character, as.numeric)

## Van Raden(2008) Genomic relationship matrix (additive)##
Genomic_add_matrix <- Gmatrix(SNPmatrix=as.matrix(SNP_data_G),
                              maf=0.05, method="VanRaden")

#Computing the dominance relationship matrix based on Vitezica 2013
Genomic_dom_matrix <- Gmatrix(SNPmatrix=as.matrix(SNP_data_G), missingValue=-9, 
                              maf=0.05, method="Vitezica")

#5 folds cross validation 20 times ####

#1. MEAN ABSOLUTE PERCENTAGE ERROR (MAPE)
MAPE = function(y_actual,y_predict){
  mean(abs((y_actual-y_predict)/y_actual))*100
}

#2. R SQUARED error metric -- Coefficient of Determination
RSQUARE = function(y_actual,y_predict){
  cor(y_actual,y_predict)^2
}

#3 ROOT MEAN SQUARE ERROR
RMSE = function(y_actual,y_predict){
  sqrt(mean((y_actual - y_predict)^2))
}

#4 Selection Efficiency 

selection_efficiency <- function(observed, predicted, top_percent = 20) {
  
  stopifnot(length(observed) == length(predicted))  # Ensure same length
  
  Total <- length(observed)  # Total number of hybrids
  top_k <- ceiling((top_percent / 100) * Total)  # Top 20% hybrids
  
  # Get indices of top observed and predicted hybrids
  top_observed_idx <- order(observed, decreasing = TRUE)[1:top_k]
  top_predicted_idx <- order(predicted, decreasing = TRUE)[1:top_k]
  
  # Compute B (number of common hybrids in top selections)
  B <- length(intersect(top_observed_idx, top_predicted_idx))
  
  # Compute R (expected number selected by chance)
  R <- ceiling((top_k^2) / Total)
  
  # Compute Coincidence Index (CI)
  CI <- (B - R) / (top_k - R)
  
  return(CI)
}

#number of codes
n.cpus <- Sys.getenv("SLURM_CPUS_PER_TASK")
# we need this to be numeric below so:
n.cpus <- as.numeric(n.cpus)
n.cpus
class(n.cpus)

# register a parallel backend specifying the number of CPUs as the number we imported via Sys.getenv()

registerDoMC(cores = n.cpus) 

###Genomic prediction with Addictive effects only.

All_TEN <- list()

for (time in 1:1){
  seed  <- 1234
  
  set.seed(seed + time)
  
  y <- as.matrix(Pheno_data[[trait]])
  
  n=length(y) #Computing the length of the vector containing the phenotypic data
  
  Z=diag(n)
  
  #x=matrix(1,n,1) ##Computing the vector of fixed effects
  Year =c(unique(Pheno_data$Year))
  folds <- Pheno_data$Year
  nfolds <- length(Year)
  
  cross_predict_add <- foreach (
    i=1:nfolds, .multicombine=T,.combine="rbind"
  ) %dopar% {
    
    Train <- which(folds != Year[i])
    Test  <- which(folds == Year[i])
    
    funout<-emmremlMultiKernel(y=y[Train], 
                               X=matrix(rep(1, n)[Train], ncol=1),
                               Z=list(Z[Train,],Z[Train,]), 
                               K=list(Genomic_add_matrix,
                                      Genomic_dom_matrix))
    
    funout$weights[1]
    funout$Vu
    funout$Ve
    funout$betahat
    
    va   <- funout$Vu * funout$weights[1]
    vd   <- funout$Vu * funout$weights[2]
    beta <- funout$betahat
    ve   <- funout$Ve
    uhatmat <-matrix(funout$uhat, ncol=2) %>%
      rowSums()
    
    yobs <- y[Test]
    x=matrix(1,length(yobs),1)
    yhat <- as.matrix(uhatmat[Test]) + (x%*%beta)
    
    data_pre_observed <- data.frame(
      Observed = yobs,
      Predicted = yhat
    )
    
    file_name <- paste0(trait, 'Add_dom_OP', Year[i], '.csv')
    write.csv(data_pre_observed, file = file_name, row.names = FALSE)
    
    list(Fold = i,
         Repeat = time,
         Year=Year[i],
         Beta =beta,
         Va = va,
         Vd = vd,
         Ve = ve,
         pearsonr=cor(yobs,yhat),
         r2_score=RSQUARE(yobs,yhat),
         rmse=RMSE(yobs,yhat),
         Select_eff=selection_efficiency(yobs,yhat),
         Trait=trait
    )
  }
  
  RESULT_DF  <- data.frame(cross_predict_add)
  RESULT_DF <- apply(RESULT_DF,2,as.character)
  
  All_TEN[[time]] <- cbind(RESULT_DF)
  
}

RESULT_add_dom <- ldply(All_TEN,rbind)

write.table(RESULT_add_dom ,file = paste0(trait,"_Adddom_Result_Year", ".csv"),col.names = T,row.names = F,sep = '\t',quote = F)

mean(as.numeric(RESULT_add_dom$COR))
