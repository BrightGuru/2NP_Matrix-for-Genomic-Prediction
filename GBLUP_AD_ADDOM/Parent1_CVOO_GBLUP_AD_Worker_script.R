# Load required libraries
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

args <- commandArgs(trailingOnly = TRUE)

trait <- args[1]
i <- args[2]

# Set seed for reproducibility
set.seed(99)

print(paste("Trait:", trait))
print(paste("Fold:", i))

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
Pheno_data <- Pheno_Blues_2019_2023  %>%
  dplyr::select(c("Hybrid",trait))

Geno_data <- Pheno_data %>% select(c("Hybrid"))
Geno_data <- left_join(Geno_data,SNP_data)

SNP_data_G <- Geno_data %>% select(-Hybrid)%>%
  mutate_if(is.character, as.numeric)

## Van Raden(2008) Genomic relationship matrix (additive)##
Genomic_add_matrix <- Gmatrix(SNPmatrix=as.matrix(SNP_data_G),
                              maf=0.05, method="VanRaden")

Genomic_add_matrix <- data.frame(Hybrid=Geno_data[,1],Genomic_add_matrix)


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

###Genomic prediction with Addictive effects only.

# Load train/test splits
file_train <- paste0("/home/uni08/osatohanmwen/2NP_final/Dataset_Parent1_CVOO/xtrain_", i, ".csv")
file_test <- paste0("/home/uni08/osatohanmwen/2NP_final/Dataset_Parent1_CVOO/xtest_", i, ".csv")

train_df <- fread(file_train)
test_df <- fread(file_test)

# Filter data based on train/test Hybrid names

y <- as.matrix(Pheno_data[[trait]])

n=length(y) #Computing the length of the vector containing the phenotypic data

Z=diag(n)


Train <- which(Genomic_add_matrix$Hybrid %in%  train_df$Hybrid)
Test  <- which(Genomic_add_matrix$Hybrid %in%  test_df$Hybrid)

Genomic_add_matrix <- as.matrix(Genomic_add_matrix %>% select(-Hybrid))

funout<-emmreml(y=y[Train], X=matrix(rep(1, n)[Train], ncol=1), 
                Z=Z[Train,], K=Genomic_add_matrix)


va   <- funout$Vu
vd   <- 0
beta <- funout$betahat
ve   <- funout$Ve

yobs <- y[Test]
x=matrix(1,length(yobs),1)

yhat <- as.matrix(funout$uhat[Test]) + (x%*%beta) 

data_pre_observed <- data.frame(
  Observed = yobs,
  Predicted = yhat
)


data_accuracy <- data.frame(Fold = i,
                            Beta =beta,
                            Va = va,
                            Vd = vd,
                            Ve = ve,
                            pearsonr=cor(yobs,yhat),
                            r2_score=RSQUARE(yobs,yhat),
                            rmse=RMSE(yobs,yhat),
                            Select_eff=selection_efficiency(yobs,yhat),
                            Trait=trait)

output_dir <- "/home/uni08/osatohanmwen/2NP_final/GBLUP_Result_Parent1_CVOO"

# Create the directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message("Directory created: ", output_dir)
} else {
  message("Directory already exists: ", output_dir)
}

# Construct file name for observed data
file_name_observed <- file.path(output_dir, paste0(trait, 'Ad_OP_Parent1_CVOO', i, ".csv"))
write.csv(data_pre_observed, file = file_name_observed, row.names = FALSE)

file_name_result <- file.path(output_dir, paste0(trait, "_Add_Result_Parent1_CVOO", i, ".csv"))
write.csv(data_accuracy, file = file_name_result, row.names = FALSE)