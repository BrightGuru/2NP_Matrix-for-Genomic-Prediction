####Packages####
library(tidyverse)
library(data.table)
library(plyr)
library(dplyr)

#### 5F ####
# Define the directory where your CSV files are located
directory_path <- "~/Downloads/2NP_final/Result_5F"

setwd(directory_path)

#Getting result of all  traits for 5F

ML_Result_5F <- list.files(pattern = "_ACC_5F_.*csv") %>% 
  lapply(
    function(file){
      data <- read.csv(file, sep = ',', stringsAsFactors = FALSE)
      data$File <-  sub(".csv.*", "", file)  # Add a new column 'File' with the file name
      #data$Trait =sub("LGBM.*", "", data$File ) 
      data$Model =paste0("2NPLGBM")
      data <- data %>%
        separate(Fold, into = c("Part1","Repeat", "Repeat_Num", "Fold", "Fold_Num"), sep = "_") %>%
        mutate(CV = paste0("5F"),
               Repeat = paste(Repeat, Repeat_Num, sep = "_"),
               Fold = paste(Fold, Fold_Num, sep = "_")) %>%
        select(r2_score,pearsonr,rmse,Trait,Select_eff,Repeat, Fold,Model,CV) 
      return(data)
    }
  ) %>% 
  bind_rows()

#ML_Result_5F$Select_eff <- (ML_Result_5F$Select_eff) /100

directory_path <- "~/Downloads/2NP_final/GBLUP_Result_5F"

setwd(directory_path)

GBLUP_Result_5F <- list.files(pattern = "_Result_.*csv") %>% 
  lapply(
    function(file){
      data <- read.csv(file, sep = ',', stringsAsFactors = FALSE)
      data$File <-  sub(".csv.*", "", file)  # Add a new column 'File' with the file name
      data$Part1 <- gsub("_Result_.*", "", data$File)  # Extract everything before "_Result_"
      data$Part2 <- gsub(".*_Result_", "", data$File)
      data$Part2 <- gsub("Repeat", "_Repeat", data$Part2)
      data$Model <- gsub(".*_", "GBLUP_", data$Part1)
      data <- data %>%
        separate(Part2, into = c("Part1", "Repeat", "Repeat_Num", "Fold", "Fold_Num"), sep = "_") %>%
        mutate(CV = paste(Part1, sep="_"),
               Repeat = paste(Repeat, Repeat_Num, sep = "_"),
               Fold = paste(Fold, Fold_Num, sep = "_")) %>%
        select(r2_score,pearsonr,rmse,Trait,Select_eff,Repeat, Fold,Model,CV)
      return(data)
    }
  ) %>% 
  bind_rows()


F5 <- rbind(ML_Result_5F,GBLUP_Result_5F) %>%
  filter(Trait%in% c("Yield.Mg.ha", "Plant.Height.cm"))

X2NPlGBM_5F <- read_csv("~/Downloads/2NP_final/New_result_DTAS_corrected/2NPlGBM_5F.csv")
X2NPlGBM_5F <- X2NPlGBM_5F %>%
  separate(Fold, into = c("Repeat", "Repeat_Num", "Fold", "Fold_Num"), sep = "_") %>%  
  separate(CV, into = c("Model", "OP", "CV"), sep = "_") %>%
  mutate(Repeat = paste(Repeat, Repeat_Num, sep = "_"),
            Fold = paste(Fold, Fold_Num, sep = "_")) %>%
  select(r2_score,pearsonr,rmse,Trait,Select_eff,Repeat, Fold,Model,CV)

X2NPlGBM_5F$Model <- recode(X2NPlGBM_5F$Model,"LGBM" = "2NPLGBM")

GBLUP_Result_5F <- read_csv("~/Downloads/2NP_final/New_result_DTAS_corrected/GBLUP_Result_5F.csv")  

GBLUP_Result_5F <- GBLUP_Result_5F %>%
  separate(Fold, into = c("Repeat", "Repeat_Num", "Fold", "Fold_Num"), sep = "_") %>%  
  separate(CV, into = c("Model1", "Model","OP", "CV"), sep = "_") %>%
  mutate(Repeat = paste(Repeat, Repeat_Num, sep = "_"),
         Fold = paste(Fold, Fold_Num, sep = "_")) %>%
  select(r2_score,pearsonr,rmse,Trait,Select_eff,Repeat, Fold,Model,CV)
  
GBLUP_Result_5F$Model <- recode(GBLUP_Result_5F$Model,"OP" = "GBLUP_ADD","dom" = "GBLUP_ADDOM")
GBLUP_Result_5F$CV <- "5F"

F5 <- rbind(GBLUP_Result_5F, X2NPlGBM_5F,F5)%>%
  filter(Trait!= "Grain.Moisture")

F5$Model <- recode(F5$Model,
                   "GBLUP_Add" = "GBLUP_ADD",
                   "GBLUP_Addom" = "GBLUP_ADDOM")


F5$Model <- factor(F5$Model, 
                   levels = c( "GBLUP_ADD", "GBLUP_ADDOM","2NPLGBM"))

F5$Trait <- recode(F5$Trait,
                   "Plant.Height.cm" = "PH",
                   "Pollen.DAP.days" = "DTA",
                   "Silk.DAP.days"   = "DTS",
                   "Yield.Mg.ha"     = "GY"
)
F5$Trait <- factor(F5$Trait,levels = c("PH","DTA","DTS","GY"))

data_summary_5F <- F5 %>%
  group_by(Trait, Model) %>%
  dplyr::summarise(
    pearsonr_mean = mean(pearsonr, na.rm = TRUE),
    pearsonr_sd   = sd(pearsonr, na.rm = TRUE),
    pearsonr_se   = sd(pearsonr, na.rm = TRUE) / sqrt(sum(!is.na(pearsonr))),
    
    rmse_mean     = mean(rmse, na.rm = TRUE),
    rmse_sd       = sd(rmse, na.rm = TRUE),
    rmse_se       = sd(rmse, na.rm = TRUE) / sqrt(sum(!is.na(rmse))),
    
    Select_eff_mean = mean(Select_eff, na.rm = TRUE),
    Select_eff_sd   = sd(Select_eff, na.rm = TRUE),
    Select_eff_se   = sd(Select_eff, na.rm = TRUE) / sqrt(sum(!is.na(Select_eff)))
  )


##Plot the results for 5F
ggplot(data_summary_5F, aes(x = Trait, y = pearsonr_mean, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = pearsonr_mean - pearsonr_sd,
                    ymax = pearsonr_mean + pearsonr_sd),
                position = position_dodge(width = 0.9),
                width = 0.2,
                linewidth = 3) +
  coord_cartesian(ylim = c(0.00, 0.9)) +
  theme_minimal() +
  labs(title = "\nPredictive ability for maize traits\n",
       y = "\nPearson Correlation\n",
       x = "\nTraits\n",
       fill = "Genomic Model") +
  theme(axis.text.x = element_text(face="bold", 
                                   size=70, hjust = 0.5,angle = 0),
        axis.text.y = element_text(face="bold", 
                                   size=70),
        axis.title = element_text(face="bold",size=50),
        plot.title = element_text(face="bold",size=84),
        legend.text = element_text(size=70,face = "bold"),
        legend.title = element_text(size=70,face = "bold"),
        legend.key.size = unit(4, 'cm'))+
  scale_fill_manual(values = c("GBLUP_ADD" = "#8ecae6",
                               "GBLUP_ADDOM" = "#711DB0",
                               "2NPLGBM"="#5F0F40"),
                    name = "Genomic Model", 
                    labels=c("GBLUP_ADD", "GBLUP_ADDOM","2NPLGBM"))+
  scale_x_discrete(expand = c(0.0, 0))

file_name_result <- file.path("~/Downloads/2NP_final/New_result_DTAS_corrected/", paste0("F5_PA",".jpg"))
ggsave(file_name_result, plot=last_plot(), dpi = 700, width = 40, height = 20, units = "in")

ggplot(data_summary_5F, aes(x = Trait, y = Select_eff_mean, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = Select_eff_mean - Select_eff_sd,
                    ymax = Select_eff_mean + Select_eff_sd),
                position = position_dodge(width = 0.9),
                width = 0.2,
                linewidth = 3) +
  coord_cartesian(ylim = c(0.00, 0.9)) +
  theme_minimal() +
  labs(title = "\nSelection Efficiency for maize traits\n",
       y = "\nEff. of Selection\n",
       x = "\nTraits\n",
       fill = "Genomic Model") +
  theme(axis.text.x = element_text(face="bold", 
                                   size=70, hjust = 0.5,angle = 0),
        axis.text.y = element_text(face="bold", 
                                   size=70),
        axis.title = element_text(face="bold",size=50),
        plot.title = element_text(face="bold",size=84),
        legend.text = element_text(size=70,face = "bold"),
        legend.title = element_text(size=70,face = "bold"),
        legend.key.size = unit(4, 'cm'))+
  scale_fill_manual(values = c("GBLUP_ADD" = "#8ecae6",
                               "GBLUP_ADDOM" = "#711DB0",
                               "2NPLGBM"="#5F0F40"),
                    name = "Genomic Model", 
                    labels=c("GBLUP_ADD", "GBLUP_ADDOM","2NPLGBM"))+
  scale_x_discrete(expand = c(0.0, 0))

file_name_result <- file.path("~/Downloads/2NP_final/New_result_DTAS_corrected/", paste0("F5_SE",".jpg"))
ggsave(file_name_result, plot=last_plot(), dpi = 700, width = 40, height = 20, units = "in")

####LOYO####
# Define the directory where your CSV files are located
directory_path <- "~/Downloads/2NP_final/LOYO_RW_result"

setwd(directory_path)

#Getting result of all  traits for LOYO

Sim_Result_LOYO <- list.files(pattern = "_ACC_LOYO.*csv") %>% 
  lapply(
    function(file){
      data <- read.csv(file, sep = ',', stringsAsFactors = FALSE)
      data$File <-  sub(".csv.*", "", file)  # Add a new column 'File' with the file name
      data$Model =paste0("2NPLGBM")
      data$Repeat <- 1
      data$CV <- "LOYO"
      data$Fold <- paste0("Fold","_",data$Year)
      data <- data %>%
        select(r2_score,pearsonr,rmse,Trait,Select_eff,Repeat, Fold,Model,CV) 
      return(data)
    }
  ) %>% 
  bind_rows()

AdddomYear <-list.files(pattern = "_Adddom_Re.*csv") %>% 
  lapply(
    function(file){
      data <- read.csv(file, sep = '\t', stringsAsFactors = FALSE)
      data$File <-  sub(".csv.*", "", file)  # Add a new column 'File' with the file name
      data$Model =paste0("GBLUP_ADDOM")
      #data$Repeat <- 1
      data$CV <- "LOYO"
      data$Fold <- paste0("Fold","_",data$Year)
      data <- data %>%
        select(r2_score,pearsonr,rmse,Trait,Select_eff,Repeat, Fold,Model,CV) 
      return(data)
    }
  ) %>% 
  bind_rows()

AddYear <- list.files(pattern = "_Add_Re.*csv") %>% 
  lapply(
    function(file){
      data <- read.csv(file, sep = '\t', stringsAsFactors = FALSE)
      data$File <-  sub(".csv.*", "", file)  # Add a new column 'File' with the file name
      data$Model =paste0("GBLUP_ADD")
      #data$Repeat <- 1
      data$CV <- "LOYO"
      data$Fold <- paste0("Fold","_",data$Year)
      data <- data %>%
        select(r2_score,pearsonr,rmse,Trait,Select_eff,Repeat, Fold,Model,CV) 
      return(data)
    }
  ) %>% 
  bind_rows() 


LOYO_model1 <- rbind(Sim_Result_LOYO,AdddomYear,AddYear) %>%
  filter(Trait%in% c("Yield.Mg.ha", "Plant.Height.cm"))

GBLUP_Result_LOYO <- read_csv("~/Downloads/2NP_final/New_result_DTAS_corrected/GBLUP_Result_LOYO.csv")
GBLUP_Result_LOYO <- GBLUP_Result_LOYO %>%
  separate(CV, into = c("Model", "Fold"), sep = "_(?=[^_]+$)") %>%
  mutate(Fold = gsub("^OP", "Fold_", Fold),
         Repeat = 1,
         CV = paste("LOYO"),
         Model = recode(Model,
                       `Add` = "GBLUP_ADD",
                       `Add_dom` = "GBLUP_ADDOM"
                       # add more specific values if needed
         )) %>%
  select(r2_score,pearsonr,rmse,Trait,Select_eff,Repeat,Fold,Model,CV)

LGBM_Result_RWLOYO <- read_csv("~/Downloads/2NP_final/New_result_DTAS_corrected/LGBM_Result_RWLOYO.csv")

LGBM_Result_RWLOYO <- LGBM_Result_RWLOYO %>%
  separate(CV, into = c("Model", "OP", "Fold", "f"), sep = "_") %>%
  filter(Fold != "RW") %>%
  mutate(
    Fold = gsub("^L0YO", "Fold_", Fold),   # replace LOYO prefix with Fold_
    Repeat = 1,                             # add Repeat column
    CV = "LOYO",                            # set CV column to "LOYO"
    Model = recode(Model,
                   `LGBM` = "2NPLGBM"       # recode specific model names
                   # add more recodes if needed
    )
  ) %>%
  select(r2_score, pearsonr, rmse, Trait, Select_eff, Repeat, Fold, Model, CV)

LOYO_model <- rbind(LOYO_model1,LGBM_Result_RWLOYO,GBLUP_Result_LOYO)%>%
  filter(Trait!="Grain.Moisture")

LOYO_model$Model <- factor(LOYO_model$Model, 
                           levels = c( "GBLUP_ADD", "GBLUP_ADDOM","2NPLGBM"))

LOYO_model$Trait <- recode(LOYO_model$Trait,
                           "Plant.Height.cm" = "PH",
                           "Pollen.DAP.days" = "DTA",
                           "Silk.DAP.days"   = "DTS",
                           "Yield.Mg.ha"     = "GY"
)
LOYO_model$Trait <- factor(LOYO_model$Trait,levels = c("PH","DTA","DTS","GY"))


data_LOYO <- LOYO_model  %>%
  select(c("Trait", "Model","Fold","pearsonr"))%>%
  group_by(Trait, Model)%>%
  pivot_wider(
    names_from = Fold,
    values_from = pearsonr,
    names_prefix = "pearsonr_"
  )

data_LOYO <- LOYO_model  %>%
  select(c("Trait", "Model","Fold","Select_eff"))%>%
  group_by(Trait, Model)%>%
  pivot_wider(
    names_from = Fold,
    values_from = Select_eff,
    names_prefix = "Select_eff_"
  )


data_summary_LOYO <- LOYO_model %>%
  group_by(Trait, Model) %>%
  dplyr::summarise(
    pearsonr_mean = mean(pearsonr, na.rm = TRUE),
    pearsonr_sd   = sd(pearsonr, na.rm = TRUE),
    pearsonr_se   = sd(pearsonr, na.rm = TRUE) / sqrt(sum(!is.na(pearsonr))),
    
    rmse_mean     = mean(rmse, na.rm = TRUE),
    rmse_sd       = sd(rmse, na.rm = TRUE),
    rmse_se       = sd(rmse, na.rm = TRUE) / sqrt(sum(!is.na(rmse))),
    
    Select_eff_mean = mean(Select_eff, na.rm = TRUE),
    Select_eff_sd   = sd(Select_eff, na.rm = TRUE),
    Select_eff_se   = sd(Select_eff, na.rm = TRUE) / sqrt(sum(!is.na(Select_eff)))
  )


##Plot the results for LOYO
ggplot(data_summary_LOYO, aes(x = Trait, y = pearsonr_mean, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = pearsonr_mean - pearsonr_se,
                    ymax = pearsonr_mean + pearsonr_se),
                position = position_dodge(width = 0.9),
                width = 0.2,
                linewidth = 3) +
  coord_cartesian(ylim = c(0.00, 0.9)) +
  theme_minimal() +
  labs(title = "\nPredictive ability for maize traits\n",
       y = "\nPearson Correlation\n",
       x = "\nTraits\n",
       fill = "Genomic Model") +
  theme(axis.text.x = element_text(face="bold", 
                                   size=70, hjust = 0.5,angle = 0),
        axis.text.y = element_text(face="bold", 
                                   size=70),
        axis.title = element_text(face="bold",size=50),
        plot.title = element_text(face="bold",size=84),
        legend.text = element_text(size=70,face = "bold"),
        legend.title = element_text(size=70,face = "bold"),
        legend.key.size = unit(4, 'cm'))+
  scale_fill_manual(values = c("GBLUP_ADD" = "#8ecae6",
                               "GBLUP_ADDOM" = "#711DB0",
                               "2NPLGBM"="#5F0F40"),
                    name = "Genomic Model", 
                    labels=c("GBLUP_ADD", "GBLUP_ADDOM","2NPLGBM"))+
  scale_x_discrete(expand = c(0.0, 0))

file_name_result <- file.path("~/Downloads/2NP_final/New_result_DTAS_corrected/", paste0("LOYO_PA",".jpg"))
ggsave(file_name_result, plot=last_plot(), dpi = 700, width = 40, height = 20, units = "in")

ggplot(data_summary_LOYO, aes(x = Trait, y = Select_eff_mean, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = Select_eff_mean - Select_eff_sd,
                    ymax = Select_eff_mean + Select_eff_sd),
                position = position_dodge(width = 0.9),
                width = 0.2,
                linewidth = 3) +
  coord_cartesian(ylim = c(0.00, 0.8)) +
  theme_minimal() +
  labs(title = "\nSelection Efficiency for maize traits\n",
       y = "\nEff. of Selection\n",
       x = "\nTraits\n",
       fill = "Genomic Model") +
  theme(axis.text.x = element_text(face="bold", 
                                   size=70, hjust = 0.5,angle = 0),
        axis.text.y = element_text(face="bold", 
                                   size=70),
        axis.title = element_text(face="bold",size=50),
        plot.title = element_text(face="bold",size=84),
        legend.text = element_text(size=70,face = "bold"),
        legend.title = element_text(size=70,face = "bold"),
        legend.key.size = unit(4, 'cm'))+
  scale_fill_manual(values = c("GBLUP_ADD" = "#8ecae6",
                               "GBLUP_ADDOM" = "#711DB0",
                               "2NPLGBM"="#5F0F40"),
                    name = "Genomic Model", 
                    labels=c("GBLUP_ADD", "GBLUP_ADDOM","2NPLGBM"))+
  scale_x_discrete(expand = c(0.0, 0))

file_name_result <- file.path("~/Downloads/2NP_final/New_result_DTAS_corrected/", paste0("LOYO_SE",".jpg"))
ggsave(file_name_result, plot=last_plot(), dpi = 700, width = 40, height = 20, units = "in")


####RW####
####Getting result of all  traits for RW ####

Sim_Result_RW <- list.files(pattern = "_ACC_RW.*csv") %>% 
  lapply(
    function(file){
      data <- read.csv(file, sep = ',', stringsAsFactors = FALSE)
      data$File <-  sub(".csv.*", "", file)  # Add a new column 'File' with the file name
      data$Model =paste0("2NPLGBM")
      data$Repeat <- 1
      data$CV <- "RW"
      data$Fold <- paste0("Fold","_",data$Year)
      data <- data %>%
        select(r2_score,pearsonr,rmse,Trait,Select_eff,Repeat, Fold,Model,CV) 
      return(data)
    }
  ) %>% 
  bind_rows()

AdddomRW <-list.files(pattern = "_Addom_Rolling_Result.*csv") %>% 
  lapply(
    function(file){
      data <- read.csv(file, sep = '\t', stringsAsFactors = FALSE)
      data$File <-  sub(".csv.*", "", file)  # Add a new column 'File' with the file name
      data$Model =paste0("GBLUP_ADDOM")
      #data$Repeat <- 1
      data$CV <- "RW"
      data$Fold <- paste0("Fold","_",data$Year)
      data <- data %>%
        select(r2_score,pearsonr,rmse,Trait,Select_eff,Repeat, Fold,Model,CV) 
      return(data)
    }
  ) %>% 
  bind_rows()

AddRW <- list.files(pattern = "_Add_Rolling_Result.*csv") %>% 
  lapply(
    function(file){
      data <- read.csv(file, sep = '\t', stringsAsFactors = FALSE)
      data$File <-  sub(".csv.*", "", file)  # Add a new column 'File' with the file name
      data$Model =paste0("GBLUP_ADD")
      #data$Repeat <- 1
      data$CV <- "RW"
      data$Fold <- paste0("Fold","_",data$Year)
      data <- data %>%
        select(r2_score,pearsonr,rmse,Trait,Select_eff,Repeat, Fold,Model,CV) 
      return(data)
    }
  ) %>% 
  bind_rows() 


RW_model2 <- rbind(Sim_Result_RW,AdddomRW,AddRW)%>%
  filter(Trait%in% c("Yield.Mg.ha", "Plant.Height.cm"))%>%
  filter(Model%in% c("2NPLGBM"))


GBLUP_Result_RW_add <- read_csv("~/Downloads/2NP_final/New_result_DTAS_corrected/GBLUP_Result_RW_add.csv")
GBLUP_Result_RW_add3 <- read_csv("~/Downloads/2NP_final/New_result_DTAS_corrected/GBLUP_Result_RW_add3.csv")
GBLUP_Result_RW_addom <- read_csv("~/Downloads/2NP_final/New_result_DTAS_corrected/GBLUP_Result_RW_addom.csv")
GBLUP_Result_RW_addom3 <- read_csv("~/Downloads/2NP_final/New_result_DTAS_corrected/GBLUP_Result_RW_addom3.csv")

RW_model1 <- rbind(GBLUP_Result_RW_add, GBLUP_Result_RW_add3, GBLUP_Result_RW_addom, GBLUP_Result_RW_addom3)%>%
  separate(CV, into = c("Model", "OP", "Fold", "f","F"), sep = "_") %>%
  mutate(
    Fold = paste0("Fold", "_", Fold),   # replace LOYO prefix with Fold_
    Repeat = 1,                             # add Repeat column
    CV = "RW",                            # set CV column to "LOYO"
    Model = recode(Model,
                   `Add` = "GBLUP_ADD",
                   `Addom` = "GBLUP_ADDOM"# recode specific model names
                   # add more recodes if needed
    )
  ) %>%
  select(r2_score, pearsonr, rmse, Trait, Select_eff, Repeat, Fold, Model, CV)

LGBM_Result_RW <- read_csv("~/Downloads/2NP_final/New_result_DTAS_corrected/LGBM_Result_RWLOYO.csv")

LGBM_Result_RW <- LGBM_Result_RW %>%
  separate(CV, into = c("Model", "OP", "Fold", "f"), sep = "_") %>%
  filter(Fold == "RW") %>%
  mutate(
    Fold = paste0("Fold","_" ,as.character(f)),   # replace LOYO prefix with Fold_
    Repeat = 1,                             # add Repeat column
    CV = "RW",                            # set CV column to "LOYO"
    Model = recode(Model,
                   `LGBM` = "2NPLGBM"       # recode specific model names
                   # add more recodes if needed
    )
  ) %>%
  select(r2_score, pearsonr, rmse, Trait, Select_eff, Repeat, Fold, Model, CV)


RW_model <- rbind(RW_model2,RW_model1,LGBM_Result_RW)%>%
  filter(Trait!="Grain.Moisture")

RW_model$Model <- factor(RW_model$Model, 
                         levels = c( "GBLUP_ADD", "GBLUP_ADDOM","2NPLGBM"))

RW_model$Trait <- recode(RW_model$Trait,
                         "Plant.Height.cm" = "PH",
                         "Pollen.DAP.days" = "DTA",
                         "Silk.DAP.days"   = "DTS",
                         "Yield.Mg.ha"     = "GY"
                         
)
RW_model$Trait <- factor(RW_model$Trait,levels = c("PH","DTA","DTS","GY"))

data_RW <- RW_model  %>%
  select(c("Trait", "Model","Fold","Select_eff"))%>%
  group_by(Trait, Model)%>%
  pivot_wider(
    names_from = Fold,
    values_from = Select_eff,
    names_prefix = "Select_eff_"
  )

data_RW <- RW_model  %>%
  select(c("Trait", "Model","Fold","Select_eff"))%>%
  group_by(Trait, Model)%>%
  pivot_wider(
    names_from = Fold,
    values_from = Select_eff,
    names_prefix = "Select_eff_"
  )
data_summary_RW <- RW_model %>%
  group_by(Model,Trait) %>%
  dplyr::summarise(
    pearsonr_mean = mean(pearsonr, na.rm = TRUE),
    pearsonr_sd   = sd(pearsonr, na.rm = TRUE),
    pearsonr_se   = sd(pearsonr, na.rm = TRUE) / sqrt(sum(!is.na(pearsonr))),
    
    rmse_mean     = mean(rmse, na.rm = TRUE),
    rmse_sd       = sd(rmse, na.rm = TRUE),
    rmse_se       = sd(rmse, na.rm = TRUE) / sqrt(sum(!is.na(rmse))),
    
    Select_eff_mean = mean(Select_eff, na.rm = TRUE),
    Select_eff_sd   = sd(Select_eff, na.rm = TRUE),
    Select_eff_se   = sd(Select_eff, na.rm = TRUE) / sqrt(sum(!is.na(Select_eff)))
  )

##Plot the results for RW
ggplot(data_summary_RW, aes(x = Trait, y = pearsonr_mean, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = pearsonr_mean - pearsonr_se,
                    ymax = pearsonr_mean + pearsonr_se),
                position = position_dodge(width = 0.9),
                width = 0.2,
                linewidth = 3) +
  coord_cartesian(ylim = c(0.00, 0.9)) +
  theme_minimal() +
  labs(title = "\nPredictive ability for maize traits\n",
       y = "\nPearson Correlation\n",
       x = "\nTraits\n",
       fill = "Genomic Model") +
  theme(axis.text.x = element_text(face="bold", 
                                   size=70, hjust = 0.5,angle = 0),
        axis.text.y = element_text(face="bold", 
                                   size=70),
        axis.title = element_text(face="bold",size=50),
        plot.title = element_text(face="bold",size=84),
        legend.text = element_text(size=70,face = "bold"),
        legend.title = element_text(size=70,face = "bold"),
        legend.key.size = unit(4, 'cm'))+
  scale_fill_manual(values = c("GBLUP_ADD" = "#8ecae6",
                               "GBLUP_ADDOM" = "#711DB0",
                               "2NPLGBM"="#5F0F40"),
                    name = "Genomic Model", 
                    labels=c("GBLUP_ADD", "GBLUP_ADDOM","2NPLGBM"))+
  scale_x_discrete(expand = c(0.0, 0))

file_name_result <- file.path("~/Downloads/2NP_final/New_result_DTAS_corrected/", paste0("RW_PA",".jpg"))
ggsave(file_name_result, plot=last_plot(), dpi = 700, width = 40, height = 20, units = "in")

ggplot(data_summary_RW, aes(x = Trait, y = Select_eff_mean, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = Select_eff_mean - Select_eff_se,
                    ymax = Select_eff_mean + Select_eff_se),
                position = position_dodge(width = 0.9),
                width = 0.2,
                linewidth = 3) +
  coord_cartesian(ylim = c(0.00, 0.8)) +
  theme_minimal() +
  labs(title = "\nSelection Efficiency for maize traits\n",
       y = "\nEff. of Selection\n",
       x = "\nTraits\n",
       fill = "Genomic Model") +
  theme(axis.text.x = element_text(face="bold", 
                                   size=70, hjust = 0.5,angle = 0),
        axis.text.y = element_text(face="bold", 
                                   size=70),
        axis.title = element_text(face="bold",size=50),
        plot.title = element_text(face="bold",size=84),
        legend.text = element_text(size=70,face = "bold"),
        legend.title = element_text(size=70,face = "bold"),
        legend.key.size = unit(4, 'cm'))+
  scale_fill_manual(values = c("GBLUP_ADD" = "#8ecae6",
                               "GBLUP_ADDOM" = "#711DB0",
                               "2NPLGBM"="#5F0F40"),
                    name = "Genomic Model", 
                    labels=c("GBLUP_ADD", "GBLUP_ADDOM","2NPLGBM"))+
  scale_x_discrete(expand = c(0.0, 0))

file_name_result <- file.path("~/Downloads/2NP_final/New_result_DTAS_corrected/", paste0("RW_SE",".jpg"))
ggsave(file_name_result, plot=last_plot(), dpi = 700, width = 40, height = 20, units = "in")

##CVOO_Parent 1

# Define the directory where your CSV files are located
directory_path <- "~/Downloads/2NP_final/Result_Parent1_CVOO"

setwd(directory_path)

#Getting result of all  traits for Parent1_CVOO

ML_Result_Parent1_CVOO <- list.files(pattern = "_ACC_Parent1_CVOO_.*csv") %>% 
  lapply(
    function(file){
      data <- read.csv(file, sep = ',', stringsAsFactors = FALSE)
      data$File <-  sub(".csv.*", "", file)  # Add a new column 'File' with the file name
      #data$Trait =sub("LGBM.*", "", data$File ) 
      data$Model =paste0("2NPLGBM")
      data <- data %>%
        separate(Fold, into = c("Part1","Part2","Repeat", "Repeat_Num", "Fold", "Fold_Num"), sep = "_") %>%
        mutate(CV = paste0("Parent1_CVOO"),
               Repeat = paste(Repeat, Repeat_Num, sep = "_"),
               Fold = paste(Fold, Fold_Num, sep = "_")) %>%
        select(r2_score,pearsonr,rmse,Trait,Select_eff,Repeat, Fold,Model,CV) 
      return(data)
    }
  ) %>% 
  bind_rows()

ML_Result_Parent1_CVOO$Select_eff <- (ML_Result_Parent1_CVOO$Select_eff) /100

directory_path <- "~/Downloads/2NP_final/GBLUP_Result_Parent1_CVOO"

setwd(directory_path)

GBLUP_Result_Parent1_CVOO <- list.files(pattern = "_Result_.*csv") %>% 
  lapply(
    function(file){
      data <- read.csv(file, sep = ',', stringsAsFactors = FALSE)
      data$File <-  sub(".csv.*", "", file)  # Add a new column 'File' with the file name
      data$Part1 <- gsub("_Result_.*", "", data$File)  # Extract everything before "_Result_"
      data$Part2 <- gsub(".*_Result_", "", data$File)
      data$Part2 <- gsub("Repeat", "_Repeat", data$Part2)
      data$Model <- gsub(".*_", "GBLUP_", data$Part1)
      data <- data %>%
        separate(Part2, into = c("Part1","Part", "Repeat", "Repeat_Num", "Fold", "Fold_Num"), sep = "_") %>%
        mutate(CV = paste(Part1, sep="_"),
               Repeat = paste(Repeat, Repeat_Num, sep = "_"),
               Fold = paste(Fold, Fold_Num, sep = "_")) %>%
        select(r2_score,pearsonr,rmse,Trait,Select_eff,Repeat, Fold,Model,CV)
      return(data)
    }
  ) %>% 
  bind_rows()


Parent1_CVOO <- rbind(ML_Result_Parent1_CVOO,GBLUP_Result_Parent1_CVOO) %>%
  filter(Trait%in% c("Yield.Mg.ha", "Plant.Height.cm"))

X2NPlGBM_Parent1_CVOO <- read_csv("~/Downloads/2NP_final/New_result_DTAS_corrected/LGBM_Result_Parent1_CVOO.csv")
X2NPlGBM_Parent1_CVOO <- X2NPlGBM_Parent1_CVOO %>%
  separate(Fold, into = c("Repeat", "Repeat_Num", "Fold", "Fold_Num"), sep = "_") %>%  
  separate(CV, into = c("Model", "OP", "CV","CV1"), sep = "_") %>%
  mutate(Repeat = paste(Repeat, Repeat_Num, sep = "_"),
         Fold = paste(Fold, Fold_Num, sep = "_"),
         CV = paste(CV ,CV1, sep = "_")) %>%
  select(r2_score,pearsonr,rmse,Trait,Select_eff,Repeat, Fold,Model,CV)

X2NPlGBM_Parent1_CVOO$Model <- recode(X2NPlGBM_Parent1_CVOO$Model,"LGBM" = "2NPLGBM")

GBLUP_Result_Parent1_CVOO <- read_csv("~/Downloads/2NP_final/New_result_DTAS_corrected/GBLUP_Result_Parent1_CVOO.csv")  

GBLUP_Result_Parent1_CVOO <- GBLUP_Result_Parent1_CVOO %>%
  separate(Fold, into = c("Repeat", "Repeat_Num", "Fold", "Fold_Num"), sep = "_") %>%  
  separate(CV, into = c("Model1", "Model","OP", "CV"), sep = "_") %>%
  mutate(Repeat = paste(Repeat, Repeat_Num, sep = "_"),
         Fold = paste(Fold, Fold_Num, sep = "_")) %>%
  select(r2_score,pearsonr,rmse,Trait,Select_eff,Repeat, Fold,Model,CV)

GBLUP_Result_Parent1_CVOO$Model <- recode(GBLUP_Result_Parent1_CVOO$Model,"OP" = "GBLUP_ADD","dom" = "GBLUP_ADDOM")
GBLUP_Result_Parent1_CVOO$CV <- "Parent1_CVOO"

Parent1_CVOO <- rbind(GBLUP_Result_Parent1_CVOO, X2NPlGBM_Parent1_CVOO,Parent1_CVOO)%>%
  filter(Trait!= "Grain.Moisture")

Parent1_CVOO$Model <- recode(Parent1_CVOO$Model,
                   "GBLUP_Add" = "GBLUP_ADD",
                   "GBLUP_Addom" = "GBLUP_ADDOM")

Parent1_CVOO$Model <- factor(Parent1_CVOO$Model, 
                   levels = c( "GBLUP_ADD", "GBLUP_ADDOM","2NPLGBM"))

Parent1_CVOO$Trait <- recode(Parent1_CVOO$Trait,
                   "Plant.Height.cm" = "PH",
                   "Pollen.DAP.days" = "DTA",
                   "Silk.DAP.days"   = "DTS",
                   "Yield.Mg.ha"     = "GY"
)
Parent1_CVOO$Trait <- factor(Parent1_CVOO$Trait,levels = c("PH","DTA","DTS","GY"))

data_summary_Parent1_CVOO <- Parent1_CVOO %>%
  group_by(Trait, Model) %>%
  dplyr::summarise(
    pearsonr_mean = mean(pearsonr, na.rm = TRUE),
    pearsonr_sd   = sd(pearsonr, na.rm = TRUE),
    pearsonr_se   = sd(pearsonr, na.rm = TRUE) / sqrt(sum(!is.na(pearsonr))),
    
    rmse_mean     = mean(rmse, na.rm = TRUE),
    rmse_sd       = sd(rmse, na.rm = TRUE),
    rmse_se       = sd(rmse, na.rm = TRUE) / sqrt(sum(!is.na(rmse))),
    
    Select_eff_mean = mean(Select_eff, na.rm = TRUE),
    Select_eff_sd   = sd(Select_eff, na.rm = TRUE),
    Select_eff_se   = sd(Select_eff, na.rm = TRUE) / sqrt(sum(!is.na(Select_eff)))
  )

##Plot the results for Parent1_CVOO
ggplot(data_summary_Parent1_CVOO, aes(x = Trait, y = pearsonr_mean, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = pearsonr_mean - pearsonr_sd,
                    ymax = pearsonr_mean + pearsonr_sd),
                position = position_dodge(width = 0.9),
                width = 0.2,
                linewidth = 3) +
  coord_cartesian(ylim = c(0.00, 0.9)) +
  theme_minimal() +
  labs(title = "\nPredictive ability for maize traits\n",
       y = "\nPearson Correlation\n",
       x = "\nTraits\n",
       fill = "Genomic Model") +
  theme(axis.text.x = element_text(face="bold", 
                                   size=70, hjust = 0.5,angle = 0),
        axis.text.y = element_text(face="bold", 
                                   size=70),
        axis.title = element_text(face="bold",size=50),
        plot.title = element_text(face="bold",size=84),
        legend.text = element_text(size=70,face = "bold"),
        legend.title = element_text(size=70,face = "bold"),
        legend.key.size = unit(4, 'cm'))+
  scale_fill_manual(values = c("GBLUP_ADD" = "#8ecae6",
                               "GBLUP_ADDOM" = "#711DB0",
                               "2NPLGBM"="#5F0F40"),
                    name = "Genomic Model", 
                    labels=c("GBLUP_ADD", "GBLUP_ADDOM","2NPLGBM"))+
  scale_x_discrete(expand = c(0.0, 0))

file_name_result <- file.path("~/Downloads/2NP_final/New_result_DTAS_corrected/", paste0("Parent1_CVOO_PA",".jpg"))
ggsave(file_name_result, plot=last_plot(), dpi = 700, width = 40, height = 20, units = "in")

ggplot(data_summary_Parent1_CVOO, aes(x = Trait, y = Select_eff_mean, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = Select_eff_mean - Select_eff_sd,
                    ymax = Select_eff_mean + Select_eff_sd),
                position = position_dodge(width = 0.9),
                width = 0.2,
                linewidth = 3) +
  coord_cartesian(ylim = c(0.00, 0.9)) +
  theme_minimal() +
  labs(title = "\nSelection Efficiency for maize traits\n",
       y = "\nEff. of Selection\n",
       x = "\nTraits\n",
       fill = "Genomic Model") +
  theme(axis.text.x = element_text(face="bold", 
                                   size=70, hjust = 0.5,angle = 0),
        axis.text.y = element_text(face="bold", 
                                   size=70),
        axis.title = element_text(face="bold",size=50),
        plot.title = element_text(face="bold",size=84),
        legend.text = element_text(size=70,face = "bold"),
        legend.title = element_text(size=70,face = "bold"),
        legend.key.size = unit(4, 'cm'))+
  scale_fill_manual(values = c("GBLUP_ADD" = "#8ecae6",
                               "GBLUP_ADDOM" = "#711DB0",
                               "2NPLGBM"="#5F0F40"),
                    name = "Genomic Model", 
                    labels=c("GBLUP_ADD", "GBLUP_ADDOM","2NPLGBM"))+
  scale_x_discrete(expand = c(0.0, 0))

file_name_result <- file.path("~/Downloads/2NP_final/New_result_DTAS_corrected/", paste0("Parent1_CVOO_SE",".jpg"))
ggsave(file_name_result, plot=last_plot(), dpi = 700, width = 40, height = 20, units = "in")


##CVO_Parent 1

# Define the directory where your CSV files are located
directory_path <- "~/Downloads/2NP_final/Result_Parent1_CVO"

setwd(directory_path)

#Getting result of all  traits for Parent1_CVO

ML_Result_Parent1_CVO <- list.files(pattern = "_ACC_Parent1_CVO_.*csv") %>% 
  lapply(
    function(file){
      data <- read.csv(file, sep = ',', stringsAsFactors = FALSE)
      data$File <-  sub(".csv.*", "", file)  # Add a new column 'File' with the file name
      #data$Trait =sub("LGBM.*", "", data$File ) 
      data$Model =paste0("2NPLGBM")
      data <- data %>%
        separate(Fold, into = c("Part1","Part2","Repeat", "Repeat_Num", "Fold", "Fold_Num"), sep = "_") %>%
        mutate(CV = paste0("Parent1_CVOO"),
               Repeat = paste(Repeat, Repeat_Num, sep = "_"),
               Fold = paste(Fold, Fold_Num, sep = "_")) %>%
        select(r2_score,pearsonr,rmse,Trait,Select_eff,Repeat, Fold,Model,CV) 
      return(data)
    }
  ) %>% 
  bind_rows()

ML_Result_Parent1_CVO$Select_eff <- (ML_Result_Parent1_CVO$Select_eff) /100

directory_path <- "~/Downloads/2NP_final/GBLUP_Result_Parent1_CVO"

setwd(directory_path)

GBLUP_Result_Parent1_CVO <- list.files(pattern = "_Result_.*csv") %>% 
  lapply(
    function(file){
      data <- read.csv(file, sep = ',', stringsAsFactors = FALSE)
      data$File <-  sub(".csv.*", "", file)  # Add a new column 'File' with the file name
      data$Part1 <- gsub("_Result_.*", "", data$File)  # Extract everything before "_Result_"
      data$Part2 <- gsub(".*_Result_", "", data$File)
      data$Part2 <- gsub("Repeat", "_Repeat", data$Part2)
      data$Model <- gsub(".*_", "GBLUP_", data$Part1)
      data <- data %>%
        separate(Part2, into = c("Part1","Part", "Repeat", "Repeat_Num", "Fold", "Fold_Num"), sep = "_") %>%
        mutate(CV = paste(Part1, sep="_"),
               Repeat = paste(Repeat, Repeat_Num, sep = "_"),
               Fold = paste(Fold, Fold_Num, sep = "_")) %>%
        select(r2_score,pearsonr,rmse,Trait,Select_eff,Repeat, Fold,Model,CV)
      return(data)
    }
  ) %>% 
  bind_rows()


Parent1_CVO <- rbind(ML_Result_Parent1_CVO,GBLUP_Result_Parent1_CVO) %>%
  filter(Trait%in% c("Yield.Mg.ha", "Plant.Height.cm"))

X2NPlGBM_Parent1_CVO <- read_csv("~/Downloads/2NP_final/New_result_DTAS_corrected/LGBM_Result_Parent1_CVO.csv")
X2NPlGBM_Parent1_CVO <- X2NPlGBM_Parent1_CVO %>%
  separate(Fold, into = c("Repeat", "Repeat_Num", "Fold", "Fold_Num"), sep = "_") %>%  
  separate(CV, into = c("Model", "OP", "CV","CV1"), sep = "_") %>%
  mutate(Repeat = paste(Repeat, Repeat_Num, sep = "_"),
         Fold = paste(Fold, Fold_Num, sep = "_"),
         CV = paste(CV ,CV1, sep = "_")) %>%
  select(r2_score,pearsonr,rmse,Trait,Select_eff,Repeat, Fold,Model,CV)

X2NPlGBM_Parent1_CVO$Model <- recode(X2NPlGBM_Parent1_CVO$Model,"LGBM" = "2NPLGBM")

GBLUP_Result_Parent1_CVO <- read_csv("~/Downloads/2NP_final/New_result_DTAS_corrected/GBLUP_Result_Parent1_CV0.csv")  

GBLUP_Result_Parent1_CVO <- GBLUP_Result_Parent1_CVO %>%
  separate(Fold, into = c("Repeat", "Repeat_Num", "Fold", "Fold_Num"), sep = "_") %>%  
  separate(CV, into = c("Model1", "Model","OP", "CV"), sep = "_") %>%
  mutate(Repeat = paste(Repeat, Repeat_Num, sep = "_"),
         Fold = paste(Fold, Fold_Num, sep = "_")) %>%
  select(r2_score,pearsonr,rmse,Trait,Select_eff,Repeat, Fold,Model,CV)

GBLUP_Result_Parent1_CVO$Model <- recode(GBLUP_Result_Parent1_CVO$Model,"OP" = "GBLUP_ADD","dom" = "GBLUP_ADDOM")
GBLUP_Result_Parent1_CVO$CV <- "Parent1_CVO"

Parent1_CVO <- rbind(GBLUP_Result_Parent1_CVO, X2NPlGBM_Parent1_CVO,Parent1_CVO)%>%
  filter(Trait!= "Grain.Moisture")

Parent1_CVO$Model <- recode(Parent1_CVO$Model,
                             "GBLUP_Add" = "GBLUP_ADD",
                             "GBLUP_Addom" = "GBLUP_ADDOM")

Parent1_CVO$Model <- factor(Parent1_CVO$Model, 
                             levels = c( "GBLUP_ADD", "GBLUP_ADDOM","2NPLGBM"))

Parent1_CVO$Trait <- recode(Parent1_CVO$Trait,
                             "Plant.Height.cm" = "PH",
                             "Pollen.DAP.days" = "DTA",
                             "Silk.DAP.days"   = "DTS",
                             "Yield.Mg.ha"     = "GY"
)
Parent1_CVO$Trait <- factor(Parent1_CVO$Trait,levels = c("PH","DTA","DTS","GY"))

data_summary_Parent1_CVO <- Parent1_CVO %>%
  group_by(Trait, Model) %>%
  dplyr::summarise(
    pearsonr_mean = mean(pearsonr, na.rm = TRUE),
    pearsonr_sd   = sd(pearsonr, na.rm = TRUE),
    pearsonr_se   = sd(pearsonr, na.rm = TRUE) / sqrt(sum(!is.na(pearsonr))),
    
    rmse_mean     = mean(rmse, na.rm = TRUE),
    rmse_sd       = sd(rmse, na.rm = TRUE),
    rmse_se       = sd(rmse, na.rm = TRUE) / sqrt(sum(!is.na(rmse))),
    
    Select_eff_mean = mean(Select_eff, na.rm = TRUE),
    Select_eff_sd   = sd(Select_eff, na.rm = TRUE),
    Select_eff_se   = sd(Select_eff, na.rm = TRUE) / sqrt(sum(!is.na(Select_eff)))
  )

##Plot the results for Parent1_CVO
ggplot(data_summary_Parent1_CVO, aes(x = Trait, y = pearsonr_mean, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = pearsonr_mean - pearsonr_sd,
                    ymax = pearsonr_mean + pearsonr_sd),
                position = position_dodge(width = 0.9),
                width = 0.2,
                linewidth = 3) +
  coord_cartesian(ylim = c(0.00, 0.9)) +
  theme_minimal() +
  labs(title = "\nPredictive ability for maize traits\n",
       y = "\nPearson Correlation\n",
       x = "\nTraits\n",
       fill = "Genomic Model") +
  theme(axis.text.x = element_text(face="bold", 
                                   size=70, hjust = 0.5,angle = 0),
        axis.text.y = element_text(face="bold", 
                                   size=70),
        axis.title = element_text(face="bold",size=50),
        plot.title = element_text(face="bold",size=84),
        legend.text = element_text(size=70,face = "bold"),
        legend.title = element_text(size=70,face = "bold"),
        legend.key.size = unit(4, 'cm'))+
  scale_fill_manual(values = c("GBLUP_ADD" = "#8ecae6",
                               "GBLUP_ADDOM" = "#711DB0",
                               "2NPLGBM"="#5F0F40"),
                    name = "Genomic Model", 
                    labels=c("GBLUP_ADD", "GBLUP_ADDOM","2NPLGBM"))+
  scale_x_discrete(expand = c(0.0, 0))

file_name_result <- file.path("~/Downloads/2NP_final/New_result_DTAS_corrected/", paste0("Parent1_CVO_PA",".jpg"))
ggsave(file_name_result, plot=last_plot(), dpi = 700, width = 40, height = 20, units = "in")

ggplot(data_summary_Parent1_CVO, aes(x = Trait, y = Select_eff_mean, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = Select_eff_mean - Select_eff_sd,
                    ymax = Select_eff_mean + Select_eff_sd),
                position = position_dodge(width = 0.9),
                width = 0.2,
                linewidth = 3) +
  coord_cartesian(ylim = c(0.00, 0.9)) +
  theme_minimal() +
  labs(title = "\nSelection Efficiency for maize traits\n",
       y = "\nEff. of Selection\n",
       x = "\nTraits\n",
       fill = "Genomic Model") +
  theme(axis.text.x = element_text(face="bold", 
                                   size=70, hjust = 0.5,angle = 0),
        axis.text.y = element_text(face="bold", 
                                   size=70),
        axis.title = element_text(face="bold",size=50),
        plot.title = element_text(face="bold",size=84),
        legend.text = element_text(size=70,face = "bold"),
        legend.title = element_text(size=70,face = "bold"),
        legend.key.size = unit(4, 'cm'))+
  scale_fill_manual(values = c("GBLUP_ADD" = "#8ecae6",
                               "GBLUP_ADDOM" = "#711DB0",
                               "2NPLGBM"="#5F0F40"),
                    name = "Genomic Model", 
                    labels=c("GBLUP_ADD", "GBLUP_ADDOM","2NPLGBM"))+
  scale_x_discrete(expand = c(0.0, 0))

file_name_result <- file.path("~/Downloads/2NP_final/New_result_DTAS_corrected/", paste0("Parent1_CVO_SE",".jpg"))
ggsave(file_name_result, plot=last_plot(), dpi = 700, width = 40, height = 20, units = "in")
