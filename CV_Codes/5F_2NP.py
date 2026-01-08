# installing basic libraries
import datatable as dt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.decomposition import TruncatedSVD
from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold, StratifiedKFold, cross_val_score
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.model_selection import GroupKFold
import sklearn.metrics as metrics

from scipy.stats import pearsonr
import lightgbm as lgb
import xgboost as xgb
import time
import pickle
import json
import random
import time
import os

# Set random seed for reproducibility
random.seed(99)
np.random.seed(99)

trait = "Yield.Mg.ha"

#Phenotype
Trait = dt.fread("/home/uni08/osatohanmwen/2NP_final/Yield.Mg.ha_BLUES_All_final.csv")
Trait=Trait.to_pandas()
# Remove "Hybrid" from the Hybrid column
Trait['Hybrid']=Trait['Hybrid'].str.replace("Hybrid", "",regex=False)
print(Trait.head())
print(Trait.shape)

#Genotype
SNP_addev = dt.fread('/home/uni08/osatohanmwen/G2F_maize/2NP_matrix_4_GP/LGBM_GBLUP_2NP/Genomic_addev_matrix.csv')
SNP_dodev = dt.fread('/home/uni08/osatohanmwen/G2F_maize/2NP_matrix_4_GP/LGBM_GBLUP_2NP/Genomic_dodev_matrix.csv')

SNP_addev=SNP_addev.to_pandas()
SNP_dodev=SNP_dodev.to_pandas()

# Filter genotype data based on Hybrids present in phenotype data
SNP_addev= SNP_addev[SNP_addev['Hybrid'].isin(Trait['Hybrid'])]
SNP_dodev= SNP_dodev[SNP_dodev['Hybrid'].isin(Trait['Hybrid'])]

# Filter phenotype data based on Hybrids present in genotype data
Trait = Trait[Trait['Hybrid'].isin(SNP_addev['Hybrid'])]

Full_data1 = SNP_addev.merge(Trait,how='inner', on="Hybrid")
print(Full_data1.shape)
Full_data2 = SNP_dodev.merge(Trait,how='inner', on="Hybrid")
print(Full_data2.shape)
###create the features and the outcome variables

features1 = Full_data1.drop(columns=['Yield.Mg.ha'])
features2 = Full_data2.drop(columns=['Hybrid', 'Yield.Mg.ha'])
features2 = features2.add_suffix("_dom") 

features = pd.concat([features1, features2], axis=1)

#features = np.concatenate([features1,features2], axis=1)

outcome =Full_data1['Yield.Mg.ha'] 

##split into training and testing test and reapet 10 times

# Define output directory
output_dir = "/home/uni08/osatohanmwen/2NP_final/Dataset_5F"
os.makedirs(output_dir, exist_ok=True)  # Ensure the folder exists

# Set the number of times you want to run k-fold cross-validation
n_repeats = 10


for repeat in range (n_repeats):

    random.seed(repeat)  # Set a different random seed for each repeat  

    kf =KFold(n_splits=5, shuffle=True, random_state=repeat)

    # split data into training and test set, then run Bayessearch to determine the best parameters and run the Random forest for prediction.

    cnt = 1

    for train_index, test_index in kf.split(features, outcome):

        t = time.time()

        xtrain = features.iloc[train_index]
        xtest  = features.iloc[test_index]
        ytrain = outcome[train_index]
        ytest  = outcome[test_index]

        # Store Hybrid IDs for Train/Test Sets
        hybrid_train = xtrain['Hybrid'] 
        hybrid_test  = xtest['Hybrid']

        hybrid_train.to_csv(os.path.join(output_dir, f"xtrain_Repeat_{repeat+1}_fold_{cnt}.csv"), index=False)
        hybrid_test.to_csv(os.path.join(output_dir, f"xtest_Repeat_{repeat+1}_fold_{cnt}.csv"), index=False)

        

        cnt += 1

