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
import random
import time
import os

# Set random seed for reproducibility
random.seed(99)
np.random.seed(99)

trait = "Yield.Mg.ha"

#Phenotype
Trait = dt.fread("/home/uni08/osatohanmwen/2NP_final/Yield.Mg.ha_BLUES_Year.csv")
Trait=Trait.to_pandas()
# Remove "Hybrid" from the Hybrid column
Trait['Hybrid']=Trait['Hybrid'].str.replace("Hybrid", "",regex=False)
Trait[['Parent_1','Parent_2']] = Trait['Hybrid'].str.split('/', expand=True)

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

###create the features, year and the outcome variables

features1 = Full_data1
features2 = Full_data2.drop(columns=['Hybrid','Year', 'Yield.Mg.ha','Parent_1','Parent_2'])
features2 = features2.add_suffix("_dom") 

features = pd.concat([features1, features2], axis=1)
years   = Full_data1['Year']
outcome =Full_data1['Yield.Mg.ha'] 

# train known hybrids, predict in unknown year

Test_set = features[features['Year'].isin([2023])].reset_index(drop=True)
Train_set= features[features['Year'].isin([2018, 2019, 2020, 2021, 2022])].reset_index(drop=True)

print(Test_set.shape)
print(Train_set.shape)

# Define GroupKFold

gkf = GroupKFold(n_splits=5)

# split data into training and test set

# Define output directory
output_dir = output_dir = "/home/uni08/osatohanmwen/2NP_final/Dataset_Parent2_CVO"
os.makedirs(output_dir, exist_ok=True)  # Ensure the folder exists


Data_OP =[]
data_accuracy     = {}
data_pre_observed = list()

n_repeats = 10

for repeat in range (n_repeats):

    print(f"\n🔄 Running Repeat {repeat + 1}/{n_repeats}...")
    random.seed(repeat)  # Set a different random seed for each repeat  
    print(f"\n🔄 Running Repeat {repeat + 1}/{n_repeats}...")

    random.seed(repeat)  # Set a different random seed for each repeat

    shuffled_indices = random.sample(range(len(Test_set)), len(Test_set))
    shuffled_indices1 = random.sample(range(len(Train_set)), len(Train_set))

    cnt = 1

    # GroupKFold on the Test Set
    for test_index, _ in gkf.split(Test_set, groups=Test_set['Hybrid']):

        t_start = time.time()

        Test_set_shuffle = Test_set.iloc[shuffled_indices]
        # Create Test Fold
        xtest = Test_set_shuffle.iloc[test_index]
        # Filter to includeentries from Train_set where either Parent_1 or Parent_2 is present in xtest
        xtrain = Train_set[
        Train_set['Parent_2'].isin(xtest['Parent_2'])]


        # Test-to-Training Hybrid Ratio (Before Adjustment)
        Test_train_ratio_before = len(set(xtest['Hybrid'])) / len(set(xtrain['Hybrid']))

        print(f"Test to train ratio (before adjustment): {Test_train_ratio_before:.3f}") 

      
        # No overlap in years between training and test sets
        assert set(xtrain['Year']) & set(xtest['Year']) == set(), "Training and test sets have overlapping years."
        #assert set(additional_hybrids) & set(xtest['Hybrid']) == set(), "Overlap detected between additional hybrids and test hybrids!"

        print(f"xtest shape: {xtest.shape}, xtrain shape: {xtrain.shape}")

        
        # Store Hybrid IDs for Train/Test Sets
        hybrid_train = xtrain[['Hybrid', 'Year','Parent_1','Parent_2']]
        hybrid_test  = xtest[['Hybrid', 'Year','Parent_1','Parent_2']]   
        
        hybrid_train.to_csv(os.path.join(output_dir, f"xtrain_Repeat_{repeat+1}_fold_{cnt}.csv"), index=False)
        hybrid_test.to_csv(os.path.join(output_dir, f"xtest_Repeat_{repeat+1}_fold_{cnt}.csv"), index=False)

        cnt += 1  # Move to next fold
    


