# installing basic libraries
import datatable as dt
import numpy as np
import pandas as pd

from sklearn.model_selection import KFold, StratifiedKFold,cross_val_score,LeaveOneGroupOut
import sklearn.metrics as metrics
from sklearn.inspection import permutation_importance
from skopt import BayesSearchCV

from scipy.stats import pearsonr
import lightgbm as lgb
import time
import random
import time
import sys

trait = sys.argv[1]

print(trait)

# Set random seed for reproducibility
random.seed(99)
np.random.seed(99)

#Phenotype
file_pheno = f"/home/uni08/osatohanmwen/2NP_final/{trait}_BLUES_Year.csv"

Trait = dt.fread(file_pheno)
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
features1 = Full_data1
features2 = Full_data2.drop(columns=['Hybrid','Year', trait])
features2 = features2.add_suffix("_dom") 

Fulldata = pd.concat([features1, features2], axis=1)

features = Fulldata.drop(columns=['Hybrid','Year',trait])
outcome = Fulldata[trait] 
years   = Fulldata['Year']

print(features.shape)

#feature importance

def feature_importances(model,colnames):
    d = dict(zip(colnames, model.feature_importances_))
    d = dict(sorted(d.items(), key=lambda x: -x[1]))
    df = pd.DataFrame(d.items(), columns=['feature', 'imp'])
    return df

def selection_efficiency(observed, predicted, top_percent=20):
   
    assert len(observed) == len(predicted), "Observed and predicted sets must have the same length."

    T = len(observed)  # Total number of hybrids
    top_k = int(np.ceil((top_percent / 100) * T))  # Top 20% hybrids
    # Get indices of top observed and top predicted hybrids
    top_observed_idx = np.argsort(observed)[-top_k:]  # Top 20% observed
    top_predicted_idx = np.argsort(predicted)[-top_k:]  # Top 20% predicted
    # Compute B (number of common hybrids)
    B = len(set(top_observed_idx) & set(top_predicted_idx))
    # Compute R (expected number selected by chance)
    R = int(np.ceil((top_k ** 2) / T))
    # Compute Coincidence Index
    CI =(B - R) / (top_k - R)
    
    return CI

##split into training and testing test using the leave one year CV

logo = LeaveOneGroupOut()
# split data into training and test set, then run Bayessearch to determine the best parameters and run the Random forest for prediction.

cnt = 1
Data_OP =[]
data_accuracy     = {}
data_pre_observed = list()

for train_index, test_index in logo.split(features,outcome,groups=years):
    
    t_start = time.time()

    xtrain = features.iloc[train_index]
    column_names = xtrain.columns.tolist()
    xtrain =xtrain.values

    xtest  = features.iloc[test_index].values
    ytrain = outcome[train_index]
    ytest  = outcome[test_index]


    LGBr_clf= lgb.LGBMRegressor(random_state = 99)

    # Bayesian optimization using an iterative Gaussian process 

    params_lg = {

        'n_estimators': (5000,15000, "log-uniform"), # No of trees# 220,520,620,
        'learning_rate' : (0.001, 0.3,"log-uniform"), 


        'num_leaves' : (10, 10000), 
        'min_child_samples' : (10, 200, "log-uniform"), 

        'max_depth':(2,50) , # maximum depth to explore (meaning the longest path between the root node and the leaf node.'max_depth': [10,15],
        'min_child_weight' : (0.01,20, "log-uniform"),
        'subsample': (0.1, 1), 
        'colsample_bytree' : (0.1, 1), 


        'reg_alpha' : (1, 200, "log-uniform"), 
        'reg_lambda' : (1, 10, "log-uniform"), 

        'n_jobs' : [-1] ,
        'seed' : [99] ,
        'verbose' : [-20]
        
            }

    Bayes_search = BayesSearchCV(
        LGBr_clf,params_lg,n_iter=25, 
        n_jobs = -1,cv=5,scoring='neg_mean_squared_error',
        random_state=99,
        verbose=0)

    print(Bayes_search.total_iterations)

    Bayes_search.fit(xtrain,ytrain)

    print(Bayes_search.best_estimator_)

    #update XGBM parameters usingb the best estimator which was obtained by using BayesSearchCV.

    LGBr_clf_best = Bayes_search.best_estimator_

    LGBr_clf_best.fit(xtrain, ytrain)

    y_pred_grid = LGBr_clf_best.predict(xtest)

        ##observed and predicted dataframe
    data_pre_observed = pd.DataFrame(
        {'Observed' : ytest,
        'Predicted': y_pred_grid
        }
    )

    ##Accuracy
    data_accuracy = pd.DataFrame(
        {'r2_score' : metrics.r2_score(ytest, y_pred_grid),
        'pearsonr'  : pearsonr(ytest, y_pred_grid)[0],
        'rmse'      : np.sqrt(metrics.mean_squared_error(ytest, y_pred_grid)),
        'Trait'     : f'{trait}',
        'Select_eff': selection_efficiency(ytest, y_pred_grid),
        'Year'     : f'{years[test_index].values[1]}'
        }, index=[0]
    )

    data_pre_observed.to_csv(f'{trait}_LGBM_OP_L0YO' + f'{years[test_index].values[1]}' + '.csv',index=False)

    data_accuracy.to_csv(f'{trait}_LGBM_ACC_LOYO' + f'{years[test_index].values[1]}' + '.csv' ,index=False)
    
    ###feature_importance
    FI_Lgbm= feature_importances(LGBr_clf_best,column_names)
    FI_Lgbm.to_csv(f'{trait}_FI_Lgbm_{years[test_index].values[1]}_LOYO.csv',index=False)


    FI_perm = permutation_importance(LGBr_clf_best,xtrain,ytrain, random_state=99, n_jobs=-1)
    forest_importances = pd.Series(FI_perm.importances_mean, index=column_names)
    forest_importances.to_csv(f'{trait}_FI_perm_{years[test_index].values[1]}_LOYO.csv',index=[1])

    print(cnt)
    print(time.time() - t_start)
    #print((ytest, y_pred_grid)[0])
    print(pearsonr(ytest, y_pred_grid)[0])
    cnt += 1
