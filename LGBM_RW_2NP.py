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

features = Fulldata.drop(columns=['Hybrid',trait])
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

# Main loop for Rolling Window CV
window_size = 3  # You can adjust the window size (e.g., 3 years of training data)
cnt = 1
Data_OP = []
data_accuracy = {}
data_pre_observed = []

# Assuming `features`, `outcome`, and `years` are defined as:
# features: DataFrame of feature columns
# outcome: Series/array of target variable (e.g., trait values)
# years: array or Series containing the years corresponding to each sample

n_years = len(sorted(years.unique()))
print(n_years)

# Rolling window iteration
for i in range(n_years - window_size):

    t_start = time.time()

    # Define the training and test sets
    train_years = years[i:i + window_size]
    test_year = years[i + window_size]

    print(f"Training years: {train_years}, Test year: {test_year}")

    # Select the training and test data based on the rolling window
    train_data = features[features['Year'].isin(train_years)]
    test_data = features[features['Year'] == test_year]

    xtrain = train_data.drop('Year', axis=1).values
    ytrain = outcome[features['Year'].isin(train_years)]
    xtest = test_data.drop('Year', axis=1).values
    ytest = outcome[features['Year'] == test_year]

    column_names = train_data.drop('Year', axis=1).columns.tolist()

    # Initialize the LGBMRegressor model
    LGBr_clf = lgb.LGBMRegressor(random_state=99)

    # Define parameter search space for Bayesian optimization
    params_lg = {
        'n_estimators': (5000, 15000, "log-uniform"),
        'learning_rate': (0.001, 0.3, "log-uniform"),
        'num_leaves': (10, 10000),
        'min_child_samples': (10, 200, "log-uniform"),
        'max_depth': (2, 50),
        'min_child_weight': (0.01, 20, "log-uniform"),
        'subsample': (0.1, 1),
        'colsample_bytree': (0.1, 1),
        'reg_alpha': (1, 200, "log-uniform"),
        'reg_lambda': (1, 10, "log-uniform"),
        'n_jobs': [-1],
        'seed': [99],
        'verbose': [-20]
    }

    # Perform Bayesian Optimization with BayesSearchCV
    Bayes_search = BayesSearchCV(LGBr_clf, params_lg, n_iter=25,
                                 n_jobs=-1, cv=5, scoring='neg_mean_squared_error',
                                 random_state=99, verbose=0)

    print(Bayes_search.total_iterations)
    # Fit Bayesian search to training data
    Bayes_search.fit(xtrain, ytrain)

    # Get best model from Bayesian search
    LGBr_clf_best = Bayes_search.best_estimator_
    print(Bayes_search.best_estimator_)
    # Fit the best model
    LGBr_clf_best.fit(xtrain, ytrain)

    # Make predictions
    y_pred_grid = LGBr_clf_best.predict(xtest)

    # Store observed vs predicted values
    data_pre_observed = pd.DataFrame({
        'Observed': ytest,
        'Predicted': y_pred_grid
    })

    # Accuracy metrics
    data_accuracy = pd.DataFrame({
        'r2_score': metrics.r2_score(ytest, y_pred_grid),
        'pearsonr': pearsonr(ytest, y_pred_grid)[0],
        'rmse': np.sqrt(metrics.mean_squared_error(ytest, y_pred_grid)),
        'Trait': f'{trait}',
        'Select_eff': selection_efficiency(ytest, y_pred_grid),
        'Year': f'{test_year}'
    }, index=[0])

    # Save results to CSV
    data_pre_observed.to_csv(f'{trait}_LGBM_OP_RW_{test_year}.csv', index=False)
    data_accuracy.to_csv(f'{trait}_LGBM_ACC_RW_{test_year}.csv', index=False)

    # Feature importance (LGBM native)
    FI_Lgbm = feature_importances(LGBr_clf_best, column_names)
    FI_Lgbm.to_csv(f'{trait}_FI_Lgbm_RW_{test_year}.csv', index=False)

    # Permutation importance (additional method)
    FI_perm = permutation_importance(LGBr_clf_best, xtrain, ytrain, random_state=99, n_jobs=-1)
    forest_importances = pd.Series(FI_perm.importances_mean, index=column_names)
    forest_importances.to_csv(f'{trait}_FI_perm_RW_{test_year}.csv', index=False)

    # Print iteration information
    print(f'Fold {cnt}')
    print(f'Time elapsed: {time.time() - t_start:.2f} seconds')
    print(f'Pearson correlation: {pearsonr(ytest, y_pred_grid)[0]:.2f}')
    
    cnt += 1
