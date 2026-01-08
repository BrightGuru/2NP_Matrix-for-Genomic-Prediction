# Load packages
import datatable as dt
import pandas as pd
import numpy as np
import lightgbm as lgb
from sklearn.model_selection import train_test_split, LeaveOneGroupOut
import sklearn.metrics as metrics
from scipy.stats import pearsonr
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix

from skopt import BayesSearchCV
import matplotlib.pyplot as plt
import seaborn as sns
import shap
import random
import time
import sys
import os


# Set random seed for reproducibility
random.seed(99)
np.random.seed(99)


trait = sys.argv[1]

#feature importance

def feature_importances(model,colnames):
    d = dict(zip(colnames, model.feature_importances_))
    d = dict(sorted(d.items(), key=lambda x: -x[1]))
    df = pd.DataFrame(d.items(), columns=['feature', 'imp'])
    return df


#Define parameter search space for Bayesian optimization
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

target = trait 

print(f"\nRunning analysis for target: {target}")

#Phenotype
Trait = dt.fread(f"/home/uni08/osatohanmwen/2NP_final/{target}_BLUES_All_final.csv")
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

features1 = Full_data1.drop(columns=['Hybrid', target])

features2 = Full_data2.drop(columns=['Hybrid', target])
features2 = features2.add_suffix("_dom") 

# Concatenate while keeping column names
features = pd.concat([features1, features2], axis=1)

# Convert to NumPy if needed
outcome =Full_data1[target] 

X_train, X_test, y_train, y_test = train_test_split(
    features, outcome, test_size=0.2, random_state=99)

####Analysis####

# Choose model type based on target variable
if y_test.nunique() > 10 and y_test.dtype in [float, int]:
    model_type = 'regression'
    LGB_clf = lgb.LGBMRegressor(random_state=99)
    scoring = 'neg_mean_squared_error'
else:
    model_type = 'classification'
    LGB_clf = lgb.LGBMClassifier(random_state=99)
    scoring = 'accuracy'

# Perform Bayesian Optimization with BayesSearchCV
Bayes_search = BayesSearchCV(
    LGB_clf,
    params_lg,
    n_iter=20,  # Number of iterations for optimization
    cv=5,  # Cross-validation
    scoring=scoring,  # Scoring metric for classification
    n_jobs=-1,
    random_state=99,
    verbose=1
)

# Fit Bayesian search to training data
Bayes_search.fit(X_train, y_train)

# Get the best model
LGB_clf_best = Bayes_search.best_estimator_

# Print best model and score
print("Best score:", Bayes_search.best_score_)
print("Best parameters:", Bayes_search.best_params_)

# Fit the best model
LGB_clf_best.fit(X_train, y_train)

# Make predictions
y_pred = LGB_clf_best.predict(X_test)

if model_type == "classification":
        
        accuracy = accuracy_score(y_test, y_pred)
        class_report = classification_report(y_test, y_pred, output_dict=True, zero_division=0)

        class_metrics = {
            "Seed_Lot": "Across",
            "Accuracy": accuracy,
            "Precision_0": class_report["0"]["precision"],
            "Recall_0": class_report["0"]["recall"],
            "F1_0": class_report["0"]["f1-score"],
            "Precision_1": class_report["1"]["precision"],
            "Recall_1": class_report["1"]["recall"],
            "F1_1": class_report["1"]["f1-score"],
            "Precision_2": class_report["2"]["precision"],
            "Recall_2": class_report["2"]["recall"],
            "F1_2": class_report["2"]["f1-score"]
        }
        acc_across_df = pd.DataFrame([class_metrics])
        acc_across_df.to_csv(f"{target}_LGBM_Evaluation_Summary.csv" ,index=False)

else:  # regression
        reg_metrics = {
            "Target_var": target,
            "Seed_Lot"  : "Across",
            'r2_score'  : metrics.r2_score(y_test, y_pred),
            "Accuracy"  : pearsonr(y_test, y_pred)[0],
            'rmse'      : np.sqrt(metrics.mean_squared_error(y_test, y_pred)),
        }

        acc_across_df = pd.DataFrame([reg_metrics])
        acc_across_df.to_csv(f"{target}_LGBM_Evaluation_Summary.csv" ,index=False)

# Create explainer
explainer = shap.TreeExplainer(LGB_clf_best)
shap_values = explainer.shap_values(X_test)

# Convert SHAP values to DataFrame
if isinstance(shap_values, list):  # classification returns a list per class
    shap_df = pd.DataFrame(np.abs(shap_values[0]), columns=X_test.columns)
else:  # regression returns a single array
    shap_df = pd.DataFrame(np.abs(shap_values), columns=X_test.columns)

# Compute mean absolute SHAP value per feature
mean_abs_shap = shap_df.mean().reset_index()
mean_abs_shap.columns = ['feature', 'mean_abs_shap']

print(mean_abs_shap.head)

 #Assign colors: blue for additive, purple for dominance
mean_abs_shap['color'] = mean_abs_shap['feature'].apply(
    lambda f: 'purple' if f.endswith('_dom') else 'blue'
)

# Create hue column for effect type
mean_abs_shap['Effect'] = mean_abs_shap['feature'].apply(
    lambda f: 'Dominance' if f.endswith('_dom') else 'Additive'
)

# Top 20
top20 = mean_abs_shap.sort_values(by='mean_abs_shap', ascending=False).head(20)

# Plot with proper hue mapping and palette
plt.figure(figsize=(10, 6))
sns.barplot(
    data=top20,
    y='feature',
    x='mean_abs_shap',
    hue='Effect',
    palette={'Additive': 'blue', 'Dominance': 'purple'},
    dodge=False,
    legend=False  # You can set True if you want to show legend
)
plt.xlabel("Mean(|SHAP value|)")
plt.title(f"SHAP Feature Importance\nTrait: {target}", fontsize=14)
plt.tight_layout()
plt.savefig(f"SHAP_{target}_Top20_Colored.png", bbox_inches='tight')
plt.close()


# Separate additive and dominance
shap_add = mean_abs_shap[~mean_abs_shap['feature'].str.endswith('_dom')]
shap_dom = mean_abs_shap[mean_abs_shap['feature'].str.endswith('_dom')]

# Sum contributions
shap_add_total = shap_add['mean_abs_shap'].sum()
shap_dom_total = shap_dom['mean_abs_shap'].sum()
shap_total = shap_add_total + shap_dom_total

# Relative (%) contribution
shap_add_pct = 100 * shap_add_total / shap_total
shap_dom_pct = 100 * shap_dom_total / shap_total

# Print
print(f"[SHAP] Additive contribution: {shap_add_total:.4f} ({shap_add_pct:.2f}%)")
print(f"[SHAP] Dominance contribution: {shap_dom_total:.4f} ({shap_dom_pct:.2f}%)")

# Plot additive vs dominance SHAP contributions
shap_contrib_df = pd.DataFrame({
    'Effect': ['Additive', 'Dominance'],
    'Mean Absolute SHAP Value': [shap_add_total, shap_dom_total],
    'Relative (%)': [shap_add_pct, shap_dom_pct]
})

shap_contrib_df.to_csv(f"SHAP_Additive_vs_Dominance_{target}.csv",index=False)

# Define custom colors
custom_palette = {'Additive': 'royalblue', 'Dominance': 'purple'}

sns.barplot(data=shap_contrib_df, x='Effect', y='Relative (%)', palette=custom_palette)
plt.title(f"SHAP-Based Contribution\n(Trait: {target})")
plt.ylabel("Relative Contribution (%)")
plt.tight_layout()
plt.savefig(f"SHAP_Additive_vs_Dominance_{target}.png")
plt.close()

# Select top 20 features based on mean absolute SHAP value
top20_features = mean_abs_shap.sort_values(by='mean_abs_shap', ascending=False).head(20)['feature'].tolist()

# Slice X and SHAP values for top 20
X_top20 = X_test[top20_features]
if isinstance(shap_values, list):  # classification
    shap_top20 = shap_values[0][:, [X_test.columns.get_loc(f) for f in top20_features]]
else:  # regression
    shap_top20 = shap_values[:, [X_test.columns.get_loc(f) for f in top20_features]]

# Beeswarm plot
shap.summary_plot(shap_top20, X_top20, show=False)
plt.title(f"SHAP Beeswarm Plot (Top 20)\nTrait: {target} | Across", fontsize=14)
plt.tight_layout()
plt.savefig(f"SHAP_Beeswarm_Top20_{target}.png", bbox_inches='tight')
plt.close()

# ---- Create explainer and compute interaction values safely ----
interaction_values = explainer.shap_interaction_values(X_test)

# ---- Use class 0 (for classification); skip if regression ----
if isinstance(interaction_values, list):
    interaction_values = interaction_values[0]

# ---- Create masks ----
features = X_test.columns
is_dom = features.str.endswith('_dom')
is_add = ~is_dom

# ---- Get interaction contributions ----
print(interaction_values.shape)
main_effects = interaction_values.diagonal(axis1=1, axis2=2)

# 1. Additive main effects (diagonal)

add_main = np.abs(main_effects[:, is_add]).sum()

# 2. Dominance main effects (diagonal)
dom_main = np.abs(main_effects[:, is_dom]).sum()

# 3. Additive × Additive (off-diagonal)
aa_mask = np.outer(is_add, is_add)
np.fill_diagonal(aa_mask, False)
aa_inter = np.abs(interaction_values[:, aa_mask]).sum()

# 4. Additive × Dominance (symmetric)
ad_mask = np.outer(is_add, is_dom) | np.outer(is_dom, is_add)
ad_inter = np.abs(interaction_values[:, ad_mask]).sum()

# 5. Dominance × Dominance (off-diagonal)
dd_mask = np.outer(is_dom, is_dom)
np.fill_diagonal(dd_mask, False)
dd_inter = np.abs(interaction_values[:, dd_mask]).sum()

# ---- Summarize ----
summary_df = pd.DataFrame({
    'Interaction Type': [
        'Additive Main',
        'Dominance Main',
        'Additive × Additive',
        'Additive × Dominance',
        'Dominance × Dominance'
    ],
    'Mean SHAP Interaction Value': [
        add_main,
        dom_main,
        aa_inter,
        ad_inter,
        dd_inter
    ]
})

# ---- Compute Relative Contribution ----
total = summary_df['Mean SHAP Interaction Value'].sum()
summary_df['Relative (%)'] = (summary_df['Mean SHAP Interaction Value'] / total) * 100
summary_df.to_csv(f"SHAP_Interaction_{target}.csv",index=False)

print(summary_df)

# Custom color for each interaction type label
custom_palette = {
    'Additive Main': 'royalblue',
    'Dominance Main': 'purple',
    'Additive × Additive': 'skyblue',
    'Additive × Dominance': 'mediumorchid',
    'Dominance × Dominance': 'indigo'
}
      
  # Plot
plt.figure(figsize=(10, 6))
ax = sns.barplot(
    data=summary_df,
    x='Interaction Type',
    y='Mean SHAP Interaction Value',
    hue='Interaction Type',        # Add this line to assign palette correctly
    palette=custom_palette,
    dodge=False                    # Prevents bars from splitting by hue
)


# Customize plot
plt.title('Mean SHAP Interaction Value by Genetic Effect Type', fontsize=14)
plt.xticks(rotation=45, ha='right')
plt.xlabel('')
plt.ylabel('Mean SHAP Interaction Value')
plt.tight_layout()
plt.savefig(f"SHAP_Interaction_{target}.png")
plt.close()

