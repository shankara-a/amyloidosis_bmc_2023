# Imports ---------------------------
import numpy as np
import pandas as pd
import numpy as np
import os
import pickle
from funcs import utils
import funcs.amyloid as amyloid
import glob

# ML imports ---------------------------
from sklearn.preprocessing import StandardScaler
import shap

# This script computes cross validaiton metrics
# and shapley values

# Vars ---------------------------
NJOBS = -1
RANDOM_STATE = 122
CV = 5
PROCESSED_DIR = "data/processed"
OUT_DIR = "data/subgroup_prediction_abbr/shap"
MODELS = "data/subgroup_prediction_abbr/*.pkl"
METRICS_FILE = "tables/ml_metrics_abbr.xlsx"

os.makedirs(OUT_DIR, exist_ok=True)

# --------------------------------
# Load Data
# --------------------------------
data_df = pd.read_csv(os.path.join(PROCESSED_DIR, "AL_with_ccp_03.tsv"), sep='\t', index_col=0).rename(columns=amyloid.ddict_unclean)

# Load top-performing imputed dataset
Xi_mice_v1 = pd.read_csv(os.path.join(PROCESSED_DIR, "Xi_abbr_echo.tsv"), sep='\t', index_col=0)
Xi_mice_v2 = pd.read_csv(os.path.join(PROCESSED_DIR, "Xi_abbr.tsv"), sep='\t', index_col=0)

Xi_mice_v1 = Xi_mice_v1.loc[Xi_mice_v1.join(data_df['fna3_cluster_n']).dropna(subset='fna3_cluster_n').index]
Xi_mice_v2 = Xi_mice_v2.loc[Xi_mice_v2.join(data_df['fna3_cluster_n']).dropna(subset='fna3_cluster_n').index]

X_v1 = pd.DataFrame(StandardScaler().fit_transform(Xi_mice_v1.values), index=Xi_mice_v1.index, columns=Xi_mice_v1.columns)
X_v2 = pd.DataFrame(StandardScaler().fit_transform(Xi_mice_v2.values), index=Xi_mice_v2.index, columns=Xi_mice_v2.columns)

Xd = {"with_echo": X_v1, "no_echo": X_v2}
y = data_df.loc[Xi_mice_v1.index, 'fna3_cluster_n'].map({'Low':0,'Intermediate':1, 'High': 2})

# --------------------------------
# Load Models
# --------------------------------
hyperopt_estimators = dict()
hyperopt_models = dict()

for filename in glob.glob(MODELS):
    with open(filename, "rb") as file:
        _name = filename.split("/")[-1].split(".pkl")[0].split("hyperopt_")[1]
        hyperopt_estimators[_name] = pickle.load(file)
        hyperopt_models[_name] = hyperopt_estimators[_name].best_model()['learner']
    
if __name__ == "__main__":
    # -----------------------------------
    # Calculate cross validation metrics
    # -----------------------------------
    metrics_echo_df = utils.cv_scorer(
        X_v1, 
        y,
        {key: hyperopt_models[key] for key in hyperopt_models.keys() if "with_echo" in key},
        random_state=RANDOM_STATE,
        folds=CV
    )

    metrics_no_echo_df = utils.cv_scorer(
        X_v2, 
        y,
        {key: hyperopt_models[key] for key in hyperopt_models.keys() if "no_echo" in key},
        random_state=RANDOM_STATE,
        folds=CV
    )

    pd.concat([metrics_echo_df, metrics_no_echo_df]).groupby(['Classifier','metric']).agg(lambda x: utils.mu_ci(x)).rename(
        columns={'0':'Low','1':'Intermediate','2':'High'}).drop(
        columns='k').sort_index(level=['Classifier','metric']).sort_values("accuracy").to_excel(METRICS_FILE)
    
    # -----------------------------------
    # Compute shapley values
    # -----------------------------------
    shap_explainers = dict()
    shap_values = dict()

    for mod in hyperopt_models.keys():
        print("------------ Running for : {} ------------".format(mod))
        if "with_echo" in mod:
            hyperopt_models[mod].fit(X_v1, y)
            shap_explainers[mod] = shap.Explainer(hyperopt_models[mod].predict, X_v1, seed=RANDOM_STATE)
            shap_values[mod] = shap_explainers[mod](X_v1)
        else:
            hyperopt_models[mod].fit(X_v2, y)
            shap_explainers[mod] = shap.Explainer(hyperopt_models[mod].predict, X_v2, seed=RANDOM_STATE)
            shap_values[mod] = shap_explainers[mod](X_v2)

        # Save files
        with open(os.path.join(OUT_DIR, "shap_values_{}.pkl".format(mod)),"wb") as f:
            pickle.dump(shap_values[mod], f)