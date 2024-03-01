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
OUT_DIR = "data/subgroup_prediction/shap"
MODELS = "data/subgroup_prediction/*.pkl"
METRICS_FILE = "tables/ml_metrics.xlsx"

os.makedirs(OUT_DIR, exist_ok=True)

# --------------------------------
# Load Data
# --------------------------------
data_df = pd.read_csv(os.path.join(PROCESSED_DIR, "AL_with_ccp_03.tsv"), sep='\t', index_col=0).rename(columns=amyloid.ddict_unclean)

# Load top-performing imputed dataset
Xi_mice = pd.read_csv("data/imputed/full_dataset/mice_qvars_02.tsv", sep="\t")
Xi_mice = Xi_mice.loc[Xi_mice.join(data_df['fna3_cluster_n']).dropna(subset='fna3_cluster_n').index]

X = pd.DataFrame(StandardScaler().fit_transform(Xi_mice.values), index=Xi_mice.index, columns=Xi_mice.columns)
y = data_df.loc[Xi_mice.index, 'fna3_cluster_n'].map({'Low':0,'Intermediate':1, 'High': 2})

# Drop due to high missingness
to_drop = ["Amyloid type","Secondary organ","Arrhythmia ","(Autonomic)",
           "(Peripheral)","SIFE M-component","UIFE M-component",
           "Education","Abdominal fat pad CR staining", "Bone marrow CR staining"]

# Create filtered dataset (for complete categoricals)
X_filt = utils.data_formatter(
    X.join(data_df.loc[Xi_mice.index, amyloid.catvars]).drop(columns=to_drop).dropna(),
    verbose=False
)
y_filt = data_df.loc[X_filt.index, 'fna3_cluster_n'].map({'Low':0,'Intermediate':1, 'High': 2})

# --------------------------------
# Load Models
# --------------------------------
hyperopt_estimators = dict()
hyperopt_models = dict()

for filename in glob.glob("data/subgroup_prediction/*.pkl"):
    with open(filename, "rb") as file:
        _name = filename.split("/")[-1].split(".pkl")[0].split("hyperopt_")[1]
        hyperopt_estimators[_name] = pickle.load(file)
        hyperopt_models[_name] = hyperopt_estimators[_name].best_model()['learner']
    
if __name__ == "__main__":
    # -----------------------------------
    # Calculate cross validation metrics
    # -----------------------------------
    metrics_df = utils.cv_scorer(
        X_filt, 
        y_filt,
        {key: hyperopt_models[key] for key in hyperopt_models.keys() if "categ" in key},
        random_state=RANDOM_STATE,
        folds=CV
    )

    metrics_filt_df = utils.cv_scorer(
        X, 
        y,
        {key: hyperopt_models[key] for key in hyperopt_models.keys() if "categ" not in key},
        random_state=RANDOM_STATE,
        folds=CV
    )

    pd.concat([metrics_df, metrics_filt_df]).groupby(['Classifier','metric']).agg(lambda x: utils.mu_ci(x)).rename(
        columns={'0':'Low','1':'Intermediate','2':'High'}).drop(
        columns='k').sort_index(level=['Classifier','metric']).sort_values("accuracy").to_excel(METRICS_FILE)
    
    # -----------------------------------
    # Compute shapley values
    # -----------------------------------
    shap_explainers = dict()
    shap_values = dict()

    for mod in hyperopt_models.keys():
        print("------------ Running for : {} ------------".format(mod))
        if "categ" in mod:
            hyperopt_models[mod].fit(X_filt, y_filt)
            shap_explainers[mod] = shap.Explainer(hyperopt_models[mod].predict, X_filt, seed=RANDOM_STATE)
            shap_values[mod] = shap_explainers[mod](X_filt)
        else:
            hyperopt_models[mod].fit(X, y)
            shap_explainers[mod] = shap.Explainer(hyperopt_models[mod].predict, X, seed=RANDOM_STATE)
            shap_values[mod] = shap_explainers[mod](X)

        # Save files
        with open(os.path.join(OUT_DIR, "shap_values_{}.pkl".format(mod)),"wb") as f:
            pickle.dump(shap_values[mod], f)