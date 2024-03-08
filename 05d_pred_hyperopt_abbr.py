# Imports ---------------------------
import numpy as np
import pandas as pd
import numpy as np
import os
import pickle
import funcs.amyloid as amyloid

# ML imports ---------------------------
from sklearn.preprocessing import StandardScaler
from hyperopt import fmin, tpe, hp, STATUS_OK, Trials
from hpsklearn import HyperoptEstimator
from hpsklearn import random_forest_classifier, gaussian_nb, \
        k_neighbors_classifier, xgboost_classification, \
        svc, decision_tree_classifier, mlp_classifier

# This script uses Hyperopt for bayesian hyperparameter optimization
# for classification of subgroups from lab values

# Vars ---------------------------
NJOBS = -1
RANDOM_STATE = 122
CV = 5
PROCESSED_DIR = "data/processed"
OUT_DIR = "data/subgroup_prediction_abbr"

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
# Run Hyperopt
# --------------------------------
if __name__ == "__main__":
    names = ["mlp","rf","nb","knn","xgb","svc","tree"]

    classifiers_to_use = [
        mlp_classifier("mlp"),
        random_forest_classifier("rf"),
        gaussian_nb("nb"),
        k_neighbors_classifier("knn"),
        xgboost_classification("xgb"),
        svc("svc"),
        decision_tree_classifier("tree")
    ]

    for idx, classifier in enumerate(classifiers_to_use):
        print("---------------- Running for {} ----------------".format(names[idx]))

        for feature_set in Xd.keys():
            print("---------------- Feature-set: {} ----------------".format(feature_set))

            estim = HyperoptEstimator(
                classifier=classifier,
                preprocessing=[],
                algo=tpe.suggest,
                max_evals=100,
                trial_timeout=120,
                seed=RANDOM_STATE,
                n_jobs=-1
            )
            
            # Search the hyperparameter space based on the data
            estim.fit(Xd[feature_set], y)
            
            # Show the results
            print("Full training accuracy: ", estim.score(Xd[feature_set], y))
            print(estim.best_model())

            # save
            with open(os.path.join(OUT_DIR, "hyperopt_{}_{}.pkl".format(names[idx], feature_set)),"wb") as f:
                pickle.dump(estim, f)