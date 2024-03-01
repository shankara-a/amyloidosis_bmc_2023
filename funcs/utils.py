# -- import packages: --------------------------------------------------------------------
import pandas as pd
import numpy as np
import scipy
import sklearn
from tqdm import tqdm
from typing import Union

# -- import packages: --------------------------------------------------------------------
import funcs.amyloid as amyloid

#----------------------------
# Dataset specific
#----------------------------
def map_encoding(x, d):
    """Map encoding."""
    try:
        return d[x]
    except:
        return None
    
def get_eras(date):
    """Get era.

    Based on Staron et al.
    Blood  Cancer  Journal          
        (2021) 11:139  ; https://doi.org/10.1038/s41408-021-00529-w
    
    """
    if date <= pd.to_datetime("01-01-2000"):
        return "Era_1-2"
    elif date <= pd.to_datetime("01/01/2010"):
        return "Era_3"
    else:
        return "Era_4"
    
def data_formatter(
        X: pd.DataFrame, 
        collapse_cats: bool = True, 
        collapse_race: bool = True,
        verbose: bool = True
        ) -> pd.DataFrame:
    """Data formatter

    Sets up amyloid data for ML methods.

    Args:
        X (pd.DataFrame): _description_
        collapse_cats (bool, optional): _description_. Defaults to True.
        collapse_race (bool, optional): _description_. Defaults to True.
        verbose (bool, optional): _description_. Defaults to True.

    Returns:
        pd.DataFrame: _description_
    """
    from sksurv.preprocessing import OneHotEncoder

    X = X.copy()
    
    # Collapse Race
    if collapse_race:
        X["Race"] = X["Race"].apply(lambda x: "Other" if x in ['American_Indian_Alaska_Native','Multiracial','Native_Hawaiian_Pacific', 'Unknown/other'] else x)

    # Add computed BU Staging
    # X["BU Stage (Computed)"] = X.apply(lambda row: assign_bu_stage(row),1)

    # Split to categorical variables
    categorical_cols = X.loc[:,X.dtypes == "object"].columns
    numeric_cols = X.loc[:,X.dtypes != "object"].columns

    # Collapse symptom and system involvement variables to whether or not patient experienced them
    if collapse_cats:
        for cat in np.intersect1d(categorical_cols, np.array(amyloid.amyloid_ros + amyloid.amyloid_symptoms)):
            X[cat] = X[cat].apply(lambda x: True if x in ['involved','yes'] else False)

    # One-hot encode all categorical variables
    X = OneHotEncoder().fit_transform(X.loc[:,categorical_cols].astype("category")).join(X.loc[:,numeric_cols])

    if verbose:
        print("Using {} quantitative variables.".format(numeric_cols.shape[0]))
        print("Using {} categorical variables.".format(categorical_cols.shape[0]))
        print("Total samples x feaures (one-hot encoded): ({} x {})".format(X.shape[0], X.shape[1]))

    return X

def create_imputed_tableone(imp:dict, outfile: str, dataset: str = "full_dataset"):
    """Create imputed tableone.

    Args:
        imp (dict): dictionary of imputed results
        dataset (str, optional): dataset to create for. Defaults to "full_dataset".
    """
    from pandas.api.types import CategoricalDtype
    from tableone import TableOne
    
    best_mice_run = imp[dataset]['mse'][[x for x in imp[dataset]['mse'] if x.startswith("mice")]].sum(0).idxmin()

    # Labels for each imputed dataset
    imp[dataset]["X"]["data"] = "Original"
    imp[dataset]["Xi_median"]["data"] = "Median"
    imp[dataset]["Xi_knn"]["data"] = "KNN (K=5)"
    imp[dataset]["Xi_mice"][int(best_mice_run.split("_")[-1])]["data"] = "MICE"

    # Add mean square errors
    imp[dataset]["X"]["MSE"] = 0
    imp[dataset]["Xi_median"]["MSE"] = imp[dataset]["mse"].sum(0)["median_imputed"]
    imp[dataset]["Xi_knn"]["MSE"] = imp[dataset]["mse"].sum(0)["knn_imputed"]
    imp[dataset]["Xi_mice"][int(best_mice_run.split("_")[-1])]["MSE"] = imp[dataset]["mse"].sum(0)[best_mice_run]
    
    # Order by missing values
    s = imp[dataset]["X"].isna().sum() / imp[dataset]["X"].shape[0]
    columns = list(s.sort_values().index)

    X_full = pd.concat([
        imp[dataset]["X"], 
        imp[dataset]["Xi_median"], 
        imp[dataset]["Xi_knn"],
        imp[dataset]["Xi_mice"][int(best_mice_run.split("_")[-1])]])
    
    X_full["data"] = X_full["data"].astype(CategoricalDtype(categories=["Original", "MICE", "KNN (K=5)", "Median"], ordered=True))

    # Create table one
    nonnormal = []
    categorical = []
    groupby = ["data"]

    mytable = TableOne(
        X_full.reset_index(), 
        columns, 
        categorical, 
        groupby,
        nonnormal,
        pval=True, 
        overall=False, 
        decimals = {'MSE':4,'WBC':2, 'Hemoglobin':2, 'Troponin': 3, 'Calcium':2, 
                    'Bone marrow plasma cells (%)':2, 'Uric acid':2, 
                    'Albumin':2, 'kappa:lambda ratio':2},
        rename=amyloid.tableone_names
    )

    mytable.to_excel(outfile)

#----------------------------
# Medical
#----------------------------
def compute_egfr(sCr: float, age: float, is_female: bool) -> float:
    """
    Compute eGFR based on CKD 2021
    
    The CKD-EPI equation, expressed as a single equation, is:

        GFR = 141 * min(Scr/κ,1)α * max(Scr/κ, 1)-1.209 * 0.993Age * 1.018 [if female] * 1.159 [if black]

    New England Journal of Medicine 2021 November 4, 385 (19): 1737-1749

    Args:
        sCr (float): input creatinine
        age (float): age
        is_female (bool): flag if patient is female

    Returns:
        (float): _description_
    """
    if is_female:
        A = 0.7

        if sCr <= 0.7:
            B = -0.241
        else:
            B = -1.2

        return 142 * np.power((sCr / A),B) * np.power(0.9938,age) * 1.012
    else:
        A = 0.9

        if sCr <= 0.9:
            B = -0.302
        else:
            B = -1.2

        return 142 * np.power((sCr / A),B) * np.power(0.9938,age)
    
def assign_bu_stage(row):
    """"
    Assign BU Staging (2019):

        + 1 point if Troponin >= 0.1 ng/ml
        + 1 point if BNP >= 81 pg/ml
        + 1 point if BNP >= 700 pg/ml

    Stage I: 0 points
    Stage II: 1 point
    Stage III: 2 points
    Stage IIIb: 3 points

    See: Blood vol. 133,3 (2019): 215-223. doi:10.1182/blood-2018-06-858951
    """
    score = 0
    if row["Troponin"] >= 0.1:
        score += 1
    if row["BNP"] >= 81:
        score +=1
    if row["BNP"] >= 700:
        score +=1 
    
    return {0:"stage I",1:"stage II",2:"stage III",3:"stage IIIb"}[score]
    
def get_palladini_renal_stage(row):
    """
    Assign Renal Staging.

    Based on Palladini et al, Bloo 2014
    10.1182/blood-2014-04-570010

    # Renal stage I: both proteinuria <= 5 g/24 h and eGFR >= 50 mL/min per 1.73 m2. 
    # Renal stage II: either proteinuria > 5 g/24 h or eGFR < 50 mL/min per 1.73 m2. 
    # Renal stage III: both proteinuria > 5 g/24 h and eGFR < 50 mL/min per 1.73 m2.

    """
    if row["24-hr UTP"] <= 5000 and row["eGFR"] >= 50:
        return "Stage I"
    elif row["24-hr UTP"] > 5000 and row["eGFR"] < 50:
        return "Stage III"
    elif row["24-hr UTP"] > 5000 or row["eGFR"] < 50:
        return "Stage II"
    else:
        return None

def get_median_os(data_df, duration="OS (yr)", event="status", groupby=None):
    """_summary_

    Args:
        data_df (_type_): _description_
    """
    from lifelines.utils import median_survival_times
    from lifelines import KaplanMeierFitter

    if groupby is None:
        kmf = KaplanMeierFitter()
        kmf.fit(durations = data_df[duration], event_observed = data_df[event])
        return kmf.median_survival_time_
    else:
        result = {}
        for group in np.unique(data_df[groupby]):
            _data_df = data_df[data_df[groupby]==group]
            kmf = KaplanMeierFitter()
            kmf.fit(durations = _data_df[duration], event_observed = _data_df[event])
            result[group] = kmf.median_survival_time_
        return result
    
def get_time_eskd(row, start_time: str = "Date of diagnosis"):
    """Return time to ESKD in years.

    Censor by last encounter (either death or last visit).

    Args:
        row (pd.Series): pd.DataFrame row with
            Date of RRT Start
            treatment_eskd (whether or not patient was placed on RRT)
            Date of death
            Date of last visit
        start_time (str): column to use for start time
    """
    import pandas as pd

    if not pd.isnull(row["Date of RRT Start"]):
        end_time = row["Date of RRT Start"]
    else:
        if row["treatment_eskd"] == 1:
            return pd.NaT
        
        if not pd.isnull(row['Date of death']):
            end_time = row["Date of death"]
        else:
            end_time = row["Date of last visit"]

    time = end_time - row[start_time] 
    return time / pd.Timedelta(days=365.25)

# Competing risk
def get_cr_event(row: pd.Series, risk_1="treatment_eskd", risk_2="status") -> int:
    """Get competing risk event.

    Args:
        row (pd.Series): Series
        risk_1 (str, optional): Risk 1 (will be assigned 1). Defaults to "treatment_eskd".
        risk_2 (str, optional): Risk 2 (will be assigned 2). Defaults to "status".

    Returns:
        int: _description_
    """
    if row[risk_1] == 1:
        return 1
    elif row[risk_2] == 1:
        return 2
    else:
        return 0

def get_cr_time(row: pd.Series, start_time="Date of diagnosis"):
    """Get competing risk time

    Args:
        row (pd.Series): _description_
        start_time (str, optional): _description_. Defaults to "Date of diagnosis".

    Returns:
        _type_: _description_
    """
    import pandas as pd
    
    if row['treatment_eskd'] == 1:
        cr_time = row['Date of RRT Start'] - row[start_time]
    elif row['status'] == 1:
        cr_time = row['Date of death'] - row[start_time]
    else:
        cr_time = row['Date of last visit'] - row[start_time]
    
    return cr_time / pd.Timedelta(days=365.25)
    
#----------------------------
# Dimensionality Reduction
#----------------------------
def get_pcs(X_df: pd.DataFrame, normalize: bool = True, n_components: int = 5, **pca_kwargs) :
    """Get principal components form an input array

    Args:
        X_df (pd.DataFrame): Input Dataframe (samples x features)
        normalize (bool, optional): Whether to standard scale input data for gaussian distribution. Defaults to True.
        n_components (int, optional): Number of principal components to use. Defaults to 5.

    Returns:
        (tuple): transformed values, PCA object, feature names
    """
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    if normalize:
        X = StandardScaler().fit_transform(X_df.values)
    else:
        X = X_df.values
    
    pca = PCA(n_components=n_components, **pca_kwargs)
    P = pca.fit_transform(X)
    P_df = pd.DataFrame(P, index=X_df.index)

    return P_df, pca, X_df.columns.values

def get_umap(X_df: pd.DataFrame, normalize: bool = True, **umap_kwargs):
    """Get UMAP projection

    Args:
        X_df (pd.DataFrame): _description_
        normalize (bool, optional): _description_. Defaults to True.

    Returns:
        _type_: _description_
    """
    import umap
    from sklearn.preprocessing import StandardScaler

    if normalize:
        X = StandardScaler().fit_transform(X_df.values)
    else:
        X = X_df.values

    embedding = umap.UMAP(**umap_kwargs).fit_transform(X)
    return pd.DataFrame(embedding, columns=["umap1","umap2"], index=X_df.index)

#----------------------------
# Stats
#----------------------------
def load_and_aggregate_imputation(original_file: str, mice_imputed_dir: str, knn_neighbors: int = 5):
    """Load & Aggregate Imputations

    Args:
        original_file (str): _description_
        mice_imputed_dir (str): _description_
        knn_neighbors (int, optional): _description_. Defaults to 5.
    """
    from sklearn.impute import KNNImputer
    import glob
    import os

    # Includes only continuous variables
    # Comparison of multiple imputation methods
    X = pd.read_csv(original_file, sep='\t', index_col=0).rename(columns=amyloid.ddict_unclean)

    # Median imputation
    Xi_median = X.fillna(X.median())

    # KNN Imputation
    imputer = KNNImputer(n_neighbors=knn_neighbors)
    Xi_knn = pd.DataFrame(
        imputer.fit_transform(X), 
        index=X.index, 
        columns=X.columns
    )

    # MICE imputation
    Xi_mice_dict = {}

    for mice_result in glob.glob(os.path.join(mice_imputed_dir, "mice_qvars*.tsv")):
        Xi_mice_dict[int(mice_result.split("/")[-1].split("_")[-1].split(".tsv")[0])] = pd.read_csv(
            mice_result, sep="\t").rename(columns={'X24_hr_UTP':'24_hr_UTP'}).rename(columns=amyloid.ddict_unclean)

    # Comparison of imputation to mean
    Xi_comparison = pd.DataFrame(X.mean(), index=X.columns, columns=['original']).join(
        pd.DataFrame(Xi_median.mean(), index=X.columns, columns=['median_imputed'])).join(
        pd.DataFrame(Xi_knn.mean(), index=X.columns, columns=['knn_imputed'])).join(
        pd.concat(
            [pd.DataFrame(Xi_mice_dict[k].mean(), index=X.columns, columns=['mice_imputed_{}'.format(k)]) for k,v in Xi_mice_dict.items()],
            axis=1
        )
    )

    # Mean Square AERror
    err = (Xi_comparison.T - Xi_comparison['original']) / Xi_comparison['original']
    mse = err.T**2

    result = {}
    result["X"] = X
    result["Xi_median"] = Xi_median
    result["Xi_knn"] = Xi_knn
    result["Xi_mice"] = Xi_mice_dict
    result["mean_comparison"] = Xi_comparison
    result["mse"] = mse

    # Save files
    Xi_median.rename(columns=amyloid.ddict_clean).to_csv(os.path.join(mice_imputed_dir, "median_qvars_01.tsv"), sep="\t")
    Xi_knn.rename(columns=amyloid.ddict_clean).to_csv(os.path.join(mice_imputed_dir, "knn_qvars_01.tsv"), sep="\t")

    return result

def compute_ranksum(in_group: pd.DataFrame, out_group: pd.DataFrame, feature: str):
    """Compute Rank-Sum p-values.

    Args:
        in_group (pd.DataFrame): _description_
        out_group (pd.DataFrame): _description_
        feature (str): _description_

    Returns:
        _type_: _description_
    """
    from scipy.stats import ranksums

    res = ranksums(
        in_group[feature].values, 
        out_group[feature].values,
        alternative='two-sided',
        nan_policy="omit"
    )
    return res.pvalue, res.statistic

def compute_ttest(in_group: pd.DataFrame, out_group: pd.DataFrame, feature: str):
    """Compute T-test.
    Defaults to Welch's T-Test.

    Args:
        in_group (pd.DataFrame): _description_
        out_group (pd.DataFrame): _description_
        feature (str): _description_

    Returns:
        _type_: _description_
    """
    from scipy.stats import ttest_rel
    from scipy.stats import ttest_ind

    # continuous
    res = ttest_ind(
        in_group[feature].values, 
        out_group[feature].values,
        alternative='two-sided',
        nan_policy="omit",
        equal_var= False
    )
    return res.pvalue, res.statistic

def get_contrasts(
        X: pd.DataFrame, 
        idx:str, 
        qvars, 
        fdr_alpha:float=0.5, 
        fdr_method:str="fdr_bh"
    ):
    """_summary_

    Args:
        X (pd.DataFrame): _description_
        idx (str): _description_
        qvars (_type_): _description_
        fdr_alpha (float, optional): _description_. Defaults to 0.5.
        fdr_method (str, optional): _description_. Defaults to "fdr_bh".

    Returns:
        _type_: _description_
    """
    from statsmodels.stats.multitest import multipletests

    results = list()

    for k in tqdm(np.unique(X[idx].dropna())):
        in_group = X[X[idx]==k]
        out_group = X[X[idx]!=k]

        _results = {}

        for m in qvars:
            d = {}
            d["p_value"], d["statistic"] = compute_ranksum(in_group, out_group, m)
            _results[m] = d

        _results_df = pd.DataFrame.from_dict(_results).T.sort_values("p_value")
        _,_results_df['fdr_bh'],_,_ = multipletests(_results_df["p_value"], alpha=fdr_alpha, method=fdr_method)
        _results_df[idx] = k
        results.append(_results_df)
    
    return pd.concat(results)

def fisher_exact(
        X: pd.DataFrame,
        meta_s: pd.Series,
        fdr_alpha: float = 0.05,
        fdr_method: str = 'fdr_bh',
        **kwargs
    ):
    """
    Fisher Exact Test.
    -------------------
    Performs fisher exact test by comparing proportions of binary features for full
    population and cluster specific populations:

        (present vs absent)
        (within cluster vs outside cluster)

    Args:
        X (pd.DataFrame): input one-hot encoded matrix
        meta_s (pd.Series): series with metadata to cluster by
        fdr_alpha (float, optional): Defaults to 0.05.
        fdr_method (str, optional): Defaults to 'fdr_bh'.
        **kwargs: for FE

    Returns:
        _type_: pd.DataFrame with pval, adj-pval, and odds ratio
    """
    from statsmodels.stats.multitest import multipletests
    import scipy.stats as stats

    # Dropps missing values
    meta_s = meta_s.dropna()

    # Filter input encodings by missing clusters
    X = X.loc[meta_s.index,:].copy()

    # Present & Absent
    X_present = X.sum(0)
    X_absent = X.shape[0] - X.sum(0)

    # Within cluster & Out cluster
    X_ci = X.join(meta_s).groupby(meta_s.name).sum()[X.columns]
    X_co = X_present - X_ci

    # Initialize results
    pval = X_ci.T.copy()
    pval_adj = X_ci.T.copy()
    odds_r = X_ci.T.copy()

    # Perform fisher exact test
    for clust in X_ci.index:
        for feat in X_ci.columns:
            odds_r.loc[feat,clust],pval.loc[feat,clust] = stats.fisher_exact([[X_present[feat], X_absent[feat]], [X_ci.loc[clust,feat], X_co.loc[clust,feat]]], alternative='less', **kwargs)

        _,pval_adj[clust],_,_ = multipletests(pval[clust], alpha=fdr_alpha, method=fdr_method)

    # Melt results
    pval_adj = pd.melt(pval_adj.reset_index(), id_vars=['index'], value_vars=pval_adj.columns).rename(
        columns={'index':'feat', 'value':'pval_adj'}).set_index(['feat',meta_s.name])

    pval = pd.melt(pval.reset_index(), id_vars=['index'], value_vars=pval.columns).rename(
        columns={'index':'feat', 'value':'pval'}).set_index(['feat',meta_s.name])

    odds_r = pd.melt(odds_r.reset_index(), id_vars=['index'], value_vars=odds_r.columns).rename(
        columns={'index':'feat', 'value':'odds_r'}).set_index(['feat',meta_s.name])

    return pval.join(pval_adj).join(odds_r).reset_index()

def run_statistical_comparisons(data_df, var, out_dir, tag=""):
    """Run statistical comparisons
        Wilcoxon Rank Sum Tests for numerical
        Fisher's exact test for categorical

    Args:
        data_df (_type_): raw data
        var (_type_): clustering to use
        out_dir (_type_): output directory
    """
    import os
    from sksurv.preprocessing import OneHotEncoder

    # Run rank-sum tests
    contrasts_df = get_contrasts(data_df, var, amyloid.qvars)
    contrasts_df.to_csv(os.path.join(out_dir, "contrasts_qvars{}.tsv".format(tag)), sep="\t")

    # Fisher exacts
    # Drop Age & vars due to high-missingness
    to_drop = [
        "Age","Amyloid type","Secondary organ","Arrhythmia ","(Autonomic)",
        "(Peripheral)","SIFE M-component","UIFE M-component",
        "Education", "Abdominal fat pad CR staining", "Bone marrow CR staining"
    ]

    catvars = list(set(amyloid.catvars)-set(to_drop))

    # If uncertain or equivocal, do not include in fisher exact
    _cat_df = data_df[catvars].replace({
        "uncertain":np.nan,
        "equivocal":np.nan,
        "involved":"yes",
        "not_involved":"no"})

    # Collapse race
    _cat_df["Race"] = _cat_df["Race"].apply(lambda x: "Other" if x in ['Multiracial','Native_Hawaiian_Pacific', 'Unknown/other'] else x)
    _cat_df = _cat_df.astype("category")

    # One hot encoding
    _cat_df = OneHotEncoder().fit_transform(_cat_df)
    contrasts_fe_df = fisher_exact(_cat_df, data_df[var])
    contrasts_fe_df['feat'] = contrasts_fe_df['feat'].str.replace("=yes","")
    contrasts_fe_df.to_csv(os.path.join(out_dir, "contrasts_fe{}.tsv".format(tag)), sep="\t")

#----------------------------
# Machine Learning
#----------------------------
def compute_permutation_importance(model, X_test, y_test, n_repeats = 15, random_state = 123):
    """_summary_

    Args:
        model (_type_): _description_
        X_test (_type_): _description_
        y_test (_type_): _description_
        n_repeats (_type_, optional): _description_. Defaults to None.
        random_state (int, optional): _description_. Defaults to 123.
    """
    from sklearn.inspection import permutation_importance

    result = permutation_importance(model, X_test, y_test, n_repeats=n_repeats, random_state=random_state)

    return pd.DataFrame(
        {
            k: result[k]
            for k in (
                "importances_mean",
                "importances_std",
            )
        },
        index=X_test.columns,
    ).sort_values(by="importances_mean", ascending=False)

def get_agg_clust(X: pd.DataFrame, n_clust: int, metric: str = "euclidean", linkage: str = "average"):
    """_summary_

    Args:
        X (pd.DataFrame): _description_
        n_clust (int): _description_
        affinity (str, optional): _description_. Defaults to "euclidean".
        linkage (str, optional): _description_. Defaults to "average".
    """
    from sklearn.cluster import AgglomerativeClustering

    # Agglomerative Clustering
    cluster = AgglomerativeClustering(n_clusters=n_clust, metric=metric, linkage=linkage)

    # Derive clusters
    clusters = cluster.fit_predict(X.values)

    return pd.DataFrame(clusters+1, index=X.index, columns=['agg_clust_{}'.format(n_clust)])

def load_ccp_result(output_dir: str) -> dict:
    """Load consensus cluster result into python

    Requires rpy2

    Args:
        output_rds (str): output of .rds file from consensus clsuter plus
        input_path (str): input clustering matrix to extract sample ids

    Returns:
        dict: dicitonary of results
    """
    import rpy2.robjects as robjects
    from rpy2.robjects import pandas2ri
    import os

    pandas2ri.activate()
    
    # Load RObject
    readRDS = robjects.r['readRDS']
    X = readRDS(os.path.join(output_dir, "ccp.rds"))
    result = {}

    # Load sample ids
    sample_id = pd.read_csv(os.path.join(output_dir, "input_matrix.tsv"), sep="\t", index_col=0).index
    
    for i in range(1,len(X)):
        results_i = {}
        
        # Consensus matrix
        cm_df = pd.DataFrame(X[i][0], index=sample_id, columns=sample_id)
        results_i['cm'] = cm_df
        
        # Annot (ie cluster)
        annot_s = pd.Series(X[i][2])
        annot_s.index = sample_id
        results_i['annot'] = annot_s
        result[i] = results_i
    
    return result

from typing import Union
def compute_ari(s1: Union[pd.DataFrame, pd.Series, list], s2: Union[pd.Series, list] = None):
    """Compute Adjusted Rand Index

    Args:
        s1 (Union[pd.DataFrame, pd.Series, list]): _description_
        s2 (Union[pd.Series, list], optional): _description_. Defaults to None.

    Returns:
        float if passed two series or lists
        pd.Dataframe if passed 
    """
    from sklearn.metrics.cluster import adjusted_rand_score
    from itertools import combinations

    def func(a,b):
        overlap = np.intersect1d(a.dropna().index, b.dropna().index)
        return adjusted_rand_score(a[overlap],b[overlap])
    
    if isinstance(s1, pd.DataFrame):
        a = np.zeros((s1.shape[1], s1.shape[1]))
        a[np.triu_indices(s1.shape[1], k=1)] = [func(s1[a],s1[b]) for a,b in combinations(s1.columns, r=2)]
        a += a.T
        np.fill_diagonal(a,1)
        return pd.DataFrame(a, index=s1.columns, columns=s1.columns)
    else:
        assert s2 is not None, "Pass two lists or series to compute adjusted rand index!"
        return func(s1,s2)
    
def cv_scorer(
    X: pd.DataFrame, 
    y: pd.DataFrame, 
    models: dict, 
    folds: int = 5, 
    random_state: int = 42):
    """Cross validation scoring.

    Args:
        X (pd.DataFrame): _description_
        y (pd.DataFrame): _description_
        models (dict): _description_
        folds (int, optional): _description_. Defaults to 5.
        random_state (int, optional): _description_. Defaults to 42.

    Returns:
        _type_: _description_
    """
    from sklearn.model_selection import KFold
    from sklearn.metrics import classification_report
    from sklearn.metrics import cohen_kappa_score

    # Initialize fold split
    kf = KFold(n_splits=folds, shuffle=True, random_state=random_state)

    # Results
    results_df = list()

    for key in models.keys():
        for i, (train_index, test_index) in enumerate(kf.split(X)):
            models[key].fit(X.iloc[train_index,:], y.iloc[train_index])
            y_pred = models[key].predict(X.iloc[test_index,:])

            # Generate classificaiton report
            cr = classification_report(y.iloc[test_index], y_pred, output_dict=True)
            cr = pd.DataFrame.from_dict(cr).drop(
                columns=['macro avg','weighted avg'], 
                index=['support']).reset_index()
            
            # Add kappa score
            cr['kappa'] = cohen_kappa_score(y.iloc[test_index], y_pred)

            # Naming
            cr['Classifier'] = key
            cr['k'] = i
            cr = cr.rename(columns={'index':'metric'})

            results_df.append(cr)
    
    return pd.concat(results_df)

def mu_ci(data, confidence:float = 0.95):
    """Create string of mean & CI for tables.

    Args:
        data (_type_): _description_
        confidence (float, optional): _description_. Defaults to 0.95.

    Returns:
        _type_: _description_
    """
    from scipy.stats import sem,t
    
    n = len(data)
    m, se = np.mean(data), sem(data)
    h = se * t.ppf((1 + confidence) / 2, n - 1)
    return "{:.3f} ({:.2f}-{:.2f})".format(m, m - h, m + h)