# -- import packages: --------------------------------------------------------------------
import pandas as pd
import numpy as np
import scipy
import sklearn
from tqdm import tqdm

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
    X["BU Stage (Computed)"] = X.apply(lambda row: assign_bu_stage(row),1)

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

def load_ccp_result(output_rds: str, input_path: str) -> dict:
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
    pandas2ri.activate()
    
    # Load RObject
    readRDS = robjects.r['readRDS']
    X = readRDS(output_rds)
    result = {}

    # Load sample ids
    sample_id = pd.read_csv(input_path, sep="\t", index_col=0).index
    
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
