# -- import packages: --------------------------------------------------------------------
import pandas as pd
import numpy as np
import scipy
import sklearn
from tqdm import tqdm

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

#----------------------------
# Functions
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