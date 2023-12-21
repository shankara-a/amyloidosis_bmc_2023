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
    
#----------------------------
# Dimensionality Reduction
#----------------------------
def get_pcs(X: np.array, normalize: bool = True, n_components: int = 5, **pca_kwargs) :
    """Get principal components form an input array

    Args:
        X (np.array): Input Array (samples x features)
        normalize (bool, optional): Whether to z-score input data for gaussian distribution. Defaults to True.
        n_components (int, optional): Number of principal components to use. Defaults to 5.

    Returns:
        (tuple): transformed values, PCA object, feature names
    """
    if normalize:
        X = scipy.stats.zscore(X, axis=0)
    
    pca = sklearn.decomposition.PCA(n_components=n_components, **pca_kwargs)
    pca.fit(X.T)
    P = pca.transform(X.T)
    P_df = pd.DataFrame(P, index=X.columns)

    return P_df, pca, X.index.values

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
