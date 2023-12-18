# -- import packages: --------------------------------------------------------------------
import pandas as pd
import numpy as np
import scipy
import sklearn

def map_encoding(x, d):
    """Map encoding."""
    try:
        return d[x]
    except:
        return None

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
        X = scipy.stats.zscore(axis=0)
    
    pca = sklearn.decomposition.PCA(n_components=n_components, **pca_kwargs)
    pca.fit(X.T)
    P = pca.transform(X.T)
    P_df = pd.DataFrame(P, index=X.columns)

    return P_df, pca, X.index.values

