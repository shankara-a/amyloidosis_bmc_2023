# -- import packages: --------------------------------------------------------------------
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Union
import matplotlib.ticker as ticker

# -- import packages: --------------------------------------------------------------------
import funcs.amyloid as amyloid

#----------------------------
# Data Exploration
#----------------------------
def missing_barplot(df: pd.DataFrame, ax: Union[plt.axes, None] = None, **barplot_kwargs):
    """Missing barplot.

    Args:
        df (pd.DataFrame): input dataframe
        ax (Union[plt.axes, None], optional): matplotlib axis to include. Defaults to None.
    """
    if ax is None:
        fig,ax = plt.subplots(figsize=(3,7))

    pd.DataFrame(
        df.isna().sum().sort_values(ascending=False) / df.shape[0],
        columns=['% Missing']
    ).plot.barh(ax=ax, width=0.8, edgecolor='k', **barplot_kwargs)

    ax.legend().remove()
    ax.set_xlabel("% Missing", fontsize=12)

def plot_clustermap(X, figsize=(8,8), xlabel=None, ylabel=None, vmin=-1, vmax=1, 
                    yticklabels=True, xticklabels=True, linecolor='k', linewidths=0.5, cmap='coolwarm', **kwargs):
    """Plot clustermp."""
    import seaborn as sns

    clustermap = sns.clustermap(
        X,
        xticklabels=xticklabels,
        yticklabels=yticklabels,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        linecolor=linecolor,
        linewidths=linewidths,
        figsize=figsize,
        **kwargs
    )

    _ = clustermap.ax_heatmap.set_xticklabels(clustermap.ax_heatmap.get_xticklabels(), fontsize=6)
    _ = clustermap.ax_heatmap.set_yticklabels(clustermap.ax_heatmap.get_yticklabels(), fontsize=6)

    if xlabel is not None:
        _ = clustermap.ax_heatmap.set_xlabel(xlabel, fontsize=8)

    if ylabel is not None:
        _ = clustermap.ax_heatmap.set_ylabel(ylabel, fontsize=8)

def plot_cmatrix(result, run_name="full_na_dataset", k=3, metas=[], **kwargs):
    """_summary_

    Args:
        result (_type_): _description_
        run_name (str, optional): _description_. Defaults to "full_na_dataset".
        k (int, optional): _description_. Defaults to 3.
    """
    from scipy.cluster import hierarchy
    from scipy.spatial import distance

    cm_df = result[run_name][k-1]['cm']

    # Compute linkages
    row_linkage = hierarchy.linkage(distance.pdist(cm_df.values), method='average')
    col_linkage = hierarchy.linkage(distance.pdist(cm_df.values.T), method='average')

    # Create clustermap
    sns.clustermap(
        cm_df,
        row_linkage=row_linkage,
        col_linkage=col_linkage,
        xticklabels=[],
        yticklabels=[],
        rasterized=True,
        figsize=(6,6),
        row_colors=[x[cm_df.index] for x in metas],
        **kwargs
        )
    
def plot_feature_importance(model, figsize=(4,5), min_importance=0.01, ax=None, color='lightblue', **kwargs):
    """Plot feature importances

    From tree-based prediction mdoels

    Args:
        model (_type_): _description_
        figsize (tuple, optional): _description_. Defaults to (4,5).
        min_importance (float, optional): _description_. Defaults to 0.01.
    """
    features_df = pd.DataFrame(
        model.feature_importances_, 
        index=model.feature_names_in_, 
        columns=['importance']
    )

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    features_df[features_df['importance']>min_importance]['importance'].sort_values(ascending=True).rename(
        index=amyloid.ddict_unclean).plot.barh(
            ax=ax, 
            width=1, 
            edgecolor='k', 
            color=color,
            **kwargs
    )

    ax.set_xlabel("Feature Importance")
    _ = ax.set_yticklabels(ax.get_yticklabels(), fontsize=7)

    return features_df

#----------------------------
# Survival
#----------------------------
def plot_cumulative_dynamic_auc(x_times, rsf_auc, rsf_mean_auc, line_color="k", ax=None):
    """Plot cumulative dynamic AUC"""
    if ax is None:
        fig,ax = plt.subplots(figsize=(4,3))

    ax.plot(x_times, rsf_auc, marker="o")
    ax.axhline(rsf_mean_auc, linestyle="--", c=line_color, alpha=0.7)
    ax.set_xlabel("Time from Diagnosis (M)")
    ax.set_ylabel("Time-dependent AUC\n($\mu$ = {:.2})".format(rsf_mean_auc))

    ax.grid(True)

#----------------------------
# Dimensionality Reduction
#----------------------------
def plot_pca(P_df, pca, c=None, cohort_s=None, cohort_colors=None, cohort_args=None, order=[1,2,3], outliers=None, title='',
    vmin=None, vmax=None, alpha=1, lw=0, s=30, cmap=plt.cm.Spectral_r, cticks=None, cticklabels=None, clabel='',
    show_legend=True, show_ax2=True, axis_length=5.5, xlim=None, ylim=None):
    """
    cohort_s: Series encoding cohorts
    cohort_colors: dict
    Modes:
    """
    if cohort_s is not None:
        cohorts = cohort_s.unique()
        nc = len(cohorts)
        if cohort_colors is None and cohort_args is None:
            # cohort_colors = {i:j for i,j in zip(cohorts, cm.get_cmap(cmap, nc)(np.arange(nc)))}
            cohort_colors = {i:j for i,j in zip(cohorts, sns.husl_palette(nc, s=1, l=0.6))}
        if cohort_args is None:
            cohort_args = {}
            for k in np.unique(cohort_s):
                cohort_args[k] = {'color': cohort_colors[k], 'marker':'o', 'edgecolor':'none', 's':s}

    if show_ax2:
        fig = plt.figure(facecolor=(1,1,1), figsize=(2*axis_length,axis_length))
        ax1 = fig.add_axes(np.array([1/2*axis_length, 0.75/axis_length, 4/2*axis_length, 4/axis_length]))
    else:
        fig = plt.figure(facecolor=(1,1,1), figsize=(axis_length,axis_length))
        ax1 = fig.add_axes(np.array([1/axis_length, 0.75/axis_length, 4/axis_length, 4/axis_length]))
    
    if cohort_s is None:  # c[P_df.index]
        sa = ax1.scatter(P_df[order[1]-1], P_df[order[0]-1], c=c, cmap=cmap, vmin=vmin, vmax=vmax, lw=lw, alpha=alpha, s=s)
    else:
        for k in np.unique(cohort_s):
        # for k in cohort_s.unique():
            i = cohort_s[cohort_s==k].index
            ax1.scatter(P_df.loc[i,order[1]-1], P_df.loc[i,order[0]-1], alpha=alpha, label=k, **cohort_args[k])
    format_plot(ax1, fontsize=10)
    ax1.set_xlabel('PC {0} ({1:.2f}%)'.format(order[1], pca.explained_variance_ratio_[order[1]-1]*100), fontsize=12)
    ax1.set_ylabel('PC {0} ({1:.2f}%)'.format(order[0], pca.explained_variance_ratio_[order[0]-1]*100), fontsize=12)

    if show_ax2:
        ax2 = fig.add_axes(np.array([6/2*axis_length, 0.75/axis_length, 4/2*axis_length, 4/axis_length]))
        if cohort_s is None:
            ax2.scatter(P_df[order[2]-1], P_df[order[0]-1], c=c, cmap=cmap, vmin=vmin, vmax=vmax, lw=lw, alpha=alpha, s=s)
        else:
            for k in np.unique(cohort_s):
                i = cohort_s[cohort_s==k].index
                ax2.scatter(P_df.loc[i,order[2]-1], P_df.loc[i,order[0]-1], alpha=alpha, label=k, **cohort_args[k])
            # ax2.legend(loc=3, fontsize=10, scatterpoints=1, handletextpad=0.1, framealpha=0.5, bbox_to_anchor=(-0.5,-0.1))

        format_plot(ax2, fontsize=10)
        ax2.set_xlabel('PC {0} ({1:.2f}%)'.format(order[2], pca.explained_variance_ratio_[order[2]-1]*100), fontsize=12)
        ax2.set_ylabel('PC {0} ({1:.2f}%)'.format(order[0], pca.explained_variance_ratio_[order[0]-1]*100), fontsize=12)

    if outliers is not None:
        ax1.scatter(P_df.loc[outliers, order[1]-1], P_df.loc[outliers, order[0]-1], c='none', edgecolors='r', marker='s', lw=1, alpha=1, s=50, label=None)
        if show_ax2:
            ax2.scatter(P_df.loc[outliers, order[2]-1], P_df.loc[outliers, order[0]-1], c='none', edgecolors='r', marker='s', lw=1, alpha=1, s=50, label=None)

    fig.suptitle(title, fontsize=12)

    if cohort_s is not None and show_legend:
        # ax2.legend(loc=0, fontsize=10, scatterpoints=1, handletextpad=0.1, framealpha=0.5, bbox_to_anchor=(-0.5,-0.1))
        leg = ax1.legend(loc=0, fontsize=9, scatterpoints=1, handletextpad=0.1, framealpha=1, labelspacing=0.35)
        for lh in leg.legendHandles:
            lh.set_alpha(1)

    # if cohort_s is None and c is not None and not isinstance(c, list) and not isinstance(c, str):
    if cohort_s is None and c is not None and len(c)==P_df.shape[0]:
        if show_ax2:
            cax = fig.add_axes(np.array([3.5/2*axis_length, 5/axis_length, 1.5/2*axis_length, 0.15/axis_length]))
        else:
            cax = fig.add_axes(np.array([3.5/axis_length, 5/axis_length, 1.5/axis_length, 0.15/axis_length]))
        # cax = fig.add_axes(np.array([3.5/2*axis_length, 4.85/axis_length, 1.5/2*axis_length, 0.15/axis_length]))
        hc = plt.colorbar(sa, cax=cax, orientation='horizontal')
        if cticks is not None:
            hc.set_ticks(cticks)
        if cticklabels is not None:
            # hc.set_ticks([0,0.5,1])
            hc.ax.tick_params(labelsize=9)
            # cax.invert_xaxis()
            cax.set_xticklabels(cticklabels, fontsize=10)

        hc.locator = ticker.MaxNLocator(integer=True, min_n_ticks=2, nbins=5)
        hc.update_ticks()

        cax.set_ylabel(clabel, rotation=0, ha='right', va='center', fontsize=12)


    if xlim is not None:
        ax1.set_xlim(xlim)
    if ylim is not None:
        ax1.set_ylim(ylim)

    return fig

def plot_pca_ax(
    P_df,
    pca,
    ax=None,
    c=None,
    cohort_s=None,
    cohort_colors=None,
    cohort_args=None,
    order=[1,2,3],
    outliers=None,
    title='',
    vmin=None,
    vmax=None,
    alpha=1,
    lw=0,
    s=30,
    cmap=plt.cm.Spectral_r,
    cticks=None,
    cticklabels=None,
    clabel='',
    show_legend=True,
    xlim=True,
    ylim=True,
    add_mean_marker=False,
    shuffle_rows=True
    ):
    """
    PCA Plot by axis.
    -------------------
    cohort_s: Series encoding cohorts
    cohort_colors: dict
    Modes:
    """
    if cohort_s is not None:
        cohort_s = cohort_s.fillna("n/a")
        cohorts = cohort_s.unique()
        nc = len(cohorts)

        if cohort_colors is None and cohort_args is None:
            cohort_colors = {i:j for i,j in zip(cohorts, sns.husl_palette(nc, s=1, l=0.6))}
            
        cohort_colors["n/a"] = "lightgrey"

        if cohort_args is None:
            cohort_args = {}
            for k in np.unique(cohort_s):
                cohort_args[k] = {'color': cohort_colors[k], 'marker':'o', 'edgecolor':'none', 's':s}
    
    if ax is None:
        fig,ax = plt.subplots(figsize=(6,6))

    if cohort_s is None:
        sa = ax.scatter(P_df[order[1]-1], P_df[order[0]-1], c=c, cmap=cmap, vmin=vmin, vmax=vmax, lw=lw, alpha=alpha, s=s, rasterized=True)
    else:
        if shuffle_rows is False:
            for k in np.unique(cohort_s):
                i = cohort_s[cohort_s==k].index
                ax.scatter(P_df.loc[i,order[1]-1], P_df.loc[i,order[0]-1], alpha=alpha, label=k, rasterized=True, **cohort_args[k])
        else:
            # Always plot missing values first
            if "n/a" in cohort_s.values:
                k = "n/a"
                i = cohort_s[cohort_s==k].index
                ax.scatter(P_df.loc[i,order[1]-1], P_df.loc[i,order[0]-1], alpha=alpha, rasterized=True, **cohort_args[k])
                
            # Randomly shuffle other rows
            for i in cohort_s[cohort_s!="n/a"].sample(frac=1).index:
                k = cohort_s[i]
                ax.scatter(P_df.loc[i,order[1]-1], P_df.loc[i,order[0]-1], alpha=alpha, rasterized=True, **cohort_args[k])

            for k in np.unique(cohort_s):
                if k!="n/a":
                    ax.scatter(0, 0, alpha=0, label=k, rasterized=True, **cohort_args[k])

        if add_mean_marker:
            for k in np.unique(cohort_s):
                if k!="n/a":
                    a,b = P_df.join(cohort_s).groupby(cohort_s.name).mean().loc[k,[order[1]-1,order[0]-1]]
                    ax.plot(a, b,'o', c=cohort_args[k]["color"], markersize=10, markeredgecolor='k')

    format_plot(ax, fontsize=10)
    ax.set_xlabel('PC {0} ({1:.2f}%)'.format(order[1], pca.explained_variance_ratio_[order[1]-1]*100), fontsize=12)
    ax.set_ylabel('PC {0} ({1:.2f}%)'.format(order[0], pca.explained_variance_ratio_[order[0]-1]*100), fontsize=12)

    if outliers is not None:
        ax.scatter(P_df.loc[outliers, order[1]-1], P_df.loc[outliers, order[0]-1], c='none', edgecolors='r', marker='s', lw=1, alpha=1, s=50, label=None, rasterized=True)

    ax.set_title(title, fontsize=12)

    if cohort_s is not None and show_legend:
        leg = ax.legend(loc=0, fontsize=9, scatterpoints=1, handletextpad=0.1, framealpha=1, labelspacing=0.35)
        for lh in leg.legendHandles:
            lh.set_alpha(1)

    if cohort_s is None and c is not None and len(c)==P_df.shape[0]:
        x1 = ax.get_position().x1
        y1 = ax.get_position().y1

        fig = plt.gcf()
        #cax = fig.add_axes(np.array([x1, y1*5/5.5, 0.15/5.5, 1/5.5]))
        cax = fig.add_axes(np.array([4.5/5.5, 5/5.5, .5/5.5, 0.15/5.5]))
        hc = plt.colorbar(sa, cax=cax, orientation='horizontal')

        if cticks is not None:
            hc.set_ticks(cticks)
        if cticklabels is not None:
            hc.ax.tick_params(labelsize=9)
            cax.set_xticklabels(cticklabels, fontsize=10)

        hc.locator = ticker.MaxNLocator(integer=True, min_n_ticks=2, nbins=5)
        hc.update_ticks()

        cax.set_ylabel(clabel, rotation=0, ha='right', va='center', fontsize=12)

    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    return ax

def format_plot(ax, tick_direction='out', tick_length=4, hide=['top', 'right'], hide_spines=True, lw=1, fontsize=9):

    for i in ['left', 'bottom', 'right', 'top']:
        ax.spines[i].set_linewidth(lw)

    # ax.axis["left"].major_ticklabels.set_ha("left")
    ax.tick_params(axis='both', which='both', direction=tick_direction, labelsize=fontsize)

    # set tick positions
    if 'top' in hide and 'bottom' in hide:
        ax.get_xaxis().set_ticks_position('none')
    elif 'top' in hide:
        ax.get_xaxis().set_ticks_position('bottom')
    elif 'bottom' in hide:
        ax.get_xaxis().set_ticks_position('top')
    else:
        ax.get_xaxis().set_ticks_position('both')

    if 'left' in hide and 'right' in hide:
        ax.get_yaxis().set_ticks_position('none')
    elif 'left' in hide:
        ax.get_yaxis().set_ticks_position('right')
    elif 'right' in hide:
        ax.get_yaxis().set_ticks_position('left')
    else:
        ax.get_yaxis().set_ticks_position('both')

    if hide_spines:
        for i in hide:
            ax.spines[i].set_visible(False)

    # adjust tick size
    for line in ax.xaxis.get_ticklines() + ax.yaxis.get_ticklines():
    #for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(tick_length) # tick length
        line.set_markeredgewidth(lw) # tick line width

    for line in (ax.xaxis.get_ticklines(minor=True) + ax.yaxis.get_ticklines(minor=True)):
        line.set_markersize(tick_length/2) # tick length
        line.set_markeredgewidth(lw/2) # tick line width