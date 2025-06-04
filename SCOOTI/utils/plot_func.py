#packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import warnings; warnings.simplefilter('ignore')
import scipy.stats as ss
from tqdm import tqdm
import sys, os
from scipy.io import loadmat
from sklearn.preprocessing import StandardScaler
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import pdist
from sklearn.preprocessing import normalize
from sklearn.decomposition import PCA
import umap, hdbscan, phate
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score, silhouette_score
from sklearn.manifold import TSNE
from scipy.spatial import distance
from scipy.cluster import hierarchy
#from clustering_toolkit import ClusteringToolkit  # Assuming you've saved it as clustering_toolkit.py
from SCOOTI.utils.MatplotProp import CanvasStyle


# Add statistical annotations
from scipy.stats import mannwhitneyu, kruskal
from matplotlib.patches import ConnectionPatch

class ComparePlotToolkit:
    def __init__(self, coef_df, labels=[], expName='', mets_category={}, save_root_path='./'):
        self.save_root_path = save_root_path
        self.input_df = coef_df
        self.labels = labels
        self.mets_category = mets_category
        self.prefix = expName

    def _set_plot_style(self):
        sns.set_context("notebook", font_scale=2.0)
        sns.set_style("whitegrid")

    def stripplot_fig(
            self, 
            col_str, # reference label
            cols, # [col1, col2, col3...]
            value_ordering=True,
            pv=0.05,
            fc=1.0,
            portion=0.0,
            norm=False,
            ylabel='Normalized coefficients',
            y_ub=1.0
            ):
        # retrieve attributes
        input_df, mets_category, labels, prefix = self.input_df, self.mets_category, self.labels, self.prefix
        self._set_plot_style() # set figure style
        dfs = input_df.loc[input_df.any(axis=1)].copy()
        dfs.columns = labels
        # using portion as the threshold to filter out objectives
        print('Remove unused rows...')
        if portion > 0:
            dfs_arr = [dfs.loc[:, labels == col] for col in cols]
            valid_mets = np.concatenate([df.index[(df > 0).mean(axis=1) >= portion].to_numpy() for df in dfs_arr])
            dfs = dfs.loc[dfs.index.isin(np.unique(valid_mets))]
        if norm: # simplex normalization
            dfs = dfs.div(dfs.max(axis=1), axis=0)
        # rearrange
        dfs['Mets_type'] = dfs.index.map(mets_category)
        dfs['Objective'] = dfs.index
        dfs['tmp'] = dfs.filter(regex=col_str).median(axis=1)
        dfs = dfs.sort_values(by=['Mets_type', 'tmp']).drop(columns=['tmp'])
        dfs.index = np.arange(len(dfs))
        
        print('Calculating pvalues and foldchanges...')
        print('Check the size of the df', dfs.shape)
        pvalues, foldChanges, xcolor = np.zeros(len(dfs)), np.zeros(len(dfs)), np.array(['k']*len(dfs), dtype='O')
        for idx, row in tqdm(dfs.iterrows()):
            tmp_fc, tmp_p, c_sel = 1, 1, 'k'
            for col in cols:
                fold = np.mean(dfs.loc[idx, dfs.iloc[:,:-2].columns[labels == col]]) / np.mean(
                    dfs.loc[idx, dfs.iloc[:,:-2].columns[labels != col]])
                _, pval = ss.kruskal(
                    dfs.loc[idx, dfs.iloc[:,:-2].columns[labels == col]],
                    dfs.loc[idx, dfs.iloc[:,:-2].columns[labels != col]])
                if fold > tmp_fc and pval < tmp_p:
                    tmp_fc, tmp_p, c_sel = fold, pval, col
            c_sel = c_sel if tmp_fc > fc and tmp_p < pv else 'k'
            foldChanges[idx], pvalues[idx], xcolor[idx] = tmp_fc, tmp_p, c_sel

        dfs['foldChanges'] = foldChanges
        dfs['pvalues'] = pvalues
        dfs['xcolor'] = xcolor
        print(dfs)
        # remove not significant rows and a fc threshold
        print('Convert df into plot-ready df...')
        dfs = dfs[(dfs['pvalues'] < pv) & ((dfs['foldChanges'] > fc) | (dfs['foldChanges'] < 1 / fc))]
        dfs = dfs.sort_values(by=['xcolor', 'foldChanges']) if value_ordering else dfs.sort_values(by='Mets_type')
        dfs['Objective'] = dfs[['Objective', 'pvalues']].apply(lambda x: f"{x[0]} \n {self._pval_to_stars(x[1])}", axis=1)
        plot_df = pd.melt(dfs, id_vars=['Objective'], value_vars=cols)
        plot_df.columns = ['Objective', 'Stage', 'Coefficient']
        plot_df['Coefficient'] = plot_df['Coefficient'].astype(float)

        # calculating width of the figure
        width = 3*(len(plot_df['Objective'].unique())<=5)+len(plot_df.Objective.unique())
        scaler1, scaler2 = 1+0.05*(len(cols)>3), 3+1*(len(cols)>3)
        # setup
        alpha = 0.2
        palette1 = {k:c for k,c in zip(cols, sns.color_palette("Pastel2")[:len(cols)])}
        palette2 = {k:c for k,c in zip(cols, sns.color_palette("Dark2")[:len(cols)])}
        dodge = .7 - .8 /scaler2
        width = width*scaler1
        HANDLES_ind = len(cols)
        hue_order = cols

        print('Generating figures...')
        # create an empty plot
        fig, ax = plt.subplots(1,1,figsize=(width+2, 8))
        #sns.stripplot(data=plot_df, x='Objective', y='Coefficient', hue='Stage', palette="Pastel2", dodge=True, alpha=0.5, ax=ax)
        #sns.pointplot(data=plot_df, x='Objective', y='Coefficient', hue='Stage', palette="Dark2", dodge=0.5, join=False, markers="X", ci=None, ax=ax)

        # stripplot
        g = sns.stripplot(y='Coefficient', x='Objective', s=10, palette=palette1, data=plot_df, hue='Stage', alpha=alpha, dodge=True, zorder=0, hue_order=hue_order, ax=ax)
        g = sns.pointplot(y='Coefficient', x='Objective', data=plot_df, dodge=dodge, join=False, hue='Stage', palette=palette2, markers="X", scale=1, ci=None, hue_order=hue_order, ax=ax)
        ax.set_ylabel(ylabel)
        ax.set_xlabel('')
        if y_ub!=1:
            ax.set_ylim([-0.001, y_ub])
        handles, labels = ax.get_legend_handles_labels()
        # enlarge dot size
        #legend = ax.legend()
        #for handle in legend.legendHandles:
        #    handle.set_markersize(handle.get_markersize() * 10)

        #for dot in handles:
        #    dot.set_size(dot.get_sizes() * 10)
        l = plt.legend(handles[-HANDLES_ind:], labels[-HANDLES_ind:], bbox_to_anchor=(1.00, 1), loc=2, borderaxespad=0., frameon=False)

        # mark yticklabels with colors by foldchange thresholds
        ticklabels = g.axes.get_xticklabels()
        for i, tick_label in enumerate(ticklabels):
            print(tick_label)
            print(dfs['xcolor'].iloc[i])
            print(len(ticklabels))
            tick_label.set_color(palette1[dfs['xcolor'].iloc[i]])
        # output figures
        CanvasStyle(g, lw=12, ticks_lw=4)
        plt.xticks(rotation=90)
        plt.tight_layout()
        plt.savefig(f"{self.save_root_path}/{prefix}_coef_compare.png")
        #plt.close()

        return dfs, plot_df

    def _pval_to_stars(self, pval):
        if pval < 0.001:
            return '***'
        elif pval < 0.01:
            return '**'
        elif pval < 0.05:
            return '*'
        else:
            return 'ns'

    def plot_allocation(self, ref_col=None, cols=[], norm=True, cutoff=0.0, mode='single', special_labels=[]):

        # retrieve attributes
        input_df, labels, prefix = self.input_df, self.labels, self.prefix
        # set plotting styke
        self._set_plot_style()

        # coefficients
        bp_df = input_df.copy()
        bp_df = bp_df[bp_df.any(axis=1)]

        mets_sort = bp_df[bp_df.columns[labels==ref_col]].mean(axis=1).sort_values(ascending=False).index
        bp_df = bp_df.reindex(mets_sort)
        # remove rows with lower values
        bp_df = bp_df[bp_df.mean(axis=1)>cutoff]
        # set up the width of the figure
        width = len(bp_df.index.unique())/2+3*(len(bp_df.index.unique())<=5)+(len(bp_df.index.unique())/2)*(len(bp_df.index.unique())<=5)
        width = 100 if width>100 else width
        # organize the dataframe
        bp_df.columns = labels
        bp_df['metabolites'] = bp_df.index
        bp_df = pd.melt(bp_df, id_vars=['metabolites'])
        bp_df.columns = ['metabolites', 'class', 'allocation']
        # create a new figure
        fig, ax = plt.subplots(1,1,figsize=(width,6))
        cols = cols if len(cols)>0 else np.unique(labels)
        colors = {k:c for k,c in zip(cols, sns.color_palette("Pastel2")[:len(cols)])}
        # Draw a nested barplot by celltype
        g = sns.barplot(
                data=bp_df, x="metabolites", y="allocation", hue="class",
                palette=colors, alpha=.5, ax=ax, hue_order=cols
        )


        # mark yticklabels with colors by foldchange thresholds
        ticklabels = g.get_xticklabels()
        for i, tick_label in enumerate(ticklabels):
            print(tick_label.get_text())
            extra_c = 'r' if tick_label.get_text() in special_labels else 'k'
            tick_label.set_color(extra_c)
        # output
        plt.xticks(rotation=90)
        CanvasStyle(ax, lw=8, ticks_lw=3)
        plt.tight_layout()
        plt.savefig(f"{self.save_root_path}/{prefix}_{ref_col}_{mode}_allocations.png")

        #plt.close()
        return bp_df



    def lollipop_fig(self, cols, cutoff=0.1, value_ordering=True, special_labels=[]):
        # retrieve attributes
        input_df, mets_category, labels, prefix = self.input_df, self.mets_category, self.labels, self.prefix
        # set plotting style
        self._set_plot_style()

        # prepare dataframe for boxplots
        dfs = input_df[input_df.any(axis=1)].copy() # remove rows with all-zero entries
        dfs.columns = labels # cell types as columns

        # celltype-specific tables
        dfs_arr = [dfs[dfs.columns[labels==col]] for col in cols]

        # remove metabolites selected only a few portion of cells
        ms = np.concatenate([df.index[(df>0).mean(axis=1)>=cutoff].to_numpy() for df in dfs_arr])
        dfs = dfs[dfs.index.isin(np.unique(ms))]

        # binarize coefficients
        dfs_arr = [dfs[dfs.columns[labels==col]] for col in cols]
        dfs_binary = pd.concat([(df>0).mean(axis=1) for df in dfs_arr], axis=1)
        dfs_binary.columns = cols
        
        # organize tables
        dfs_binary['Mets_type'] = dfs_binary.index.map(mets_category)
        dfs_binary['Objective'] = dfs_binary.index
        dfs_binary['tmp'] = dfs_binary[dfs_binary.columns[dfs_binary.columns.str.contains(cols[0])]].median(axis=1)
        dfs_binary = dfs_binary.sort_values(by=['Mets_type','tmp'])
        dfs_binary = dfs_binary.drop(columns=['tmp'])
        dfs_binary.index = np.arange(len(dfs_binary))

        # make the dataframe for plots
        plot_df = pd.melt(dfs_binary, id_vars=['Objective'], value_vars=[ele for ele in dfs_binary.columns[:-2]])
        plot_df.columns = ['Objective', 'Class', 'Proportion']
        plot_df['Mets_type'] = plot_df['Objective'].apply(lambda x: x.split(' ')[0]).map(mets_category)
        width = len(plot_df.Objective.unique())+4*(len(plot_df.Objective.unique())<3)
        width = 100 if width>100 else width
        # for vertical plot
        tmp_df = plot_df.copy()
        tmp_df.index = np.arange(len(tmp_df))
        tmp_max = tmp_df.groupby('Objective').max()
        tmp_min = tmp_df.groupby('Objective').min()
        # magnitude ordering
        if value_ordering==True:
            reind = pd.DataFrame({
                'max':plot_df.groupby('Objective').max()['Proportion'].to_numpy(),
                'var':plot_df[['Objective', 'Proportion']].groupby('Objective').var().to_numpy().flatten(),
                'obj':plot_df[['Objective', 'Proportion']].groupby('Objective').var().index}).sort_values(
                        by=['max', 'var'])['obj'].to_numpy()
            tmp_max = tmp_max.T[reind].T
            tmp_min = tmp_min.T[reind].T
            reind = np.repeat(reind, len(cols))
            plot_df.index = plot_df['Objective']
            plot_df = plot_df.T[reind].T
            extra_c = ['r' if ele else 'k' for ele in tmp_max.index.isin(special_labels)]
        
        # create an empty plot
        print('width:', width)
        fig, ax = plt.subplots(1,1,figsize=(width//2, 6))
        colors = {k:c for k,c in zip(cols, sns.color_palette("Dark2")[:len(cols)])}
        print('Ready to plot...')
        ax.vlines(
                x=tmp_max.index.to_numpy(),
                ymin=tmp_min['Proportion'].to_numpy().astype(float),
                ymax=tmp_max['Proportion'].to_numpy().astype(float),
                linewidth=3,
                color='grey'
                )
        g = sns.scatterplot(
                y='Proportion', x='Objective', hue_order=cols,
                data=plot_df, hue='Class', palette=colors,
                markers="o", s=200, zorder=7, ax=ax
                      )

        handles, labels_ = ax.get_legend_handles_labels()
        plt.legend(frameon=False, borderaxespad=0.)
        ticklabels = ax.get_xticklabels()
        for tick in ticklabels:
            tick.set_color('red' if tick.get_text() in special_labels else 'black')
        ax.set_ylim([-0.05, 1.05])
        ax.set_xlabel('')
        plt.xticks(rotation=90)
        CanvasStyle(ax, lw=12, ticks_lw=4)
        plt.tight_layout()
        plt.savefig(f"{self.save_root_path}/{prefix}_lollipop_feature_compare.png")
        
        return plot_df















        #dfs = input_df.loc[input_df.any(axis=1)].copy()
        #dfs.columns = labels
        #dfs_arr = [dfs.loc[:, labels==col] for col in cols]
        #selected_mets = np.unique(np.concatenate([df.index[(df > 0).mean(axis=1) >= cutoff].to_numpy() for df in dfs_arr]))
        #dfs = dfs.loc[dfs.index.isin(selected_mets)]
        #dfs_arr = [dfs.loc[:, labels==col] for col in cols]
        #dfs_binary = pd.concat([(df > 0).mean(axis=1) for df in dfs_arr], axis=1)
        #dfs_binary.columns = cols
        #dfs_binary['Mets_type'] = dfs_binary.index.map(mets_category)
        #dfs_binary['Objective'] = dfs_binary.index
        #dfs_binary['tmp'] = dfs_binary[cols[0]]
        #dfs_binary = dfs_binary.sort_values(by=['Mets_type', 'tmp']).drop(columns=['tmp'])
        #plot_df = pd.melt(dfs_binary, id_vars=['Objective'], value_vars=cols, var_name='Class', value_name='Proportion')
        #tmp_max = plot_df.groupby('Objective')['Proportion'].max()
        #tmp_min = plot_df.groupby('Objective')['Proportion'].min()
        #n_objectives = plot_df['Objective'].nunique()
        #width = max(6, n_objectives // 2)
        #fig, ax = plt.subplots(figsize=(width, 6))
        #ax.vlines(x=tmp_max.index, ymin=tmp_min.values, ymax=tmp_max.values, color='grey', linewidth=3)
        #sns.scatterplot(data=plot_df, x='Objective', y='Proportion', hue='Class', palette='Dark2', s=200, zorder=5, ax=ax)
        #handles, labels_ = ax.get_legend_handles_labels()
        ##for dot in handles:
        ##    dot.set_sizes(dot.get_sizes() * 10)
        ##legend = ax.legend()
        ##for handle in legend.legendHandles:
        ##    handle.set_markersize(handle.get_markersize() * 10)
        #plt.legend(frameon=False, borderaxespad=0.)
        #ticklabels = ax.get_xticklabels()
        #for tick in ticklabels:
        #    tick.set_color('red' if tick.get_text() in special_labels else 'black')
        #ax.set_ylim([-0.05, 1.05])
        #ax.set_xlabel('')
        #plt.xticks(rotation=90)
        #plt.tight_layout()
        #plt.savefig(f"{self.save_root_path}/{prefix}_lollipop_feature_compare.png")
        ##plt.close()
        
        return plot_df





class ClusteringToolkit:
    """
    Example
    -------
    # --- Step 1: Simulate input data ---
    # Create a fake dataframe (100 metabolites × 50 cells)
    np.random.seed(42)
    fake_data = pd.DataFrame(
        np.random.rand(100, 50),
        index=[f"met_{i}" for i in range(100)],
        columns=[f"cell_{i}" for i in range(50)]
    )
    
    # Simulated cell type labels (2 groups)
    labels = ['TypeA'] * 25 + ['TypeB'] * 25
    
    # Simulated metabolite categories
    mets_category = {f"met_{i}": "Group1" if i < 50 else "Group2" for i in range(100)}
    
    # --- Step 2: Initialize toolkit ---
    toolkit = ClusteringToolkit(
        data=fake_data,
        save_path='./results/',
        prefix='test_',
        mets_category=mets_category
    )
    
    # --- Step 3: Generate correlation heatmap ---
    toolkit.corrmap()
    
    # --- Step 4: Generate clustered heatmap ---
    toolkit.corr_clustermap(labels=labels)
    
    # --- Step 5: Run dimensionality reduction and plot ---
    reduced_df = toolkit.reduction_scatter(
        labels=labels,
        method='UMAP',        # Could also be 'PCA', 'tSNE', 'PHATE'
        title_suffix='demo',
        continuous=False,
        alpha=0.7
    )
    
    # --- Step 6: Perform reclustering using HDBSCAN ---
    reclustered_df = toolkit.reclustering(
        cluster_res=reduced_df,
        cols=['UMAP1', 'UMAP2'],
        min_size=5,
        method='UMAP',
        title_suffix='demo'
    )
    
    # --- Step 7 (Optional): Evaluate clustering quality ---
    score = toolkit.clustering_evaluation(reclustered_df, true_col='label', pred_col='cluster', method='MI')
    print(f"Adjusted Mutual Info Score: {score:.3f}")

    """
    def __init__(self, data: pd.DataFrame, save_path: str, prefix: str, mets_category: dict):
        self.data = data.replace([np.inf, -np.inf, np.nan], 0)
        self.save_path = save_path
        self.prefix = prefix
        self.mets_category = mets_category
        self.cluster_res = pd.DataFrame()

    def _get_palette(self, labels):
        unique_labels = pd.Series(labels).unique()
        if len(unique_labels) > 12:
            colormap = plt.cm.gist_ncar
            norm = plt.Normalize(vmin=0, vmax=len(unique_labels) - 1)
            return [colormap(norm(i)) for i in range(len(unique_labels))], dict(zip(unique_labels, [colormap(norm(i)) for i in range(len(unique_labels))]))
        else:
            colors = sns.color_palette("Set2", len(unique_labels))
            return colors, dict(zip(unique_labels, colors))

    def clustering_evaluation(self, df, true_col, pred_col, method='MI'):
        if method == 'MI':
            return adjusted_mutual_info_score(df[true_col], df[pred_col])
        elif method == 'RI':
            return adjusted_rand_score(df[true_col], df[pred_col])
        else:
            X = df[["UMAP1", "UMAP2"]].values
            labels = df[true_col].values
            return silhouette_score(X, labels, metric="sqeuclidean")


    def extract_cluster_centroids(self, X, method='ward', metric='euclidean', n_clusters=5):
        """
        Perform hierarchical clustering, label samples, and extract cluster centroids.
        Returns:
            labels: (n_samples,) cluster assignments
            centroids: (n_clusters, n_features) mean vectors per cluster
        """
        distance_matrix = pdist(X, metric=metric)
        linkage_matrix = sch.linkage(distance_matrix, method=method)
        if n_clusters==0:
            # Extract linkage distances
            linkage_dists = linkage_matrix[:, 2]   
            # Compute the differences between sorted distances
            dist_diffs = np.diff(linkage_dists)
            # Find the largest gap (can indicate optimal cut)
            cut_idx = np.argmax(dist_diffs)
            optimal_threshold = linkage_dists[cut_idx] # estimate
            
            print("Auto-selected threshold:", optimal_threshold)
        
            labels = sch.fcluster(linkage_matrix, t=optimal_threshold, criterion='distance')
        else:
            # Assign cluster labels
            labels = sch.fcluster(linkage_matrix, t=n_clusters, criterion='maxclust')
    
        # Compute cluster centroids
        centroids = np.array([X[labels == i].mean(axis=0) for i in range(1, len(np.unique(labels)) + 1)])
    
    
        return labels, centroids


    def corrmap(self, method='spearman'):
        corr_map = self.data.corr(method)
        plt.figure(figsize=(20, 20))
        sns.heatmap(corr_map, cmap='viridis')
        plt.savefig(f"{self.save_path}/{self.prefix}_corrmap.png")

    def corr_clustermap(self, labels=[], xlabel=False, show_cluster=True, method='average', return_clusters=[False, 0]):
        
        # get correlation map
        corr_df = self.data.corr().replace([np.nan, np.inf, -np.inf], 0)
        labels = labels if len(labels)>0 else self.data.columns
        # Create colors for labels
        sample_pal = sns.color_palette('Set3', n_colors=pd.Series(labels).unique().size)
        # map the colors to cell types
        sample_lut = dict(zip(map(str, pd.Series(labels).unique()), sample_pal))
        sample_colors = pd.Series(labels).map(sample_lut)
        sample_colors.index = corr_df.columns
        # Create a custom colormap for the heatmap values
        cmap = sns.diverging_palette(h_neg=210, h_pos=350, s=90, l=30, as_cmap=True)
        #colors, lut = self._get_palette(labels)
        #sample_colors = pd.Series(labels).map(lut)
        
        if show_cluster:
            row_linkage = hierarchy.linkage(distance.pdist(corr_df), method=method)
            col_linkage = hierarchy.linkage(distance.pdist(corr_df.T), method=method)
            g = sns.clustermap(corr_df, row_linkage=row_linkage, col_linkage=col_linkage,
                               cmap=cmap, col_colors=sample_colors, row_colors=sample_colors,
                               figsize=(20, 20))
        else:
            g = sns.clustermap(corr_df, cmap=cmap, row_colors=sample_colors, col_colors=sample_colors,
                               figsize=(20, 20))

            # edit legend boxes for the additioanl colors
            from matplotlib.pyplot import gcf
            
            for label in pd.Series(labels).unique():
                g.ax_row_dendrogram.bar(0, 0, color=sample_lut[label], label=label, linewidth=0)
            
            l2 = g.ax_row_dendrogram.legend(title='Class', loc="upper right", bbox_to_anchor=(1., 1.),
                    ncol=4, bbox_transform=gcf().transFigure,frameon=False, fontsize=14)

        plt.savefig(f"{self.save_path}/{self.prefix}_corr_clustermap.png")

        if return_clusters[0]:
            return self.extract_cluster_centroids(corr_df, n_clusters=return_clusters[1])

    def dimension_reduction(self, method='PCA', para=[30, 2, 0.1], seed=8):
        data_T = self.data.T
        if method == 'UMAP':
            reducer = umap.UMAP(n_neighbors=para[0], min_dist=para[2], n_components=para[1], random_state=seed)
        elif method == 'tSNE':
            reducer = TSNE(n_components=para[1], perplexity=30, random_state=seed)
        elif method == 'PHATE':
            reducer = phate.PHATE(k=15, t=100)
        else:
            reducer = PCA(n_components=para[1])
        return reducer.fit_transform(data_T)

    def reduction_scatter(self, labels, method='UMAP', title_suffix='', continuous=False, alpha=0.5, para=[30, 2, 0.1]):

        sns.set_context("notebook", font_scale=2.0)
        sns.set_style("whitegrid")
        pc = self.dimension_reduction(method, para)
        df = pd.DataFrame(pc, columns=[f"{method}1", f"{method}2"])
        df['label'] = labels
        self.cluster_res = df

        colors, lut = self._get_palette(labels)
        palette = dict(zip(np.unique(labels), colors))

        #plt.figure(figsize=(12, 8))
        fig, ax = plt.subplots(1,1,figsize=(12, 8))
        ax = sns.scatterplot(data=df, x=f"{method}1", y=f"{method}2", hue='label', palette=palette, alpha=alpha, ax=ax)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.title(f"{method} projection")
        plt.tight_layout()
        CanvasStyle(ax, lw=12, ticks_lw=4)
        plt.savefig(f"{self.save_path}/{self.prefix}_{method}_{title_suffix}_reduction.png")
        return df

    def reclustering(self, cluster_res: pd.DataFrame, cols: list, min_size=10, method='UMAP', title_suffix=''):

        sns.set_context("notebook", font_scale=2.0)
        sns.set_style("whitegrid")
        clusterer = hdbscan.HDBSCAN(min_cluster_size=min_size)
        cluster_res['cluster'] = clusterer.fit_predict(cluster_res[cols])
        self.cluster_res = cluster_res

        #plt.figure(figsize=(12, 8))
        fig, ax = plt.subplots(1,1,figsize=(12, 8))
        ax = sns.scatterplot(data=cluster_res, x=cols[0], y=cols[1], hue='cluster', palette='Set2', alpha=0.6, ax=ax)
        plt.title("HDBSCAN Reclustering")
        CanvasStyle(ax, lw=12, ticks_lw=4)
        plt.tight_layout()
        plt.savefig(f"{self.save_path}/{self.prefix}_{method}_{title_suffix}_recluster.png")
        return cluster_res



def boxplot_fig(
    input_df, mets_category, labels, col_str, col1, col2,
    prefix, value_ordering=True, pv=0.05, fc=1.0, col3='',
    portion=0.0, norm=False, plottype='boxplot', xlabel='Normalized coefficients',
    save_root_path='./', y_ub=1.0, verbose=False
):
    """
    Visualize metabolite comparisons across samples using boxplot or stripplot.

    Parameters
    ----------
    input_df : pd.DataFrame
        Metabolite coefficients (rows = metabolites, cols = cells).

    mets_category : dict
        Mapping of metabolite name to category.

    labels : array-like
        Sample/cell labels (len = number of columns in input_df).

    col_str : str
        Reference column keyword used for median sorting.

    col1, col2 : str
        Group names to compare.

    prefix : str
        Filename prefix for saving plots.

    value_ordering : bool, optional
        If True, order by effect size (fold-change), default True.

    pv : float
        P-value threshold.

    fc : float
        Fold change threshold.

    col3 : str, optional
        Optional third group for 3-way comparison.

    portion : float
        Minimum cell proportion to retain a metabolite.

    norm : bool
        If True, normalize values to [0, 1].

    plottype : str
        One of 'boxplot' or 'stripplot'.

    xlabel : str
        Label for y-axis.

    save_root_path : str
        Folder to save the plot.

    y_ub : float
        Upper bound for y-axis.

    Returns
    -------
    plot_df : pd.DataFrame
        DataFrame used for plotting.

    means : list
        Assigned colors for fold-change direction.

    tick_labels : list
        Y-tick labels used in the plot.
    """
    labels = pd.Series(labels)
    dfs = input_df.loc[input_df.any(axis=1)].copy()
    dfs.columns = labels

    # Filter metabolites by portion
    if portion > 0:
        dfs_arr = [dfs.loc[:, labels == g] for g in [col1, col2, col3] if g]
        valid = np.concatenate([
            df.index[(df > 0).mean(axis=1) >= portion].to_numpy() for df in dfs_arr
        ])
        dfs = dfs.loc[dfs.index.isin(np.unique(valid))]

    # Normalize if required
    if norm:
        dfs = dfs.div(dfs.max(axis=1), axis=0)

    dfs['Mets_type'] = dfs.index.map(mets_category)
    dfs['Objective'] = dfs.index
    dfs['tmp'] = dfs.filter(regex=col_str).median(axis=1)
    dfs = dfs.sort_values(by=['Mets_type', 'tmp']).drop(columns='tmp')

    # Preallocate results
    pvalues, fold_changes, colors_assigned = [], [], []

    group_list = [col1, col2] + ([col3] if col3 else [])
    color_palette = sns.color_palette("Dark2")[:len(group_list)]

    for idx, row in dfs.iterrows():
        row_vals = row[group_list]

        if col3:
            tmp_fc, tmp_p, selected_color = 1, 1, 'k'
            for i, group in enumerate(group_list):
                vals = dfs[group].loc[idx]
                other_groups = [g for g in group_list if g != group]
                other_vals = dfs[other_groups].loc[idx].values.ravel()
                fold = np.mean(vals) / max(np.mean(other_vals), 1e-6)
                _, pval = ss.mannwhitneyu(vals, other_vals)
                if fold > tmp_fc and pval < tmp_p:
                    tmp_fc, tmp_p, selected_color = fold, pval, color_palette[i]
            fold_changes.append(tmp_fc)
            pvalues.append(tmp_p)
            colors_assigned.append(selected_color if tmp_fc > fc and tmp_p < pv else 'k')
        else:
            vals1, vals2 = dfs[col1].loc[idx], dfs[col2].loc[idx]
            fold = np.mean(vals2) / max(np.mean(vals1), 1e-6)
            _, pval = ss.mannwhitneyu(vals1, vals2)
            pvalues.append(pval)
            fold_changes.append(fold)
            if fold > fc and pval < pv:
                colors_assigned.append(color_palette[1])
            elif fold < 1/fc and pval < pv:
                colors_assigned.append(color_palette[0])
            else:
                colors_assigned.append('k')

    # Add stats to DataFrame
    dfs['pvalues'] = pvalues
    dfs['foldChanges'] = fold_changes
    dfs['means'] = colors_assigned

    # Filter by significance
    sig = (dfs['pvalues'] < pv) & ((dfs['foldChanges'] > fc) | (dfs['foldChanges'] < 1/fc))
    dfs = dfs[sig]

    if value_ordering:
        dfs = dfs.sort_values(by=['means', 'foldChanges'])
    else:
        dfs = dfs.sort_values(by='Mets_type')

    # Convert Objective label with stars
    dfs['Objective'] = dfs.apply(lambda x: f"{x['Objective']}\n{stars(x['pvalues'])}", axis=1)

    # Prepare for plotting
    plot_df = dfs.melt(id_vars='Objective', value_vars=[g for g in group_list])
    plot_df.columns = ['Objective', 'Stage', 'Coefficient']
    plot_df['Coefficient'] = plot_df['Coefficient'].astype(float)

    width = max(8, len(dfs) + 2)
    sns.set_context("notebook", font_scale=2.0)

    fig, ax = plt.subplots(figsize=(width, 8))

    palette_light = dict(zip(group_list, sns.color_palette("Pastel2")[:len(group_list)]))
    palette_dark = dict(zip(group_list, color_palette))
    dodge = 0.5

    if plottype == 'boxplot':
        sns.boxplot(data=plot_df, x='Objective', y='Coefficient', hue='Stage',
                    palette=palette_light, ax=ax, showmeans=False)
    else:
        sns.stripplot(data=plot_df, x='Objective', y='Coefficient', hue='Stage',
                      palette=palette_light, dodge=True, alpha=0.4, ax=ax)
        sns.pointplot(data=plot_df, x='Objective', y='Coefficient', hue='Stage',
                      palette=palette_dark, dodge=dodge, join=False, markers='X', ci=None, ax=ax)

    ax.set_ylabel(xlabel)
    ax.set_xlabel('')
    if y_ub != 1.0:
        ax.set_ylim([-0.001, y_ub])

    # Highlight tick labels by effect
    tick_labels = ax.get_xticklabels() if plottype == 'stripplot' else ax.get_yticklabels()
    for i, tick in enumerate(tick_labels):
        tick.set_color(colors_assigned[i] if i < len(colors_assigned) else 'k')

    # Tweak and save
    plt.xticks(rotation=90)
    plt.tight_layout()
    CanvasStyle(ax, lw=12, ticks_lw=4)
    plt.savefig(f'{save_root_path}/{prefix}_coef_compare.png')

    return plot_df, colors_assigned, tick_labels



def plot_allocation(coef_df, labels, ref_celltype=None, prefix='', norm=True, cutoff=0.0,
                    mode='single', special_labels=[], save_dir='/path/to/save'):
    """
    Plot allocation of metabolite coefficients for different cell types.

    Parameters
    ----------
    coef_df : pd.DataFrame
        Shape (n_metabolites, n_samples), with coefficient values.

    labels : array-like
        List or array of labels for each sample (columns in coef_df).

    ref_celltype : str or None
        Reference cell type to order metabolites. Required for 'overlap' and 'groupbar'.

    prefix : str
        Prefix for saving the plot.

    norm : bool
        If True, normalize each column by its sum.

    cutoff : float
        Show only metabolites with mean allocation above this threshold.

    mode : {'single', 'overlap', 'groupbar'}
        Plot mode. 'single' for one cell type, 'overlap' for stacked,
        'groupbar' for grouped by cell types.

    special_labels : list of str
        Metabolites to color red on x-axis.

    save_dir : str
        Directory to save the figure.
    """
    labels = pd.Series(labels)
    bp_df = coef_df.loc[coef_df.any(axis=1)]

    if norm:
        bp_df = bp_df.div(bp_df.sum(axis=0), axis=1)

    if mode in ['overlap', 'groupbar']:
        if ref_celltype is None:
            raise ValueError("`ref_celltype` must be specified for 'overlap' or 'groupbar' mode")
        ref_cols = bp_df.columns[labels == ref_celltype]
        sort_index = bp_df[ref_cols].mean(axis=1).sort_values(ascending=False).index
    else:  # mode == 'single'
        ref_cols = bp_df.columns[labels == ref_celltype]
        sort_index = bp_df[ref_cols].mean(axis=1).sort_values(ascending=False).index

    bp_df = bp_df.loc[sort_index]
    mean_vals = bp_df.mean(axis=1)
    bp_df = bp_df[mean_vals > cutoff]

    width = max(6, len(bp_df.index) / 2)
    sns.set_context("notebook", font_scale=2.)

    if mode == 'single':
        melted = bp_df[ref_cols].copy()
        melted['Metabolite'] = melted.index
        melted = melted.melt(id_vars='Metabolite', value_name='Allocation')
        fig, ax = plt.subplots(figsize=(width, 6))
        sns.barplot(data=melted, x='Metabolite', y='Allocation', color='k', ax=ax)

    elif mode == 'groupbar':
        plot_df = bp_df.copy()
        plot_df.columns = labels
        plot_df['Metabolite'] = plot_df.index
        plot_df = plot_df.melt(id_vars='Metabolite', var_name='CellType', value_name='Allocation')
        fig, ax = plt.subplots(figsize=(width, 6))
        unique_labels = labels.unique()
        cp = sns.color_palette('Set2', len(unique_labels))
        sns.barplot(data=plot_df, x='Metabolite', y='Allocation', hue='CellType', ax=ax,
                    palette=cp, alpha=0.7, hue_order=unique_labels)

    elif mode == 'overlap':
        fig, ax = plt.subplots(figsize=(width, 6))
        unique_labels = labels.unique()
        colors = sns.color_palette('Set2', len(unique_labels))
        patches = []
        for i, ct in enumerate(unique_labels):
            cols = bp_df.columns[labels == ct]
            plot_df = bp_df[cols].copy()
            plot_df['Metabolite'] = plot_df.index
            melted = plot_df.melt(id_vars='Metabolite', value_name='Allocation')
            sns.barplot(data=melted, x='Metabolite', y='Allocation', color=colors[i], alpha=0.5, ax=ax)
            patches.append(mpatches.Patch(color=colors[i], label=ct))
        ax.legend(handles=patches)

    # Highlight special metabolites
    ticklabels = ax.get_xticklabels()
    for tick in ticklabels:
        if tick.get_text() in special_labels:
            tick.set_color('red')

    ax.set_xlabel('')
    ax.set_ylabel('Allocation')
    plt.xticks(rotation=90)
    CanvasStyle(ax, lw=8, ticks_lw=3)
    plt.tight_layout()

    filename = f"{prefix}_{ref_celltype}_{mode}_allocations.png" if ref_celltype else f"{prefix}_{mode}_allocations.png"
    plt.savefig(f"{save_dir}/{filename}")
    plt.close()

# ---------------------------
# Visualize Latent Space
# ---------------------------

def cluster_with_umap_hdbscan(X, n_neighbors=15, min_dist=0.1, min_cluster_size=5, n_components=2, figsize=(8,6), title='HDBSCAN', prefix='', save_dir='./'):
    """
    Reduce high-dimensional data using UMAP and cluster using HDBSCAN.
    Returns:
        embedding: UMAP-reduced data (n_samples x n_components)
        labels: cluster labels from HDBSCAN (-1 = noise)
    """
    # Step 1: Scale the data
    #X_scaled = StandardScaler().fit_transform(X)

    # Step 2: Reduce to low dimensions
    reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist,
                        n_components=n_components, random_state=42)
    embedding = reducer.fit_transform(X)

    # Step 3: Cluster using HDBSCAN
    clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size)
    labels = clusterer.fit_predict(embedding)


    # Plot
    plt.figure(figsize=figsize)
    if labels is not None:
        labels = np.array(labels)
        for label in np.unique(labels):
            idx = labels == label
            plt.scatter(
                    embedding[idx, 0],
                    embedding[idx, 1],
                    label=str(label),
                    alpha=0.7
                    )
        if len(np.unique(labels))<10:
            plt.legend(title="Class", fontsize=20)
    else:
        plt.scatter(embedding[:, 0], embedding[:, 1], alpha=0.7)
    
    plt.title(title, fontsize=28)
    plt.xlabel("PC1", fontsize=20)
    plt.ylabel("PC2", fontsize=20)
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(f"{save_dir}/{prefix}_pca.png")

    return embedding, labels



def plot_pca(z, labels=None, title="PCA", figsize=(8, 6), save_dir='./', prefix=''):
    """
    Visualize 2D PCA of latent space (e.g., output z from DAAE encoder)
    
    Parameters:
        z (ndarray): Latent matrix (n_samples x latent_dim), usually from model.encoder(x)
        labels (array-like, optional): Categorical labels for coloring (e.g., cell types)
        title (str): Plot title
        figsize (tuple): Size of the plot
    """
    z = z.detach().cpu().numpy() if hasattr(z, "detach") else np.asarray(z)

    # Standardize before PCA
    z_scaled = z#StandardScaler().fit_transform(z)
    pca = PCA(n_components=2)
    pca.fit(z_scaled)
    # Get explained variance
    explained = pca.explained_variance_ratio_ * 100  # convert to percent
    z_pca = pca.transform(z_scaled) 
    # Plot
    plt.figure(figsize=figsize)
    if labels is not None:
        labels = np.array(labels)
        for label in np.unique(labels):
            idx = labels == label
            plt.scatter(
                    z_pca[idx, 0],
                    z_pca[idx, 1],
                    label=str(label),
                    alpha=0.7
                    )
        if len(np.unique(labels))<10:
            plt.legend(title="Class", fontsize=20)
    else:
        plt.scatter(z_pca[:, 0], z_pca[:, 1], alpha=0.7)
    
    plt.title(title, fontsize=28)
    plt.xlabel(f"PC1 ({explained[0]:.1f}% variance)", fontsize=20)
    plt.ylabel(f"PC2 ({explained[1]:.1f}% variance)", fontsize=20)
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(f"{save_dir}/{prefix}_pca.png")


def plot_umap(
        X,
        labels=None,
        archetypes=None,
        title="UMAP Projection",
        random_state=42,
        figsize=(8, 6),
        save_dir='./',
        prefix='',
        nN=15,
        nC=2,
        ):
    """
    Visualize data and optional archetypes using UMAP.
    
    Parameters:
        X : np.ndarray
            Input data of shape (n_samples, n_features)
        labels : np.ndarray or list (optional)
            Cluster or group labels for each sample
        archetypes : np.ndarray (optional)
            Archetypes to overlay (shape: n_archetypes x n_features)
        title : str
            Title for the plot
    """
    reducer = umap.UMAP(random_state=random_state, n_components=nC, n_neighbors=nN)
    X_umap = reducer.fit_transform(X)

    plt.figure(figsize=(8, 6))
    if labels is not None:
        plt.scatter(
                X_umap[:, 0],
                X_umap[:, 1],
                c=labels,
                cmap='tab10',
                s=20,
                alpha=0.7,
                label='Samples'
                )
    else:
        plt.scatter(X_umap[:, 0], X_umap[:, 1], s=20, alpha=0.7, label='Samples')

    if archetypes is not None:
        A_umap = reducer.transform(archetypes)
        plt.scatter(A_umap[:, 0], A_umap[:, 1], color='black', marker='X', s=120, label='Archetypes')

    plt.title(title, fontsize=28)
    plt.xlabel("UMAP1", fontsize=20)
    plt.ylabel("UMAP2", fontsize=20)
    plt.legend()
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(f"{save_dir}/{prefix}_umap.png")



def extract_cluster_centroids(X, method='ward', metric='euclidean', n_clusters=5):
    """
    Perform hierarchical clustering, label samples, and extract cluster centroids.
    Returns:
        labels: (n_samples,) cluster assignments
        centroids: (n_clusters, n_features) mean vectors per cluster
    """
    distance_matrix = pdist(X, metric=metric)
    linkage_matrix = sch.linkage(distance_matrix, method=method)

    # Assign cluster labels
    labels = sch.fcluster(linkage_matrix, t=n_clusters, criterion='maxclust')

    # Compute cluster centroids
    centroids = np.array([X[labels == i].mean(axis=0) for i in range(1, n_clusters + 1)])


    return labels, centroids



# *--------------------*
# SECTION 2
# *--------------------*
def plot_clustermap(
        df, method='euclidean',
        metric='ward',
        cmap='viridis',
        save_dir='./',
        prefix='',
        figsize=[10, 10]
        ):
    """
    Plot a clustermap of a non-negative DataFrame using seaborn.
    """
    #assert (df >= 0).all().all(), "DataFrame must have non-negative values"
    sns.set(style="white")
    clustermap = sns.clustermap(
            df, cmap=cmap, method=metric,
            metric=method, col_cluster=True,
            row_cluster=True,
            xticklabels=True,
            yticklabels=False,
            figsize=(figsize[0], figsize[1])
            )
    plt.title("Clustermap of DataFrame")
    plt.savefig(f"{save_dir}/{prefix}_clustermap.png")
