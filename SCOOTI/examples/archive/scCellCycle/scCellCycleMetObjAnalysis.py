"""Analysis sinlge-cell cell-cycle datasets

This script is to analyze the metabolic objectives of single-cell datasets in G1, S, and G2 phases.
"""

# packages
import pandas as pd
import pickle
import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import cobra
import openpyxl
import sys
sys.path.append('./GeneralMethods/')
from regressionAnalyzer import *
from AnalysisKits import *
from pvalue_for_bars import significance
from MatplotProp import CanvasStyle, PltProps
PltProps()
import warnings; warnings.simplefilter('ignore')
from statsmodels.stats.multitest import fdrcorrection
import scipy.stats as ss
import os
from tqdm.notebook import tqdm, trange
import networkx as nx
import pingouin as pg

# Regression models
from regressorCollection import *

# Set cobra solver to glpk in order to avoid err msg
config = cobra.Configuration()
config.solver = "glpk"



# +++++++++++++++++++
# + Data processing +
# +++++++++++++++++++
# path to regression results
coef_paths = {
        'SuperLearner':'./regression_models/scCellCycle/flux_sl_scscCellCycle_input_norm_outcome_nonorm_Top40_k0.01_r0.01.csv',
        }


# collect models
fdfs = {}
for k in coef_paths.keys():
    print(k)
    res = pd.read_csv(coef_paths[k], index_col=0)
    res.columns = pd.Series(res.columns).apply(lambda x: f'{k}-{x}')
    fdfs[k] = res

# get fluxes of constrained models without objectives
coef_df = pd.concat((fdfs[k] for k in fdfs.keys()), axis=1)

# rename the columns
plot_df = coef_df.copy()
#plot_df = plot_df.div(plot_df.sum(axis=0), axis=1)
plot_df.columns = pd.Series(plot_df.columns).apply(
        lambda x: '_'.join(x.split('_')[2:]).split('_k')[0]
        )

# get pca results and labels
labels = pd.Series(plot_df.columns).apply(lambda x: x.split('_')[1])
cellCycle_map = {'G1':'G1', 'G2':'G2', 'S':'S'}
labels = labels.replace(cellCycle_map)


# +++++++++++++++++++++++
# + Metabolic Tradeoffs +
# +++++++++++++++++++++++
# Trade-off of objective coefficients
tradeoff_df = plot_df.copy().T
tradeoff_df['cellType'] = labels.to_numpy()
sampling_objFlux = read_sampling_objFlux(
        path='./fluxPrediction/RandomObjCoef/cellcycle_dmemf12_50percent/',
        medium='DMEMF12'
        )

# Tradeoff of objective coef
tradeoff_df_norm = tradeoff_df.copy()
tradeoff_df_norm.iloc[:,:-1] = tradeoff_df.iloc[:, :-1].div(
        tradeoff_df.iloc[:, :-1].sum(axis=1), axis=0
        ).fillna(0)

# Relationship of coefficients
ObjectiveTradeOffs(
        tradeoff_df,
        f'scCellCycle',
        'cellType', sampling_objFlux=[], corrplot=True, corr_th=0.25,
        save_root_path='./result_figures/',
        compare_mets=['gthrd'],# 'xolest_hs', 'dcmp'],
        theory_ref=False,
        cellType_curvefit=True,
        pareto_line='connected',
        hue_order=['G1', 'S', 'G2']
        )

# 2D Pareto fronts based on the total ratio
ObjectiveTradeOffs(
        tradeoff_df,
        f'scCellCycle',
        'cellType', sampling_objFlux=sampling_objFlux,
        corrplot=False,
        save_root_path='./result_figures/',
        compare_mets=['gthrd'],# 'gthrd', 'chsterol'],
        theory_ref=True,
        cellType_curvefit=False,
        pareto_line='connected',
        hue_order=['G1', 'S', 'G2']
        )


# distance to the 3d surface
dist3d = pareto_surface_to_points(
        tradeoff_df_norm,
        sampling_objFlux,
        'gthox', 'chsterol', 'gthrd', prefix='scCellCycle',
        )
# pvalues
dist_p = []
for phase in ['G1', 'S', 'G2']:
    dist_p.append(ss.ttest_ind(
           dist3d[dist3d['cellTypes']=='Random']['Distances'],
            dist3d[dist3d['cellTypes']==phase]['Distances']
            )[1])

# Triangle plots among three metabolites
#tradeoff_df_norm['chsterol+dcmp+cys]
triangle_plot(
        tradeoff_df_norm,
        [],
        'dcmp', 'gthox', 'gthrd', 
        'scCellCycle',
        ['G1', 'S', 'G2']
        )

triangle_plot(
        tradeoff_df_norm,
        [],
        'h2o', 'ala-L', 'mag_hs', 'scCellCycle'
        )
# 3D pareto plots
pareto_surface_3d(
        tradeoff_df_norm,
        sampling_objFlux,
        'gthox', 'dcmp', 'chsterol', 'scCellCycle_display',
        rot=[35, 60]# 35, 220
        )
 
# +++++++++++++++++++++++++++++++++++
# + PCA identifies metabolic traits +
# +++++++++++++++++++++++++++++++++++
# PCA of randomly simulated fluxes       
archetype_df = sampling_objFlux.copy().T
archetype_df = archetype_df.iloc[:-1, :]

# Merge ratio of fluxes and ratio of coefficients
sel_df = plot_df.T[archetype_df.index].T
sel_tmp_df = sel_df.div(sel_df.sum(axis=0), axis=1)
arch_tmp_df = archetype_df.div(archetype_df.sum(axis=0), axis=1)
merge_tmp = pd.concat((arch_tmp_df, sel_tmp_df,), axis=1)
labels_tmp = ['Control']*arch_tmp_df.shape[1]+labels.to_list()
# PCA
cf = clustering_func(
        merge_tmp,
            './result_figures/',
            f'scCellCycle_50percent_mets',
            mets_category(),
        )

# clustering
# 3 traits
cf.reduction_scatter(
        labels_tmp, continuous=False, func='PCA', para=[2,2]
        )
# 4 traits
cf.reduction_scatter3D(
        labels_tmp, continuous=False, func='PCA',
        para=[3,3], plot_order=['Control', 'G1', 'S', 'G2'],
        rot=[30, 120], alpha=[0.3, 0.5, 0.5, 0.5]
        )
 

# ++++++++++++++++++++++++++++
# + Metabolic trait analysis +
# ++++++++++++++++++++++++++++
# find the archetype points
def sample_trait_func(pc_coord_df):

    # sample traits
    Strait1_ind = np.arange(len(labels))[(pc_coord_df['PC1']==pc_coord_df['PC1'].min())]
    Strait2_ind = np.arange(len(labels))[(pc_coord_df['PC2']==pc_coord_df['PC2'].max())]
    Strait3_ind = np.arange(len(labels))[(pc_coord_df['PC2']==pc_coord_df['PC2'].min())]
    Strait4_ind = np.arange(len(labels))[(pc_coord_df['PC3']==pc_coord_df['PC3'].max())]

    return np.array([Strait1_ind, Strait2_ind, Strait3_ind, Strait4_ind])

def control_trait_func(pc_coord_df, arch_tmp_df):

    # mechanistic traits
    trait1_ind = np.arange(len(pc_coord_df))[(pc_coord_df['PC1']==pc_coord_df['PC1'].min())]
    trait2_ind = np.arange(arch_tmp_df.shape[1])[(pc_coord_df['PC2']==pc_coord_df['PC2'].max())]
    trait3_ind = np.arange(arch_tmp_df.shape[1])[(pc_coord_df['PC2']==pc_coord_df['PC2'].min())]
    trait4_ind = np.arange(arch_tmp_df.shape[1])[(pc_coord_df['PC3']==pc_coord_df['PC3'].max())]
    
    return np.array([trait1_ind, trait2_ind, trait3_ind, trait4_ind])


# trait analysis
merge_tmp, labels_tmp, arch_df = trait_analysis(
        sampling_objFlux,
        plot_df,
        labels,
        sample_trait_func,
        control_trait_func,
        n_pc=3,
        sample_name='CellCycle',
        plot=True
        )

# Overlay coefficients of each metabolite onto the pca plot
for met in sampling_objFlux.columns[:-1]:
    # visualize with PCA plot
    cf = clustering_func(
            merge_tmp,
                './result_figures/',
                f'scCellCycle_50percent_{met}',#_woH2O',
                mets_category(),
            )
    coef = merge_tmp[merge_tmp.index==met].to_numpy()[0].astype(float)
    min_v = np.min(coef[coef>0])
    met_sum = merge_tmp[merge_tmp.index==met].to_numpy()[0].astype(float)
    #met_sum = np.log10(
    #        merge_tmp[merge_tmp.index==met].to_numpy()[0].astype(float)+min_v/10
    #        )
    # 3 traits
    pca3D = cf.reduction_scatter3D(
            met_sum, continuous=True, func='PCA',
            save_animate=False, projection=True, high_dimension=False,
            para=[3,3], plot_order=['Random', 'G1', 'S', 'G2'],
            rot=[30, 120+90+30], alpha=[0.2, 0.5, 0.5, 0.5]
            )



# ++++++++++++++++++++++++++
# + Entropy and cell types +
# ++++++++++++++++++++++++++
# calculate entropy
Enz = logVariance_allocation(
        plot_df, labels, 'scCellCycle', cellCycle_map, dataType='Objective'
        )

# Overlay the entropy with the UMAP plot
# initiate the objective for clustering and analysis
plot_df = plot_df[Enz.index]
cf = clustering_func(
            plot_df.div(plot_df.sum(axis=0), axis=1),
            './result_figures/',
            f'scCellCycle_Top40_Enz',
            mets_category(),
        )
        
# clustering
cf.reduction_scatter(
        Enz['Entropy of Objective'].to_numpy(),
        continuous=True, func='UMAP', para=[50, 100]
        )

# Cell type separation
# initiate clustering methods
cf = clustering_func(
            plot_df.div(plot_df.sum(axis=0), axis=1),
            './result_figures/',
            f'scCellCycle_Top40_2stages',
            mets_category(),
        )

# clustering and labeling with cell types
#cf.clustermap(labels)
#cf.corr_clustermap(labels)
#cf.reduction_scatter(labels, continuous=False, func='PCA')
cf.reduction_scatter(labels, continuous=False, func='UMAP', para=[50,100])


# comparisons of coefficients
target_pairs = [
        ('G1','S'),
        ('S', 'G2'),
        ('G2', 'G1'),
        ]
cols = ['G1', 'S', 'G2']
lollipop_fig(
    plot_df,
    mets_category(),
    np.array(labels),
    f'scCellCycle_coef_Top40_{col1}vs{col2}vs{col3}',
    cols, cutoff=0.1,
    save_root_path='./result_figures/'
    )
# three sample plot
col1, col2, col3 = 'G1', 'S', 'G2'
boxplot_fig(
    plot_df, mets_category(),
    np.array(labels), col2, col1, col2,
    f'scCellCycle_coef_Top40_{col1}vs{col2}vs{col3}',
    col3=col3, norm=True, fc=1, portion=0.1,
    plottype='stripplot', value_ordering=True,
    xlabel='Normalized coefficient',
    save_root_path='./result_figures/'
        )
# ratio
boxplot_fig(
    plot_df.div(plot_df.sum(axis=0), axis=1), mets_category(),
    np.array(labels), col2, col1, col2,
    f'scCellCycle_coef_Top40_{col1}vs{col2}vs{col3}_ratio',
    col3=col3, norm=False, fc=1, portion=0.1,
    plottype='stripplot', value_ordering=True, xlabel='Ratio',
    save_root_path='./result_figures/'
        )




# ++++++++++++++++++++++++++++++++++++++++++++++
# + Distances to the default biomass objective +
# ++++++++++++++++++++++++++++++++++++++++++++++
# distance to biomass
coef_dist = coef_distance_to_biomassObj(
        plot_df, labels, 'scCellCycleHela_noH2O',
        norm=True, rank=False, func='euclidean',
        save_root_path='./result_figures/',
        histplot=True, boxplot_cols=[], boxplot_order=['S', 'G2', 'G1']
        )

