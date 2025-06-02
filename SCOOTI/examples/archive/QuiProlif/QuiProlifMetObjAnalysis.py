"""Constrained model analysis for bulk transcriptomics of Quiescent-Proliferative states

This script is to analyze the metabolic objectives
- what metabolites are significant
- how many proportion of cells/conditions are metabolites used
- allocation of metabolites
- distance from the biomass objective to inferred metabolic objectives
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
from stat_tests import *
from pvalue_for_bars import significance
from MatplotProp import CanvasStyle, PltProps
PltProps()
import warnings; warnings.simplefilter('ignore')
from statsmodels.stats.multitest import fdrcorrection
import statsmodels.api as sm
import scipy.stats as ss
import os
from tqdm.notebook import tqdm, trange
import networkx as nx
import pingouin as pg
from sklearn.metrics import matthews_corrcoef
# Regression models
from regressorCollection import *

# Set cobra solver to glpk in order to avoid err msg
config = cobra.Configuration()
config.solver = "glpk"


"""Supportive function

Adjust and fix the column names
"""

def columns_name_process(df):
    
    labels = []
    df_labels = [] 
    for name in df.columns:
        if 'Sharma_21' in name:
            if '_P_' in name:
                labels.append('Quiescent')
                df_labels.append('_'.join(name.split('_P_')))
            else:
                labels.append('Proliferative')
                df_labels.append('_'.join(name.split('_')[:2])+'_P_'+'_'.join(name.split('_')[2:]))
        else:
            if '_P_' in name:
                labels.append('Proliferative')
            else:
                labels.append('Quiescent')
            df_labels.append('_'.join(name.split('_')[:2])+'_'+'_'.join(name.split('_')[3:]))
    
    return labels, df_labels





# +++++++++++++++++++
# + Data processing +
# +++++++++++++++++++
## path to regression results modeled with recon2.2
#coef_paths = {
#        'SuperLearner':'./regression_models/QuiProlif/flux_sl_BulkRNAseq_QuiProlif_recon2.2_recon2_2_input_norm_outcome_nonorm_k0.1_r0.01.csv',
#        }
#
#para = 'k0.01_r0.01'
#
#
## path to regression results modeled with recon3d
#coef_paths = {
#        'SuperLearner':'./regression_models/QuiProlif/flux_sl_BulkRNAseq_QuiProlif_recon3D_recon3d_input_norm_outcome_nonorm_paraScan.csv',
#        }
#para = 'k0.1_r0.01'


# path to regression results modeled with recon1
coef_paths = {
        'SuperLearner':'./regression_models/QuiProlif/flux_sl_BulkRNAseq_QuiProlif_input_norm_outcome_nonorm_k0.1_r10.csv',
        }
para = 'k0.1_r10'

# selected model
mod = 'SuperLearner'
# collect models
fdfs = {}
for k in coef_paths.keys():
    print(k)
    res = pd.read_csv(coef_paths[k], index_col=0)
    res.columns = pd.Series(res.columns).apply(lambda x: f'{k}-{x}')
    fdfs[k] = res

# get fluxes of constrained models without objectives
coef_df = pd.concat((fdfs[k] for k in fdfs.keys()), axis=1)
paras = pd.Series(coef_df.columns).apply(
        lambda x: '_'.join(x.split('_')[-2:])
        )
# rename the columns
plot_df = coef_df[coef_df.columns[paras==para]].copy()
plot_df.columns = pd.Series(plot_df.columns).apply(
        lambda x: '_'.join(x.split('-')[1:]).split('_k')[0]
        )
dataset_cols = plot_df.columns#[pd.Series(plot_df.columns).str.contains('min_18')]
plot_df = plot_df[dataset_cols]
plot_df = plot_df[plot_df.columns[plot_df.any(axis=0)]]

# Get labels
labels, tmp = columns_name_process(plot_df) 
# replace names
qp_map = {
        'Proliferative':'Proliferative',
        'Quiescent':'Quiescent'
        }
# initiate clustering methods
cf = clustering_func(
            plot_df,
            './result_figures/',
            f'QuiProlif_norm_nonorm_{para}_{mod}',
            mets_category(),
        )

# clustering and labeling with cell types
cf.corr_clustermap(labels)
#cf.reduction_scatter(labels, continuous=False, func='PCA')
cf.reduction_scatter(labels, continuous=False, func='UMAP', para=[5,10])


# +++++++++++++++++++
# + Data processing +
# +++++++++++++++++++
# comparisons of coefficients
col1, col2 = 'Proliferative','Quiescent'

# proportion plot
portion_df = lollipop_fig(
    plot_df, mets_category(), np.array(labels),
    f'QuiProlif_coef_{col1}vs{col2}', [col1, col2],
    cutoff=0.1,
    save_root_path='./result_figures/'
    )
# unadjusted coefficients
boxplot_fig(
        plot_df, mets_category(),
        np.array(labels), col2, col1, col2,
        f'QuiProlif_coef_{col1}vs{col2}',
        fc=1, portion=0., norm=True, plottype='stripplot',
        value_ordering=True,
        save_root_path='./result_figures/'
        )


# +++++++++++++++++++++++++++++++++++
# + Distances to biomass objectives +
# +++++++++++++++++++++++++++++++++++
# distance to biomass
biomass_model_coef = getRxnCoefficients(
        model_path='./models/Shen2019.mat',
        model_name='Recon1',
        )
# Get labels
labels, tmp = columns_name_process(plot_df) 
coef_dist = coef_distance_to_biomassObj(
        plot_df, labels, biomass_model_coef, 'QuiProlif_ratio',
        norm=True, rank=False, func='euclidean', 
        save_root_path='./result_figures/',
        histplot=True, boxplot_cols=[], boxplot_order=['Quiescent', 'Proliferative']
        )


"""Coefficients of metabolic objectives

Visualize the objectives by each dataset
"""


# path to regression results
coef_paths = {
        'SuperLearner':'./regression_models/QuiProlif/flux_sl_BulkRNAseq_QuiProlif_input_norm_outcome_nonorm_k0.1_r10.csv',
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

# Separate dataframes by parameters and models
paras = pd.Series(coef_df.columns).apply(lambda x: '_'.join(x.split('_')[-2:]))
mods = pd.Series(coef_df.columns).apply(lambda x: x.split('-')[0])


#
# collect data modeled by different parameters
para_dict = {}
for para in ['k0.1_r10']:#paras.unique():
    for mod in coef_paths.keys():
        # rename the columns
        plot_df = coef_df.copy()
        plot_df = plot_df[plot_df.columns[paras==para] & plot_df.columns[mods==mod]]
        plot_df.columns = pd.Series(plot_df.columns).apply(
                lambda x: '_'.join(x.split('-')[1:]).split('_k')[0]
                )
        plot_df = plot_df[plot_df.columns[plot_df.any(axis=0)]]
        # get pca results and labels
        data_sources = pd.Series(plot_df.columns).apply(
                lambda x: '_'.join(x.split('_')[:2])
                )
        # get labels
        cellTypes, labels = columns_name_process(plot_df) 
        plot_df.columns = labels
        plot_tmp = plot_df.T
        plot_tmp['cellTypes'] = cellTypes
        plot_tmp['dataSources'] = data_sources.to_numpy()
        plot_tmp = plot_tmp.sort_values(by=['cellTypes', 'dataSources'])
        #labels = plot_tmp['cellTypes']
        # drop additional columns
        plot_tmp = plot_tmp.drop(columns=['cellTypes', 'dataSources'])
        plot_tmp = plot_tmp.T

        for source in data_sources.unique():
            # make plots
            pca_df = plot_tmp[
                    plot_tmp.columns[pd.Series(plot_tmp.columns).str.contains(source)]
                    ]
            sns.set(font_scale=2)
            g = sns.clustermap(
                        pca_df[pca_df.any(axis=1)].rank().T,
                        cmap='Blues',
                        z_score=0,
                        row_cluster=False,
                        col_cluster=True,
                        yticklabels=True,
                        figsize=(len(pca_df[pca_df.any(axis=1)].index)+2,
                            len(pca_df[pca_df.any(axis=1)].columns)+4)
                    )

            # save the figures
            plt.savefig(
                    f'./result_figures/QuiProlif_{para}_{mod}_{source}_clustermap.png'
                    )



