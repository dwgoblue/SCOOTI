"""Constrained model analysis for bulk transcriptomics of proliferative-quiescent states

Cell-specific constrained models without objective functions
- The model predicts the fluxes of ideal optimality
- We aim to model fluxes that can separate cell types
- Clustering performance was applied to the clusters of UMAP of flux predictions
- Best parameter combination is firstly determined by silhouette score
- If the score is all negative, mutual info and rand index were considered
"""


# for ipython
#%load_ext autoreload
#%autoreload 2

# packages
import pandas as pd
import pickle
import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import cobra
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


def columns_name_process(df):

    """Supportive function
    
    Adjust and fix the column names
    """
    
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
# + Data Processing +
# +++++++++++++++++++
# path to regression results
flux_paths = {
        'QP':'./fluxPrediction/prolif_qui/'
        }
# initiate clustering score dictionary
scores = {'Coverage':[], 'MI':[], 'RI':[], 'SI':[], 'kappa':[], 'rho':[]}
for kappa in [10, 1, 0.1, 0.01, 0.001]:
    for rho in [10, 1, 0.1, 0.01, 0.001]:
        # collect models
        fdfs = {}
        for k in flux_paths.keys():
            print(k)
            res = load_multiObj_models(
                    flux_paths[k], medium='DMEMF12',
                    return_variables=True, norm=False,
                    CFR_paraScan=True, CFR_k=[kappa], CFR_r=[rho],
                    file_suffix='_fluxes.csv.gz'
                    )
            fdfs[k] = res
        
        # get fluxes of constrained models without objectives
        flux_df = pd.concat((fdfs[k] for k in fdfs.keys()), axis=1)
        flux_df = flux_df[flux_df.columns[flux_df.any(axis=0)]].replace(
                [np.inf, -np.inf], [0, 0]
                )
        dataset_cols = flux_df.columns#[pd.Series(flux_df.columns).str.contains('min_18')]
        flux_df = flux_df[dataset_cols]
        # remove parameter suffix
        flux_df.columns = pd.Series(flux_df.columns).apply(
                lambda x: x.split('_k')[0]
                )
        # Get labels
        labels, tmp = columns_name_process(flux_df) 
        #labels = pd.Series(flux_df.columns).apply(lambda x: x.split('_k')[0])
        cf = clustering_func(
                    flux_df,
                    './result_figures/',
                    f'QP_flux_no_obj_{kappa}_{rho}',
                    mets_category(),
                )
        
        # clustering
        cf.corr_clustermap(labels)
        #cf.reduction_scatter(labels, continuous=False, func='PCA')
        # clustering with umap
        try:
            umaps = cf.reduction_scatter(labels, continuous=False, func='UMAP', para=[3,10])
            # cluster the umap results
            recluster = cf.reclustering(umaps, ['UMAP1', 'UMAP2'], min_size=3)
            # remove noises
            clustered_res = recluster[recluster['cluster']>=0]
            if len(clustered_res):
                # evaluation 
                scoreSI = cf.clustering_evaluation(
                        clustered_res, 'Cell type', 'cluster', method='SI'
                        )
                scoreRI = cf.clustering_evaluation(
                        clustered_res, 'Cell type', 'cluster', method='RI'
                        )
                scoreMI = cf.clustering_evaluation(
                        clustered_res, 'Cell type', 'cluster', method='MI'
                        )
                scores['Coverage'].append(len(clustered_res)/len(recluster))
                scores['SI'].append(scoreSI)
                scores['RI'].append(scoreRI)
                scores['MI'].append(scoreMI)
                scores['kappa'].append(kappa)
                scores['rho'].append(rho)
        except:
            print('Error', kappa, rho)


# +++++++++++++++++
# + Output report +
# +++++++++++++++++
# clustering report
clustering_scores = pd.DataFrame(scores)
clustering_scores.to_csv('./QP_clustering_scores.csv')
