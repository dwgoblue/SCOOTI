"""Constrained model analysis for single-cell cell-cycle phases in HeLa

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



# path to regression results
flux_paths = {
        'scCellCycle':'/nfs/turbo/umms-csriram/daweilin/fluxPrediction/scCellCycleHela/Top40Scan/',
        }

# initiate clustering score dictionary
scores = {'MI':[], 'RI':[], 'SI':[], 'kappa':[], 'rho':[]}
for kappa in [0.01,]:
    for rho in [0.01,]:
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
        # Get labels
        labels = pd.Series(flux_df.columns).apply(lambda x: ''.join(x.split('_')[3]))
        cf = clustering_func(
                    flux_df,
                    '/home/daweilin/StemCell/Project_mESC_JinZhang/regressor_results_new/',
                    f'scCellCycle_Top40_flux_no_obj_{kappa}_{rho}',
                    mets_category(),
                )
        
        # clustering
        cf.corr_clustermap(labels)
        #cf.reduction_scatter(labels, continuous=False, func='PCA')
        umaps = cf.reduction_scatter(labels, continuous=False, func='UMAP', para=[5,500])
        recluster = cf.reclustering(umaps, ['UMAP1', 'UMAP2'])
        scoreSI = cf.clustering_evaluation(
                recluster, 'Cell type', 'cluster', method='SI'
                )
        scoreRI = cf.clustering_evaluation(
                recluster, 'Cell type', 'cluster', method='RI'
                )
        scoreMI = cf.clustering_evaluation(
                recluster, 'Cell type', 'cluster', method='MI'
                )
        scores['SI'].append(scoreSI)
        scores['RI'].append(scoreRI)
        scores['MI'].append(scoreMI)
        scores['kappa'].append(kappa)
        scores['rho'].append(rho)

# clustering report
clustering_scores = pd.DataFrame(scores)
clustering_scores.to_csv('/home/daweilin/StemCell/Project_mESC_JinZhang/scCellCycle_clustering_scores.csv')

