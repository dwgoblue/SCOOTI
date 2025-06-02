"""
scEmbryoMetObjAnalysis.py
==============================================================
Metabolic objectives analysis for single-cell mouse embryogenesis
"""

#packages
import pandas as pd
import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import cobra
import warnings; warnings.simplefilter('ignore')
import scipy.stats as ss
import os
from tqdm.notebook import tqdm, trange
from adjustText import adjust_text
import sys, os
from SCOOTI.metObjAnalyzer import metObjAnalyzer

# function to change labels
def label_func(df):
    return pd.Series(df.columns).apply(
        lambda x: x.split('_')[0]
        ).replace({'1C':'sc1C2C', '2C':'sc2CBC'})

# initiate the object
method = 'sl' # superLearner
moa = metObjAnalyzer(
        flux_paths={
        'sc1C2C':'/nfs/turbo/umms-csriram/daweilin/fluxPrediction/sc1C2C/paraScan/',
        'sc2CBC':'/nfs/turbo/umms-csriram/daweilin/fluxPrediction/sc2CBC/paraScan/',
            },
        coef_paths={
            'sc1C2C':f'/nfs/turbo/umms-csriram/daweilin/regression_models/scEmbryo_paraScan/flux_{method}_sc1C2C_input_norm_outcome_nonorm_k0.1_r0.01.csv',
            'sc2CBC':f'/nfs/turbo/umms-csriram/daweilin/regression_models/scEmbryo_paraScan/flux_{method}_sc2CBC_input_norm_outcome_nonorm_k0.1_r0.01.csv',
            },
        save_root_path='/home/daweilin/StemCell/Project_mESC_JinZhang/test_output/',
        GEM_path='/home/daweilin/StemCell/cancer_model.mat',
        uncon_model_path='/nfs/turbo/umms-csriram/daweilin/fluxPrediction/unconstrained_models/pfba/KSOM/',
        col_map={'sc1C2C':'sc1C2C', 'sc2CBC':'sc2CBC'},
        label_func=label_func,
        samplingFlux_path='/nfs/turbo/umms-csriram/daweilin/fluxPrediction/RandomObjCoef/embryo_ksom_50percent/',
        sel_para='k0.1_r0.01',
        prefix='scEmbryo',
        medium='KSOM',
      )

# flux analysis
moa.fluxAnalysis(kappa_arr=[0.1], rho_arr=[0.01])
# load objective coefficients
moa.get_coef()
# analysis of metabolic objectives
moa.coefAnalysis(clustering=True, entropy=True, distance=True, umap_para=[5, 50])

# ++++++++++++++++++++++++++++
# + Metabolic trait analysis +
# ++++++++++++++++++++++++++++
# find the archetype points
def sample_trait_func(pc_coord_df):

    # sample traits
    trait1_ind = np.arange(len(labels))[(pc_coord_df['PC2']==pc_coord_df['PC2'].max())]
    trait2_ind = np.arange(len(labels))[(pc_coord_df['PC1']==pc_coord_df['PC1'].max())]
    trait3_ind = np.arange(len(labels))[(pc_coord_df['PC1']==pc_coord_df['PC1'].min())]

    return np.array([trait1_ind, trait2_ind, trait3_ind])

def control_trait_func(pc_coord_df, arch_tmp_df):

    # mechanistic traits
    trait1_ind = np.arange(len(pc_coord_df))[(pc_coord_df['PC2']==pc_coord_df['PC2'].max())]
    trait2_ind = np.arange(arch_tmp_df.shape[1])[(pc_coord_df['PC1']==pc_coord_df['PC1'].max())]
    trait3_ind = np.arange(arch_tmp_df.shape[1])[(pc_coord_df['PC1']==pc_coord_df['PC1'].min())]
        
    return np.array([trait1_ind, trait2_ind, trait3_ind])

# tradeoff analysis
moa.tradeoff_analysis(
        tri_mets=['gh', 'chsterol', 'glutathione'],
        pareto_mets=['gh'],
        triangle=False,
        pareto3D=False,
        traitAnalysis=True,
        sample_trait_func=sample_trait_func,
        control_trait_func=control_trait_func,
        )



