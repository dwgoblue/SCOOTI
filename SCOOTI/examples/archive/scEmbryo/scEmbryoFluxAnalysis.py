"""
scEmbryoFluxAnalysis.py
==============================================================
Constrained model analysis for single-cell mouse embryogenesis
- The model predicts the fluxes of ideal optimality
- We aim to model fluxes that can separate cell types
- Clustering performance was applied to the clusters of UMAP of flux predictions
- Best parameter combination is firstly determined by silhouette score
- If the score is all negative, mutual info and rand index were considered
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
from SCOOTI.GeneralMethods.findSigGenes import findSigGenes

def scEmbryoTransitions_sigGenes():
    
    # single cell embryo path
    sc1C2C = '/nfs/turbo/umms-csriram/daweilin/data/scEmbryo/GSE136714/old_sigGenes/1C2C/'
    sc2CBC = '/nfs/turbo/umms-csriram/daweilin/data/scEmbryo/GSE136714/old_sigGenes/2CBC/'
    # function to make tables into binary map
    def read_data(root, fn):
        path = root+'/'+fn
        # get sheets
        exc = pd.ExcelFile(path)
        sn = exc.sheet_names
        # read excel
        udf = pd.read_excel(
                path,
                index_col=0,
                sheet_name=sn[0]
                )
        ddf = pd.read_excel(
                path,
                index_col=0,
                sheet_name=sn[1]
                )
        udf.index = udf.iloc[:,0].to_numpy()
        udf.loc[:, udf.columns[0]] = 1
        ddf.index = ddf.iloc[:,0].to_numpy()
        ddf.loc[:, ddf.columns[0]] = -1
        udf.columns = [fn.split('.xlsx')[0]]
        ddf.columns = [fn.split('.xlsx')[0]]
        df = pd.concat((udf, ddf), axis=0)

        return df


    # read cells
    f1C2C = os.listdir(sc1C2C)
    df1C2C = pd.concat([
        read_data(sc1C2C, f) for f in f1C2C
        ], axis=1).fillna(0)
    f2CBC = os.listdir(sc2CBC)
    df2CBC = pd.concat([
        read_data(sc2CBC, f) for f in f2CBC
        ], axis=1).fillna(0)

    # save data
    df1C2C.to_csv('/nfs/turbo/umms-csriram/daweilin/result_files/embryoProject/gene1C2C.csv')
    df2CBC.to_csv('/nfs/turbo/umms-csriram/daweilin/result_files/embryoProject/gene2CBC.csv')

    
    


    # GEM
    glist = pd.read_csv('/home/daweilin/StemCell/glist.csv')
    cols = dw1.columns#[pd.Series(dw1.columns).isin(glist['Var1'])]
    # boolean gene states
    plot_dw = []
    for dw in [dw1, dw2]:
        plot_dw.append(dw[cols].copy())
    plot_up = []
    for up in [up1, up2]:
        plot_up.append(up[cols].copy())
    ups = pd.concat(plot_up, axis=0)
    dws = pd.concat(plot_dw, axis=0)

    labels = pd.Series(ups.index).apply(lambda x: x.split('_')[0])

    cf = clustering_func(
            ups.T,
            '/home/daweilin/StemCell/Project_mESC_JinZhang/regressor_results_new/',
            f'scEmbryo_upGene',
            mets_category(),
        )

    # show correlations
    cf.corr_clustermap(labels)

    cf = clustering_func(
            dws.T,
            '/home/daweilin/StemCell/Project_mESC_JinZhang/regressor_results_new/',
            f'scEmbryo_dwGene',
            mets_category(),
        )

    # show correlations
    cf.corr_clustermap(labels)

    # +1, 0, -1 for significantly up-regulated, neutral, and down-regulated genes
    merge_df = ups.astype(int)-dws.astype(int)
    cf = clustering_func(
            merge_df.T,
            '/home/daweilin/StemCell/Project_mESC_JinZhang/regressor_results_new/',
            f'scEmbryo_sigGene',
            mets_category(),
        )
    # show correlations
    labels = labels.replace(['2cell', '32cell'], ['sc1C2C', 'sc2CBC'])
    cf.corr_clustermap(labels.to_numpy())
    cf.reduction_scatter(labels, continuous=False, func='UMAP', para=[5,50])


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
        coef_paths={},
        save_root_path='/home/daweilin/StemCell/Project_mESC_JinZhang/test_output/',
        GEM_path='/nfs/turbo/umms-csriram/daweilin/data/models/Shen2019.mat',
        uncon_model_path='/nfs/turbo/umms-csriram/daweilin/fluxPrediction/unconstrained_models/pfba/KSOM/',
        col_map={'sc1C2C':'sc1C2C', 'sc2CBC':'sc2CBC'},
        label_func=label_func,
        samplingFlux_path='/nfs/turbo/umms-csriram/daweilin/fluxPrediction/RandomObjCoef/embryo_ksom_50percent/',
        sel_para='',
        prefix='scEmbryo',
        medium='KSOM',
      )

# flux analysis
moa.fluxAnalysis(kappa_arr=[10, 1, 0.1, 0.01, 0.001], rho_arr=[10, 1, 0.1, 0.01, 0.001])


# get flux
moa.get_flux(
        kappa=kappa,#[10, 1, 0.1, 0.01, 0.001],
        rho=rho,#[10, 1, 0.1, 0.01, 0.001],
        rank=False,
        )



for kappa in [10, 1, 0.1, 0.01, 0.001]:
    for rho in [10, 1, 0.1, 0.01, 0.001]:
        # initiate an analyzer
        moa = metObjAnalyzer(
                    flux_paths = {
    'sc1C2C':'/nfs/turbo/umms-csriram/daweilin/fluxPrediction/sc1C2C/paraScan/',
    'sc2CBC':'/nfs/turbo/umms-csriram/daweilin/fluxPrediction/sc2CBC/paraScan/',
                },
        
                    coef_paths={},
                    save_root_path='/home/daweilin/StemCell/Project_mESC_JinZhang/regressor_results_new/',
                    GEM_path=GEM_path,
                    uncon_model_path=uncon_model_path,
                    col_map={},
                    label_func=label_func,
                    samplingFlux_path='',
                    sel_para=f'k{kappa}_r{rho}',
                    prefix='DiseasePertubations',
                    medium='DMEMF12',
                    )
        # get flux data
        moa.get_flux(
                kappa=kappa,#[10, 1, 0.1, 0.01, 0.001],
                rho=rho,#[10, 1, 0.1, 0.01, 0.001],
                rank=False,
                )
