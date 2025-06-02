
#packages
import pandas as pd
import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import cobra
import sys
sys.path.append('/home/daweilin/StemCell/GeneralMethods/')
from regressionAnalyzer import *
from stat_tests import *
from AnalysisKits import *
from MatplotProp import CanvasStyle, PltProps
PltProps()
import warnings; warnings.simplefilter('ignore')
import scipy.stats as ss
import os
from tqdm.notebook import tqdm, trange
from adjustText import adjust_text
# Regression models
from regressorCollection import *

# Set cobra solver to glpk in order to avoid err msg
config = cobra.Configuration()
config.solver = "glpk"


from scipy.stats import mannwhitneyu



"""

  ██████╗ ██████╗ ███╗   ██╗███████╗████████╗
██╔════╝██╔═══██╗████╗  ██║██╔════╝╚══██╔══╝
██║     ██║   ██║██╔██╗ ██║███████╗   ██║   
██║     ██║   ██║██║╚██╗██║╚════██║   ██║   
╚██████╗╚██████╔╝██║ ╚████║███████║   ██║   
 ╚═════╝ ╚═════╝ ╚═╝  ╚═══╝╚══════╝   ╚═╝   
                                            
███╗   ███╗ ██████╗ ██████╗ ███████╗██╗     
████╗ ████║██╔═══██╗██╔══██╗██╔════╝██║     
██╔████╔██║██║   ██║██║  ██║█████╗  ██║     
██║╚██╔╝██║██║   ██║██║  ██║██╔══╝  ██║     
██║ ╚═╝ ██║╚██████╔╝██████╔╝███████╗███████╗
╚═╝     ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝╚══════╝
                                            
---
Bulk multi-omic analysis
Analysis of constrained models without objectives
 
>> The goal is to find a parameter combination that makes DFA correlate with CFR


"""
# path to regression results
CFR_paths = {
        'CFR1C2C':'./fluxPrediction/1C2C/CFR/paraScan/',
        'CFR2CBC':'./fluxPrediction/2CBC/CFR/paraScan/',
        }

CFR_map = {
    'Israel1C2C_proteomics':'CFR Isreal, et. al. 1C2C',
    'Gao_1C_2C':'CFR Gao, et. al. 1C2C',
    'Yu_1C_2C':'CFR Yu, et. al. 1C2C',
    'Gao_2C_BC':'CFR Gao, et. al. 2CBC',
    'Zhang_2C_BC':'CFR Zhao, et. al. 2CBC',
    'Israel2CBC_proteomics':'CFR Isreal, et. al. 2CBC',
    'Sharpley_2C_BC':'CFR Sharpley, et. al. 2CBC'
    }
# initiate clustering score dictionary
scores = {'MI':[], 'RI':[], 'SI':[], 'kappa':[], 'rho':[]}
CFR_dict = []
for kappa in [10,1,0.1,0.01,0.001]:
    for rho in [10,1,0.1,0.01,0.001]:
        # collect models
        fdfs = {}
        for k in CFR_paths.keys():
            print(k)
            res = load_multiObj_models(
                    CFR_paths[k], medium='KSOM',
                    return_variables=True, norm=False,
                    CFR_paraScan=True, CFR_k=[kappa], CFR_r=[rho],
                    file_suffix='_fluxes.csv.gz'
                    )
            fdfs[k] = res
        # get fluxes of constrained models without objectives
        CFR_df = pd.concat((fdfs[k] for k in fdfs.keys()), axis=1)
        CFR_df = CFR_df[
                CFR_df.columns[CFR_df.any(axis=0)]
                ].replace(
                        [np.inf, -np.inf], [0, 0]
                        )
        CFR_dict.append(CFR_df)

# merge CFR data
CFR_merge = pd.concat(CFR_dict, axis=1)
CFR_merge.columns = pd.Series(CFR_merge.columns).apply(
        lambda x: '_'.join(x.split('_')[:-2]).replace(
            '_'.join(x.split('_')[:-2]),
            CFR_map['_'.join(x.split('_')[:-2])]
            )+'_'+'_'.join(x.split('_')[-2:])
        )
# path to access DFA data
DFA_paths = {
        'DFA1C2C':'./fluxPrediction/1C2C/DFA/paraScan/',
        'DFA2CBC':'./fluxPrediction/2CBC/DFA/paraScan/',
        }
# dictionary for replaceing column names
DFA_map = {
    'DFA1C2C_Jin_Time_Course_Metabolomics':'DFA Zhang, et. al. 1C2C',
    'DFA1C2C_Jin_Metabolomics_DFA':'DFA Zhao, et. al. ESC2CL',
    'DFA2CBC_Jin_Metabolomics_DFA':'DFA Zhao, et. al. 2CBC',
    'DFA2CBC_GSE159484_HPLC_Metabolomics_DFA':'DFA Sharpley, et. al. 2CBC',
    }

# DFA
DFA_dict = []
for kappa in [10,1,0.1,0.01,0.001]:
    # collect models
    fdfs = {}
    for k in DFA_paths.keys():
        print(k)
        res = load_multiObj_models(
                DFA_paths[k], medium='KSOM',
                return_variables=True, norm=False,
                DFA_paraScan=True, DFA_k=[kappa],
                file_suffix='_fluxes.csv.gz'
                )
        if k=='DFA1C2C':
            res2 = load_multiObj_models(
                    DFA_paths[k], medium='DMEMF12',
                    return_variables=True, norm=False,
                    DFA_paraScan=True, DFA_k=[kappa],
                    file_suffix='_fluxes.csv.gz'
                    )
            res = pd.concat((res, res2), axis=1)
        res.columns = pd.Series(res.columns).apply(
                lambda x: k+'_'+x
                )
        fdfs[k] = res
    
    # get fluxes of constrained models without objectives
    DFA_df = pd.concat((fdfs[k] for k in fdfs.keys()), axis=1)
    DFA_df = DFA_df[DFA_df.columns[DFA_df.any(axis=0)]].replace(
            [np.inf, -np.inf], [0, 0]
            )
    DFA_dict.append(DFA_df)

# merge CFR data
DFA_merge = pd.concat(DFA_dict, axis=1)
DFA_merge.columns = pd.Series(DFA_merge.columns).apply(
        lambda x: '_'.join(x.split('_')[:-1]).replace(
            '_'.join(x.split('_')[:-1]),
            DFA_map['_'.join(x.split('_')[:-1])]
            )+'_'+x.split('_')[-1]
        )

# get parameters
CFR_paras = pd.Series(CFR_merge.columns).apply(
        lambda x: '_'.join(x.split('_')[-2:])
        )
DFA_paras = pd.Series(DFA_merge.columns).apply(
        lambda x: x.split('_')[-1]
        )

# iterate thru the parameters
for CFR_para in CFR_paras.unique():
    for DFA_para in DFA_paras.unique():
        # get columns with specific parameters
        CFR_tmp = CFR_merge[
                CFR_merge.columns[CFR_paras==CFR_para]
                ]
        DFA_tmp = DFA_merge[
                DFA_merge.columns[DFA_paras==DFA_para]
                ]
        # merge dfs for correlations
        merge_tmp = pd.concat((DFA_tmp, CFR_tmp), axis=1).corr(
                method='spearman'
                )
        # clustermap
        pca_df = merge_tmp.copy()

        # Create a custom colormap for the heatmap values
        cmap = sns.diverging_palette(220, 20, as_cmap=True)
        # make plots
        sns.set(font_scale=2.5)
        print('Correlation size')
        g = sns.clustermap(
                    pca_df.replace(
                        [np.nan, np.inf, -np.inf], [0,0,0]
                        ),#.reset_index(drop=True),
                    cmap=cmap,
                    center=0,
                    row_cluster=True,
                    col_cluster=True,
                    #yticklabels=False,
                    #xticklabels=False,
                    figsize=(20,20)
                )
        
        # save the figures
        plt.savefig(
                f'./result_figures/embryoFlux_corr_clustermap_batch_{CFR_para}_{DFA_para}.png'
                )
    



# change labels
CFR_merge.columns = pd.Series(CFR_merge.columns).apply(
        lambda x: 'CFR_'+x
        )
DFA_merge.columns = pd.Series(DFA_merge.columns).apply(
        lambda x: 'DFA_'+x
        )
# merge two dataframes
bulk_merge = pd.concat((CFR_merge, DFA_merge), axis=1)
bulk_merge.corr()

# Get labels
DFA_labels = pd.Series(DFA_merge.columns).apply(
        lambda x: 'DFA_'+''.join(x.split('_')[-1])
        )
CFR_labels = pd.Series(CFR_merge.columns).apply(
        lambda x: 'CFR_'+'_'.join(x.split('_')[-2:])
        )
CFRcopy.columns = CFR_labels
DFA_merge.columns = DFA_labels


for i in range(4):
    print(25*i, 25*i+25)
    tmp_labels = CFR_labels.iloc[25*i:25*i+25]
    CFR_merge = CFRcopy.iloc[:, 25*i:25*i+25]

    # merge labels
    labels = pd.concat((DFA_labels, tmp_labels))

    # calculate correlations
    flux_merge = pd.concat((DFA_merge, CFR_merge), axis=1)
    
    # clustermap
    pca_df = flux_merge.copy()
    
    # get sample names without IDs
    cols = labels
    # get colors for type of cells/samples
    cellType_pal = sns.color_palette(
            'cubehelix',#'Set3',
            n_colors=pd.Series(cols).unique().size,
            )
    # map the colors to cell types
    cellType_lut = dict(zip(map(str, pd.Series(cols).unique()), cellType_pal))
    cellType_colors = pd.Series(cols).map(cellType_lut)
    cellType_colors.index = pca_df.columns
    # Create a custom colormap for the heatmap values
    cmap = sns.diverging_palette(h_neg=210, h_pos=350, s=90, l=30, as_cmap=True)
    # make plots
    sns.set(font_scale=2.5)
    print('Correlation size')
    #print(pca_df[pca_df.any(axis=1)].corr().shape)
    #print(pca_df.sum())
    #print(pca_df[pca_df.any(axis=1)].corr())
    g = sns.clustermap(
                pca_df[pca_df.any(axis=1)].corr().replace([np.nan, np.inf, -np.inf], [0,0,0]),#.reset_index(drop=True),
                cmap='viridis',
                #z_score=0,
                row_colors=cellType_colors,
                col_colors=cellType_colors,
                row_cluster=False,
                col_cluster=False,
                yticklabels=False,
                xticklabels=False,
                figsize=(20,20)
            )
    
    # edit legend boxes for the additioanl colors
    from matplotlib.pyplot import gcf
    
    for label in pd.Series(cols).unique():
        g.ax_row_dendrogram.bar(0, 0, color=cellType_lut[label], label=label, linewidth=0)
    
    l2 = g.ax_row_dendrogram.legend(title='Cell Type', loc="upper right", bbox_to_anchor=(1., 1.),
            ncol=3, bbox_transform=gcf().transFigure,frameon=False, fontsize=24)
    
    # save the figures
    plt.savefig(f'./result_figures/embryoFlux_corr_clustermap_batch{i}_2CBC.png')
    
