"""

██████╗ ██╗   ██╗██╗     ██╗  ██╗                    
██╔══██╗██║   ██║██║     ██║ ██╔╝                    
██████╔╝██║   ██║██║     █████╔╝                     
██╔══██╗██║   ██║██║     ██╔═██╗                     
██████╔╝╚██████╔╝███████╗██║  ██╗                    
╚═════╝  ╚═════╝ ╚══════╝╚═╝  ╚═╝                    
                                                     
███████╗███╗   ███╗██████╗ ██████╗ ██╗   ██╗ ██████╗ 
██╔════╝████╗ ████║██╔══██╗██╔══██╗╚██╗ ██╔╝██╔═══██╗
█████╗  ██╔████╔██║██████╔╝██████╔╝ ╚████╔╝ ██║   ██║
██╔══╝  ██║╚██╔╝██║██╔══██╗██╔══██╗  ╚██╔╝  ██║   ██║
███████╗██║ ╚═╝ ██║██████╔╝██║  ██║   ██║   ╚██████╔╝
╚══════╝╚═╝     ╚═╝╚═════╝ ╚═╝  ╚═╝   ╚═╝    ╚═════╝ 
                                                     

-----------
Description
-----------
- Infer metabolic objectives in 1C2C/2CBC embryogenesis
- Compare the objectives of transcriptome/proteome/metabolome
- Look for parameter combination with consensus

--------
Overview
--------
- Constrained models without objectives
- Analysis of metabolic objectives



"""
# for ipython
#%load_ext autoreload
#%autoreload 2

#packages
import pandas as pd
import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import cobra
import sys
sys.path.append('./GeneralMethods/')
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

 ██████╗ ██████╗      ██╗        
██╔═══██╗██╔══██╗     ██║        
██║   ██║██████╔╝     ██║        
██║   ██║██╔══██╗██   ██║        
╚██████╔╝██████╔╝╚█████╔╝        
 ╚═════╝ ╚═════╝  ╚════╝         
                                 
 ██████╗ ██████╗ ███████╗███████╗
██╔════╝██╔═══██╗██╔════╝██╔════╝
██║     ██║   ██║█████╗  █████╗  
██║     ██║   ██║██╔══╝  ██╔══╝  
╚██████╗╚██████╔╝███████╗██║     
 ╚═════╝ ╚═════╝ ╚══════╝╚═╝     

---
Coefficients of metabolic objective based on 



"""

# Step1:
# iterate thru all different parameters and models
CFR_parameters = ['obj52_paraScan_norm_nonorm']
DFA_parameters = ['obj52_paraScan_norm_nonorm']
reg_mods = ['sl']
save_root_path = './regressor_results_new/'
# Step2:
# iterations for plots
for reg_mod in reg_mods:
    for CFR_parameter, DFA_parameter in zip(CFR_parameters, DFA_parameters):
        
        # Step2.1: data path
        # get regression models 1C2C
        EmbryoFiles1C2C = {
                        'CFR1C2C':'./regression_models/flux_{1}_CFR1C2C_{0}.csv'.format(CFR_parameter, reg_mod),
                        'DFA1C2C':'./regression_models/flux_{1}_DFA1C2C_{0}.csv'.format(DFA_parameter, reg_mod),
                        }
        # get regression models for 2CBC
        EmbryoFiles2CBC = {
                        'CFR2CBC':'./regression_models/flux_{1}_CFR2CBC_{0}.csv'.format(CFR_parameter, reg_mod),
                        'DFA2CBC':'./regression_models/flux_{1}_DFA2CBC_{0}.csv'.format(DFA_parameter, reg_mod),
                        }

        # calculate 1C2C
        EmbryoFiles = EmbryoFiles1C2C
        # read files and put into a dict
        EMdfs = {k:pd.read_csv(EmbryoFiles[k], index_col=0) for k in EmbryoFiles.keys()}
        # select parameter
        para_sel = pd.Series(EMdfs['DFA1C2C'].columns).str.contains('dk10')
        EMdfs['DFA1C2C'] = EMdfs['DFA1C2C'][EMdfs['DFA1C2C'].columns[para_sel]]
        # replace the column names for 1C2C
        EMdfs['DFA1C2C'].columns = [
                'DFA Zhang, et. al. 1C2C',
                'DFA Zhao, et. al. ESC2CL',
                ]

        # select parameter
        #para_sel = pd.Series(EMdfs['CFR1C2C'].columns).str.contains('k0.1_r1')
        para_sel = EMdfs['CFR1C2C'].columns[pd.Series(EMdfs['CFR1C2C'].columns).apply(lambda x: '_'.join(x.split('_')[-2:]))=='k10_r0.001']
        EMdfs['CFR1C2C'] = EMdfs['CFR1C2C'][para_sel]
        EMdfs['CFR1C2C'].columns = [
                'CFR Isreal, et. al. 1C2C',
                'CFR Gao, et. al. 1C2C',
                'CFR Yu, et. al. 1C2C'
                ]
        # merge and normalize the table
        pca1C2C, labels_1C2C = matrixProcessing(EMdfs, norm='None')

        # calculate 2CBC
        EmbryoFiles = EmbryoFiles2CBC
        # read files and put into a dict
        EMdfs = {k:pd.read_csv(EmbryoFiles[k], index_col=0) for k in EmbryoFiles.keys()}

        # select parameter
        para_sel = pd.Series(EMdfs['DFA2CBC'].columns).str.contains('dk10')
        EMdfs['DFA2CBC'] = EMdfs['DFA2CBC'][EMdfs['DFA2CBC'].columns[para_sel]]
        # replace the column names for 1C2C
        EMdfs['DFA2CBC'].columns = [
                'DFA Zhao, et. al. 2CBC',
                'DFA Sharpley, et. al. 2CBC',
                ]
        # select parameter
        #para_sel = pd.Series(EMdfs['CFR2CBC'].columns).str.contains('k0.1_r1|dk10')
        para_sel = EMdfs['CFR2CBC'].columns[pd.Series(EMdfs['CFR2CBC'].columns).apply(lambda x: '_'.join(x.split('_')[-2:]))=='k10_r0.001']
        EMdfs['CFR2CBC'] = EMdfs['CFR2CBC'][para_sel]
        EMdfs['CFR2CBC'].columns = [
                'CFR Gao, et. al. 2CBC',
                'CFR Zhao, et. al. 2CBC',
                'CFR Isreal, et. al. 2CBC',
                'CFR Sharpley, et. al. 2CBC'
                ]
        # merge and normalize the table
        pca2CBC, labels_2CBC = matrixProcessing(EMdfs, norm='None')

        # remove metabolites did not show in DFA res
        #pca1C2CDFA = pca1C2C[pca1C2C.columns[
        #        pd.Series(pca1C2C.columns).str.contains('DFA')
        #        ]]
        #pca1C2Cind = pca1C2CDFA[pca1C2CDFA.any(axis=1)].index

        ## remove metabolites did not show in DFA res
        #pca2CBCDFA = pca2CBC[pca2CBC.columns[
        #        pd.Series(pca2CBC.columns).str.contains('DFA')
        #        ]]
        #pca2CBCind = pca2CBCDFA[pca2CBCDFA.any(axis=1)].index
        #
        ## rearrange tables
        #pca1C2C = pca1C2C[pca1C2C.index.isin(pca1C2Cind)]
        #pca2CBC = pca2CBC[pca2CBC.index.isin(pca2CBCind)]
        #pca1C2C = pca1C2C[(pca1C2C>0).sum(axis=1)>1]
        #pca2CBC = pca2CBC[(pca2CBC>0).sum(axis=1)>1]
        # Merge
        coefs_df = pd.concat((pca1C2C, pca2CBC), axis=1)
        labels = np.append(labels_1C2C, labels_2CBC)
        
        # make heatmap for 1C2C and 2CBC respectively
        pca_df = pd.concat((pca1C2C, pca2CBC), axis=1).fillna(0)
        rows = pca_df.any(axis=1)


        # normalize the heatmaps with the maximum value of the plots
        from matplotlib.colors import LogNorm, Normalize
        vmax = max(pca1C2C.values.max(), pca2CBC.values.max())
        vmin = min(
                pca1C2C[pca1C2C>0].fillna(1).values.min(),
                pca2CBC[pca2CBC>0].fillna(1).values.min()
                )
        log_norm = LogNorm(vmin=vmin, vmax=vmax)
        
        # merge 1C2C and 2CBC
        pca_merge = pd.concat(
                    (pca1C2C[rows].T,
                    pca2CBC[rows].T),axis=0
                    )
        fig, ax = plt.subplots(1,1,figsize=(12, 6))
        mask = np.zeros_like(pca_merge)
        mask[pca_merge==0] = True
        sns.set(font_scale=2.)
        g= sns.heatmap(
                pca_merge,
                mask=mask,
                cmap='Blues',
                annot=False,
                linewidth=1,
                square=True,
                vmin=0,
                vmax=vmax,
                fmt='.2g',
                cbar_kws={"shrink": 0.5},
                linecolor='k',
                norm=log_norm,
                ax=ax
                )
        g.set_facecolor('xkcd:white')
        plt.tight_layout()
        # save the figures
        CanvasStyle(ax, square=True)
        plt.savefig(f'{save_root_path}/BulkEmbryo_{CFR_parameter+DFA_parameter}_{reg_mod}_merge_heatmap_dfafocus.png')




