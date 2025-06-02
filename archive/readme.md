![ICON](./icon.png)
# SCOOTI
`SCOOTI: Single-Cell Optimization Objective and Tradeoff Inference` is a data-driven computational framework to identify cell-specific objectives and trade-offs from transcriptomics data.

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Preparation of Unconstrained Models](#preparation-of-constrained-models)
- [Omic-based Constrained Models](#omic-based-constrained-models)
- [Inference of Metabolic Objectives](#identify-metabolic-objectives)
- [Metabolic Objective Analysis](#metabolic-objective-analysis)
- [License](#license)

# Overview
`SCOOTI` is able to
1. Formulate metabolic objectives for various cell types, such as cancer cell lines and embryos, in different conditions.
2. Incorporate only transcriptomics data to identify metabolic objectives that is sufficient, but other omics data also works.
3. Extract information of metabolism with cell-specific objectives to analyze the difference/similarity of each cell fate.
4. Analyze the relationships among metabolites and integrate mechanistic trade-offs to identify metabolic traits.

> What is objective function in metabolic modeling?

Biologically, metabolic objective is the metabolic goal a cell wants to achieve. Take bacterial cells for example, scientists working on metabolic modeling usually assumed that growth/doubling is the most important goal for the cells. They built a objective function called `biomass` consisted of nucleotides, amino acides, lipids and so on to represent the growth/doubling of cells which sometimes considered as the phenotypes of the cells.

Mathematically, objective function is a linear combination of several different metabolites. For example,
```
Objective = 3 ATP + 2 NADPH + GSH + Ala
```
where `ATP`, `NADPH`, `GSH`, and `Ala` are metabolites contributing to the metabolic objective. In other words, the phenotypes of cells or the "activities" of cells are determined by the magnitude of cells.

---

Assumes
1. Every cell can optimize their functions to achieve certain metabolic goals depending on different conditions.
2. Biomass, namely cell growth/doubling, is not always the goal for different cells, especially for mammalian cells.
3. Different metabolic distributions and expression profiles can lead to the same metabolic objectives.

# System Requirements

## OS Requirements

The framework has been tested on the following systems:

- Linux: Centos Linux 7 (core) and Ubuntu 22.04

## Python dependencies
```
numpy
scipy
pandas
matplotlib
seaborn
```

## MATLAB dependencies
```
gurobi
cobratoolbox
```

# Installation

To use SCOOTI, first install it using pip:

```
   (venv) $ pip install git+https://github.com/dwgoblue/SCOOTI.git --upgrade
```

> Cobratoolbox is required for SCOOTI to model contrained and unconstrained models. The instruction of cobratoolbox and its installation can be found in the [OpenCobra](https://opencobra.github.io/cobratoolbox/stable/installation.html). In addition, optimization solvers are required for cobratoolbox. We recommend Gurobi solver for the linear version of iMAT and MOOMIN, and CPLEX 20.0.1 for COMPASS. 



# Preparation of Unconstrained Models

Unconstrained models optimized the demand reactions of single metabolites across different compartments. For example, ``atp[c]+atp[n]+atp[m]-->`` is one of the single objective functions. Optimizing this function is equivalent to maximizing the production of ATP. The bash script below shows the example to generate 52 unconstrained models which optimize 52 different single metabolites recorded in the csv file.


```
#!/bin/bash
# input settings
# path to access your matlab-version cobratoolbox
COBRA_path='./cobratoolbox/'
# path to access the metabolic model
GEM_path='./GEMs/Shen2019.mat'
# name of the model
model_name='Recon1'
# leave it blank if no user-defined objectives
obj_candidate_list_file='./obj52_metabolites_recon1.csv'

# path to access the significant genes data
data_dir='./sigGenes/prolif_qui/'
prefix_name='model' # name of the experiment pls set to 'model' for unconstraint models
medium='DMEMF12' # KSOM for embryos and DMEMF12 for cell culture
save_root_path='./fluxPrediction/unconstrained_models/pfba/' # path to save predicted fluxes

# start the simulation of flux predictions
matlab -nosplash -noFigureWindows -r "multiObj_CBM(~, $COBRA_path, $GEM_path, $model_name, $obj_candidate_list_file, $data_dir, $prefix_name, $medium, $save_root_path)"
```

# Omics-based Constrained Models
Constrained models are predicted metabolic fluxes with respect to transcriptomics data. This reflects the control of the metabolic genes. Here, we show the example to generate constrained models for proliferative and quiescent states in a bash file.
```
#!/bin/bash
# input settings
# path to access your matlab-version cobratoolbox
COBRA_path='./cobratoolbox/'
# path to access the metabolic model
GEM_path='./GEMs/Shen2019.mat'
# name of the model
model_name='Recon1'
# leave it blank if no user-defined objectives
obj_candidate_list_file='./obj52_metabolites_recon1.csv'
# objective values
input_obj_tb=''

# parameter settings
DFA_kappa=-1
CFR_kappa=0.1
CFR_rho=10
paraLen=1 # how many kappa/rho used for scanning
random_para=0 # bool, 1 to enable random sampling
init_objective=1 # 1 for none, 2 for biomass objective
genekoflag=0 # bool
rxnkoflag=0 # bool
FVAflag=0 # bool
pfba=1 # 0 for fba and 1 for pfba (minimize sum of fluxes)
medium_perturbation=0 # 1 for depletion or excess of metabolites in medium
pairwise_CFR_model=0
algorithm='iMAT'

# path to access the significant genes data
data_dir='./sigGenes/prolif_qui/'
prefix_name='QuiProlif' # name of the experiment pls set to 'model' for unconstraint models
medium='DMEMF12' # KSOM for embryos and DMEMF12 for cell culture
late_stage='upgenes' # suffix of the file names of significant up-genes
early_stage='dwgenes' # suffix of the file names of significant down-genes
simulation='CFR' # CFR for transcriptomics and proteomics; DFA for metabolomics
constraint=1 # apply constraints to the model
save_root_path='./fluxPrediction/prolif_qui/' # path to save predicted fluxes
CFR_model_path=''
extraWeight=0

# for-loop settings
# if scanning kappa from 1E-3 to 1E1 and rho from 1E-3 to 1E1,
# there will be 25 combinations in total. 
START_NUM=1
END_NUM=$paraLen*$paraLen

# execute parameter scanning
for (( run=$START_NUM; run<=END_NUM; run++ )); do
  matlab -nosplash -noFigureWindows -r "multiObj_CBM($run, $DFA_kappa, $CFR_kappa, $CFR_rho, $COBRA_path,$GEM_path, $model_name, $obj_candidate_list_file, $input_obj_tb, $paraLen, $random_para, $init_objective, $genekoflag, $rxnkoflag, $FVAflag, $pfba, $medium_perturbation, $data_dir, $prefix_name, $medium, $late_stage, $early_stage, $simulation, $constraint, $save_root_path, $CFR_model_path, $pairwise_CFR_model, $extraWeight, $algorithm)"
done
```

# Inference of Metabolic Objectives
After we prepared the constrained models and unconstrained models, we are able to infer metabolic objectives. Here, we show the example to infer metabolic objectives with meta-learner models for proliferative and quiescent states in a bash file.
```
# General settings
unconModel=./fluxPrediction/unconstrained_models/pfba/DMEMF12/
conModel=./fluxPrediction/prolif_qui/
savePath=./regression_models/QuiProlif/
geneListPath=./unique_gene_list.mat
kappa=0.1
rho=10
# run model
python SCOOTI_trainer.py --unconModel $unconModel --conModel $conModel --savePath $savePath --kappaArr $kappa --rhoArr $rho --expName QuiProlif --unconNorm T --conNorm F --medium DMEMF12 --method CFR --model recon1 --inputType flux --rank F --stackModel F --sampling F --learner L --geneKO F --geneListPath $geneListPath
```

# Metabolic Objective Analysis
Once we got the metabolic objectives, we can rely on a Python class `metObjAnalyzer` to perform statistical analysis, phenotype analysis, trade-off analysis, and so on. Here, we show the example to analyze the inferred metabolic objectives for proliferative and quiescent states in a python script.

```
#!/usr/bin/env python
# -*- coding: utf-8 -*-

# initiate a instance with a python object
moa = metObjAnalyzer(
            flux_paths = {
    'QP':'./fluxPrediction/prolif_qui/',
                    },
            coef_paths={
    'QP':f'./regression_models/QuiProlif_paraScan/flux_sl_BulkRNAseq_QuiProlif_input_norm_outcome_nonorm_k0.1_r10.csv',
                },
            save_root_path='./result_files/',
            GEM_path='./metabolicModel/GEMs/Shen2019.mat',
            uncon_model_path='./fluxPrediction/unconstrained_models/pfba/DMEMF12/',
            col_map={},
            label_func=label_func,
            samplingFlux_path='',
            sel_para='k0.1_r10',
            prefix='QuiProlif',
            medium='DMEMF12',
            )
# get coefficients of metabolic objectives
moa.get_coef(metType_cluster=False)
# analyze the coefficients with distance to biomass and statistical comparisons between groups
moa.coefAnalysis(
            norm=False,
            unknown_clustering=False,
            clustering=False,
            entropy=False,
            distance=True,
            compare=True,
            umap_para=[5, 50],
            recluster_min_size=10,
            method='average',
            ref_col='Proliferative',
            cutoff=0.01
            )
```







