"""
SCOOTI_trainer.py
===============================================================
executable regressorTraining.py to initiate regression training
"""
# packages
import numpy as np
import pandas as pd
import os, sys
from tqdm import tqdm
from SCOOTI.regressorTraining import *
import argparse



# parser
parser = argparse.ArgumentParser()
parser.add_argument("--unconModel", help="Fluxes of unconstrained models")
parser.add_argument("--conModel", help="Fluxes of constrained models")
parser.add_argument("--savePath", help="Root path to save regression models")

# optional
parser.add_argument(
        "--kappaArr",
        default='10,1,0.1,0.01,0.001',
        help="Parameter kappa for CFR; 0 for non-CFR models"
        )
parser.add_argument(
        "--rhoArr",
        default='10,1,0.1,0.01,0.001',
        help="Parameter rho for CFR; 0 for non-CFR models"
        )
parser.add_argument(
        "--dkappaArr",
        default='10,1,0.1,0.01,0.001',
        help="Parameter kappa for DFA; 0 for non-CFR models"
        )
parser.add_argument("--expName", default='regression', help="Name of the experiment")
parser.add_argument(
        "--unconNorm", choices=['T', 'F'], default='T', help="Normalization"
        )
parser.add_argument(
        "--conNorm", choices=['T', 'F'], default='F', help="Normalization"
        )
parser.add_argument("--medium", default='', help="Medium name, '', KSOM, or DMEMF12")
parser.add_argument("--method", default='cfr', help="Method name; model, CFR, DFA, COMPASS, or INIT")
parser.add_argument("--model", default='recon1', help="Name of the metabolic model")
parser.add_argument(
        "--inputType", default='flux', help="Type of inputs; flux, pageRank, or stack"
        )
parser.add_argument(
        "--clusterPath", default='', help="path to access the .csv file that clusters metabolites"
        )
parser.add_argument(
        "--objListPath", default='', help="path to access the .csv file that select metabolites as objectives"
        )
parser.add_argument(
        "--rank",
        choices=['T', 'F'],
        default='F',
        help="Rank data if set true"
        )

parser.add_argument(
        "--stackModel",
        choices=['T', 'F'],
        default='F',
        help="stacked models used if set true"
        )

parser.add_argument(
        "--sampling",
        choices=['U', 'C', 'F'], # U: sampling unconstrained models, C: sampling constrained models, F: no sampling
        default='F',
        help="Iterate through all different sets of unconstrained/constrained models to infer objectives"
        )

parser.add_argument(
        "--learner",
        choices=['L', 'A', 'Lasso', 'MLP'], # L: linear regression, A: ADALINE single-layered neural network
        default='L',
        help="choose different meta-learner to stack the results of base models"
        )

parser.add_argument(
        "--geneKO",
        choices=['T', 'F'],
        default='F',
        help="Single gene knockout models used in the training or not"
        )

parser.add_argument(
        "--geneListPath",
        default='',
        help="List of unique genes from genome scale metabolic models (.mat)"
        )


parser.add_argument(
        "--learningRate",
        default=0.001,
        help="Learning rate of ADALINE model"
        )
parser.add_argument(
        "--epo",
        default=10000,
        help="Limitation of learning cycles for ADALINE model"
        )


# get arguments
args = parser.parse_args()
if args.kappaArr:
    kappaArr = [float(s) for s in args.kappaArr.split(',')]

if args.rhoArr:
    rhoArr = [float(s) for s in args.rhoArr.split(',')]

if args.dkappaArr:
    dkappaArr = [float(s) for s in args.dkappaArr.split(',')]
# adjust the inputs
args.unconNorm = True if args.unconNorm=='T' else False
args.conNorm = True if args.conNorm=='T' else False
args.rank = True if args.rank=='T' else False
args.stackModel = True if args.stackModel=='T' else False
args.geneKO = True if args.geneKO=='T' else False


# run regression
if args.sampling=='U':
    # samples of unconstrained models
    regression_exe = unconstrained_models_sampling_regressorTraining(
        args.unconModel,
        args.conModel,
        args.savePath,
        kappa_arr=kappaArr,
        rho_arr=rhoArr,
        dkappa_arr=dkappaArr,
        expName=args.expName,
        uncon_norm=args.unconNorm,
        con_norm=args.conNorm,
        medium=args.medium,
        method=args.method,
        model=args.model,
        input_type=args.inputType,
        cluster_path=args.clusterPath,
        objList_path=args.objListPath,
        rank=args.rank,
        stack_model=args.stackModel,
        learner=args.learner,
        geneKO=args.geneKO,
        geneList_path=args.geneListPath,
        learning_rate=args.learningRate,
        epo=args.epo
    )

elif args.sampling=='C':
    # samples of constrained models
    regression_exe = constrained_model_sampling_regressorTraining(
        args.unconModel,
        args.conModel,
        args.savePath,
        kappa_arr=kappaArr,
        rho_arr=rhoArr,
        dkappa_arr=dkappaArr,
        expName=args.expName,
        uncon_norm=args.unconNorm,
        con_norm=args.conNorm,
        medium=args.medium,
        method=args.method,
        model=args.model,
        input_type=args.inputType,
        cluster_path=args.clusterPath,
        objList_path=args.objListPath,
        rank=args.rank,
        stack_model=args.stackModel,
        learner=args.learner,
        geneKO=args.geneKO,
        geneList_path=args.geneListPath,
        learning_rate=args.learningRate,
        epo=args.epo
    )
else:
    regression_exe = regressorTraining(
        args.unconModel,
        args.conModel,
        args.savePath,
        kappa_arr=kappaArr,
        rho_arr=rhoArr,
        dkappa_arr=dkappaArr,
        expName=args.expName,
        uncon_norm=args.unconNorm,
        con_norm=args.conNorm,
        medium=args.medium,
        method=args.method,
        model=args.model,
        input_type=args.inputType,
        cluster_path=args.clusterPath,
        objList_path=args.objListPath,
        rank=args.rank,
        stack_model=args.stackModel,
        learner=args.learner,
        geneKO=args.geneKO,
        geneList_path=args.geneListPath,
        learning_rate=args.learningRate,
        epo=args.epo
    )






