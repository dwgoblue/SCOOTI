"""
SCOOTI_sampler.py
==============================================================
executable fluxCoefSampler.py to sample fluxes or coefficients
"""
# packages
import numpy as np
import pandas as pd
import os, sys
from tqdm import tqdm
from SCOOTI.pyfluxModel.fluxModel import modelSetter, coefsampler
import argparse



# parser
parser = argparse.ArgumentParser()
parser.add_argument("--datatype", help="Type of output data", choices=['flux', 'coef'])
parser.add_argument("--GEM_path", help="path to access GEM", default="")
parser.add_argument("--objective_path", help="path to access the CSV of objective items")
parser.add_argument("--sample_num", help="Sample number", default=100)
parser.add_argument("--rootpath", help="path to save data", default='./')
parser.add_argument("--medium_path", help="path to access the medium setups", default='')
parser.add_argument("--medium_name", help="Name of the medium you choose", default='DMEMF12', choices=['DMEMF12', 'KSOM'])
parser.add_argument("--suffix", help="Name of experiment", default='recon1')
parser.add_argument("--func", help="Function used to generate coefficients", default='random', choices=['random', 'single'])
parser.add_argument("--coef_path", help="path to access the coefficients of objectives", default='')


# parse arguments
args = parser.parse_args()
if args.datatype=='flux':
    # load flux sampler
    sampler = modelSetter(
        GEM_path=args.GEM_path,
        objective_path=args.objective_path,
        medium_path=args.medium_path,
        medium_name=args.medium_name
    )
    # add demand reactions
    gem_tmp = sampler.build_objective_candidates()
    if len(args.coef_path)>0:
        obj_coef = pd.read_csv(args.coef_path, index_col=0)
        sampler.assign_multi_objectives(obj_coef, gem_tmp, sample_num=args.sample_num,
            rootpath=args.rootpath)
    else:
        sampler.assign_single_objectives(gem_tmp, sample_num=args.sample_num,
            rootpath=args.rootpath)
else: # args.datatype=='coef'
    # load coefficient sampler
    coef = coefSampler(
            single_obj=args.objective_path,
            sample_num=self.sample_num,
            save_path=args.rootpath,
            suffix=args.suffix,
            func=args.func
            )


