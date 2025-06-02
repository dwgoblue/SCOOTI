"""
fluxSampler.py
=======================================================
Analysis of metabolic objectives and fluxes in diseases
"""


import os
import json
import time
import warnings
from datetime import datetime

import pandas as pd
import numpy as np
import cobra
from cobra.sampling import sample
import cobra.flux_analysis.parsimonious as pFBA

warnings.simplefilter('ignore')
cobra.Configuration().solver = "glpk"


class FluxSampler:
    """
    A class to perform analysis of metabolic objectives and fluxes using a COBRA model.

    Attributes:
        GEM_path (str): Path to the GEM MATLAB file.
        objective_path (str): Path to the CSV file listing objective metabolites.
        medium_path (str): Path to the Excel file defining the culture medium.
        gem (cobra.Model): Loaded GEM model.
        objectives (pd.DataFrame): Objective metabolites.
        objective_candidates (list): List of reactions used as metabolic objectives.
    """
    def __init__(self, GEM_path, objective_path, medium_path, medium_name):
        self.GEM_path = GEM_path
        self.objective_path = objective_path
        self.medium_path = medium_path
        self.medium_name = medium_name
        self.gem = self.load_model()
        self.objectives = self.load_objectives()
        self.objective_candidates = []

    def load_model(self):
        """
        Load and configure the GEM model by applying the culture medium constraints.

        Returns:
            cobra.Model: Configured GEM model with updated medium.
        """
        gem = cobra.io.load_matlab_model(self.GEM_path)
        media = pd.read_excel(self.medium_path, sheet_name=self.medium_name)
        print('Loading GEM and loading the enviromental setup...')
        with gem:
            medium = gem.medium
            ex_mets = [k for k in medium.keys()]
            for EX_rxn, new_lb in zip(media.iloc[:, 5], media.iloc[:, 2]):
                if EX_rxn in ex_mets:
                    medium[EX_rxn] = new_lb
                    gem.medium = medium
                    print('Add:', EX_rxn)
                #except Exception as e:
                else:
                    print('Skipping')

        return gem

    def load_objectives(self):
        """
        Load the list of objective metabolites from a CSV file.

        Returns:
            pd.DataFrame: Filtered objectives table.
        """
        return pd.read_csv(self.objective_path, index_col=0).iloc[1:]

    def build_objective_candidates(self):
        """
        Build demand reactions for all candidate metabolites and append them to the GEM.

        Returns:
            cobra.Model: Modified GEM with added demand reactions.
        """
        compartments = ['c', 'm', 'n', 'x', 'r', 'g', 'l']
        gem_tmp = self.gem.copy()

        print('Adding demand reactions of single objectives...')
        for obj in self.objectives['metabolites']:
            if obj == 'gh':
                self.objective_candidates.append('biomass_objective')
            else:
                for c in compartments:
                    met_id = f'{obj}[{c}]'
                    try:
                        print('Add demand reaction:', met_id)
                        gem_tmp.metabolites.get_by_id(met_id)
                        reaction = cobra.Reaction(f'{obj}_demand')
                        reaction.name = f'Objective candidate {obj}'
                        reaction.add_metabolites({met_id: -1})
                        gem_tmp.add_reactions([reaction])
                        self.objective_candidates.append(f'{obj}_demand')
                        break
                    except:
                        continue
        return gem_tmp


    def assign_single_objectives(self, gem_tmp, sample_num, rootpath):
        """
        Sample flux distributions from the GEM for each individual objective.

        Args:
            gem_tmp (cobra.Model): The GEM with added demand reactions.
            sample_num (int): Number of flux samples to draw.
            rootpath (str): Output directory path.
        """
        for i, candidate in enumerate(self.objective_candidates):
            gem_tmp2 = gem_tmp.copy()
            pFBA.add_pfba(gem_tmp2, candidate)

            obj_c = np.zeros(len(self.objectives['metabolites']))
            obj_c[i] = 1

            samples = pd.DataFrame(sample(gem_tmp2, sample_num, processes=20))
            samples = samples.sample(frac=1)
            samples['Obj'] = samples[candidate].to_numpy()
            samples.index = np.arange(len(samples))
            samples = samples.T

            for col in samples.columns:
                folder = os.path.join(rootpath, f'fs_{col}')
                os.makedirs(folder, exist_ok=True)
                out_name = f'model_ct1_obj{i+1}_data1'
                excelname = self.save_metadata(
                    candidate,
                    self.objectives['metabolites'].to_numpy(),
                    obj_c,
                    folder,
                    self.GEM_path,
                    out_name
                )
                df = pd.DataFrame(samples[col], columns=['upgene'])
                df.to_csv(f'{excelname}_fluxes.csv.gz', compression='gzip')


    def assign_multi_objectives(self, obj_coef, gem_tmp, sample_num, rootpath):
        """
        Sample flux distributions based on linear combinations of multiple objectives.

        Args:
            obj_coef (pd.DataFrame): Objective coefficient matrix.
            gem_tmp (cobra.Model): The GEM with added demand reactions.
            sample_num (int): Number of flux samples to draw.
            rootpath (str): Output directory path.
        """
        obj_df = obj_coef.copy()
        obj_df.index = obj_df.index.to_series().apply(lambda x: f'{x}_demand' if x != 'gh' else 'biomass_objective')

        for col in obj_df.columns:
            gem_tmp2 = gem_tmp.copy()
            obj_dict = {gem_tmp2.reactions.get_by_id(k): v for k, v in obj_df[col].items()}
            pFBA.add_pfba(gem_tmp2, obj_dict)

            obj_c = self.objectives['metabolites'].apply(
                lambda x: obj_df[col].get(x, 0)
            ).to_numpy()

            samples = pd.DataFrame(sample(gem_tmp2, sample_num, processes=20))
            samples = samples.sample(frac=1)
            samples['Obj'] = samples[obj_df[col][obj_df[col] > 0].index].sum(axis=1).to_numpy()
            samples.index = np.arange(len(samples))
            samples = samples.T

            for s in samples.columns:
                folder = os.path.join(rootpath, f'fs_{col}')
                os.makedirs(folder, exist_ok=True)
                out_name = f'model_ct1_obj{s}_data1'
                excelname = self.save_metadata(
                    obj_df[col][obj_df[col] > 0].index[0],
                    self.objectives['metabolites'].to_numpy(),
                    obj_c,
                    folder,
                    self.GEM_path,
                    out_name
                )
                df = pd.DataFrame(samples[s], columns=['upgene'])
                df.to_csv(f'{excelname}_fluxes.csv.gz', compression='gzip')


    def save_metadata(
            self,
            candidate,
            objectives,
            obj_c,
            root_path,
            model_path, 
            out_name='',
            data_path='',
            sample_name='',
            upsheet='',
            dwsheet='',
            ctrl=0,
            kappa=1,
            rho=1,
            medium='DMEMF12',
            genekoflag=False,
            rxnkoflag=False,
            media_perturbation=False
            ):
        """
        Save sampling configuration and metadata as a JSON file.

        Returns:
            str: Base filename prefix (excluding file extension).
        """
        # initiate a metadata dictionary
        metadata = {
            'obj': list(objectives),
            'obj_type': 'demand',
            'obj_c': list(obj_c),
            'output_path': root_path,
            'input_path': data_path if ctrl else sample_name,
            'file_name': out_name,
            'with_constraint': ctrl,
            'CFR_kappa': kappa,
            'CFR_rho': rho,
            'medium': self.medium_name,
            'genekoflag': genekoflag,
            'rxnkoflag': rxnkoflag,
            'media_perturbation': media_perturbation,
            'objWeights': 1,
            'objRxns': candidate,
            'model_path': model_path,
            'upStage': upsheet,
            'dwStage': dwsheet
        }

        file_prefix = time.strftime("%b%d%Y%H%M%S")
        filename = os.path.join(root_path, f'[{file_prefix}]{out_name}')
        with open(f'{filename}_metadata.json', 'w') as f:
            json.dump(metadata, f)
        print('Metadata saved at:', filename)
        return filename

class coefsampler:
        """sample coefficients or single-objective coefficients

        the function will output and save the table of coefficients for 
        a list of metabolites (objectives).

        attributes
        ----------
        single_obj : list,
            a list of metabolites. the name is better to match the name used in modeling.
        sample_num : int,
            the number of samples if running coefficient sampling.
        save_path : str,
            path to save the table
        suffix : str, default='general'
            name of the experiment
        func : str, default='random'
            choose a function to generate coefficient;
            "random" for random_objective_coefficients;
            otherwise, single_objective_coefficients

        returns
        -------
        df : pandas.dataframe,
            coefficients table with metabolites as index and samples as columns

        """
        def __init__(self, single_obj, sample_num, save_path, suffix, func):
            if type(single_obj)==str:
                self.single_obj = pd.read_csv(single_obj, index_col=0).iloc[1:].values.flatten()

            self.sample_num = sample_num
            self.save_path = save_path
            self.suffix = suffix

            coef = self.random_objective_coefficients() if func=='random' else self.single_objective_coefficients()
            self.coef = coef

        def random_objective_coefficients(self):
            # sampling and save
            df = pd.dataframe(np.random.rand(len(self.single_obj), self.sample_num))
            df.index = single_obj
            df.columns = [f'sample_{ind}' for ind in np.arange(len(df.columns))]
            df.to_csv(self.save_path+'/'+'samplingobjcoef_{self.suffix}.csv')
            return df

        def single_objective_coefficients(self):
            # sampling and save
            df = pd.dataframe(
                    np.eye(len(self.single_obj)),
                    columns=self.single_obj,
                    index=self.single_obj)
            df.to_csv(self.save_path+'/'+f'singleobj_{self.suffix}.csv')
            return df


if __name__ == "__main__":
    
    # load flux sampler
    sampler = fluxsampler(
        gem_path="./scooti/scooti/metabolicmodel/gems/shen2019.mat",
        objective_path="./scooti/scooti/metabolicmodel/obj52_metabolites_shen2019.csv",
        medium_path="./scooti/scooti/metabolicmodel/final_medium_map_recon1.xlsx",
        medium_name='dmemf12'
    )
    gem_tmp = sampler.build_objective_candidates()
    
    #
    # load coefficient sampler
    coef = coefsampler(
            single_obj="./scooti/scooti/metabolicmodel/allmets_metabolites_shen2019.csv",
            sample_num=100,
            save_path="/nfs/turbo/umms-csriram/daweilin/fluxprediction/randomobjcoef/synthetic_data/",
            suffix="allmets_recon1",
            func="single"
            )
    coef_path = pd.read_csv("/nfs/turbo/umms-csriram/daweilin/fluxprediction/randomobjcoef/synthetic_data/singleobjcoef_singleobj.csv")

    if coef_path:
        obj_coef = pd.read_csv(coef_path, index_col=0)
        sampler.assign_multi_objectives(obj_coef, gem_tmp, sample_num=5,
            rootpath="/nfs/turbo/umms-csriram/daweilin/fluxprediction/unconstrained_models/fluxsampling_test/")
    else:
        sampler.assign_single_objectives(gem_tmp, sample_num=5,
            rootpath="/nfs/turbo/umms-csriram/daweilin/fluxprediction/unconstrained_models/fluxsampling_test/")

