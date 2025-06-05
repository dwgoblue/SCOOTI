"""
fluxModel.py
=======================================================
Analysis of metabolic objectives and fluxes in diseases
"""
import os
import json
import time
import warnings
from datetime import datetime
from tqdm import tqdm
import pandas as pd
import numpy as np
import cobra
from cobra.sampling import sample
import cobra.flux_analysis.parsimonious as pFBA
from cobra.manipulation.delete import knock_out_model_genes
#from metabolicmodelpipeline.utils.cfr_suite import cfr_optimize, apply_cfr, summarize, get_fluxes, process_flux

warnings.simplefilter('ignore')
cobra.Configuration().solver = "glpk"


class modelSetter:
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
    def __init__(self, GEM_path, objective_path, medium_path, medium_name, mapping_file='./SCOOTI/SCOOTI/metabolicModel/GEMs/recon1_genes.json'):
        self.GEM_path = GEM_path
        self.objective_path = objective_path
        self.medium_path = medium_path
        self.medium_name = medium_name
        self.mapping_file = mapping_file
        # load metabolic network model
        self.model_loader = {
                'mat':cobra.io.load_matlab_model,
                'xml':cobra.io.read_sbml_model
                }
        self.gem = self.load_model()
        # load objective metabolites
        self.objectives = self.load_objectives()

        # Initialize the list attributes
        self.objective_candidates = []
        self._cfr_on = []
        self._cfr_off = []
        self._cfr_on_params = -1
        self._cfr_off_params = -1

    def load_model(self):
        """
        Load and configure the GEM model by applying the culture medium constraints.

        Returns:
            cobra.Model: Configured GEM model with updated medium.
        """
        print(self.GEM_path.split('.')[-1])
        gem = self.model_loader[self.GEM_path.split('.')[-1]](self.GEM_path)
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

    def build_objective_candidates(self, update_model=True):
        """
        Build demand reactions for all candidate metabolites and append them to the GEM.

        Returns:
            cobra.Model: Modified GEM with added demand reactions.
        """
        compartments = ['c', 'm', 'n', 'x', 'r', 'g', 'l']
        gem_tmp = self.gem.copy()
        print('Number of reactions before updates:', len([r for r in gem_tmp.reactions]))

        print('Adding demand reactions of single objectives...')
        for obj in self.objectives['metabolites']:
            if obj == 'gh':
                self.objective_candidates.append('biomass_objective')
            else:
                gem_mets = [met.id for met in gem_tmp.metabolites]
                met_ids = [
                        f'{obj}[{c}]' for c in compartments if f'{obj}[{c}]' in gem_mets
                        ]
                met_objs = [gem_tmp.metabolites.get_by_id(met) for met in met_ids]
                print('Add demand reaction:', met_ids)
                reaction = cobra.Reaction(f'{obj}_demand')
                reaction.name = f'Objective candidate {obj}'
                reaction.add_metabolites({
                    met_obj: -1 for met_obj in met_objs
                    })
                gem_tmp.add_reactions([reaction])
                self.objective_candidates.append(f'{obj}_demand')

        print('Number of reactions after updates:', len([r for r in gem_tmp.reactions]))
        if update_model:
            print('Updating the model with new demand reactions...')
            print('[WARNING] This will not update objective functions in the model.')
            self.gem = gem_tmp.copy()
        return gem_tmp


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
            'CFR_kappa': self._cfr_off_params[0],
            'CFR_rho': self._cfr_on_params[0],
            'CFR_epsilon1': self._cfr_on_params[1],
            'CFR_epsilon2': self._cfr_off_params[1],
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
        filename = os.path.join(
                root_path,
                f'[{file_prefix}]{out_name}'
                )
        with open(f'{filename}_metadata.json', 'w') as f:
            json.dump(metadata, f)
        print('Metadata saved at:', filename)
        return filename

    def apply_constraints(
            self,
            on_list=[],
            off_list=[], 
            on_params=(0.01, 0.001),
            off_params=(0.01, 0.001)
            ):
        """
        Apply structural constraints to self.gem based on active/inactive reactions or genes.
        No objective function is changed or solved here.
    
        Args:
            on_list (list): Genes or reactions to encourage (keep active).
            off_list (list): Genes or reactions to suppress (shut down).
            on_params (tuple): (rho, epsilon1) for 'on' reactions.
            off_params (tuple): (kappa, epsilon2) for 'off' reactions.
        """
        model = self.gem
    
        def parse_rxns(rxn_list):
            if not rxn_list:
                return []
            elif any(item in model.reactions for item in rxn_list):
                return [rxn.id if hasattr(rxn, 'id') else rxn for rxn in rxn_list]
            elif any(item in model.genes for item in rxn_list):
                with model:
                    return [rxn.id for rxn in knock_out_model_genes(model, rxn_list)]
            else:
                return []
        if self.mapping_file is not None:
            on_list, off_list = self.convert_gene_list(on_list, off_list, self.mapping_file)
        on_rxns = parse_rxns(on_list)
        off_rxns = parse_rxns(off_list)
    
        for rxn_id in on_rxns:
            if rxn_id in model.reactions:
                rxn = model.reactions.get_by_id(rxn_id)
                if abs(rxn.lower_bound) < on_params[1] and rxn.lower_bound < 0:
                    rxn.lower_bound = -on_params[1]
                if rxn.upper_bound < on_params[1]:
                    rxn.upper_bound = on_params[1]
    
        for rxn_id in off_rxns:
            if rxn_id in model.reactions:
                rxn = model.reactions.get_by_id(rxn_id)
                rxn.lower_bound = min(rxn.lower_bound, off_params[1])
                rxn.upper_bound = max(rxn.upper_bound, -off_params[1])
                if abs(rxn.lower_bound) < off_params[1]:
                    rxn.lower_bound = 0.0
                if abs(rxn.upper_bound) < off_params[1]:
                    rxn.upper_bound = 0.0
    
        self._cfr_on = on_rxns
        self._cfr_off = off_rxns
        self._cfr_on_params = on_params
        self._cfr_off_params = off_params
    
        print(f"CFR constraints applied: {len(on_rxns)} ON, {len(off_rxns)} OFF reactions.")
    
    @staticmethod
    def add_soft_flux_penalty(model, excluded_rxns, penalty=1e-6):
        """
        Apply soft penalty to all reactions not in excluded_rxns.
        """
        penalty_targets = [
            (rxn.forward_variable, rxn.reverse_variable)
            for rxn in model.reactions
            if rxn.id not in excluded_rxns
        ]
        linear_coeffs = {v: -penalty for pair in penalty_targets for v in pair}
        model.objective.set_linear_coefficients(linear_coeffs)
    
    def assign_single_objectives(self, sample_num=20, rootpath='./', sample_flux=True, scan_id=None, data_id=None):
        """
        Assign each candidate objective individually and sample or optimize.
    
        Args:
            sample_num (int): Number of samples (or runs).
            rootpath (str): Directory to save flux data.
            sample_flux (bool): If True, sample fluxes; otherwise optimize.
        """
        for i, candidate in enumerate(self.objective_candidates):
            gem_tmp = self.gem.copy()
            gem_tmp.objective = candidate
            self.add_soft_flux_penalty(gem_tmp, excluded_rxns=[candidate])
    
            print('testing', gem_tmp.optimize().fluxes.sum())
            if sample_flux:
                samples = pd.DataFrame(sample(gem_tmp, sample_num, processes=sample_num))
                samples = samples.sample(frac=1)
                samples['Obj'] = samples[candidate].to_numpy()
                samples.index = np.arange(len(samples))
                samples = samples.T
                samples.index = [r.id for r in gem_tmp.reactions]
            else:
                sol = gem_tmp.optimize()
                samples = pd.DataFrame(sol.fluxes)
                samples.columns = ['flux']
            print(candidate, samples.head())
            obj_c = np.zeros(len(self.objectives['metabolites']))
            obj_c[self.objectives['metabolites'].to_numpy() == candidate] = 1.0
            for col in samples.columns:
                out_path = self.save_flux_sample(
                    sample_vector=samples[col],
                    candidate_name=candidate,
                    obj_c=obj_c,
                    rootpath=rootpath,
                    sample_id=col,
                    obj_id=i,
                    sample_flux=sample_flux,
                    scan_id=scan_id,
                    data_id=data_id
                )


    def assign_multi_objectives(self, obj_coef, sample_num=20, rootpath='./', sample_flux=True, scan_id=None, data_id=None):
        """
        Assign linear combinations of objectives and sample or optimize.
    
        Args:
            obj_coef (pd.DataFrame): Reaction weights (columns = samples).
            sample_num (int): Number of samples (or runs).
            rootpath (str): Output folder.
            sample_flux (bool): Whether to sample (True) or optimize (False).
        """
        obj_df = obj_coef.copy()
        obj_df.index = obj_df.index.to_series().apply(lambda x: f'{x}_demand' if x != 'gh' else 'biomass_objective')
    
        print('Multiobjective coefficient assignment...')
        for obj_id, col in tqdm(enumerate(obj_df.columns)):
            gem_tmp = self.gem.copy()
            obj_dict = {gem_tmp.reactions.get_by_id(k): v for k, v in obj_df[col].items()}
            for r, c in obj_dict.items():
                r.objective_coefficient = c
    
            self.add_soft_flux_penalty(gem_tmp, excluded_rxns=obj_df[col][obj_df[col] != 0].index)
    
            if sample_flux:
                samples = pd.DataFrame(sample(gem_tmp, sample_num, processes=sample_num))
                samples = samples.sample(frac=1)
                samples['Obj'] = samples[obj_df[col][obj_df[col] > 0].index].sum(axis=1).to_numpy()
                samples.index = np.arange(len(samples))
                samples = samples.T
            else:
                sol = gem_tmp.optimize()
                samples = pd.DataFrame(sol.fluxes)
                samples.columns = ['flux']
    
            print(samples.head())
            for s in tqdm(samples.columns):
                obj_c = self.objectives['metabolites'].apply(lambda x: obj_df[col].get(x, 0)).to_numpy()
                out_path = self.save_flux_sample(
                    sample_vector=samples[s],
                    candidate_name=obj_df[col][obj_df[col] > 0].index[0],
                    obj_c=obj_c,
                    rootpath=rootpath,
                    sample_id=s,
                    obj_id=obj_id,
                    scan_id=scan_id,
                    data_id=data_id
                )
                
    def convert_gene_list(self, on_list, off_list, mapping_file='./SCOOTI/SCOOTI/metabolicModel/GEMs/recon1_genes.json'):
        # CRITICAL: Make sure the genes in the on_list and off_list are present in the model.
        # Most of the time, the genes in the model are in the form of bigg ids, e.g. 'b0001'.
        # If you have a list of gene names, you need to convert them to bigg ids.
        # get dictionary of mapping bigg ids to genes
        with open(mapping_file) as f:
            j = json.load(f)
        
        ## load up/down example gene lists
        #on_list = pd.read_csv(
        #        './examples/example_sigGenes/johnson_18_GSE117444_upgenes.csv', index_col=0
        #        )
        #off_list = pd.read_csv(
        #        './examples/example_sigGenes/johnson_18_GSE117444_dwgenes.csv', index_col=0
        #        )
        #print('List the first 5 genes in the on_list:', on_list.head())
        # convert name into id
        converted_on_list = []
        converted_off_list = []
        for k, v in j.items():
            if v in on_list:
                converted_on_list.append(k)
            if v in off_list:
                converted_off_list.append(k)
        return converted_on_list, converted_off_list



    def save_flux_sample(
            self,
            sample_vector,
            candidate_name,
            obj_c,
            rootpath,
            sample_id=None, 
            obj_id=None,
            scan_id=None,
            data_id=None,
            sample_flux=True
            ):
        """
        Save a single flux sample or optimized flux and its metadata.
    
        Args:
            sample_vector (pd.Series): The flux vector.
            candidate_name (str): Objective reaction name.
            obj_c (np.array): Objective weight vector.
            rootpath (str): Output directory.
            sample_id (str): Sample ID or group name.
            obj_id (str): objective ID or group name.
            scan_id (str): parascan ID or group name.
            data_id (str/int): Optional ID to include in the filename (for dataset source).
            sample_flux (bool): Whether this is a flux sample (True) or an optimal flux (False).
        
        Returns:
            str: Full path to saved flux CSV.
        """
        if sample_flux:
            folder = os.path.join(rootpath, f'fs_{sample_id}')
        else:
            folder = rootpath  # Save directly into rootpath
    
        os.makedirs(folder, exist_ok=True)
    
        file_prefix = time.strftime("%b%d%Y%H%M%S")
        data_suffix = f'_data{data_id}' if data_id is not None else 1
        parascan_prefix = f'_ct{scan_id}' if scan_id is not None else 1
        obj_middle = f'_obj{obj_id}' if obj_id is not None else 1
        out_name = f'model_{parascan_prefix}{obj_middle}{data_suffix}'
    
        excelname = self.save_metadata(
            candidate_name,
            self.objectives['metabolites'].to_numpy(),
            obj_c.astype(float),
            folder,
            self.GEM_path,
            out_name
        )
    
        df = pd.DataFrame(sample_vector, columns=['flux'])
        flux_path = f'{excelname}_fluxes.csv.gz'
        df.to_csv(flux_path, compression='gzip')
        return flux_path


class MultiConstraintModelSetter(modelSetter):

    def read_csv_batch(parent_path):
        """
        Read a CSV file and return a DataFrame.
        
        Args:
            parent_path (str): Path to the directory containing the CSV file.
        
        Returns:
            pd.DataFrame: DataFrame containing the CSV data.
        """
        fs = os.listdir(parent_path)
        fs = [f for f in fs if f.endswith('.csv')]
        file_prefix_array = np.unique([f.split('_upgenes')[0].split('_dwgenes')[0] for f in fs])
        return fs, file_prefix_array

    def apply_constraint_batch(
            self,
            file_prefix_array,
            upfile_suffix='upgenes',
            dwfile_suffix='dwgenes',
            rootpath='./',
            sample_flux=False,
            obj_coef=None,
            sample_num=20,
            para_config={0:{'rho':0.01, 'epsilon1':0.001, 'kappa':0.01, 'epsilon2':0.001}},
            ):
        """
        Apply multiple sets of CFR constraints defined in up/down gene CSVs.

        Args:
            upgene_file (str): CSV file with up-regulated genes (rows = genes, cols = samples).
            downgene_file (str): CSV file with down-regulated genes (same format).
            rootpath (str): Output folder for saving.
            sample_flux (bool): Whether to sample (True) or solve (False).
        """

        gem_tmp = self.gem.copy() # placehold clean model
        for scan_id in tqdm(para_config.keys(), desc='Apply different parameters'):
            for data_id, file_prefix in tqdm(enumerate(file_prefix_array), desc='Applying constraints'):
                # assign clean model to self.gem
                self.gem = gem_tmp.copy()
                # Read up/down gene files
                df_up = pd.read_csv(f'{file_prefix}_{upfile_suffix}.csv', index_col=0)
                df_down = pd.read_csv(f'{file_prefix}_{dwfile_suffix}.csv', index_col=0)

                # Ensure the data_id is in the DataFrame columns
                up_genes = df_up.to_numpy()
                down_genes = df_down.to_numpy()

                # Clone model and apply constraints
                self.apply_constraints(
                        on_list=up_genes,
                        off_list=down_genes,
                        on_params=(para_config[scan_id]['rho'], para_config[scan_id]['epsilon1']),
                        off_params=(para_config[scan_id]['kappa'], para_config[scan_id]['epsilon2'])
                        )

                # Assign objective(s) and run (e.g. biomass)
                gem_tmp = self.gem.copy()
                # run simulation
                if obj_coef is not None:
                    self.assign_multi_objectives(
                            obj_coef,
                            sample_num=sample_num,
                            rootpath=rootpath,
                            sample_flux=sample_flux,
                            scan_id=scan_id,
                            data_id=data_id
                            )
                else:
                    self.assign_single_objectives(
                            sample_num=sample_num,
                            rootpath=rootpath,
                            sample_flux=sample_flux,
                            scan_id=scan_id,
                            data_id=data_id
                            )





class coefSampler:
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
            self.single_obj = pd.read_csv(
                    single_obj, index_col=0
                    ).iloc[1:].values.flatten()
        else:
            self.single_obj = single_obj
        self.sample_num = sample_num
        self.save_path = save_path
        self.suffix = suffix
    
        coef = self.random_objective_coefficients() if func=='random' else self.single_objective_coefficients()
        self.coef = coef
    
    def random_objective_coefficients(self):
        # sampling and save
        df = pd.DataFrame(np.random.rand(len(self.single_obj), self.sample_num))
        df.index = self.single_obj
        df.columns = [f'sample_{ind}' for ind in np.arange(len(df.columns))]
        df.to_csv(self.save_path+'/'+f'samplingObjCoef_{self.suffix}.csv')
        return df
    
    def single_objective_coefficients(self):
        # sampling and save
        df = pd.DataFrame(
                np.eye(len(self.single_obj)),
                columns=self.single_obj,
                index=self.single_obj)
        df.to_csv(self.save_path+'/'+f'singleObj_{self.suffix}.csv')
        return df




## Example usage:
#if __name__ == "__main__":
##    
#    # load flux sampler
#    modeler = modelSetter(
#        GEM_path="./SCOOTI/SCOOTI/metabolicModel/GEMs/Shen2019.xml",
#        objective_path="./SCOOTI/SCOOTI/metabolicModel/GEMs/obj52_metabolites_shen2019.csv",
#        medium_path="./SCOOTI/SCOOTI/metabolicModel/FINAL_MEDIUM_MAP_RECON1.xlsx",
#        medium_name='DMEMF12'
#    )
#    modeler.build_objective_candidates(update_model=True)
#
#    # add testing constraints
#    on_list = modeler.gem.genes[:200]
#    off_list = modeler.gem.genes[400:620]
#    modeler.apply_constraints(
#            on_list=on_list,
#            off_list=off_list,
#            on_params=(0.01, 0.001),
#            off_params=(0.01, 0.001)
#            )
#    # load coefficient sampler
#    coef = coefSampler(
#            single_obj="./SCOOTI/SCOOTI/metabolicModel/GEMs/obj52_metabolites_shen2019.csv",
#            sample_num=3,
#            save_path="./SCOOTI/SCOOTI/examples/Sampling/example_output/",
#            suffix="obj52_recon1",
#            func="random"
#            )
#    coef_path = "./SCOOTI/SCOOTI/examples/Sampling/example_output/samplingObjCoef_obj52_recon1.csv"
#
#    #if coef_path:
#    # multiple objectives
#    obj_coef = pd.read_csv(coef_path, index_col=0)
#    modeler.assign_multi_objectives(obj_coef, sample_num=3,
#        rootpath="./SCOOTI/SCOOTI/examples/Sampling/example_output/")
#
#    # single objectives
#    modeler.assign_single_objectives(
#            sample_num=3,
#            rootpath="./SCOOTI/SCOOTI/examples/Sampling/example_output/"
#            )
#
