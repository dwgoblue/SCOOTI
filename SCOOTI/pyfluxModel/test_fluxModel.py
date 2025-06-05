import os
import pandas as pd
from SCOOTI.pyfluxModel.fluxModel import MultiConstraintModelSetter

def test_multi_constraint_batch():
    # Define input paths
    example_dir = "./SCOOTI/SCOOTI/examples/example_sigGenes/"
    GEM_path = "./SCOOTI/SCOOTI/metabolicModel/GEMs/Shen2019.xml"
    objective_path = "./SCOOTI/SCOOTI/metabolicModel/GEMs/obj52_metabolites_shen2019.csv"
    medium_path = "./SCOOTI/SCOOTI/metabolicModel/FINAL_MEDIUM_MAP_RECON1.xlsx"
    medium_name = "DMEMF12"
    output_dir = "./SCOOTI/SCOOTI/examples/Sampling/test_batch_output/"

    os.makedirs(output_dir, exist_ok=True)

    # Instantiate the modeler
    modeler = MultiConstraintModelSetter(
        GEM_path=GEM_path,
        objective_path=objective_path,
        medium_path=medium_path,
        medium_name=medium_name,
        mapping_file='./SCOOTI/SCOOTI/metabolicModel/GEMs/recon1_genes.json'
    )

    # Load demand reactions
    modeler.build_objective_candidates(update_model=True)

    # Detect file prefixes like "experimentA"
    fs = os.listdir(example_dir)
    prefixes = sorted(set(
        f.split("_upgenes")[0].split("_dwgenes")[0]
        for f in fs if f.endswith(".csv")
    ))

    # Build full path to prefix files (e.g., "SCOOTI/examples/example_sigGenes/experimentA")
    file_prefix_array = [os.path.join(example_dir, p) for p in prefixes]

    # Apply constraints and run
    modeler.apply_constraint_batch(
        file_prefix_array=file_prefix_array,
        rootpath=output_dir,
        sample_flux=False,
        obj_coef=None,
        sample_num=3,
        para_config={1: {'rho': 1, 'epsilon1': 0.001, 'kappa': 1, 'epsilon2': 0.001}}
    )

if __name__ == "__main__":
    test_multi_constraint_batch()

