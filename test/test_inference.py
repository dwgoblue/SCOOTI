# -*- coding: utf-8 -*-

def test_regression_training():

    # packages
    from SCOOTI.regressorTraining import *
    regression_exe = regressorTraining(
        '../SCOOTI/examples/example_fluxPreduction/example_unconstrained_moodels/',
        '../SCOOTI/examples/example_fluxPreduction/example_constrained_models/',
        './',
        kappa_arr='10',
        rho_arr='0.1'
        dkappa_arr='-1',
        expName='regression_test',
        uncon_norm=True,
        con_norm=False,
        medium='DMEMF12',
        method='cfr',
        model='recon1',
        input_type='flux',
        cluster_path='',
        objList_path='../',
        rank=False,
        stack_model=False,
        learner='L',
        geneKO=False,
        geneList_path='',
        learning_rate=0.001,
        epo=10
    )





import subprocess
def test_cli_with_subprocess():
    result = subprocess.run(
        [
            "python",
            "../SCOOTI/SCOOTI_trainer.py",
            "--unconModel",
            "./examples/unconstrained_models/",
            "--conModel",
            "./examples/example_fluxPrediction/",
            "--savePath",
            "./examples/example_regressionModels/",
            "--kappaArr",
            "10,1,0.1",
            "--rhoArr",
            "10,1,0.1",
            "--dkappaArr",
            "10,1,0.1",
            "--expName",
            "test_run",
            "--unconNorm",
            "T",
            "--conNorm",
            "F",
            "--medium",
            "DMEMF12",
            "--method",
            "cfr",
            "--model",
            "recon1",
            "--inputType",
            "flux",
            "--clusterPath",
            "",
            "--objListPath",
            "",
            "--rank",
            "F",
            "--stackModel",
            "F",
            "--sampling",
            "F",
            "--learner",
            "L",
            "--geneKO",
            "F",
            "--geneListPath",
            "",
            "--learningRate",
            0.001,
            "--epo",
            5000
            ],
        capture_output=True,
        text=True
    )
    assert result.returncode == 0
    assert "completed" in result.stdout




