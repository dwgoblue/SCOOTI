function init_constrain(model, )

%% Predict metabolic changes using INIT method
% This script is adapted from Alec's code:
% https://github.com/sriram-lab/quiescence-depth-modulation/blob/main/flux_sampling/quiesc_depth_rat_GEM_init_chrr.m
% transcriptomics-constrained model is created with INIT
%% Initializing COBRA 

clear;

%% compute fluxes for each stage of quiescence
% creates empty cell array for constrained models
quiescence_models = {};

% runs INIT on each quiescence depth to generate series of constrained
% models
for i = 1:length(idx_qds)

    disp(strcat("Processing ", string(i), " of ", string(length(idx_qds)), " samples"))

    norm_data_subset = norm_data(:, idx_qds{i});
    
    alpha = 1.1; beta = 0.6;

    % runs INIT to obtain metabolic model constrained by transcriptomics
    quiescence_models{i} = run_init2(model, alpha, beta, norm_data_subset, fujimaki_19.Gene, 32);

end

flux_avg_array = {};
flux_all_array = {};

% creates parallel pool
pool = parpool('Processes', 9)

environment = getEnvironment();

% runs CHRR on each quiescence depth 
parfor i = 1:length(quiescence_models)

    disp(i)
    num_samples = 1000

    restoreEnvironment(environment)

    % runs CHRR
    [flux_avg_array{i}, flux_all_array{i}] = run_chrr2(model, quiescence_models{i}, num_samples);

end

delete(pool)
