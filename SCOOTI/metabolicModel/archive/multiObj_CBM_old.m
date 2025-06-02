function multiObj_CBM(jj, DFA_kappa, CFR_kappa, CFR_rho, COBRA_path, GEM_path, model_name, obj_candidate_list_file, input_obj_tb, paraLen, random_para, init_objective, genekoflag, rxnkoflag, FSflag, pfba, medium_perturbation, data_dir, prefix_name, medium, late_stage, early_stage, simulation, constraint, save_root_path, CFR_model_path, pairwise_CFR_model, extraWeight, algorithm, data_series, prefix_series, medium_series)
  % multiObj_CBM.m calls CFRinterface/DFAinterface to run CFR/DFA modeling and generate flux prediction based on transcriptomics/metabolomics data.
  %   This code is capable of predicting fluxes with GEMs Recon1/Recon2.2/Recon3D, multi-objectives, gene knockouts, and with/without pFBA.
  % Inputs:
  %   jj : integer, index number that is used to scan parameters. Default value is set to 1
  %   DFA_kappa : float, parameter to minimize flux activity coefficients
  %   CFR_kappa : float, parameter to reduce the upper bound of reactions associated with downregulated genes
  %   CFR_rho : float, parameter to increase the lower bound of reactions associated with upregulated genes
  %   COBRA_path : string, path to access the openCOBRA toolbox
  %   GEM_path : string, path to access the genome-scale metabolic models like Recon1
  %   model_name : string, name of the GEM, options including Recon1, Recon2.2, and Recon3D
  %   obj_candidate_list_file : string, path to access the candidate of multi-objectives, only applicable for single-objective models
  %   input_obj_tb : string, path to access the file of objective coefficients, only applicable for multiobjective models
  %   paraLen : integer, number to indicate how how many parameters that will be selected for scanning
  %   random_para : bool, enable random sampling of parameters in defined ranges
  %   init_objective : integer, the number is the index of objective list. In the default settings, 1 is no objective and 2 is biomass
  %   genekoflag : bool, enable in silico single-gene knock
  %   rxnkoflag : bool, enable in silico reaction knockout
  %   FSflag : bool, output flux sampling from solution space
  %   pfba : bool, minimizing sum of reaction fluxes or not
  %   medium_perturbation : bool, enable in silico single-metabolite depletion
  %   data_dir : string, path to access significant genes/proteins/metabolites as constraints
  %   prefix_name : string, name of experiment. it will be used for 
  %   medium : string, name of sources of extracellular metabolites. options including DMEMF12 and KSOM
  %   late_stage : string, name of late stage which also indicates the suffix of constraint files
  %   early_stage : string, name of early stage which also indicates the suffix of constraint files
  %   simulation : string, indicate which type of modeling method. options including CFR and DFA
  %   constraint : bool, enable adding constraints or not
  %   save_root_path : string, path to save the results of predicted fluxes
  %   CFR_model_path : string. path to access the constrained CFR models. This is only for stacking/multiomics integration
  %   pairwise_CFR_model : bool, match the name of the contrained CFR models or not
  %   extraWeight : float, apply to scale the importance of contrained CFR models for stacking different constraints
  %   algorithm : string, default value is iMAT but options including iMAT, MOOMIN, INIT (currently not working)
  %   data_series : string, default value is an empty string and the input paths should be separated by ","
  %   prefix_series : string, default value is an empty string and the input paths should be separated by ","
  %   medium_series : string, default value is an empty string and the input paths should be separated by ","
  % Output: 
  %   None; however, the predicted fluxes will be saved
 
%% Part 0. default settings for DFA and CFR parameters
  if (~exist('jj','var')) || (isempty(jj))
      jj = 1;
  end
  if (~exist('input_obj_tb','var')) || (isempty(input_obj_tb))
      input_obj_tb = '';
  end
  if (~exist('DFA_kappa','var')) || (isempty(DFA_kappa))
      dkappa = 1;
  else,
      dkappa = DFA_kappa;
  end
  if (~exist('CFR_kappa','var')) || (isempty(CFR_kappa))
      ckappa = 1;
  else,
      ckappa = CFR_kappa;
  end
  if (~exist('CFR_rho','var')) || (isempty(CFR_rho))
      crho = 1;
  else,
      crho = CFR_rho;
  end
  if (~exist('extraWeight','var')) || (isempty(extraWeight))
      extraWeight = 1;
  end
  if (~exist('CFR_model_path','var')) || (isempty(CFR_model_path))
      CFR_model_path = {};
  end
  if (~exist('pairwise_CFR_model','var')) || (isempty(pairwise_CFR_model))
      pairwise_CFR_model = 0;
  end
  if (~exist('init_objective','var')) || (isempty(init_objective))
      init_objective = 1;
  end
  if (~exist('algorithm','var')) || (isempty(algorithm))
      algorithm = 'iMAT';
  end
  if strcmp(prefix_name, 'model')
    simulation = 'CFR'
    dkappa = -1
    ckappa = 1
    crho = 1
    constraint = 0
    paraLen=1
    random_para=0
  end
  if (~exist('genekoflag','var')) || (isempty(genekoflag))
      genekoflag = 0;
  end
  if (~exist('rxnkoflag','var')) || (isempty(rxnkoflag))
      rxnkoflag = 0;
  end
  if (~exist('FSflag','var')) || (isempty(FSflag))
      FSflag = 0;
  end
  if (~exist('medium_perturbation','var')) || (isempty(medium_perturbation))
      medium_perturbation = 0;
  end
  if (~exist('pfba','var')) || (isempty(pfba))
      pfba = 1;
  end
  if (~exist('late_stage','var')) || (isempty(late_stage))
      late_stage = 'upgenes';
  end
  if (~exist('early_stage','var')) || (isempty(early_stage))
      early_stage = 'dwgenes';
  end
  if (~exist('data_series','var')) || (isempty(data_series))
    data_series = {};
    prefix_series = {};
    medium_series = {};
  else,
    data_series_tmp = strsplit(data_series, ',')
    prefix_series_tmp = strsplit(prefix_series, ',');
    medium_series_tmp = strsplit(medium_series, ',');
    data_series = {}
    prefix_series = {}
    medium_series = {}
    for kk=1:length(data_series_tmp),
      data_series(kk, 1) = data_series_tmp(1, kk);
      prefix_series(kk, 1) = prefix_series_tmp(1, kk);
      medium_series(kk, 1) = medium_series_tmp(1, kk);
    end
  end
%% Part 1: metabolic model
  % initialize COBRA toolbox
  addpath(COBRA_path)
  run initCobraToolbox(false)
  changeCobraSolver('gurobi')
  % INPUT files
  % Load metabolic network model
  model = load(GEM_path);
  fn = fieldnames(model);
  model = getfield(model, fn{1});
  % fix the model for the default objective
  if strcmp(model_name, 'Recon1'),
    model.c(find(contains(model.rxns, 'biomass_objective'))) = 1;
  end
  % objective candidates
  if length(obj_candidate_list_file)==0,
    obj_candidate_list_file = './obj52_metabolites_shen2019.csv';
  end
% Part 2: Input data
  % random sampling
  reind = randsample(paraLen*paraLen, paraLen);
  input_data_choose = prefix_name;

  % additional input for objectives and corresponding weights
  if length(input_obj_tb)==0,
    input_obj_file = './buffer_objective_coefficients.csv';
    input_objective_weights = 0; % set to 1 if giving existing objective weights
  else,
    input_obj_file = input_obj_tb;
    input_objective_weights = 1;
  end
  
  % read coefficients for objectives
  input_obj_tb = readtable(input_obj_file);

  %% Significant genes/proteins/metabolites
  if constraint==1,
    % data input setting
    if length(data_series)==0,
      [data_series, prefix_series, medium_series] = batch_input_preprocess(data_dir, prefix_name, medium);
    end

    ctrl = 1 %set to 1 if we want to apply constraints
  else, % CFR without constraint
    data_series = {'./metabolicModel/buffer_constraint.xlsx'};
    prefix_name = 'unconstraint';
    prefix_series = {'2CBC_Israel'};
    ctrl = 0; % remove constraint
    medium_series = {medium};
  end


  % reuse models with previous constraints
  CFR_models = {}
  if length(CFR_model_path)>0,
    CFR_model_dir = dir(CFR_model_path)
    for n=3:length(CFR_model_dir),
      if contains(CFR_model_dir(n).name, '.mat'),
        CFR_models{size(CFR_models, 1)+1, 1} = sprintf('%s/%s', CFR_model_path, CFR_model_dir(n).name);
      end
    end
  end
  disp(CFR_models)

%% Part 3: Parameter settings for CBM 
  if strcmp(simulation, 'DFA')==1, % optimize reactions/new demand reactions
%% settings for parameter space for DFA kappa
    if dkappa==-1,
      DFA_para_space = logspace(1, -3, paraLen);
      disp('Parameter settings for DFA...')
      if random_para==1,
        DFA_para_space = logspace(1, -3, paraLen*paraLen);
        DFA_para_space = DFA_para_space(reind);
      end
      kappa = DFA_para_space(jj)
      rho = 0; % for bugs in DFA (unnecessary parameter)
    else,
      kappa = dkappa;
      rho = 0;
    end
  else, % CFR or model
%% settings for parameter space of CFR kappa and rho
    if ckappa==-1,
      CFR_para_space = logspace(1, -3, paraLen);
      [Y, Z] = meshgrid(CFR_para_space, CFR_para_space);
      CFR_para_pack = [Y(:), Z(:)];
      disp('Parameter settings for CFR...')
      if random_para==1,
        CFR_para_pack = CFR_para_pack(reind, :);
      end
      kappa = CFR_para_pack(jj, 1);
      rho = CFR_para_pack(jj, 2);
    else,
      kappa = ckappa;
      rho = crho;
    end
  end

%% Part4: Iterate through all different input datasets
  % setups for parallel computing
  %delete(gcp('nocreate'));
  %parpool(32);
  %disp(length(data_series))
  poolobj = gcp("nocreate"); % If no pool, do not create new one.
  if isempty(poolobj)
      poolsize = 0;
  else
      poolsize = poolobj.NumWorkers
  end
%% Parfor...
  %for data_index=1:length(data_series),
  for data_index=1:length(data_series),
    %data_index = 1;
    %size(data_series)
    %size(medium_series)
    data_path = data_series{data_index};
    %disp('Read file...')
    %disp(data_path)
    % info of the data/simulation
    prefix = prefix_series{data_index};
    % get cell-specific objective with corresponding cell names
    %if input_objective_weights & strcmp(input_data_choose, 'model')==0,
    %  prefix_str = prefix_pattern_func(prefix);
    %  input_obj_tb_ind = find([logical(1) contains(input_obj_tb.Properties.VariableNames(2:end), prefix_str)]);
    %end
    % change culture medium KSOM_AA_Jin or DMEMF12_MGSA
    medium = medium_series{data_index};
    
    % candidate objectives
    single_obj = {};
    single_obj{1, 1} = '';
    candidates = readtable(obj_candidate_list_file);
    candidate_size = size(candidates);
    for l=1:candidate_size,
      single_obj(l+1, 1) = candidates{l, 2};
    end
    
    % objective weights for single objectives
    para_pack = eye(length(single_obj));
   
%% CFR simulations with different objective weights
    % parallel computing settings
    %init_for = 1;
    %disp('input_obj_check')
    %disp(input_objective_weights)
    if input_objective_weights==0, % single objective
      init_for = init_objective;
      if strcmp(input_data_choose, 'model')==1, % iterate thru all obj
        init_for = 2; % change the initial index to start with 'gh' instead of '' (unconstrained models)
        end_for = length(single_obj);
      else, % without obj for constraint models
        end_for = init_for;%length(single_obj);%1
      end
    elseif strcmp(input_data_choose, 'model')==1,
      % use new objectives without constraint
      init_for = 2;
      ss = size(input_obj_tb);
      end_for = ss(2);
    else,
      % use new objective WITH constraint
      init_for = 2;%init_for;
      %ss = size(input_obj_tb(1:end, input_obj_tb_ind));
      ss = size(input_obj_tb(1:end, 1:end));
      end_for = ss(2);
    end
    %delete(gcp('nocreate'))
    %parpool(36);
    %for ii=init_for:end_for,
    for ii=init_for:end_for,

      % use single objective functions
      if input_objective_weights==0,
        % for single objective
        obj = single_obj; % comment out this one if you to model
        % objective weights
        %disp('Weights for objectives:')
        obj_c = [];
        for d=1:length(obj),
          obj_c = [obj_c, para_pack(ii, d)];
        end

      % use multiobjective functions from files without constraint
      elseif strcmp(input_data_choose, 'model')==1,
        
        obj = input_obj_tb{2:ss(1), 1};
        obj_c = input_obj_tb{2:ss(1), ii}; 
        get_late_stage = input_obj_tb.Properties.VariableNames{ii};

      % use multiobjective function from files WITH corresponding constraint
      else,
        obj = input_obj_tb{2:ss(1), 1};
        obj_c = input_obj_tb{2:ss(1), ii};
        get_late_stage = late_stage;

      end
      % DO NOT CHANGE ANYTHING BELOW
      % Running exps for all conditions
      %disp('Optimizing the demand reaction:')
      %disp(obj)
      if strcmp(obj, ''),
          obj_type = '';
      else,
          obj_type = 'Demand';
      end
      
      
      %%%%%%%%%%%%%%%%%%
      % start simulation
      %%%%%%%%%%%%%%%%%%
      % later stages
      
      % filename logics: 
      % 1: '[datetime]Project_simulation_fluxes.csv'
      % 2: '[datetime]Project_simulation_metadata.json'
      % 3: Other outputs
      out_name = sprintf('%s_ct%d_obj%d_data%d', simulation, jj, ii, data_index);
      disp('Results will be saved in:')
      disp(out_name)
      
      if strcmp(simulation, 'DFA')==1, 

        late_stage = prefix
        disp('DFA...')
        %rho = 0
        % change culture medium KSOM_AA_Jin or DMEMF12_MGSA
        if length(CFR_models)==0,
          DFAinterface(GEM_path, obj, obj_type, obj_c, save_root_path, data_path, late_stage, early_stage, out_name, ctrl, kappa, genekoflag, rxnkoflag, medium, medium_perturbation, FSflag, model_name, '', 0, algorithm);
        else,
          for iii=1:length(CFR_models),
            out_name_mod = sprintf('%s_recon%d', out_name, iii);
            DFAinterface(GEM_path, obj, obj_type, obj_c, save_root_path, data_path, late_stage, early_stage, out_name, ctrl, kappa, genekoflag, rxnkoflag, medium, medium_perturbation, FSflag, model_name, CFR_models{iii, 1}, extraWeight, algorithm);
          end
        end
      else,
      %simulation=='CFR' | simulation=='model',
        % CFR with fluxes
        if input_objective_weights==0,
          if length(CFR_models)==0,
            CFRinterface(GEM_path, pfba, obj, obj_type, obj_c, save_root_path, data_path, out_name, late_stage, early_stage, ctrl, kappa, rho, medium, genekoflag, rxnkoflag, medium_perturbation, FSflag, model_name, '', 0, algorithm);
          else, % multiple CFR_models
            %for iii=1:length(CFR_models),
            %delete(gcp('nocreate'));
            %parpool(36);
            %parfor iii=1:length(CFR_models),
            for iii=1:length(CFR_models),
              disp(iii)
              out_name_mod = sprintf('%s_recon%d', out_name, iii);
              if pairwise_CFR_model==1,
                % match sample names
                cfr_model = CFR_models{iii, 1};
                fprefix = strsplit(cfr_model, '_model_CFR');
                fname = sprintf('%s_metadata.json', fprefix{1});
                fid = fopen(fname); 
                raw = fread(fid, inf); 
                str = char(raw'); 
                fclose(fid); 
                val = jsondecode(str);
                ftar = strsplit(val.input_path, '/');
                fsample = strsplit(data_path, '/');
                % run CFR only when the names are matched
                if strcmp(ftar{end}, fsample{end}),
                  CFRinterface(GEM_path, pfba, obj, obj_type, obj_c, save_root_path, data_path, out_name_mod, late_stage, early_stage, ctrl, kappa, rho, medium, genekoflag, rxnkoflag, medium_perturbation, FSflag, model_name, CFR_models{iii, 1}, extraWeight, algorithm);
                end % else skip the model
              else, % run all the model without pairing the samples
                CFRinterface(GEM_path, pfba, obj, obj_type, obj_c, save_root_path, data_path, out_name_mod, late_stage, early_stage, ctrl, kappa, rho, medium, genekoflag, rxnkoflag, medium_perturbation, FSflag, model_name, CFR_models{iii, 1}, extraWeight, algorithm);
              end % end if pairwise_CFR_model
            end % end for the for loop of CFR models
            end % end for using CFR_models or not
        else, % with input obj weight
          CFRinterface(GEM_path, pfba, obj, obj_type, obj_c, save_root_path, data_path, out_name, get_late_stage, early_stage, ctrl, kappa, rho, medium, genekoflag, rxnkoflag, medium_perturbation, FSflag, model_name, '', 0, algorithm);
        end % end for using input obj weight or not
      end % end for DFA or CFR
    end % end for the for-loop of data_series
  end % end for the function





