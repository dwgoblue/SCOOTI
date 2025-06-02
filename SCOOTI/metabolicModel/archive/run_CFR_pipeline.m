addpath(genpath('utils'));

function run_CFR_pipeline(config)
  % Unpack config
  arguments
    config.model_path
    config.pfba logical = true
    config.obj
    config.obj_type
    config.obj_c
    config.root_path
    config.data_path
    config.out_name
    config.upsheet
    config.dwsheet
    config.ctrl logical = true
    config.kappa double = 0.3
    config.rho double = 1
    config.medium string = ""
    config.genekoflag logical = false
    config.rxnkoflag logical = false
    config.media_perturbation logical = false
    config.FSflag logical = false
    config.model_name
    config.CFR_model string = ""
    config.extra_weight double = 0.0
    config.algorithm
  end

  % Set up output name and file prefix
  sample_name = config.upsheet + config.dwsheet;
  file_prefix = string(datetime('now','TimeZone','local','Format','MMMdyHHmm'));
  excelname = sprintf('%s/[%s]%s', config.root_path, file_prefix, config.out_name);
  fprintf('[Main] Output will be saved to: %s\n', excelname);

  % Load model and configure medium
  model = load_and_configure_model(config.model_path, config.model_name, config.medium);

  % Load and process gene expression data
  if config.ctrl
    [ups, dws] = load_expression_tables(config.data_path, config.upsheet, config.dwsheet);
    [ups, dws] = process_expression_data(ups, dws, config.model_name);
  else
    ups = {}; dws = {};
  end

  % Save metadata
  save_metadata_json(excelname, model, config.obj, config.obj_type, config.obj_c, ...
    config.root_path, config.data_path, config.out_name, config.ctrl, config.kappa, ...
    config.rho, config.medium, config.genekoflag, config.rxnkoflag, ...
    config.media_perturbation, config.upsheet, config.dwsheet, ...
    config.model_path, config.CFR_model, config.extra_weight, config.algorithm);

  % Load CFR model if provided
  if strlength(config.CFR_model) > 0
    raw = load(config.CFR_model);
    fn = fieldnames(raw);
    recon_model = raw.(fn{1});
  else
    recon_model = '';
  end

  % Run main simulation
  metadata = loadjson(sprintf('%s_metadata.json', excelname)); % or pass metadata directly
  run_simulation(model, excelname, config.algorithm, ups, dws, ...
                 config.ctrl, config.FSflag, config.kappa, config.rho, ...
                 recon_model, config.extra_weight, ...
                 config.genekoflag, config.rxnkoflag, file_prefix, ...
                 config.out_name, config.root_path, metadata);
end

