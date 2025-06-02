function CFRinterface(model_path, pfba, obj, obj_type, obj_c, root_path, data_path, out_name, upsheet, dwsheet, ctrl, kappa, rho, medium, genekoflag, rxnkoflag, media_perturbation, FSflag, model_name, CFR_model, extra_weight, algorithm)

  sample_name = sprintf('%s%s', upsheet, dwsheet);

%% Load metabolic network model
  model = load(model_path); % or model_human_duarte.mat
  fn = fieldnames(model);
  model = getfield(model, fn{1});
  if strcmp(model_name, 'Recon1'),
    model.c(find(contains(model.rxns, 'biomass_objective'))) = 1;
  end

%% manage file names
  % File name with datetime as prefix
  file_prefix = string(datetime('now','TimeZone','local','Format','MMMdyHHmm'));
  % Output file
  filename = sprintf('%s/[%s]%s', root_path, file_prefix, out_name);
  disp(sprintf('%s %s', 'The CFR result has been saved in', filename));
  % Output file
  excelname = filename;
  
%% change culture medium
  if strcmp(medium, '')==0,
    % media = readtable('FINAL_MEDIUM_MAP.xlsx','Sheet','KSOM_AA_Jin');
    media = readtable('FINAL_MEDIUM_MAP_RECON1.xlsx','Sheet',medium);
    EX_mets = media(:,3);
    EX_rxns = media(:,6);
    [len, ~] = size(EX_rxns);
    m2 = model;
    
    for i=1:len,
        ex_rxn_ind = findRxnIDs(model, EX_rxns{i,:});
        if ex_rxn_ind>0;
            m2 = changeRxnBounds(m2, model.rxns{ex_rxn_ind}, EX_mets{i,:}, 'l');
        end
    end
    model = m2;
  end
  
%% Set objectives for multi- or single-objective problems
  model = objective_setting_function(model, obj, obj_c, obj_type, model_name, algorithm);
  
%% Save settings in metadata
  metadata.obj = obj;
  metadata.obj_type = obj_type;
  metadata.obj_c = obj_c;
  metadata.output_path = root_path;
  if ctrl==1,
    metadata.input_path = data_path;
  else,
    metadata.input_path = sample_name;
  end
  metadata.file_name = out_name;
  metadata.with_constraint = ctrl;
  metadata.CFR_kappa = kappa;
  metadata.CFR_rho = rho;
  metadata.medium = medium;
  metadata.genekoflag = genekoflag;
  metadata.rxnkoflag = rxnkoflag;
  metadata.media_perturbation = media_perturbation;
  obj_check = find(model.c~=0); % find objective functions
  if length(obj_check)~=0,
    metadata.objWeights = model.c(obj_check);
    metadata.objRxns = model.rxns{obj_check};
  else
    metadata.objWeights = 0;
    metadata.objRxns = '';
  end
  metadata.model_path = model_path;
  metadata.upStage = upsheet;
  metadata.dwStage = dwsheet;

  % Load reconstructed models
  if length(CFR_model)>0,
    recon_model = load(CFR_model);
    fn = fieldnames(recon_model);
    recon_model = getfield(recon_model, fn{1});
  else,
    recon_model = '';
  end
  metadata.CFRModel = CFR_model;
  metadata.extraWeight = extra_weight;
  metadata.algorithm = algorithm;
  
%% convert structure to json files
  encodedJSON = jsonencode(metadata);
  JSONFILE_name= sprintf('%s_metadata.json', excelname);
  fid = fopen(JSONFILE_name,'w');
  fprintf(fid, encodedJSON);
  fclose('all')
  
  % Separate values and labels
  if ctrl==0,
    ups = {};
    dws = {};
  else,
    if contains(data_path, '.csv'),
      % make file names with .csv as suffix
      data_path = strsplit(data_path, '.csv');
      data_path = data_path{1};
      upfile = sprintf('%s_%s.csv', data_path, upsheet);
      dwfile = sprintf('%s_%s.csv', data_path, dwsheet);
      % read .csv
      ups = readtable(upfile);
      dws = readtable(dwfile);
    elseif contains(data_path, '.mat'),
        % make file names with .mat as suffix
        data_path = strsplit(data_path, '.mat');
        data_path = data_path{1};
        upfile = sprintf('%s_%s', data_path, '.mat');
        ups = load(upfile);
        dws = ups;
    else,
      % read .xlsx
      ups = readtable(data_path,'Sheet', upsheet);
      dws = readtable(data_path,'Sheet', dwsheet);
    end
  end
 
%% manage model parameters
  size_ups = size(ups); size_dws = size(dws);
  size_ups = size_ups(1); size_dws = size_dws(1); 
  if size_ups>0 & size_dws>0,
    % match the gene name from the model
    if strcmp(model_name, 'Recon3D'),

      % convert BiGGIDs to annotations
      BiGG_map_path = './SCOOTI/SCOOTI/metabolicModel/GEMs/Recon3D_genes.json';
      % upregulated genes
      tmp_ups = {};
      ups{:,2} = upper(ups{:,2})
      upgenes = BiGG_to_annotation(ups{:,2}, BiGG_map_path, 'to_BiGG');
      for jj=1:length(upgenes),
        tmp_ups{jj, 1} = jj;
        tmp_ups{jj, 2} = upgenes{jj};
      end

      % downregulated genes
      dws{:,2} = upper(dws{:,2})
      tmp_dws = {};
      dwgenes = BiGG_to_annotation(dws{:,2}, BiGG_map_path, 'to_BiGG');
      for jj=1:length(dwgenes),
        tmp_dws{jj, 1} = jj;
        tmp_dws{jj, 2} = dwgenes{jj};
      end

      % convert them into tables
      ups = cell2table(tmp_ups); dws = cell2table(tmp_dws);

    % Recon2.2
    elseif strcmp(model_name, 'Recon2.2'),

      % convert BiGGIDs to annotations
      BiGG_map_path = './SCOOTI/SCOOTI/metabolicModel/GEMs/Recon2.2_symbol_to_hgnc.json';
      % upregulated genes
      tmp_ups = {};
      ups{:,2} = upper(ups{:,2})
      upgenes = BiGG_to_annotation(ups{:,2}, BiGG_map_path, 'to_BiGG');
      for jj=1:length(upgenes),
        tmp_ups{jj, 1} = jj;
        tmp_ups{jj, 2} = upgenes{jj};
      end
      % downregulated genes
      tmp_dws = {};
      dws{:,2} = upper(dws{:,2})
      dwgenes = BiGG_to_annotation(dws{:,2}, BiGG_map_path, 'to_BiGG');
      for jj=1:length(dwgenes),
        tmp_dws{jj, 1} = jj;
        tmp_dws{jj, 2} = dwgenes{jj};
      end
      ups = cell2table(tmp_ups); dws = cell2table(tmp_dws);

    end
  end

  %% main function for CFR
  function run_simulation(model, excelname),
  
    % Empty container for saving output
    % plus one for WT
    %%%%%%%%%%%
    T = {};
    n = 2;
    % get gene names from the model
    genes = model.genes;
    uk = unique(genes);
    % Create a table for saving results.
    tp1 = "string"; tp2 = "double";
    %cols = ups.Properties.VariableNames;
    sample_name = 'var';
    vn = {}; vn{1} = 'rxns'; vn{2} = sample_name;
%vn(2:n) = cols(2:end);
    tp = [repelem([tp1, tp2], [1 n-1])];
    Size = [length(model.rxns)+1 n];
    T = table('Size',Size,'VariableTypes',tp,'VariableNames',vn);
    rxn_labels = model.rxns;
    rxn_labels{length(model.rxns)+1, 1} = 'Obj';
    T.rxns = rxn_labels;
    
    %%
    % Initialize progress bar
    %progressbarText(0);
    %disp('Initiating constraint-based modeling...')
    ind = 2  
    if ctrl==0,
      dwlist = {};
      uplist = {};
    elseif strcmp(algorithm, 'MOOMIN'),
      %disp('remove gene list')
      dwlist = {};
      uplist = {};
    else,
      % go through all the columns
      %disp(ups)
      up = table2cell(ups(:, ind));
      dw = table2cell(dws(:, ind));
      %disp('debug')
      %disp(up)

      % find upregulation genes that match the genes in the model
      uplist = {};
      count = 0;
      for i=2:length(uk);
          for j=1:length(up);
              if isempty(up{j})==0 && upper(string(up(j)))==upper(string(uk(i)));
                  count=count+1;
                  uplist{count, 1} = uk{i};
              end
          end
      end
      uplist = uplist(~cellfun('isempty',uplist));
      
      % find downregulation genes that match the genes in the model
      dwlist = {};
      count = 0;
      for i=1:length(uk);
          for j=1:length(dw);
              if isempty(dw{j})==0 && upper(string(dw(j)))==upper(string(uk(i)));
                  count=count+1;
                  dwlist{count, 1} = uk{i};
              end
          end
      end
      dwlist = dwlist(~cellfun('isempty',dwlist));
    
    end

    if FSflag==1,
      % flux sampling
      %achr_name = sprintf('[%s]%s', file_prefix, out_name);
      %%achr_name = 'test'
      %[modelSampling, fluxSamples] = sampleCbModel(model, [], 'ACHR');
      %%disp(modelSampling)
      %%disp(size(fluxSamples))
      %%fluxSampleInfo = load('./sampleCbModelTmp.mat');
      %%rows = fluxSampleInfo.modelSampling.rxns;
      %rows = modelSampling.rxns;
      %posArr = zeros(length(rows), 1);
      %for rxn_ind=1:length(rows),
      %  pos = find(strcmp(model.rxns, rows{rxn_ind}));
      %  posArr(rxn_ind, 1) = pos;
      %end
      %%fluxSamples = load(sprintf('%s_10.mat', achr_name)).points;
      %%fluxSamples = chrrSampler(model, 1, 1000, model.c)
      %%sample_cols = ;
      [fluxSamples] = run_CHRR2(model, 100);
      objFlux = optimizeCbModel(model).f;
      for col_num=1:size(fluxSamples, 2),
        % get a column
        gh = zeros(length(model.rxns), 1);
        gh(posArr, 1) = fluxSamples(1:end, col_num);
        % Append objective fluxes.
        gh(end+1, 1) = objFlux;
        % Save result in tb.
        %disp('bugs')
        %disp(sprintf('%d', col_num))
        T.(sprintf('%d', col_num)) = gh;
        %disp('==================================================')
        %progressbarText(ind/n); 
        % File name with datetime as prefixexcelname
        % excelsheet = 'KSOM_Medium';
        excelsheet = 'fluxes';
        subfolder = sprintf('%s/fva_%d/', root_path, col_num);
        % create a subfolder
        mkdir(subfolder);
        % edit path
        excelname = sprintf('%s/fva_%d/[%s]%s', root_path, col_num, file_prefix, out_name);
        % save to a json file
        JSONFILE_name= sprintf('%s_metadata.json', excelname);
        fid = fopen(JSONFILE_name,'w');
        fprintf(fid, encodedJSON);
        fclose('all')
        % save sampled flux
        disp(sprintf('%s %s', 'The CFR result has been saved in', excelname));
        tmp_filename = sprintf('%s_%s.csv', excelname, excelsheet);
        writetable(T, tmp_filename);
      end
    else,

      if strcmp(algorithm, 'INIT'),
        disp('INIT')
        disp('fix model GPR rules...')
        % fix model rule for mapExpressionToReactions
        initial_model = model;
        initial_model.b = model.b(:, 1);
        initial_model = fixModelRules(initial_model);
        initial_model.c = initial_model.c*0
        model = initial_model;
        % run model
        [constrained_model] = run_init2(model, uplist, dwlist, 36, 1)
        [flux_avg, flux_all, rxn_inactive] = run_CHRR2(model, constrained_model, 100);
        for col_num=1:size(flux_all, 2),
          % Reaction fluxes.
          gh = flux_all(1:end, col_num);
          % Append objective fluxes.
          gh(end+1, 1) = 0;
          % Save result in tb.
          T.(sprintf('%s', sample_name)) = gh;
          %disp('==================================================')
          %progressbarText(ind/n);
          
          % make new sub-directories to save sampled fluxes
          excelsheet = 'fluxes';
          subfolder = sprintf('%s/sample_%d/', root_path, col_num);
          % create a subfolder
          mkdir(subfolder);
          % edit path
          excelname = sprintf('%s/sample_%d/[%s]%s', root_path, col_num, file_prefix, out_name);
          % save to a json file
          JSONFILE_name= sprintf('%s_metadata.json', excelname);
          fid = fopen(JSONFILE_name,'w');
          fprintf(fid, encodedJSON);
          fclose('all')
          % save sampled flux
          disp(sprintf('%s %s', 'The CFR result has been saved in', excelname));
          tmp_filename = sprintf('%s_%s.csv', excelname, excelsheet);
          writetable(T, tmp_filename);
        end % end the for-loop of sampling 

      elseif strcmp(algorithm, 'MOOMIN'),
        disp('MOOMIN')
        % optional, do gene ko experiment
        if genekoflag, % [WARNING] there is a bug when it outputs gene names 
          disp('Start computing KO fluxes via MOOMIN...')
          % WT fluxes
          % fix model rule for mapExpressionToReactions
          model3 = model;
          model3.b = model.b(:, 1);
          model3 = fixModelRules(model3);
          [fluxstate, objFlux] = run_moomin(model3, ups);
          geneko_flux(1:length(fluxstate),1) = fluxstate;
          % Append objective fluxes.
          %geneko_flux(end+1, 1) = objFlux;
          unqgenes = unique(model.genes);
          save_genes(1) = {'WT'};
          %kcount = 0;
          for kk = 1:length(unqgenes),
              %try
              %kcount = kcount+1;
              model2 = deleteModelGenes(model,unqgenes(kk));
              % fix model rule for mapExpressionToReactions
              model3 = model2;
              model3.b = model2.b(:, 1);
              model3 = fixModelRules(model3);
              [fluxstate, objFlux] = run_moomin(model3, ups);
              save_genes(kk+1) = unqgenes(kk);
              %fluxstate_ko(1:length(fluxstate),1) = fluxstate;
              %fluxstate_ko(end+1,1) = objFlux;
              geneko_flux(1:length(fluxstate),kk+1) = fluxstate;
              %catch
              %  disp('Could not find the gene:')
              %  disp(unqgenes{kk})
              %  disp(find(contains(model.genes, unqgenes{kk})))
              %end
          end
          % save results of deletion tests
          rxn_labels = model.rxns;
          %rxn_labels{length(model.rxns)+1, 1} = 'Obj';
          excelsheet = 'CFR-geneDel';
          emptyCells = cellfun(@isempty, save_genes);
          emptyCells_ind = find(emptyCells);
          for ind=1:length(emptyCells_ind),
            empind = emptyCells_ind(ind);
            save_genes{empind} = 'unknown';
          end
          %genes{1} = 'rxns';
          %genes{1} = 'WT';
          %genes(2:1+length(save_genes)) = save_genes;
          genes = save_genes;
          gene_deletions(:, 1) = rxn_labels;
          gene_deletions(:, 2:1+length(save_genes)) = num2cell(geneko_flux);
          % export to .mat
          mat_filename = sprintf('%s_%s.mat', excelname, excelsheet);
          save(mat_filename, "rxn_labels", "genes", "geneko_flux");
        else,
          % fix model rule for mapExpressionToReactions
          model2 = model;
          model2.b = model.b(:, 1);
          model2 = fixModelRules(model2);
          [fluxstate, objFlux] = run_moomin(model2, ups);
          % Reaction fluxes.
          gh = fluxstate;
          % Append objective fluxes.
          gh(end+1, 1) = objFlux;
          % Save result in tb.
          T.(sprintf('%s', sample_name)) = gh;
          %disp('==================================================')
          %progressbarText(ind/n);
          
          % File name with datetime as prefixexcelname
          % excelsheet = 'KSOM_Medium';
          excelsheet = 'fluxes';
          disp(sprintf('%s %s', 'The MOOMIN result has been saved in', excelname));
          tmp_filename = sprintf('%s_%s.csv', excelname, excelsheet);
          writetable(T, tmp_filename);
        end % end for MOOMIN
      else, % CFR
        % with constraint
        [uncon_flux, fluxstate, grate_naive, geneko_flux, rxnko_growthrate, solverobj_naive, model_out]=constrain_flux_regulation(model, uplist, dwlist, kappa, rho, 1E-3, 0, genekoflag, rxnkoflag, [], [], recon_model, extra_weight);

        % save context-specific models
        %if length(CFR_model)==0,
        model_filename = sprintf('%s_%s.mat', excelname, 'model_CFR');
        save(model_filename, "model_out");
        %end
        % no constraint (control experiment)
        %[fluxstate_nc,grate_naive_nc,geneko_growthrate_nc, rxnko_growthrate_nc,solverobj_naive_nc]=constrain_flux_regulation(model, uplist, dwlist, kappa, rho, 1E-3, 1, genekoflag, rxnkoflag);
      %disp('---------------------------------------------------')
          
          
        % Reaction fluxes.
        if pfba==1,
          gh = fluxstate;
        else,
          gh = uncon_flux;
        end
        % Append objective fluxes.
        gh(end+1, 1) = grate_naive;
        % Save result in tb.
        T.(sprintf('%s', sample_name)) = gh;
        %disp('==================================================')
        %progressbarText(ind/n);
        
        % File name with datetime as prefixexcelname
        % excelsheet = 'KSOM_Medium';
        excelsheet = 'fluxes';
        disp(sprintf('%s %s', 'The CFR result has been saved in', excelname));
        tmp_filename = sprintf('%s_%s.csv', excelname, excelsheet);
        writetable(T, tmp_filename);
        % writetable(T,excelname,'FileType','spreadsheet','Sheet', excelsheet);
        
        if genekoflag, % [WARNING] there is a bug when it outputs gene names 
          % save results of deletion tests
          rxn_labels = model.rxns;
          %rxn_labels{length(model.rxns)+1, 1} = 'Obj';
          excelsheet = 'CFR-geneDel';
          unigenes = unique(model.genes);
          emptyCells = cellfun(@isempty, unigenes);
          emptyCells_ind = find(emptyCells);
          for ind=1:length(emptyCells_ind),
            empind = emptyCells_ind(ind);
            unigenes{empind} = 'unknown';
          end
          %genes{1} = 'rxns';
          genes{1} = 'WT';
          genes(2:1+length(unigenes)) = unigenes;%reshape(unigenes, [1, length(unigenes)]);
          %genes{1, 3} = 'unknown'; % fix
          %gdrate = geneko_growthrate;
          %gdrate(end+1, 1) = grate_naive;
          gene_deletions(:, 1) = rxn_labels;
          gene_deletions(:, 2:1+length(unigenes)+1) = num2cell(geneko_flux);
          % export to .mat
          mat_filename = sprintf('%s_%s.mat', excelname, excelsheet)
          save(mat_filename, "rxn_labels", "genes", "geneko_flux")
          %T = cell2table(gene_deletions);
          %disp('Debug')
          %emptyCells = cellfun(@isempty, genes);
          %disp(sum(emptyCells))
          %genes{find(emptyCells)} = 'unknown';
          %genes(emptyCells) = {'Test'}
          %genes(1:10)
          %emptyCells = cellfun(@isempty, genes);
          %disp(sum(emptyCells))
          %T = array2table(geneko_flux, 'RowNames', reshape(rxn_labels, [1, length(rxn_labels)]), 'VariableNames',reshape(genes, [1, length(genes)]))
          %T.Properties.VariableNames = genes;
          %disp(sprintf('%s %s', 'The CFR result has been saved in', excelname));
          %tmp_filename = sprintf('%s_%s.csv', excelname, excelsheet)
          %writetable(T, tmp_filename);
          % writetable(gene_deletions,excelname,'FileType','spreadsheet','Sheet', excelsheet);
        end
        % end if genekoflag
        
        %if rxnkoflag,
        %  excelsheet = 'CFR-rxnDel';
        %  rxns = model.rxns;
        %  rxns{end+1, 1} = 'WT';
        %  rdrate = rxnko_growthrate;
        %  rdrate(end+1, 1) = grate_naive;
        %  rxn_deletions = table(rxns, rdrate);
        %  disp(sprintf('%s %s', 'The CFR result has been saved in', excelname));
        %  tmp_filename = sprintf('%s_%s.csv', excelname, excelsheet)
        %  writetable(rxn_deletions, tmp_filename);
        %  % writetable(rxn_deletions,excelname,'FileType','spreadsheet','Sheet', excelsheet);
        %  % save the results of the model without constraints
        %  excelsheet = 'fluxes-nc';
        %  rxns = model.rxns;
        %  rxns{end+1, 1} = 'Obj';
        %  %rfrate = fluxstate_nc;
        %  %rfrate(end+1, 1) = grate_naive_nc;
        %  %rxn_deletions = table(rxns, rfrate);
        %  %disp(sprintf('%s %s', 'The CFR result has been saved in', excelname));
        %  %tmp_filename = sprintf('%s_%s.csv', excelname, excelsheet)
        %  %writetable(T, tmp_filename);
        %  % writetable(T,excelname,'FileType','spreadsheet','Sheet', excelsheet);
        %end
        %% end if rxnkoflag
        %
        %if genekoflag,
        %  % save results of deletion tests
        %  excelsheet = 'CFR-geneDel-nc';
        %  genes = unique(model.genes);
        %  genes{end+1, 1} = 'WT';
        %  gdrate = geneko_growthrate_nc;
        %  gdrate(end+1, 1) = grate_naive_nc;
        %  gene_deletions = table(genes, gdrate);
        %  disp(sprintf('%s %s', 'The CFR result has been saved in', excelname));
        %  tmp_filename = sprintf('%s_%s.csv', excelname, excelsheet)
        %  writetable(gene_deletions, tmp_filename);
        %  % writetable(gene_deletions,excelname,'FileType','spreadsheet','Sheet', excelsheet);
        %end
        %% end if genekoflag
        %
        %if rxnkoflag,
        %  excelsheet = 'CFR-rxnDel-nc';
        %  rxns = model.rxns;
        %  rxns{end+1, 1} = 'WT';
        %  rdrate = rxnko_growthrate_nc;
        %  rdrate(end+1, 1) = grate_naive_nc;
        %  rxn_deletions = table(rxns, rdrate);
        %  disp(sprintf('%s %s', 'The CFR result has been saved in', excelname));
        %  tmp_filename = sprintf('%s_%s.csv', excelname, excelsheet)
        %  writetable(rxn_deletions, tmp_filename);
        %% writetable(rxn_deletions,excelname,'FileType','spreadsheet','Sheet', excelsheet);
        %end
        % end if rxnkoflag
      end % end if INIT or iMAT (CFR)
    end % end if sampling
  end
  % end if run simulation for CFR
  
  % initialize function to run simulations
  fcn = @run_simulation;
  
  % perturbation of extracellular environment
  if media_perturbation==1,
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Perturbation of media condition following with the ideas used in egem algorithm
      model2 = model;
      egem = load('./egem/supplementary_software_code.mat');
      mediareactions1 = egem.mediareactions1;
      tmpname = excelname;
      
      delete(gcp('nocreate'))
      parpool(32);
      parfor i = 1:50
          model = model2;
          for kappatype = 1:2
              if kappatype == 1, kappa  = 10; else kappa = 0.01;end 
          % delete(gcp('nocreate'))
          % parpool(32);
          % parfor i = 1:50
              excelname = tmpname
              kappa1 = kappa;
              if (kappatype == 2) & (ismember(i,[2,3,5:19])) % trace elements
                  kappa1 = kappa/100;
              elseif (kappatype == 1) & (ismember(i,[1;4])) % glucose or glutamine
                  kappa1 = 3;
              end
  
              % new suffix of the file names
              if kappatype==1,
                  excelname = sprintf('%s_%s_%s', excelname, 'EXC', mediareactions1{i});
              else,
                  excelname = sprintf('%s_%s_%s', excelname, 'DEP', mediareactions1{i});
              end
              disp(sprintf('%s %s', 'The DFA result has been saved in', excelname));
  
              % change media..
              [ix pos]  = ismember(mediareactions1(i), model.rxns);
              model.lb(pos) = model.lb(pos)*kappa1;
              % run_simulation(model, excelname);
              feval(fcn, model, excelname);
          end
          
  
          disp(kappatype)
      end
      
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  else,
      % run_simulation(model, excelname)
      feval(fcn, model, excelname);
  
  end
  % end if media_perturbation

end
% end if CFRinterface
