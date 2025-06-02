function [data_series, prefix_name, prefix_series, root_path, late_stage, early_stage, simulation, ctrl, medium_series, prefix_pattern_func]=CBM_dataInput(input_data_choose)

  % function to get the name of each column for inferred objective functions
  prefix_pattern_func = @(x) disp(x)

  switch input_data_choose

    %% Input data (significant genes) of CCLE datasets
    case 'CFR_BulkRNAseq'
      % data input setting for bulk datasets
      data_dir = '/nfs/turbo/umms-csriram/daweilin/data/BulkRNAseq/sigGenes/BulkRNAseq/';
      prefix_name = 'BulkRNAseq';
      medium = 'DMEMF12';
      [data_series, prefix_series, medium_series] = batch_input_preprocess(data_dir, prefix_name, medium);
      disp(data_series)
      % 
      late_stage = 'upgenes'
      early_stage = 'dwgenes'
      simulation = 'CFR';
      ctrl = 1 %set to 1 if we want to apply constraints
      % path to save CFR results based on humanOrgan datasets
      root_path = '/nfs/turbo/umms-csriram/daweilin/fluxPrediction/BulkRNAseq/top20/';


    case 'DFA_cellCycleValid'
      data_series = {'/nfs/turbo/umms-csriram/daweilin/data/cellCycleValidation/nameMatch_metabolomics.xlsx'};%'/nfs/turbo/umms-csriram/daweilin/data/cellCycleValidation/nameMatch_metabolomics.xlsx'}
      prefix_name = 'cellCycleValid'
    %sheet_names = {'E1C-L1C'; 'E2C-L2C'; 'E1C-L2C'};
    %sample_names = {'E1CL1C'; 'E2CL2C'; 'E1CL2C'};
      prefix_series = {'G1G2'};%'SG2'};
    %sample_names = {'E1CL2C'};
    %root_path = '/home/daweilin/StemCell/Project_mESC_JinZhang/validation/multiobj/DFA/1C2C/Zhang1C2C'
      root_path = '/home/daweilin/StemCell/Project_mESC_JinZhang/validation/scCellCycleValidation/DFA/G1G2/'
      late_stage = 'Late'
      early_stage = 'Early'
      simulation = 'DFA';
      ctrl = 1 %set to 1 if we want to apply constraints
      medium_series = {'DMEMF12'};%'DMEMF12'};

  %% Metabolomics datasets for cell cycle DFA
    case 'CFR_cellCycleValid'
      % data input setting for single-cell human organ datasets
      data_dir = '/nfs/turbo/umms-csriram/daweilin/data/cellCycleValidation/G1_protein/';
      data_dir_files = dir(data_dir);
      data_series= {};
      prefix_series = {};
      for kk=3:length(data_dir_files),
        data_name = strsplit(data_dir_files(kk).name, '_')
        dependence = strsplit(data_name{length(data_name)}, 'genes')
        data_name = join(data_name(1:length(data_name)-1), '_')
        data_series{kk-2, 1} = sprintf('%s/%s%s', data_dir_files(kk).folder, data_name{:}, dependence{length(dependence)});
        
        %data_series{kk-2, 1} = sprintf('%s/%s', data_dir_files(kk).folder, data_dir_files(kk).name);
         
        name_arr = strsplit(data_dir_files(kk).name);
        prefix_name = 'protein_cellCycleValid';
        prefix_str = prefix_name;
        for ll=1:length(name_arr),
          if ll==length(name_arr),
            append_name = strsplit(name_arr{ll}, '.');
            if strcmp(append_name(length(append_name)), 'csv'),
              append_name = join(append_name(1:length(append_name)-1), '_');
              append_name = strsplit(append_name{:}, '_');
              append_name = join(append_name(1:length(append_name)-1), '_');
            else,
              append_name = join(append_name(1:length(append_name)-1), '_');
            end
          else,
            append_name = name_arr{ll};
          end
          prefix_str = sprintf('%s_%s', prefix_str, append_name{:});
        end
        prefix_series{kk-2, 1} = prefix_str;
      end
      data_series = unique(data_series);
      prefix_series = unique(prefix_series);
      medium_series = {};
      medium_series(1:length(prefix_series), 1) = {'DMEMF12'};
      late_stage = 'upgenes'
      early_stage = 'dwgenes'
      simulation = 'CFR';
      ctrl = 0 %set to 1 if we want to apply constraints
      % path to save CFR results based on humanOrgan datasets
      root_path = '/home/daweilin/StemCell/Project_mESC_JinZhang/validation/scCellCycle_multiObj/unconstraint/G1_protein/';

  %% Significant genes of single cell human cell datasets
    case 'CFR_scTcell'
      % data input setting for single-cell human organ datasets
      data_dir_humanOrgan = '/nfs/turbo/umms-csriram/daweilin/data/Tcell/GSE126030/lymphNode-1/sigGenes/GSM3589411/';
      data_dir_files = dir(data_dir_humanOrgan);
      data_series= {};
      prefix_series = {};
      for kk=3:length(data_dir_files),
        data_name = strsplit(data_dir_files(kk).name, '_')
        dependence = strsplit(data_name{length(data_name)}, 'genes')
        data_name = join(data_name(1:length(data_name)-1), '_')
        data_series{kk-2, 1} = sprintf('%s/%s%s', data_dir_files(kk).folder, data_name{:}, dependence{length(dependence)});
        
        %data_series{kk-2, 1} = sprintf('%s/%s', data_dir_files(kk).folder, data_dir_files(kk).name);
         
        name_arr = strsplit(data_dir_files(kk).name);
        prefix_name = 'activated-1-lymphNode';
        prefix_str = prefix_name;
        for ll=1:length(name_arr),
          if ll==length(name_arr),
            append_name = strsplit(name_arr{ll}, '.');
            if strcmp(append_name(length(append_name)), 'csv'),
              append_name = join(append_name(1:length(append_name)-1), '_');
              append_name = strsplit(append_name{:}, '_');
              append_name = join(append_name(1:length(append_name)-1), '_');
            else,
              append_name = join(append_name(1:length(append_name)-1), '_');
            end
          else,
            append_name = name_arr{ll};
          end
          prefix_str = sprintf('%s_%s', prefix_str, append_name{:});
        end
        prefix_series{kk-2, 1} = prefix_str;
      end
      data_series = unique(data_series);
      prefix_series = unique(prefix_series);
      medium_series = {};
      medium_series(1:length(prefix_series), 1) = {'DMEMF12'};
      late_stage = 'upgenes'
      early_stage = 'dwgenes'
      simulation = 'CFR';
      ctrl = 1 %set to 1 if we want to apply constraints
      % path to save CFR results based on humanOrgan datasets
      root_path = '/nfs/turbo/umms-csriram/daweilin/fluxPrediction/scTcell/activated-1-lymphNode/';
   
  %% Significant genes of single cell human cell datasets
    case 'CFR_scEmbryo'
      % data input setting for single-cell human organ datasets

      data_dir = '/nfs/turbo/umms-csriram/daweilin/data/scEmbryo/GSE136714/single_cell/sigGenes/32cell/';
      prefix_name = '32cell';
      medium = 'KSOM';
      [data_series, prefix_series, medium_series] = batch_input_preprocess(data_dir, prefix_name, medium);

      late_stage = 'upgenes'
      early_stage = 'dwgenes'
      simulation = 'CFR';
      ctrl = 1 %set to 1 if we want to apply constraints
      % path to save CFR results based on humanOrgan datasets
      root_path = '/nfs/turbo/umms-csriram/daweilin/fluxPrediction/scEmbryo/32cell/';

  %% Significant genes of single cell human cell datasets
    case 'CFR_PosNeg2C'
      % data input setting for single-cell human organ datasets

      data_dir = '/nfs/turbo/umms-csriram/daweilin/data/ESC_2Cpos_2Cneg/sigGenes/';
      prefix_name = 'PosNeg2C';
      medium = 'DMEMF12';
      [data_series, prefix_series, medium_series] = batch_input_preprocess(data_dir, prefix_name, medium);

      late_stage = 'upgenes'
      early_stage = 'dwgenes'
      simulation = 'CFR';
      ctrl = 1 %set to 1 if we want to apply constraints
      % path to save CFR results based on humanOrgan datasets
      root_path = '/nfs/turbo/umms-csriram/daweilin/fluxPrediction/PosNeg2C/biomassObj/';


  %% Significant genes of single cell human cell datasets
    case 'CFR_naivePluripotency'
      % data input setting for single-cell human organ datasets

      data_dir = '/nfs/turbo/umms-csriram/daweilin/data/GSE107060_naivePluripotency/';
      prefix_name = 'PosNeg2C';
      medium = 'DMEMF12';
      [data_series, prefix_series, medium_series] = batch_input_preprocess(data_dir, prefix_name, medium);

      late_stage = 'upgenes'
      early_stage = 'dwgenes'
      simulation = 'CFR';
      ctrl = 1 %set to 1 if we want to apply constraints
      % path to save CFR results based on humanOrgan datasets
      root_path = '/nfs/turbo/umms-csriram/daweilin/fluxPrediction/PosNeg2C/biomassObj/';

  %% Significant genes of single cell human cell datasets
    case 'CFR_ESC'
      % data input setting for single-cell human organ datasets

      data_dir = '/nfs/turbo/umms-csriram/daweilin/data/Tzelepis2016/sigGenes/Zhang2013/';
      prefix_name = 'ESC';
      medium = 'DMEMF12';
      [data_series, prefix_series, medium_series] = batch_input_preprocess(data_dir, prefix_name, medium);

      late_stage = 'upgenes'
      early_stage = 'dwgenes'
      simulation = 'CFR';
      ctrl = 1 %set to 1 if we want to apply constraints
      % path to save CFR results based on humanOrgan datasets
      root_path = '/nfs/turbo/umms-csriram/daweilin/fluxPrediction/Tzelepis2016/Recon3D_biomassGeneKO/current_version/';
   
  %% Significant genes of single cell human cell datasets
    case 'CFR_scEMT'
      % data input setting for single-cell human organ datasets

      data_dir = '/nfs/turbo/umms-csriram/daweilin/data/scEMT/sigGenes/A549-TGFB1/';
      prefix_name = 'A549-TGFB1';
      medium = 'DMEMF12';
      [data_series, prefix_series, medium_series] = batch_input_preprocess(data_dir, prefix_name, medium);

      late_stage = 'upgenes'
      early_stage = 'dwgenes'
      simulation = 'CFR';
      ctrl = 1 %set to 1 if we want to apply constraints
      % path to save CFR results based on humanOrgan datasets
      root_path = '/nfs/turbo/umms-csriram/daweilin/fluxPrediction/scEMT/A549-TGFB1/CI068_0d7d/';
   
  %% Quiesence vs proliferation
    case 'CFR_Aging'
      data_dir = '/nfs/turbo/umms-csriram/daweilin/data/Aging/GSE53960_yu_bulk_14/sigGenes/';
      prefix_name = 'Aging';
      medium = 'DMEMF12';
      [data_series, prefix_series, medium_series] = batch_input_preprocess(data_dir, prefix_name, medium);

      late_stage = 'upgenes'
      early_stage = 'dwgenes'
      simulation = 'CFR';
      ctrl = 1 %set to 1 if we want to apply constraints
      % path to save CFR results based on humanOrgan datasets
      root_path = '/nfs/turbo/umms-csriram/daweilin/fluxPrediction/Aging/';

  %% Quiesence vs proliferation
    case 'CFR_QuiProlif'
      data_dir = '/nfs/turbo/umms-csriram/daweilin/data/dawei_prolif_qui_datasets/sigGenes/prolif_qui/';
      prefix_name = 'prolif-qui';
      medium = 'DMEMF12';
      [data_series, prefix_series, medium_series] = batch_input_preprocess(data_dir, prefix_name, medium);

      late_stage = 'upgenes'
      early_stage = 'dwgenes'
      simulation = 'CFR';
      ctrl = 1 %set to 1 if we want to apply constraints
      % path to save CFR results based on humanOrgan datasets
      root_path = '/nfs/turbo/umms-csriram/daweilin/fluxPrediction/prolif_qui_Recon2.2/';
   
  %% Significant genes of single cell human cell datasets
    case 'CFR_scCellCycle'
      % input data
      data_dir = '/nfs/turbo/umms-csriram/daweilin/data/cellCycleMetabolism/scMatrix/sigGenes/scCellCycle/';
      prefix_name = 'scCellCycle';
      medium = 'DMEMF12';
      [data_series, prefix_series, medium_series] = batch_input_preprocess(data_dir, prefix_name, medium);
      late_stage = 'upgenes'
      early_stage = 'dwgenes'
      simulation = 'CFR';
      ctrl = 1 %set to 1 if we want to apply constraints
      % path to save CFR results based on humanOrgan datasets
      root_path = '/nfs/turbo/umms-csriram/daweilin/fluxPrediction/scCellCycleHela_biomass/';
      % function to get cell-specific objectives
      prefix_pattern_func = @prefix_pattern_scCellCycleHela

  %% Input data (significant genes) of CCLE datasets
    case 'CFR_CCLE'
      % input data
      data_dir = '/nfs/turbo/umms-csriram/daweilin/data/BulkRNAseq/sigGenes/ccleRand/';
      prefix_name = 'CCLE';
      medium = 'DMEMF12';
      [data_series, prefix_series, medium_series] = batch_input_preprocess(data_dir, prefix_name, medium);

      late_stage = 'upgenes'
      early_stage = 'dwgenes'
      simulation = 'CFR';
      ctrl = 1 %set to 1 if we want to apply constraints
      % path to save CFR results based on humanOrgan datasets
      root_path = '/nfs/turbo/umms-csriram/daweilin/fluxPrediction/BulkRNAseq/CCLERand/k10_r0.1/';

  
  %% Input data (significant genes) of nci60 datasets
    case 'CFR_NCI60'
      % input data
      data_dir = '/nfs/turbo/umms-csriram/daweilin/data/BulkRNAseq/sigGenes/NCI60Top40/';
      prefix_name = 'NCI60';
      medium = 'DMEMF12';
      [data_series, prefix_series, medium_series] = batch_input_preprocess(data_dir, prefix_name, medium);
      
      late_stage = 'upgenes'
      early_stage = 'dwgenes'
      simulation = 'CFR';
      ctrl = 1 % set to 1 if we want to apply constraints
      % path to save CFR results based on humanOrgan datasets
      root_path = '/nfs/turbo/umms-csriram/daweilin/fluxPrediction/BulkRNAseq/NCI60Top40_biomass/';


  %% data input setting for human protein atlas
    case 'CFR_HumanHPA'
      % input data
      data_dir = '/nfs/turbo/umms-csriram/daweilin/data/BulkRNAseq/sigGenes/HumanHPANeg/';
      prefix_name = 'HumanHPA';
      medium = 'DMEMF12';
      [data_series, prefix_series, medium_series] = batch_input_preprocess(data_dir, prefix_name, medium);

      late_stage = 'upgenes'
      early_stage = 'dwgenes'
      simulation = 'CFR';
      ctrl = 1 %set to 1 if we want to apply constraints
      % path to save CFR results based on humanOrgan datasets
      root_path = '/nfs/turbo/umms-csriram/daweilin/fluxPrediction/BulkRNAseq/HumanHPANeg/k10_r0.1/';

  %% single-cell 1C2C embryo
    case 'CFR_sc1C2C'

      %data_dir = '/nfs/turbo/umms-csriram/daweilin/data/scEmbryo/GSE136714/old_sigGenes/1C2C/';
      %prefix_name = 'sc1C2C';
      %medium = 'KSOM';
      %[data_series, prefix_series, medium_series] = batch_input_preprocess(data_dir, prefix_name, medium);

      data_dir_sc1C2C = '/nfs/turbo/umms-csriram/daweilin/data/scEmbryo/GSE136714/old_sigGenes/1C2C/';
      data_dir_files = dir(data_dir_sc1C2C);
      data_series = {};
      prefix_series = {};
      for kk=3:length(data_dir_files),
        data_series{kk-2, 1} = sprintf('%s/%s', data_dir_files(kk).folder, data_dir_files(kk).name);
        name_arr = strsplit(data_dir_files(kk).name);
        prefix_name = 'sc1C2C';
        prefix_str = prefix_name;
        for ll=1:length(name_arr),
          if ll==length(name_arr),
            append_name = strsplit(name_arr{ll}, '.');
            append_name = append_name{1};
          else,
            append_name = name_arr{ll};
          end
          prefix_str = sprintf('%s_%s', prefix_str, append_name);
        end
        prefix_series{kk-2, 1} = prefix_str;
        medium_series{kk-2, 1} = 'KSOM';
      end

      root_path = '/nfs/turbo/umms-csriram/daweilin/fluxPrediction/sc1C2C_multiObj/norm_nonorm/'
  
      late_stage = '2C'
      early_stage = '1C'
      simulation = 'CFR';
      ctrl = 0 %set to 1 if we want to apply constraints
      
      % function to get the name of each column for inferred objective functions
      prefix_pattern_func = @prefix_pattern_scEmbryo

  %% single-cell 2CBC embryo
    case 'CFR_sc2CBC'
      data_dir = '/nfs/turbo/umms-csriram/daweilin/data/scEmbryo/GSE136714/old_sigGenes/2CBC/';
      %prefix_name = 'sc2CBC';
      %medium = 'KSOM';
      %[data_series, prefix_series, medium_series] = batch_input_preprocess(data_dir, prefix_name, medium);
      data_dir_files = dir(data_dir);
      data_series = {};
      prefix_series = {};
      for kk=3:length(data_dir_files),
        data_series{kk-2, 1} = sprintf('%s/%s', data_dir_files(kk).folder, data_dir_files(kk).name);
        name_arr = strsplit(data_dir_files(kk).name);
        prefix_name = 'sc2CBC';
        prefix_str = prefix_name;
        for ll=1:length(name_arr),
          if ll==length(name_arr),
            append_name = strsplit(name_arr{ll}, '.');
            append_name = append_name{1};
          else,
            append_name = name_arr{ll};
          end
          prefix_str = sprintf('%s_%s', prefix_str, append_name);
        end
        prefix_series{kk-2, 1} = prefix_str;
        medium_series{kk-2, 1} = 'KSOM';
      end
      root_path = '/nfs/turbo/umms-csriram/daweilin/fluxPrediction/sc2CBC_multiObj/norm_nonorm/'
      
      late_stage = 'BC'
      early_stage = '2C'
      simulation = 'CFR';
      ctrl = 0 %set to 1 if we want to apply constraints

      % function to get the name of each column for inferred objective functions
      prefix_pattern_func = @prefix_pattern_scEmbryo
  
  %% Metabolomics datasets for DFA
  %% 1C2C & ESC2CL
    
  %% 1C2C embryo
    case 'CFR_1C2C'
      data_series = {'/nfs/turbo/umms-csriram/daweilin/data/Jin2021_proteomics/Israel1C2C_proteomics.xlsx';'/nfs/turbo/umms-csriram/daweilin/data/Jin2021_transcriptomics/Gao_1C_2C.xlsx';'/nfs/turbo/umms-csriram/daweilin/data/Jin2021_transcriptomics/Yu_1C_2C.xlsx'}
      prefix_series = {'1C2C_Israel';'1C2C_Gao';'1C2C_Yu'}
      root_path = '/nfs/turbo/umms-csriram/daweilin/fluxPrediction/1C2C/CFR/paraScan/'
      prefix_name = '1C2C'
      late_stage = '2C'
      early_stage = '1C'
      simulation = 'CFR';
      ctrl = 1 %set to 1 if we want to apply constraints
      for kk=1:length(data_series),
        medium_series{kk,1} = 'KSOM';
      end
    %root_path = '/home/daweilin/stemcell/project_mesc_jinzhang/validation/multiobj/model'
    %root_path = '/home/daweilin/StemCell/Project_mESC_JinZhang/validation/cancer_model/model/KSOM/'
    %root_path = '/home/daweilin/StemCell/Project_mESC_JinZhang/validation/media_perturbation/CFR/1C2C/Transcriptomics/Yu1C2C';
    %root_path = '/home/daweilin/StemCell/Project_mESC_JinZhang/validation/multiobj/CFR/2CBC/Sharpley2CBC'
  %% 2CBC embryo
    case 'CFR_2CBC'
      data_series = {'/nfs/turbo/umms-csriram/daweilin/data/Jin2021_proteomics/Israel2CBC_proteomics.xlsx';'/nfs/turbo/umms-csriram/daweilin/data/Jin2021_transcriptomics/Gao_2C_BC.xlsx';'/nfs/turbo/umms-csriram/daweilin/data/Jin2021_transcriptomics/Zhang_2C_BC.xlsx';'/nfs/turbo/umms-csriram/daweilin/data/Jin2021_transcriptomics/Sharpley_2C_BC.xlsx'}
      prefix_series = {'2CBC_Israel';'2CBC_Gao';'2CBC_Zhang';'2CBC_Sharpley'}
      root_path = '/nfs/turbo/umms-csriram/daweilin/fluxPrediction/2CBC/CFR/paraScan/'
      prefix_name = '2CBC'
      late_stage = 'BC'
      early_stage = '2C'
      simulation = 'CFR';
      ctrl = 1 %set to 1 if we want to apply constraints
      for kk=1:length(data_series),
        medium_series{kk,1} = 'KSOM';
      end
  
  %% Metabolomics datasets for DFA
  %% 1C2C & ESC2CL
    case 'DFA_1C2C'
      data_series = {'/nfs/turbo/umms-csriram/daweilin/data/Jin2021_metabolomics/Jin_Time_Course_Metabolomics.xlsx';'/nfs/turbo/umms-csriram/daweilin/data/Jin2021_metabolomics/Jin_Metabolomics_DFA.xlsx'}
      prefix_name = '1C2C'
    %sheet_names = {'E1C-L1C'; 'E2C-L2C'; 'E1C-L2C'};
    %sample_names = {'E1CL1C'; 'E2CL2C'; 'E1CL2C'};
      prefix_series = {'E1C-L2C';'ES-ES'};
    %sample_names = {'E1CL2C'};
    %root_path = '/home/daweilin/StemCell/Project_mESC_JinZhang/validation/multiobj/DFA/1C2C/Zhang1C2C'
      root_path = '/nfs/turbo/umms-csriram/daweilin/fluxPrediction/1C2C/DFA/paraScan/'
      late_stage = '2C'
      early_stage = '1C'
      simulation = 'DFA';
      ctrl = 1 %set to 1 if we want to apply constraints
      medium_series = {'KSOM';'DMEMF12'};

  %% 2CBC DFA
    case 'DFA_2CBC'
      data_series = {'/nfs/turbo/umms-csriram/daweilin/data/Jin2021_metabolomics/Jin_Metabolomics_DFA.xlsx';'/nfs/turbo/umms-csriram/daweilin/data/Jin2021_metabolomics/GSE159484_HPLC_Metabolomics_DFA.xlsx'};
      prefix_name = '2CBC'
      prefix_series = {'ES-2C';'2C_M_B'};
    %sample_names = {'EM'};
    % root_path = '/home/daweilin/StemCell/Project_mESC_JinZhang/validation/multiobj/DFA/2CBC/Zhang2CBC'
      root_path = '/nfs/turbo/umms-csriram/daweilin/fluxPrediction/2CBC/DFA/paraScan/'
      late_stage = 'BC'
      early_stage = '2C'
      simulation = 'DFA'
      ctrl = 1 %set to 1 if we want to apply constraints
      medium_series = {'KSOM';'KSOM'};
  % root_path = '/home/daweilin/StemCell/Project_mESC_JinZhang/validation/media_perturbation/DFA/1C2C/ZhangESC2CL'
    otherwise % CFR without constraint
      data_series = {'/home/daweilin/StemCell/Project_mESC_JinZhang/Jin2021_proteomics/Israel2CBC_proteomics.xlsx'}%;'/home/daweilin/StemCell/Project_mESC_JinZhang/Jin2021_transcriptomics/Gao_2C_BC.xlsx'}
      prefix_name = 'unconstraint'
      prefix_series = {'2CBC_Israel'}%;'2CBC_Gao'}
      %root_path = '/home/daweilin/StemCell/Project_mESC_JinZhang/validation/cancer_model/allmets_models/'
      root_path = '/nfs/turbo/umms-csriram/daweilin/fluxPrediction/unconstrained_models/pfba/Recon1/KSOM/'

      late_stage = 'BC'
      early_stage = '2C'
      simulation = 'model'
      ctrl = 0; % remove constraint
      medium_series = {'KSOM'}%;'KSOM'};
  end
end



function prefix_pattern_name=prefix_pattern_scEmbryo(prefix)
  prefix_pattern_name = strsplit(prefix, '_');
  prefix_pattern_name = join(prefix_pattern_name(2:end), '_');
end


function prefix_pattern_name=prefix_pattern_scCellCycleHela(prefix)
  prefix_pattern_name = strsplit(prefix, '_');
  prefix_pattern_name = join(prefix_pattern_name(4:end), '_');
end
