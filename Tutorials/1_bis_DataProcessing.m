%% Tutorial 1 — Data Processing Pipeline (bulk)
%
% Covers the full pipeline from raw Kilosort output to a FeatureStore
% ready for analysis: SpikeData → RecordingProcessor → FeatureStore.
%
% Prerequisites: DeePhys on the MATLAB path (run startup.m from repo root).
%
% Fill in the placeholder paths in Section 1 before running.

%% 1 Find relevant sorted paths
root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Maxtwo/phornauer/EI_iNeurons/"; %Root path
path_logic = {'2*','*0*','Network','w*','sorter_output','qc_output'}; %Variable parts
sorting_path_list = generate_sorting_path_list(root_path, path_logic);
fprintf("Generated %i sorting paths\n",length(sorting_path_list))

%% 2 Load Metadata and Kilosort output into SpikeData
%
% One struct per recording. Any scalar fields are stored in the FeatureStore
% and become available for subsetting and ML labels.

metadata_filepath = "/links/groups/hierlemann/Projects/lododi/EI_AdvancedSciences_Philipp/Recordings/DeePhys_Original_EI_241119.xlsx";

allSD = generate_spikedata_from_sorting_path_list(sorting_path_list, metadata_filepath);

% % % % If split from a parent recording, specify the parent path here.
% % % % If omitted, SpikeData tries to auto-detect a sibling 'qc_output' folder.
% % % if ~isempty(parent_ks_path)
% % %     metadata.ParentInputPath = parent_ks_path;
% % % end

%% 3  Run all processing steps in sequence
failed = {};
parfor iPath = 1:length(allSD)
    try
        proc = RecordingProcessor(allSD{iPath});
        proc.Parameters.QC.Amplitude = [0, 1000];
        proc.runAll(); % runQC + computeUnitFeatures + computeParentFeatures
                       %       + computeNetworkFeatures + computeConnectivity
                       %       + computeCellTypeFeatures + computeSpatialAnalysis
        proc_file = fullfile(allSD{iPath}.InputPath, 'RecordingProcessor.mat');
        proc.save(proc_file);
    catch ME
        failed = [failed; {allSD{iPath}.InputPath, ME.message, ME.stack(1).name, ME.stack(1).line}];
        warning('Failed: %s\n  %s (in %s line %d)', ...
            allSD{iPath}.InputPath, ME.message, ME.stack(1).name, ME.stack(1).line);
    end
end

% After the loop:
if ~isempty(failed)
    failedT = cell2table(failed, 'VariableNames', {'Path','Error','Function','Line'});
    disp(failedT);
end

%% 13  Save and load

proc_file = fullfile(save_dir, 'RecordingProcessor.mat');
proc.save(proc_file);

proc_loaded = RecordingProcessor.load(proc_file);
fprintf('Loaded %d units from disk.\n', numel(proc_loaded.Units));

% Inspect status flags
disp(proc_loaded.Status);

%% 14  Batch loading with parallel workers
%
% Use the same path pattern to discover saved RecordingProcessor files.

load_root  = root_path;
load_logic = {'2*', '*0*', 'Network', 'w*', 'sorter_output', 'qc_output'};
load_paths = generate_sorting_path_list(load_root, load_logic);

proc_paths = fullfile(string(load_paths), 'RecordingProcessor.mat');
proc_paths = proc_paths(isfile(proc_paths));

% Only building a FeatureStore below, so skip loading full processors
% (Connectivity/Bursts/raw SpikeData) entirely — see §15.

%% 15  Assemble FeatureStore from multiple processors
%
% FeatureStore.fromProcessorPaths reads each recording's lightweight
% feature sidecar (written automatically by RecordingProcessor.save())
% instead of the full processor — typically an order of magnitude less
% data moved, and no chunking needed even for hundreds of recordings.
% Need the actual processors afterward instead? Use
% RecordingProcessor.loadMany(proc_paths) + FeatureStore.fromProcessors(procs).

fs = FeatureStore.fromProcessorPaths(proc_files);

% Three tables
fprintf('UnitTable     : %d rows × %d cols\n', height(fs.UnitTable),      width(fs.UnitTable));
fprintf('RecordingTable: %d rows × %d cols\n', height(fs.RecordingTable), width(fs.RecordingTable));
fprintf('MetadataTable : %d rows × %d cols\n', height(fs.MetadataTable),  width(fs.MetadataTable));

%% 16  Partial loading (RAM-friendly)
%
% Load only feature tables without heavy SpikeData/Connectivity:

[uft, nft, status] = RecordingProcessor.loadFeatureTables( ...
    fullfile(string(load_paths{1}), 'RecordingProcessor.mat'));
fprintf('Loaded feature tables: %d unit rows, %d network cols\n', height(uft), width(nft));

%% 17  Inspect table structure

% Column names in UnitTable
disp(string(fs.UnitTable.Properties.VariableNames)');

% First few rows
disp(fs.UnitTable(1:min(5, height(fs.UnitTable)), 1:8));

%% 18  Save and load FeatureStore

fs_file = fullfile(save_dir, 'FeatureStore.mat');
fs.save(fs_file);

fs2 = FeatureStore.load(fs_file);
fprintf('Loaded FeatureStore: %d units, %d recordings.\n', ...
    height(fs2.UnitTable), height(fs2.RecordingTable));

%% 19  Experiment — filterable analysis container

exp = Experiment.fromFeatureStore(fs);
exp.listMetadata();

% Filter by metadata
exp_subset = exp.filter('Mutation', 'WT').filter('DIV', [14, 21, 28]);
fprintf('Filtered: %d recordings, %d units\n', ...
    height(exp_subset.FeatureStore.MetadataTable), ...
    height(exp_subset.FeatureStore.UnitTable));
