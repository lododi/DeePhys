%% Tutorial 1 — Data Processing Pipeline
%
% Covers the full pipeline from raw Kilosort output to a FeatureStore
% ready for analysis: SpikeData → RecordingProcessor → FeatureStore.
%
% Prerequisites: DeePhys on the MATLAB path (run startup.m from repo root).
%
% Fill in the placeholder paths in Section 1 before running.

%% 1  Setup — fill in your paths

% Path to one Kilosort output directory (must contain spike_times.npy etc.)
ks_path = '/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Maxtwo/phornauer/EI_iNeurons/250212/T002523/Network/well012/sorter_output'; %/links/groups/hierlemann/Projects/ADHD/AnalyzedData/datastore/.cache/70000140/data/well_20/KS-2-251030-0148-70000140-0.01-5-75';

% Optional: path to a parent (concatenated) recording.
% Set to '' if your recording was not split from a longer session.
parent_ks_path = '';

% Where to save the processed RecordingProcessor and FeatureStore files
save_dir = '/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Maxtwo/phornauer/EI_iNeurons/250212/T002523/Network/well012/sorter_output';

%% 2  Metadata struct
%
% One struct per recording. Any scalar fields are stored in the FeatureStore
% and become available for subsetting and ML labels.

metadata = struct();
metadata.ChipID         = 'T002523';
metadata.PlatingDate    = '2024-11-10';
metadata.RecordingDate  = '2025-02-12';  % DIV computed automatically if absent
metadata.DIV            = 90;
metadata.Mutation       = 'WT';
metadata.Concentration  = 0;
metadata.EI_Ratio       = '75:25';

% If split from a parent recording, specify the parent path here.
% If omitted, SpikeData tries to auto-detect a sibling 'qc_output' folder.
if ~isempty(parent_ks_path)
    metadata.ParentInputPath = parent_ks_path;
end

%% 3  Load Kilosort output into SpikeData

sd = SpikeData.fromKilosort(ks_path, metadata);

% Inspect
fprintf('RecordingID : %s\n', sd.RecordingID);
fprintf('Duration    : %.1f s\n', sd.Duration);
fprintf('N spikes    : %d\n', numel(sd.SpikeTimes));
fprintf('N templates : %d\n', size(sd.TemplateWaveforms, 1));
fprintf('Parent path : %s\n', sd.ParentPath);

%% 4  Create RecordingProcessor (no analysis yet)
%
% Pass optional parameter overrides in the second argument (struct).
% All unspecified fields fall back to RecordingProcessor.returnDefaultParams().

proc = RecordingProcessor(sd);
proc.Parameters.QC.Amplitude = [0,1000];
% Alternatively, load directly from disk in one step:
%   proc = RecordingProcessor.fromKilosort(ks_path, metadata);

%% 5  Quality control — filter templates into proc.Units

proc.runQC();

fprintf('Units after QC: %d\n', numel(proc.Units));

% % Inspect individual units
% if ~isempty(proc.Units)
%     u = proc.Units(1);
%     fprintf('  Unit 1 — TemplateID: %d, N spikes: %d, FR: %.2f Hz\n', ...
%         u.TemplateID, numel(u.SpikeTimes), numel(u.SpikeTimes) / sd.Duration);
% end

%% 6  Compute per-unit features
%
% Computes ActivityFeatures, WaveformFeatures, RegularityFeatures, Catch22.
% Results stored in proc.UnitFeatureTable.

proc.computeUnitFeatures();

% Peek at the table
disp(proc.UnitFeatureTable(1:min(5, height(proc.UnitFeatureTable)), :));

%% 7  Compute parent features (full-recording ACGs)
%
% Reads spike data from proc.SpikeData.ParentPath, computes CCG once for ALL
% templates in the parent recording (shared cache across siblings), then
% stores the autocorrelogram diagonal as Parent_ACG1..N columns.
%
% If no parent path is available, a warning is printed and the step is skipped
% silently — child ACG features serve as automatic fallback.

proc.computeParentFeatures();

% Parent ACG columns are now in proc.UnitFeatureTable with 'Parent_' prefix
parent_cols = startsWith(string(proc.UnitFeatureTable.Properties.VariableNames), 'Parent_');
fprintf('Parent feature columns added: %d\n', sum(parent_cols));

%% 8  Compute network features

proc.computeNetworkFeatures();

disp(proc.NetworkFeatureTable);

%% 9  Compute connectivity (CCG / STTC)
%
% Monosynaptic connection windows are now configurable via the CCG sub-struct
% of proc.Parameters.Connectivity. The defaults match the values used in
% previous DeePhys releases; override only when needed for a different
% preparation (e.g. iPSC cultures vs. acute slices):
%
%   proc.Parameters.Connectivity.CCG.PostStart  = 0.0008;  % post-spike excit. window start (s) — default 0.8 ms
%   proc.Parameters.Connectivity.CCG.PostEnd    = 0.004;   % post-spike excit. window end   (s) — default 4.0 ms
%   proc.Parameters.Connectivity.CCG.PreStart   = 0.0032;  % pre-spike baseline window start (s) — default 3.2 ms
%   proc.Parameters.Connectivity.CCG.BonfWindow = 0.005;   % Bonferroni-correction window    (s) — default 5.0 ms

proc.computeConnectivity();

%% 10  Compute cell-type features (E/I-stratified)
%
% Produces E/I-stratified graph, balance, activity, burst, and correlation
% features — columns with suffixes like _EE, _II, _EI, or prefixes like
% ExcitatoryFraction, MeanFiringRate_E/I.
%
% Requires CellTypeLabels to be set (populated by Tutorial 4 after running
% CellTypeClassifier.classify and RecordingProcessor.applyLabelsFromClassifier).
% If CellTypeLabels is empty or all-NaN, this step returns immediately with no
% new columns added.
%
% Feature groups produced (accessible via fs.recordingMatrix("CellTypeBalance") etc.):
%   CellTypeGraph       — graph metrics on E-only and I-only sub-networks
%   CellTypeBalance     — ExcitatoryFraction, MeanFiringRate_E/I, FiringRateRatio_EI
%   CellTypeActivity    — MeanCV2_E/I, MeanFanoFactor_E/I, MeanLvR_E/I
%   CellTypeBurst       — burst lead/participation per cell type
%   CellTypeCorrelation — MeanSTTC_EE/II/EI, STTCSynchronyIndex

proc.computeCellTypeFeatures();
fprintf('CellTypeFeatures status: %s\n', proc.Status.CellTypeFeatures);

%% 11  Compute spatial analysis
%
% Computes spatial distribution features from electrode coordinates stored in
% proc.SpikeData.ElectrodeCoordinates. Skipped with a warning if coordinates
% are absent. Enriched by CellTypeLabels when available.
%
% Sub-analyses run automatically:
%   A. Spatial spread    — ConvexHullArea, ChipCoverage, MeanPairwiseDistance, CentroidSpread
%   B. E/I mixing        — SpatialMixingIndex, MeanFractionExcNN_E/I (requires labels)
%   D. FC decay          — SpatialDecayTau, DistanceFCCorrelation (requires Connectivity)
%   E. FR gradient       — SpatialFRMoransI, CenterPeripheryFR_ratio
%   F. E/I balance map   — SpatialEIVariability, SpatialEIMoransI (requires labels)
%   G. Ripley's K/L      — RipleysL_max, ClusterScale (optional)
%   H. Burst propagation — BurstOriginDispersion, MeanBurstPropagationSpeed (requires Bursts)
%
% All spatial columns are accessible via fs.recordingMatrix("SpatialFeatures")
% or fs.unitMatrix("SpatialFeatures") after FeatureStore assembly.

proc.computeSpatialAnalysis();
fprintf('SpatialFeatures status: %s\n', proc.Status.SpatialFeatures);

% Inspect spatial columns added to the network feature table
if ~isempty(proc.NetworkFeatureTable)
    sp_names = string(proc.NetworkFeatureTable.Properties.VariableNames);
    sp_prefixes = ["ConvexHull","ChipCoverage","MeanPairwise","CentroidSpread", ...
                   "SpatialMixing","MeanFractionExc","SpatialFR","CenterPeriphery", ...
                   "SpatialDecay","DistanceFC","SpatialEI","BurstOrigin", ...
                   "MeanBurstProp","Ripley","ClusterScale"];
    sp_mask = false(size(sp_names));
    for p = sp_prefixes, sp_mask = sp_mask | startsWith(sp_names, p); end
    fprintf('Spatial network features added: %d columns\n', sum(sp_mask));
end

% --- Visualization ---
%
% Two standalone functions are available for spatial inspection:
%
  % plotSpatialUnitMap(proc.SpikeData.ElectrodeCoordinates, proc.Units, ...
  %     proc.CellTypeLabels);
%   % → Scatter map of units at reference electrode positions, colored by cell type.
%   %   Pass [] as third argument to skip cell-type coloring.
%
  % [~, viz_data] = computeSpatialEIBalance( ...
  %     proc.SpikeData.ElectrodeCoordinates, proc.Units, proc.CellTypeLabels);
  % plotSpatialEIBalance(viz_data, proc.SpikeData.ElectrodeCoordinates);
%   % → Kernel-smoothed heatmap of local excitatory fraction across the MEA chip.

%% 12  Shortcut: run all steps in sequence

proc2 = RecordingProcessor(sd);
proc2.Parameters.QC.Amplitude = [0,1000];
proc2.runAll();   % runQC + computeUnitFeatures + computeParentFeatures
                  %       + computeNetworkFeatures + computeConnectivity
                  %       + computeCellTypeFeatures + computeSpatialAnalysis

%% 13  Save and load
save_dir = '/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Maxtwo/phornauer/EI_iNeurons/250212/T002523/Network/well008/sorter_output';
proc_file_1 = fullfile(save_dir, 'RecordingProcessor.mat');
proc.save(proc_file);

proc_loaded = RecordingProcessor.load(proc_file);
fprintf('Loaded %d units from disk.\n', numel(proc_loaded.Units));

% Inspect status flags
disp(proc_loaded.Status);

%% 14  Batch loading with parallel workers
%
% RecordingProcessor.loadMany() loads multiple saved .mat files in parallel
% using parfor. It auto-detects legacy MEArecording format and migrates on
% the fly, so you can point it at a mix of old and new .mat files.
%
% Use generate_sorting_path_list to discover recording directories via a
% path pattern. Each element of path_logic is a glob that matches one
% level of the directory tree under root_path.

root_path  = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Maxtwo/phornauer";
path_logic = {'C*', '*', 'w*', 'sorter_output', 'segment_*', 'test*','*'};

proc_paths = generate_sorting_path_list(root_path, path_logic);
fprintf('Discovered %d sorting paths\n', numel(proc_paths));

% Build .mat paths from the discovered directories
proc_files = fullfile(string(proc_paths), 'RecordingProcessor.mat');
proc_files = proc_files(isfile(proc_files));   % keep only existing files
%%
procs = RecordingProcessor.loadMany(proc_files);
fprintf('Batch-loaded %d processors.\n', numel(procs));

%% 15  Assemble FeatureStore from multiple processors
%
% FeatureStore.fromProcessors concatenates all unit/recording/metadata tables
% from the loaded processors into a single FeatureStore. Processors with
% empty SpikeData or no units are skipped automatically.
%
% This is the right call when you already have procs loaded for other
% reasons (as here). If building the FeatureStore is your ONLY goal, skip
% §14's loadMany entirely and call FeatureStore.fromProcessorPaths(proc_files)
% instead — it reads each recording's lightweight feature sidecar rather
% than the full processor (Connectivity/Bursts/raw SpikeData are commonly
% 10-20x the size of what's actually used here), typically an order of
% magnitude faster:
%
%   fs = FeatureStore.fromProcessorPaths(proc_files);

fs = FeatureStore.fromProcessors(procs_converted);

fprintf('UnitTable     : %d rows x %d cols\n', height(fs.UnitTable),      width(fs.UnitTable));
fprintf('RecordingTable: %d rows x %d cols\n', height(fs.RecordingTable), width(fs.RecordingTable));
fprintf('MetadataTable : %d rows x %d cols\n', height(fs.MetadataTable),  width(fs.MetadataTable));

%% 16  Batch legacy conversion
%
% RecordingProcessor.convertMany() converts legacy MEArecording.mat files
% in parallel and saves each result to disk, keeping memory bounded.
% Returns string array of output paths (empty string for failed conversions).
%
% Use generate_sorting_path_list to discover legacy recordings by path pattern.


legacy_root = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Maxtwo/phornauer/EI_iNeurons/"; %Root path
legacy_logic = {'2*','*0*','Network','w*','sorter_output','qc_output'}; %Variable parts

legacy_dirs = generate_sorting_path_list(legacy_root, legacy_logic);
fprintf('Found %d legacy directories\n',length(legacy_dirs))

%%
legacy_root  = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Maxtwo/phornauer";
legacy_logic = {'C*', '*', 'w*', 'sorter_output', 'segment_*'};
legacy_dirs  = generate_sorting_path_list(legacy_root, legacy_logic);
fprintf('Found %d legacy directories\n', numel(legacy_dirs));

%%
legacy_mats   = fullfile(string(legacy_dirs), 'MEArecording.mat');
legacy_files   = legacy_mats(isfile(legacy_mats));
converted_dirs = fullfile(legacy_dirs,"test_proc");

converted_paths = RecordingProcessor.convertMany(legacy_files, converted_dirs(isfile(legacy_mats)));

n_ok     = sum(converted_paths ~= "");
n_failed = sum(converted_paths == "");
fprintf('Conversion complete: %d succeeded, %d failed.\n', n_ok, n_failed);

% Load all converted processors
good = converted_paths(converted_paths ~= "");
procs_converted = RecordingProcessor.loadMany(good);
fprintf('Converted and loaded %d processors.\n', numel(procs_converted));

%% 17  Chunked FeatureStore assembly for large datasets
%
% Prefer FeatureStore.fromProcessorPaths(proc_files) (§15) — it never loads
% full processors in the first place, so there's usually nothing to chunk.
% Reach for the chunked version below only when you specifically need actual
% RecordingProcessor objects in memory afterward (not just the FeatureStore),
% since loading all of those simultaneously for hundreds of recordings can
% exceed available memory. FeatureStore.fromProcessorsChunked loads them in
% batches, builds a FeatureStore per chunk, saves to disk, then combines the
% lightweight table-only results — same output as FeatureStore.fromProcessors
% (§15), just with bounded peak memory.
%
% Adjust chunk_size based on available RAM (lower = less memory, slower).

chunk_size = 20;
fs_large = FeatureStore.fromProcessorsChunked(proc_files, save_dir, chunk_size);
fs_large.save(fullfile(save_dir, 'FeatureStore.mat'));

%% 18  Partial loading (RAM-friendly)
%
% When you only need feature tables or unit metadata — not the heavy SpikeData
% or connectivity matrices — use loadUnits or loadFeatureTables.
% These use MATLAB's matfile for selective variable access without loading
% the full RecordingProcessor into memory.

proc_file = fullfile(save_dir, 'RecordingProcessor.mat');

% Load only units and status
[units, status] = RecordingProcessor.loadUnits(proc_file);
fprintf('Loaded %d units (status: QC=%s, UnitFeatures=%s)\n', ...
    numel(units), status.QC, status.UnitFeatures);

% Load only feature tables
[uft, nft, status] = RecordingProcessor.loadFeatureTables(proc_file);
fprintf('UnitFeatureTable: %d rows x %d cols\n', height(uft), width(uft));
fprintf('NetworkFeatureTable: %d rows x %d cols\n', height(nft), width(nft));

%% 19  Inspect table structure

disp(string(fs.UnitTable.Properties.VariableNames)');
disp(fs.UnitTable(1:min(5, height(fs.UnitTable)), 1:8));

%% 20  Save and load FeatureStore

save_dir = '/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Maxtwo/phornauer/Chemogenetics/';

fs_file = fullfile(save_dir, 'FeatureStore.mat');
fs.save(fs_file);

fs2 = FeatureStore.load(fs_file);
fprintf('Loaded FeatureStore: %d units, %d recordings.\n', ...
    height(fs2.UnitTable), height(fs2.RecordingTable));

%% 21  Experiment — filterable analysis container
%
% Experiment wraps a FeatureStore with convenient metadata-based filtering
% and analysis methods (classify, regress, reduce). It replaces the old
% RecordingGroup for subsetting and grouping workflows.

% From loaded processors (full pipeline)
exp = Experiment.fromProcessors(procs);

% Or from a saved FeatureStore (lightweight — no spike data in memory)
exp = Experiment.fromFeatureStore(fs);

% Inspect metadata fields and unique values
exp.listMetadata();

%% 22  Filtering and subsetting
%
% filter() and exclude() return a new Experiment with a subset of the data.
% Filters can be chained.

exp_wt  = exp.filter('Mutation', 'WT');
exp_het = exp.filter('Mutation', {'WT','HET'});
exp_d14 = exp.filter('DIV', [14, 21, 28]);
exp_conc = exp.filter('Concentration', [0,0.01,0.1,1,10,100,Inf]).filter('Region','Cortex');

% Exclude specific values
exp_no_ko = exp.exclude('Mutation', 'KO');

% Chain filters
exp_wt_d14 = exp.filter('Mutation', 'WT').filter('DIV', 14);

% Check what's left
fprintf('After filtering: %d recordings, %d units\n', ...
    height(exp_wt_d14.FeatureStore.MetadataTable), ...
    height(exp_wt_d14.FeatureStore.UnitTable));

% Unique values of a field
disp(exp.uniqueValues('Concentration'));
