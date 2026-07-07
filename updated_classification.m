% run_cell_type_classification.m
%
% Full cell-type classification and E/I analysis pipeline.
%
% Pipeline:
%   1. Load processors → assemble FeatureStore
%   2. Initialise CellTypeClassifier
%   3. Bootstrap firing-rate test (inhibitory candidates)
%   4. Generate training labels (UMAP + outlier detection)
%   5. Classify all units (supervised UMAP)
%   6. E/I analysis (EIAnalyzer)
%   7. Save results

%% 0 — Setup

deephys_root = '/home/lododi/git/DeePhysNewPhil/DeePhys';
addpath(genpath(deephys_root));

%% Generate sorting_path_list 
root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Maxtwo/phornauer/";
path_logic = {'Chemogenetics','Large*','well0*','sorter_output'};
split_prefix= 'segment_*';

sorting_path_list{1} = generate_sorting_path_list(root_path, path_logic, split_prefix);

path_logic = {'Chemogenetics','Low*','well0*','sorter_output'};
sorting_path_list{2} = generate_sorting_path_list(root_path, path_logic, split_prefix);

path_logic = {'Chemogenetics_2','Week_1','well0*','sorter_output'};
sorting_path_list{3} = generate_sorting_path_list(root_path, path_logic, split_prefix);

path_logic = {'Chemogenetics_2','Week_2','well0*','sorter_output'};
sorting_path_list{4} = generate_sorting_path_list(root_path, path_logic, split_prefix);

full_sorting_path_list = vertcat(sorting_path_list{:});

root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Maxtwo/lododi2/Neuropixel_Dataset/251121/";
path_logic = {'T*','well*'};
split_prefix= 'well*';
target_sorting_path_list = generate_sorting_path_list(root_path, path_logic, split_prefix);
% /net/bs-filesvr02/export/group/hierlemann/intermediate_data/Maxtwo/phornauer/Chemogenetics/Cellexp_Mouse_Dataset/CellExp_Dataset
%% 1 — Load data
save_dir = '/home/lododi/git/DeePhysNewPhil';
% Paths to saved RecordingProcessor .mat files (one per recording)
proc_paths = full_sorting_path_list + "/MEArecording.mat";
target_paths = target_sorting_path_list + "/MEArecording.mat";


% Only FeatureStores + Units (waveform/ACG) are needed below, not
% Connectivity/Bursts/raw SpikeData — FeatureStore.fromProcessorPaths and
% RecordingProcessor.loadForFeatureStore read each recording's lightweight
% feature sidecar instead of the full processor (commonly 10-20x smaller).
% Legacy MEArecording.mat files without a sidecar fall back to a full load
% automatically (same result, just without the speedup).
fs        = FeatureStore.fromProcessorPaths(proc_paths);
target_fs = FeatureStore.fromProcessorPaths(target_paths);
%%
cortex_idx        = fs.MetadataTable.Region == "Cortex" & fs.MetadataTable.Concentration ~= Inf;
cortex_proc_paths = proc_paths(cortex_idx);
cortex_paths      = [cortex_proc_paths(:); target_paths(:)];  % (:) avoids row/column mismatch
cortex_fs         = FeatureStore.fromProcessorPaths(cortex_paths);

% Build UnitData array in FeatureStore row order
ud_cell = cell(1, numel(cortex_paths));
parfor i = 1:numel(cortex_paths)
    s = RecordingProcessor.loadForFeatureStore(cortex_paths(i));
    ud_cell{i} = s.Units;
end
ud = [ud_cell{:}];

fprintf('FeatureStore: %d units across %d recordings\n', ...
    height(cortex_fs.UnitTable), height(cortex_fs.RecordingTable));

%% 2 — Initialise CellTypeClassifier
nan_idx = isnan(cortex_fs.MetadataTable.Concentration);
cortex_fs.MetadataTable.Concentration(nan_idx) = 0;
params = CellTypeClassifier.returnDefaultParams();

% Bootstrap firing-rate test
params.Bootstrap.PreCutout  = [0, 1200];
params.Bootstrap.PostCutout = [6000, 7200];
params.Bootstrap.BinSize    = 20;
params.Bootstrap.NIter      = 1000;
params.Bootstrap.Alpha      = 0.0001;
params.Bootstrap.GroundTruthMethod = "full_curve"; %"two_side" or "full_curve"
params.Bootstrap.Direction  = 'increase';  % 'increase' | 'decrease' | 'both'

% UMAP embedding
params.UMAP.NDims         = 2;
params.UMAP.AutoNNeighbors  = true;
params.UMAP.MinNNeighbors   = 15;
params.UMAP.MinDist         = 0.1;
params.UMAP.Spread          = 1.0;
params.UMAP.ACGWeight       = 1.0;
params.UMAP.WaveformWeight  = 1.0;
params.UMAP.NormalizationVar = 'RecordingID'; %RecordingID / ChipID

% ── Louvain community detection ────────────────────────────────────────────────
params.Community.LouvainResolution            = 1.0;
params.Community.InhibitoryCommunityRelThresh = 0.3;
params.Community.EnrichmentFactor            = 1;
params.Community.PuritySigmaThreshold        = 2.5;
params.Community.CommunityFallbackThreshold  = 0.2;
params.Community.LouvainRestarts             = 1;

% ACG source: prefer full-recording (Parent_ACG*) from FeatureStore
params.Harmonization.ACGSource  = "FullACG";
params.Harmonization.ACGBinSize = 0.001;
params.Harmonization.ACGLag     = 0.1;
params.Harmonization.WaveformPostTrough = 1;

% Outlier detection (dip test + Mahalanobis/GMM — no isolation forest)
% params.OutlierDetection.OutlierAlpha            = 0.01;
% params.OutlierDetection.DipTestAlpha            = 0.05;
% params.OutlierDetection.MaxResponsiveComponents = 3;
params.OutlierDetection.CounterexampleRatio     = 1;
params.OutlierDetection.Domain                = "umap";
params.OutlierDetection.ContaminationFraction = 0.1;
params.OutlierDetection.AutoGMMSeparation     = 1.5;
params.OutlierDetection.NPCAComponents        = 15;
params.OutlierDetection.CounterexampleDistancePercentile  = 95;
params.OutlierDetection.Method = "community";

% Optional: data-driven adaptive params (all off by default)
params.UMAP.AutoNNeighbors              = true;   % n_neighbors = max(15, sqrt(N))
params.UMAP.AutoConfidenceK             = true;   % kNN k = max(5, sqrt(N_train))
params.UMAP.FeatureSelection            = true;   % remove low-var / correlated features
params.Bootstrap.UseFDR                 = false;   % BH correction instead of fixed alpha

%
params.Diagnostics.Enable      = true;    % master toggle (default false)
% params.Diagnostics.SaveDir     = '/path/to/output';   % save PNGs here
params.Diagnostics.ShowFigures = true;    % false suppresses display (batch mode)
params.CultureKeys = ["ChipID", "PlatingDate", "RecordingDate"];

ctc = CellTypeClassifier(cortex_fs, ud, params);

%% 3 — Identify inhibitory candidates

% metadata_filter restricts which cultures contribute candidates
% e.g. only cultures where AAV == 128 (DREADD-inhibitory)
ctc.identifyResponsiveUnits({'AAV', 128});
ctc.discardSegmentUnits();
ctc.classifyUnits();

%%
results = ctc.optimizeUnsupervisedUMAP(); % Not working, probably because of the mixed datasets 

%%
ctc.generateTrainLabels();

tl = ctc.TrainLabels;
fprintf('Train set: %d excitatory, %d inhibitory (Q=%.3f)\n', ...
    sum(tl.sorted_y_train == 1), sum(tl.sorted_y_train == 2), tl.Q_modularity);

% If the community path was used, inspect community structure:
if tl.use_community_path
    fprintf('Inhibitory communities: %s\n', num2str(tl.inh_comm_ids));
end
%%
results = ctc.optimizeHyperparams();
fprintf('Optimized %d variables, best objective: %.4f\n', ...
    results.nVars, results.bestObjective);

% Apply best topology params, then run the rest of the pipeline
ctc.params = parseStructparams(ctc.params, results.bestParams);
ctc.generateTrainLabels();
ctc.classifyUnits();

% Check stability with the optimized params
stability = ctc.assessStability('NRuns', 5);

%% 4 — Generate training labels

ctc.generateTrainLabels();

%% 5 — Classify all units

ctc.classifyUnits();

%% 6 — Inspect classification quality

plotCellTypeFeatures(ctc);
% sortACGsByPeak(ctc);

% I/E fraction per chip
chip_col = string(ctc.FeatureStore.UnitTable.RecordingDate);
chips    = unique(chip_col, 'stable');
ratios = zeros(size(chips));
for c = 1:numel(chips)
    mask = chip_col == chips(c);
    ie   = sum(ctc.UnitLabels(mask) == 2) / sum(mask);
    ratios(c) = ie;
    % fprintf('  ChipID %s: I/E fraction = %.2f\n', chips(c), ie);
end
figure("Color","w"); bar(ratios); xticklabels(chips); xlabel ("Experiment"); 
xticklabels(["Exp1W3", "Exp1W2", "Exp2W1", "Exp2W2", "Rat in vivo"])
ylabel(" Inhibitory ratio"); box off
% Optional: assess label stability across seeds (requires ~5x pipeline time)
%   stability = ctc.assessStability('NRuns', 5);
%   fprintf('Stability: ARI = %.3f +/- %.3f  (stable=%d)\n', ...
%       stability.meanARI, stability.stdARI, stability.isStable);

%% 7 — E/I analysis

eia_params = struct();
eia_params.Activity.BinSize              = 0.01;
eia_params.Activity.SecCutout            = [0, 1200];
eia_params.BurstDetection.Threshold      = 2;
eia_params.BurstDetection.SmoothWindow   = 9;
eia_params.BurstDetection.PeakCutout     = 150;
eia_params.BurstDetection.PeakHeight     = 0.1;
eia_params.BurstDetection.PeakProminence = 0.1;
eia_params.BurstDetection.MinPeakDistance = 6;
eia_params.BurstDetection.SelectedVariant = 4;

eia = EIAnalyzer(ctc, eia_params);
eia.computeActivity();
eia.detectBursts();
eia.extractBurstCutouts();
eia.computeCorrelations();
eia.normalizeBurstCutouts();

%% 8 — Visualise burst dynamics

eia.PlotNetworkActivity(1);

norm_bursts = eia.NormalizedCutouts.bursts;
norm_ie     = eia.NormalizedCutouts.inh_frac;

figure('Color', 'w');
tiledlayout(2, 2, 'TileSpacing', 'compact');
nexttile; imagesc(norm_ie);     title('I/E ratio (normalised)',      'FontWeight', 'normal'); colorbar
nexttile; imagesc(norm_bursts); title('Total activity (normalised)', 'FontWeight', 'normal'); colorbar
nexttile; plot(nanmean(norm_ie), 'r'); hold on; plot(nanmean(norm_bursts), 'k');
          legend('I/E', 'Total'); xlabel('Bin'); box off
nexttile; plot(gradient(nanmean(norm_bursts)));
          xlabel('Bin'); title('d/dt total activity', 'FontWeight', 'normal'); box off

%% 9 — Save

if ~isfolder(save_dir); mkdir(save_dir); end
save_path = fullfile(save_dir, 'cell_type_classification_results');
save(save_path, 'ctc', 'eia');
fprintf('Saved to: %s.mat\n', save_path);
