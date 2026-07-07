% updated_classification_chemo_only.m
%
% Controlled variant of updated_classification.m: identical chemogenetics
% recording selection (same cortex_idx filter), but WITHOUT concatenating
% in the Neuropixels target_files. Exists to isolate one question raised
% while debugging updated_classification.m's cell-type classification
% results: does including the Neuropixels dataset (~1500/8847 units)
% genuinely improve UMAP/Louvain E/I separability, or does it just look
% that way because Louvain can trivially partition on recording modality
% first (batch effect) rather than genuine cell-type structure?
%
% Differences from updated_classification.m:
%   1. cortex_files = cortex_proc_files only -- target_files is loaded
%      for reference but never concatenated in.
%   2. Harmonization.WaveformPreTrough/PostTrough narrowed to [-0.5, +1] ms
%      (was implicitly wider via WaveformPostTrough=1 alone + default
%      PreTrough). This was validated separately during debugging: it
%      raised mean training silhouette from 0.09 to 0.17 and responsive-
%      in-inhibitory-community coherence from 64.5% to 70.9% on this same
%      chemogenetics-only recording set, by keeping the waveform window
%      focused on the trough region instead of picking up noisy tail
%      samples as spuriously "discriminative" features.
%
% Everything else (Bootstrap/UMAP/Community/OutlierDetection params) is
% kept identical to updated_classification.m for a fair, single/double
% -variable comparison. Compare this run's coherence/silhouette against
% updated_classification.m's to answer the question above.
%
% Pipeline:
%   1. Load processors → assemble FeatureStore (chemogenetics only)
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

% Neuropixels target dataset is discovered for reference only -- NOT used
% below. Kept here so this script's path-discovery section stays a direct
% diff against updated_classification.m.
root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Maxtwo/lododi2/Neuropixel_Dataset/251121/";
path_logic = {'T*','well*'};
split_prefix= 'well*';
target_sorting_path_list = generate_sorting_path_list(root_path, path_logic, split_prefix);

%% 1 — Load data (chemogenetics only -- Neuropixels excluded)

save_dir = '/home/lododi/git/DeePhysNewPhil';
proc_files   = full_sorting_path_list + "/MEArecording.mat";
target_files = target_sorting_path_list + "/MEArecording.mat";   % unused below, loaded for parity only

fs = FeatureStore.fromProcessorPaths(proc_files);

cortex_idx        = fs.MetadataTable.Region == "Cortex" & fs.MetadataTable.Concentration ~= Inf;
cortex_proc_files = proc_files(cortex_idx);
cortex_files      = cortex_proc_files(:);   % chemogenetics only -- no target_files concatenated
cortex_fs         = FeatureStore.fromProcessorPaths(cortex_files);

% Build UnitData array in FeatureStore row order
ud_cell = cell(1, numel(cortex_files));
parfor i = 1:numel(cortex_files)
    s = RecordingProcessor.loadForFeatureStore(cortex_files(i));
    ud_cell{i} = s.Units;
end
ud = [ud_cell{:}];

fprintf('FeatureStore (chemogenetics only): %d units across %d recordings\n', ...
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

% Waveform window narrowed to the trough-informative region (validated
% during debugging -- see header comment). updated_classification.m only
% sets WaveformPostTrough=1 with PreTrough left at its wider default;
% here both are pinned explicitly.
params.Harmonization.WaveformPreTrough  = 0.5;
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

% Enables detecting (and by default, automatically recomputing + re-saving)
% Parent_ACG* whose actual stored (BinSize, Lag) doesn't match
% params.Harmonization (ACGBinSize/ACGLag) above -- bin count alone can't
% reliably tell (e.g. Lag=1/BinSize=0.01 and Lag=2/BinSize=0.02 both give
% 201 bins). The recompute only touches each recording's lightweight
% feature sidecar, not the full RecordingProcessor.mat.
%   ctc.attachProcPaths(cortex_files)        % auto-recompute on mismatch (default)
%   ctc.attachProcPaths(cortex_files, false)  % detect + warn only, no disk writes
ctc.attachProcPaths(cortex_files);

%% 3 — Identify inhibitory candidates

% metadata_filter restricts which cultures contribute candidates
% e.g. only cultures where AAV == 128 (DREADD-inhibitory)
ctc.identifyResponsiveUnits({'AAV', 128});
ctc.discardSegmentUnits();
ctc.classifyUnits();

%%
results = ctc.optimizeUnsupervisedUMAP(); % Chemogenetics-only control: compare its coherence
                                           % against updated_classification.m's mixed-dataset run.

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
results = ctc.optimizeHyperparameters();
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
% No "Rat in vivo" entry here -- chemogenetics-only, 4 experiments expected.
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
save_path = fullfile(save_dir, 'cell_type_classification_results_chemo_only');
save(save_path, 'ctc', 'eia');
fprintf('Saved to: %s.mat\n', save_path);
