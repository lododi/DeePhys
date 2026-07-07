%% Tutorial 4 — Cell Type Classification
%
% Classifies neurons as excitatory (1) or inhibitory (2) using a transductive
% graph-based pipeline. All units are embedded together in a single unsupervised
% UMAP; Louvain community detection identifies inhibitory communities from their
% responsive-unit enrichment; label propagation on the UMAP graph assigns cell
% types to every unit without a separate test-set projection.
%
% Two ground-truth strategies are supported, covered in separate parts:
%
%   PART A — Drug-Response Classification (sections 1–10)
%     Bootstrap firing-rate test identifies units that increase firing after
%     stimulus — these are the inhibitory ground-truth set. Louvain communities
%     are evaluated by their fraction of responsive units.
%
%   PART B — Metadata-Based Classification (sections 11–14)
%     Labels read directly from a UnitTable column (e.g. optogenetics tag,
%     genetic marker, patch-clamp ground truth). Skips bootstrap entirely.
%     Both classes have explicit ground truth.
%
% Classification methods:
%   classify()              — top-level entry: routes to single-seed or ensemble
%   classifyUnits()         — single-seed graph label propagation
%   classifyUnitsEnsemble() — majority vote across multiple RNG seeds (default)
%
% Pipeline:
%   identifyResponsiveUnits()       <- flag drug-responsive (inhibitory) candidates
%   [optimizeUnsupervisedUMAP()]    <- optional Phase 1 BayOpt (recommended)
%   generateTrainLabels()           <- Louvain on UMAP graph -> inh/CE training set
%   classify()                      <- label propagation on UMAP graph (transductive)
%
% Prerequisites:
%   - Tutorial 1 completed (FeatureStore and RecordingProcessors saved)
%   - DeePhys on the MATLAB path
%   - Brain Connectivity Toolbox (BCT) on the MATLAB path (community_louvain)


%% ====================================================================
%% PART A — Drug-Response Classification
%% ====================================================================

%% 1  Load data
%
% Load a pre-built FeatureStore and build the UnitData array from saved
% RecordingProcessors. Both are needed by CellTypeClassifier.
%
% Use generate_sorting_path_list to discover processor directories via a
% path pattern, then build the .mat paths from the discovered directories.

root_path  = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Maxtwo/phornauer";
path_logic = {'C*', '*', 'w*', 'sorter_output', 'segment_*', 'test*','*'};

sorting_paths = generate_sorting_path_list(root_path, path_logic);
fprintf('Discovered %d sorting paths\n', numel(sorting_paths));

proc_paths = fullfile(string(sorting_paths), 'RecordingProcessor.mat');
proc_paths = proc_paths(isfile(proc_paths));


save_dir = '/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Maxtwo/phornauer/Chemogenetics/';
fs_file = fullfile(save_dir, 'FeatureStore.mat');

fs = FeatureStore.load(fs_file);

%%
sorting_paths = string(sorting_paths);  % ensure string array, 1x595
n = numel(sorting_paths);

well_id    = nan(1, n);
segment_id = nan(1, n);

for s = 1:n
    p = sorting_paths(s);

    well_tok = regexp(p, 'well(\d+)', 'tokens', 'once');
    seg_tok  = regexp(p, 'segment_(\d+)', 'tokens', 'once');

    if ~isempty(well_tok)
        well_id(s) = str2double(well_tok{1});
    end
    if ~isempty(seg_tok)
        segment_id(s) = str2double(seg_tok{1});
    end
end

keep_idx = find(well_id > 11 & segment_id < 6);

%%
% Only Units (with waveform/ACG data) are needed here, not Connectivity/
% Bursts/raw SpikeData — RecordingProcessor.loadForFeatureStore reads the
% lightweight feature sidecar (written automatically by save()) instead of
% the full processor, typically an order of magnitude less data per file.
good_proc_paths = sorting_paths(keep_idx);
proc_files = fullfile(good_proc_paths,'RecordingProcessor.mat');

ud_cell = cell(1, numel(proc_files));
parfor i = 1:numel(proc_files)
    s = RecordingProcessor.loadForFeatureStore(proc_files(i));
    ud_cell{i} = s.Units;
end
ud = [ud_cell{:}];
clear ud_cell

%% 2  Parameter setup (drug-response)
%
% Parameters below reflect defaults for the typical MEA drug-response
% experiment. A minimal call needs only GroupingVar / GroupingValues.
% Start with defaults; tune only what diagnostics flag.

params = CellTypeClassifier.returnDefaultParams();
%params = struct();

% ── Harmonization ────────────────────────────────────────────────────────
params.Harmonization.ACGBinSize = 0.0001;
params.Harmonization.ACGLag     = 0.1;
params.Harmonization.ACGSource  = 'FullACG';

% ── Bootstrap firing-rate test ───────────────────────────────────────────
%
% GroundTruthMethod controls how inhibitory candidates are identified:
%   'two_window': bootstrap permutation test comparing pre vs post FR.
%   'full_curve': per-unit Spearman rank correlation across dose levels.
%     Falls back to 'two_window' if fewer than MinRecordings are available.
params.Bootstrap.GroundTruthMethod = 'full_curve';
params.Bootstrap.Alpha          = 0.0001;
params.Bootstrap.NIter          = 1000;
params.Bootstrap.Direction      = 'increase';
params.Bootstrap.PreCutout  = [0, 1200];
params.Bootstrap.PostCutout = [6000, 7200];
params.Bootstrap.BinSize    = 20;

% ── Normalization and recording selection ────────────────────────────────
params.UMAP.NormalizationVar = 'RecordingID';
params.UMAP.GroupingVar    = 'Concentration';
%params.UMAP.GroupingValues = 0;

% ── Unsupervised UMAP geometry ──────────────────────────────────────────
params.UMAP.NDims           = 2;
params.UMAP.AutoNNeighbors  = true;
params.UMAP.MinNNeighbors   = 15;
params.UMAP.MinDist         = 0.1;
params.UMAP.Spread          = 1.0;
params.UMAP.ACGWeight       = 1.0;
params.UMAP.WaveformWeight  = 1.0;

% ── Louvain community detection ─────────────────────────────────────────
params.Community.LouvainResolution            = 1.0;
params.Community.InhibitoryCommunityRelThresh = 0.3;
params.Community.EnrichmentFactor             = 1;
params.Community.PuritySigmaThreshold         = 2.5;
params.Community.CommunityFallbackThreshold   = 0.2;
params.Community.LouvainRestarts              = 1;

% ── Classification method ───────────────────────────────────────────────
%   "graph" (default): label propagation on the UMAP NxN graph.
%   "knn": distance-weighted kNN in feature space (sanity check / small data).
params.Classification.Method               = "graph";
params.Classification.GraphMaxIter         = 100;
params.Classification.GraphConvergenceTol  = 1e-4;

% ── Ensemble (on by default) ────────────────────────────────────────────
%   classify() routes to classifyUnitsEnsemble() when Enabled=true.
%   Set Enabled=false for fast iteration during development (~5x faster).
params.Ensemble.Enabled      = true;
params.Ensemble.Seeds        = [42, 1042, 2042, 3042, 4042];
params.Ensemble.MinAgreement = 0.6;

% ── Counterexample selection ─────────────────────────────────────────────
params.OutlierDetection.CounterexampleRatio = 1;
params.OutlierDetection.CounterexampleDistancePercentile  = 95;
params.OutlierDetection.AutoGMMSeparation     = 1.5;
params.OutlierDetection.ContaminationFraction = 0.1;

% ── Reproducibility ─────────────────────────────────────────────────────
params.RNGSeed = 42;

% Optional: data-driven adaptive params (all off by default)
params.UMAP.AutoNNeighbors              = true;   % n_neighbors = max(15, sqrt(N))
params.UMAP.AutoConfidenceK             = true;   % kNN k = max(5, sqrt(N_train))
params.UMAP.FeatureSelection            = true;   % remove low-var / correlated features
params.Bootstrap.UseFDR                 = false;   % BH correction instead of fixed alpha

% ── Diagnostics ─────────────────────────────────────────────────────────
params.Diagnostics.Enable      = true;
% params.Diagnostics.SaveDir     = '/path/to/output';
params.Diagnostics.ShowFigures = true;
params.CultureKeys = ["ChipID", "PlatingDate", "RecordingDate"];

% Construct classifier
ctc = CellTypeClassifier(fs, ud, params);

%% 3  Identify responsive units (drug-response)
%
% Bootstrap permutation test: compares pre-stimulus vs post-stimulus firing
% rate across cultures. Units with a significant rate increase become the
% inhibitory candidate set (positive class).

ctc.identifyResponsiveUnits();

fprintf('Inhibitory candidates: %d / %d total units (%.1f%%)\n', ...
    sum(ctc.ResponsiveUnitIdx), numel(ctc.ResponsiveUnitIdx), ...
    100 * mean(ctc.ResponsiveUnitIdx));

%% 4  Generate training labels (Louvain community detection)
%
% Internally: global normalization -> unsupervised UMAP -> Louvain community
% detection -> inhibitory training set + excitatory counterexamples.

ctc.generateTrainLabels();

tl = ctc.TrainLabels;
fprintf('Train set: %d excitatory, %d inhibitory (Q=%.3f)\n', ...
    sum(tl.sorted_y_train == 1), sum(tl.sorted_y_train == 2), tl.Q_modularity);

if tl.use_community_path
    fprintf('Inhibitory communities: %s\n', num2str(tl.inh_comm_ids));
end

%% 5  Classify all units
%
% classify() is the recommended entry point. It routes to:
%   - classifyUnitsEnsemble() when Ensemble.Enabled = true (default)
%     Runs generateTrainLabels + classifyUnits N times with different seeds,
%     then takes the majority vote.
%   - classifyUnits() when Ensemble.Enabled = false
%     Single-seed graph label propagation on the UMAP graph.
%
% You can also call classifyUnits() or classifyUnitsEnsemble() directly
% if you want to bypass the routing.

ctc.classify();

labels = ctc.UnitLabels;
n_exc  = sum(labels == 1, 'omitnan');
n_inh  = sum(labels == 2, 'omitnan');
n_nan  = sum(isnan(labels));
fprintf('Excitatory: %d  Inhibitory: %d  Unclassified: %d\n', n_exc, n_inh, n_nan);

%% 6  Bayesian optimization of UMAP + community parameters (optional)
%
% Run this when default parameters produce unsatisfying community structure.

% ctc_opt = CellTypeClassifier(fs, ud, params);
% ctc_opt.identifyResponsiveUnits();
% results = ctc_opt.optimizeUnsupervisedUMAP();
% fprintf('Phase 1 BayOpt: best coherence = %.3f\n', -results.bestObjective);
% ctc_opt.generateTrainLabels();
% ctc_opt.classify();

%% 7  Validate training labels
%
% Leave-one-culture-out cross-validation on the training set.

% template_dir = fullfile(tempdir, 'ctc_umap_templates');
% if ~isfolder(template_dir), mkdir(template_dir); end
% ctc.Parameters.UMAP.TemplateDir = template_dir;
% ctc.validateTrainingLabels();

%% 8  Inspect UMAP embedding and classification

unsup = ctc.Reduction.Unsupervised;
labels_unique = ctc.UnitLabels(ctc.NormalizedFeatures.unique_to_rep);
conf_unique   = ctc.UnitConfidence(ctc.NormalizedFeatures.unique_to_rep);

exc_mask = labels_unique == 1;
inh_mask = labels_unique == 2;
unc_mask = isnan(labels_unique);

figure;
subplot(1, 2, 1);
hold on;
if any(unc_mask)
    scatter(unsup(unc_mask,1), unsup(unc_mask,2), 8, [.6 .6 .6], 'filled', ...
        'MarkerFaceAlpha', 0.4, 'DisplayName', 'Unclassified');
end
scatter(unsup(exc_mask,1), unsup(exc_mask,2), 15, [0.1 0.3 0.8], 'filled', ...
    'MarkerFaceAlpha', 0.7, 'DisplayName', 'Excitatory');
scatter(unsup(inh_mask,1), unsup(inh_mask,2), 15, [0.8 0.1 0.1], 'filled', ...
    'MarkerFaceAlpha', 0.7, 'DisplayName', 'Inhibitory');
hold off;
legend('Box', 'off');
title('Classification labels');
xlabel('UMAP 1'); ylabel('UMAP 2');

subplot(1, 2, 2);
msz_exc = max(5, 5 + 40 * conf_unique(exc_mask));
msz_inh = max(5, 5 + 40 * conf_unique(inh_mask));
hold on;
scatter(unsup(exc_mask,1), unsup(exc_mask,2), msz_exc(:), [0.1 0.3 0.8], 'filled', ...
    'MarkerFaceAlpha', 0.7, 'DisplayName', 'Excitatory');
scatter(unsup(inh_mask,1), unsup(inh_mask,2), msz_inh(:), [0.8 0.1 0.1], 'filled', ...
    'MarkerFaceAlpha', 0.7, 'DisplayName', 'Inhibitory');
hold off;
title('Classification confidence (size \propto conf)');
xlabel('UMAP 1'); ylabel('UMAP 2');

%% 9  Inspect harmonized features

wf  = ctc.HarmonizedWaveforms;
acg = ctc.HarmonizedACGs;
sr  = ctc.HarmonizedSR;
fprintf('Waveform : %d x %d at %.0f Hz\n', size(wf,1), size(wf,2), sr);
fprintf('ACG      : %d x %d\n', size(acg,1), size(acg,2));

%% 10  Troubleshooting (drug-response)
%
% -- Problem: fewer than ~10 responsive units per culture --
%   1. Relax Bootstrap.Alpha to 1e-6.
%   2. Switch to Bootstrap.GroundTruthMethod = 'full_curve' for dose-response.
%   3. Use 'metadata' if external labels are available (see Part B).
%
% -- Problem: inhibitory fraction < 10% or > 30% --
%   1. Check diagnosticTrainLabels community fractions.
%   2. Run optimizeUnsupervisedUMAP to improve community structure.
%   3. Increase CounterexampleRatio (e.g. to 2).
%
% -- Problem: unstable results across RNG seeds --
%   1. Use classify() with Ensemble.Enabled = true (default).
%   2. Run ctc.assessStability('NRuns', 5). Target ARI > 0.90.
%   3. Run optimizeUnsupervisedUMAP if ARI < 0.85.
%
% -- Problem: Louvain fallback triggered --
%   1. Run optimizeUnsupervisedUMAP (community coherence is the objective).
%   2. Lower Community.CommunityFallbackThreshold to 0.3.


%% ====================================================================
%% PART B — Metadata-Based Classification
%% ====================================================================
%
% Use this path when ground truth cell type labels are available from an
% external source: patch-clamp recordings, optogenetic tagging, genetic
% markers (e.g. GAD67-GFP), or pure E/I co-culture ratios.
%
% Key differences from drug-response (Part A):
%   - No bootstrap firing-rate test is run.
%   - Both classes (excitatory and inhibitory) have explicit ground truth.
%   - identifyResponsiveUnits() reads labels from a UnitTable column instead
%     of performing statistical tests.
%   - The Louvain community detection still runs for outlier filtering, but
%     counterexamples come from explicit labels, not distance-based selection.

%% 11  Parameter setup (metadata method)
%
% The key parameter is Bootstrap.GroundTruthMethod = 'metadata'.
% You must also specify:
%   LabelField:                column name in FeatureStore.UnitTable
%   ResponsiveClassValue:      value in that column for the "responsive" class
%   CounterexampleClassValue:  value for the counterexample class (optional;
%                              if empty, any non-responsive non-empty value is used)
%
% Example: if UnitTable has a column "CellType" with values "excitatory"
% and "inhibitory", set LabelField = "CellType", ResponsiveClassValue =
% "inhibitory", CounterexampleClassValue = "excitatory".

params_meta = struct();

% Same harmonization settings as Part A
params_meta.Harmonization.ACGBinSize = 0.0005;
params_meta.Harmonization.ACGLag     = 0.1;
params_meta.Harmonization.ACGSource  = 'FullACG';

% Metadata-specific parameters
params_meta.Bootstrap.GroundTruthMethod       = 'metadata';
params_meta.Bootstrap.LabelField              = 'EI_Ratio';       % column in UnitTable
params_meta.Bootstrap.ResponsiveClassValue    = '0:100';           % value -> inhibitory
params_meta.Bootstrap.CounterexampleClassValue = '100:0';          % value -> excitatory

% ResponsiveClassLabel controls which numeric label (1 or 2) the responsive
% class maps to. Default 2 = responsive -> inhibitory (standard convention).
params_meta.TrainLabels.ResponsiveClassLabel = 2;

% UMAP settings (same as Part A, but GroupingVar/GroupingValues may differ)
params_meta.UMAP.NDims           = 5;
params_meta.UMAP.AutoNNeighbors  = true;
params_meta.UMAP.MinNNeighbors   = 15;
params_meta.UMAP.NormalizationVar = 'ChipID';
params_meta.UMAP.GroupingVar     = 'Concentration';
params_meta.UMAP.GroupingValues  = 0;

% Classification
params_meta.Classification.Method = "graph";

% Ensemble
params_meta.Ensemble.Enabled      = true;
params_meta.Ensemble.Seeds        = [42, 1042, 2042, 3042, 4042];
params_meta.Ensemble.MinAgreement = 0.6;

params_meta.RNGSeed = 42;

%% 12  Construct classifier and identify labeled units

ctc_meta = CellTypeClassifier(fs, ud, params_meta);

% identifyResponsiveUnits reads labels from UnitTable — no FR test is run.
ctc_meta.identifyResponsiveUnits();

n_resp = sum(ctc_meta.ResponsiveUnitIdx);
n_ce   = sum(ctc_meta.CounterexampleUnitIdx);
n_total = numel(ctc_meta.ResponsiveUnitIdx);
fprintf('Metadata labels: %d responsive, %d counterexamples, %d unlabeled\n', ...
    n_resp, n_ce, n_total - n_resp - n_ce);

%% 13  Generate training labels and classify

ctc_meta.generateTrainLabels();
ctc_meta.classify();

labels_meta = ctc_meta.UnitLabels;
fprintf('Metadata path: %d excitatory, %d inhibitory, %d unclassified\n', ...
    sum(labels_meta == 1, 'omitnan'), ...
    sum(labels_meta == 2, 'omitnan'), ...
    sum(isnan(labels_meta)));

%% 14  Inspect metadata classification results

% Harmonized features
wf_m  = ctc_meta.HarmonizedWaveforms;
acg_m = ctc_meta.HarmonizedACGs;
sr_m  = ctc_meta.HarmonizedSR;
fprintf('Waveform : %d x %d at %.0f Hz\n', size(wf_m,1), size(wf_m,2), sr_m);
fprintf('ACG      : %d x %d\n', size(acg_m,1), size(acg_m,2));

% Optional visualization
% plotCellTypeFeatures(ctc_meta);
% sortACGsByPeak(ctc_meta.HarmonizedACGs');

% If ground-truth labels are available for all units, evaluate accuracy:
% gt_labels = ... ;   % (1 x N) ground truth: 1=exc, 2=inh
% stats = CellTypeClassifier.evaluateLabels(ctc_meta.UnitLabels, gt_labels);


%% ====================================================================
%% SHARED — Additional tools
%% ====================================================================

%% 15  Stability assessment
%
% Verify that classification is stable across RNG seeds. Target ARI > 0.90.
% Works for both drug-response and metadata-based classifiers.

% stability = ctc.assessStability('NRuns', 5);
% fprintf('Stability: ARI = %.3f +/- %.3f\n', stability.meanARI, stability.stdARI);
