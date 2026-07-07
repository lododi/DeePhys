%% Tutorial 3 — Phenotype Analysis
%
% Where Tutorial 2 is about *looking at* features, Tutorial 3 is about
% *using* them to answer a question: does this recording/culture/unit
% belong to genotype A or B? Does firing behaviour predict drug
% concentration? Phenotype analysis takes the FeatureStore built in
% Tutorial 1 and asks whether a metadata label (Mutation, EI_Ratio,
% Concentration, ...) can be predicted from the extracted features,
% using ML models with cross-validation that respects the data hierarchy
% (units nest in recordings, recordings nest in cultures/chips).
%
% Roadmap of this file:
%   1        Build an Experiment from saved RecordingProcessors — the
%            orchestration object all analysis below is called on.
%   2 - 8    Core prediction workflows at each level of the hierarchy:
%            unit, recording, and culture. Includes classification,
%            regression, and UMAP/PCA dimensionality reduction for
%            visualizing separability before/instead of classifying.
%   9        How to retrieve results already stored on the Experiment.
%   10 - 11  Direct API — call Classifier/DimReducer yourself, bypassing
%            Experiment, when you need a custom feature matrix or CV
%            scheme not covered by exp.classify/exp.reduce.
%   12 - 17  A troubleshooting/validation toolbox for when raw accuracy
%            numbers alone aren't trustworthy: batch-effect correction
%            (z-score or ComBat) when chips/days confound the label of
%            interest, feature-group ablation, confusion matrices,
%            per-feature importance, firing-rate filtering, and
%            permutation testing to check accuracy is above chance.
%
% In short: reach for this tutorial once you have a FeatureStore and a
% hypothesis of the form "does X predict Y" — not for exploring what
% features look like (that's Tutorial 2) or for the unsupervised cell
% type pipeline (that's Tutorial 4).
%
% Prerequisites:
%   - Tutorial 1 completed (FeatureStore and RecordingProcessors saved)
%   - DeePhys on the MATLAB path

%% 1  Build Experiment from saved processors
%
% Discover RecordingProcessor.mat files on disk, keep a subset matching
% this dataset's layout (the well/segment filter below is specific to this
% example run — adapt or drop it for your own data), then wrap them in an
% Experiment. Experiment is the object every section below is called on;
% it owns the merged FeatureStore plus a Results struct that accumulates
% everything computed as you go.
%
% Everything below only needs FeatureStore-level data (classify/reduce/
% regress never touch exp.Processors), so FeatureStore.fromProcessorPaths
% + Experiment.fromFeatureStore is used instead of RecordingProcessor.loadMany
% + Experiment.fromProcessors — it skips loading each recording's
% Connectivity/Bursts/raw SpikeData (commonly 10-20x the size of what's
% actually used here), typically an order of magnitude faster for large
% recording counts. If you need exp.Processors itself (e.g. to re-run raw
% analyses afterward), use Experiment.fromProcessors(RecordingProcessor.loadMany(...))
% instead.

root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Maxtwo/phornauer"; %Root path
path_logic = {'C*','*','w*','sorter_output','segment_*','test*','*'}; %Variable parts
ei_path_list = generate_sorting_path_list(root_path, path_logic);
fprintf("Generated %i sorting paths\n",length(ei_path_list))

%%
ei_path_list = string(ei_path_list);  % ensure string array, 1x595
n = numel(ei_path_list);

well_id    = nan(1, n);
segment_id = nan(1, n);

for s = 1:n
    p = ei_path_list(s);

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
good_proc_paths = ei_path_list(keep_idx);
proc_files = fullfile(good_proc_paths,'RecordingProcessor.mat');

fs  = FeatureStore.fromProcessorPaths(proc_files);
exp = Experiment.fromFeatureStore(fs);

% Need the actual processors afterward (e.g. to re-run raw analyses)?
%   procs = RecordingProcessor.loadMany(proc_files);
%   exp   = Experiment.fromProcessors(procs);
%
% Or from a saved FeatureStore directly (skip the processors):
%   exp = Experiment.fromLegacyGroup(rg);   % from old RecordingGroup

%% 2  Unit-level classification
%
% Classifies each unit using cross-validation grouped by recording,
% so units from the same recording are never split across train/test.
% This is the finest-grained level: one training example per unit, which
% gives the most data but also the most within-recording correlation —
% hence the recording-grouped CV so that correlation can't leak into the
% accuracy estimate.
%
% Unit-level X can be orders of magnitude larger than Recording/Culture
% level (thousands of units vs. tens/hundreds of recordings/cultures) —
% if this runs out of memory, the RF's OOB permuted predictor importance
% is usually the culprit (it re-predicts every OOB sample once per
% feature, per tree, per fold). Set ComputeImportance = false and/or
% lower NumTrees first, before trimming FeatureGroups.

opts = struct();

opts.Algorithm     = 'rf';      % 'rf' (random forest) or 'svm'
opts.KFold         = 5;
opts.FeatureGroups = 'all';     % or ["ActivityFeatures","WaveformFeatures"]
opts.CVLevel       = 'recording';  % group CV at recording level

% Memory knobs for large unit-level runs (uncomment if classify() runs
% out of memory — see Classifier.classify for details):
% opts.ComputeImportance = false;  % skip OOB permutation importance (biggest win)
% opts.NumTrees          = 100;    % fewer trees than the default 500
% opts.Surrogate         = 'off';  % only needed if X has missing values

result_unit = exp.classify('Unit', 'Concentration', opts);

% Inspect result — classify() returns a (1xK) ClassificationResult array, one per fold.
summary_unit = ClassificationResult.summarizeFolds(result_unit);
fprintf('Unit-level accuracy: %.2f ± %.2f\n', summary_unit.mean_accuracy, summary_unit.std_accuracy);

%% 3  Unit-level classification with parent ACGs
%
% Parent ACGs (from concatenated recording) improve classification by
% providing a fuller spike history. Pass ParentFeatures to use them.

opts_parent = struct();
opts_parent.Algorithm     = 'rf';
opts_parent.FeatureGroups = 'all';
opts_parent.ParentFeatures = 'ACG';   % prefer Parent_ACG* from FeatureStore

result_parent = exp.classify('Unit', 'EI_Ratio', opts_parent);

%% 4  Unit-level dimensionality reduction
%
% Before trusting a classification accuracy, it helps to see whether the
% label of interest is visually separable in feature space at all. UMAP
% embeds every unit in 2D; coloring by the label gives a quick sanity
% check — clear clusters suggest classification should work, a uniform
% blob suggests it won't (or that batch effects dominate, see §12/12b).

opts_umap = struct();
opts_umap.Method       = 'UMAP';
opts_umap.FeatureGroups = 'all';

exp.reduce('Unit', opts_umap);

reduction = exp.Results.DimReduction.Unit.UMAP;
disp(reduction);

% Scatter plot coloured by Mutation
embedding = reduction.Reduction;   % DimReductionResult stores coords in .Reduction
labels    = string(exp.FeatureStore.UnitTable.EI_Ratio);
if ~isempty(embedding) && size(embedding, 2) >= 2
    figure;
    gscatter(embedding(:,1), embedding(:,2), labels);
    xlabel('UMAP 1'); ylabel('UMAP 2');
    title('Unit UMAP — coloured by EI_Ratio');
end

%% 5  Recording-level classification
%
% Same classify() call as §2, but each row is now one whole recording's
% aggregated features instead of one unit's — coarser-grained, far fewer
% samples, but a natural fit when the label of interest is a property of
% the recording/culture rather than of individual units (e.g. genotype).
% No CVLevel needed here since there's nothing finer than a recording to
% group by.

opts_rec = struct();
opts_rec.Algorithm = 'rf';
opts_rec.KFold     = 5;

result_rec = exp.classify('Recording', 'EI_Ratio', opts_rec);
summary_rec = ClassificationResult.summarizeFolds(result_rec);
fprintf('Recording-level accuracy: %.2f ± %.2f\n', summary_rec.mean_accuracy, summary_rec.std_accuracy);

%% 6  Recording-level dimensionality reduction
%
% The recording-level analogue of §4 — same idea (visualize separability
% before/instead of classifying), demonstrated here with PCA instead of
% UMAP to show that exp.reduce accepts either Method interchangeably.

exp.reduce('Recording', struct('Method', 'PCA'));
pca_result = exp.Results.DimReduction.Recording.PCA;

%% 7  Culture-level classification
%
% Aggregates recordings into culture-level feature vectors (one row per culture).
% Each culture is represented by features at the specified DIV time points,
% concatenated into a wide vector.

opts_cult = struct();
opts_cult.IdentityKeys   = ["ChipID", "EI_Ratio"];
opts_cult.GroupingVar    = 'DIV';
opts_cult.GroupingValues = [12, 17, 25, 32];
opts_cult.Algorithm      = 'rf';
opts_cult.KFold          = 5;

result_cult = exp.classify('Culture', 'EI_Ratio', opts_cult);
summary_cult = ClassificationResult.summarizeFolds(result_cult);
fprintf('Culture-level accuracy: %.2f ± %.2f\n', summary_cult.mean_accuracy, summary_cult.std_accuracy);

%% 8  Culture-level regression
%
% Regress a numeric metadata variable (e.g. drug concentration) from
% culture-level features.

opts_reg = struct();
opts_reg.IdentityKeys   = ["ChipID", "PlatingDate"];
opts_reg.GroupingVar    = 'DIV';
opts_reg.GroupingValues = [7, 14, 21, 28];
opts_reg.KFold          = 5;

result_conc = exp.regress('Culture', 'Concentration', opts_reg);
summary_conc = RegressionResult.summarizeFolds(result_conc);
fprintf('Concentration regression R2: %.2f\n', summary_conc.mean_R2);

%% 9  Access stored results
%
% Every exp.classify/exp.reduce/exp.regress call above also stashed its
% output on exp.Results, keyed by level and target/method — so you don't
% need to keep the local variables around to revisit a result later.

% All classification results
disp(fieldnames(exp.Results.Classification));
disp(fieldnames(exp.Results.DimReduction));

% Classification result fields
r = exp.Results.Classification.Mutation;
disp(r);

%% 10  Direct API — bypass Experiment
%
% exp.classify covers the common cases, but sometimes you need a feature
% matrix or CV grouping it doesn't expose — e.g. a hand-picked subset of
% units, or CV groups that aren't RecordingID/ChipID. Classifier.classify
% is the lower-level function exp.classify calls internally; using it
% directly gives full control at the cost of doing the bookkeeping yourself.

% Extract feature matrix manually
[X, unit_ids] = exp.FeatureStore.unitMatrix('all');
Y = exp.FeatureStore.UnitTable.Mutation;

% Hierarchy-aware CV: units from the same recording stay together
rec_ids  = exp.FeatureStore.UnitTable.RecordingID;
cv_folds = GroupedCV.byGroups(rec_ids, 5, Y);

% Standalone classification
opts_direct = struct('Algorithm', 'rf', 'KFold', 5);
opts_direct.CVGroups = rec_ids;
result_direct = Classifier.classify(X, Y, opts_direct);

%% 11  Direct API — dimensionality reduction
%
% Same idea as §10 but for reduction: call DimReducer.reduce directly on
% a feature matrix you assembled yourself, bypassing exp.reduce.

[X_unit, ~] = exp.FeatureStore.unitMatrix('all');
opts_dr = struct('Method', 'UMAP', 'NDims', 2);
dim_result = DimReducer.reduce(X_unit, opts_dr);
disp(dim_result);

%% 12  Within-group normalization
%
% Normalize features per recording (or per ChipID) before ML to remove
% batch effects. Specify the metadata field to group by.
%
% Singleton groups (recordings with only one unit) cannot be z-scored and
% automatically fall back to pooled statistics across all units with a
% 'NormalizationPipeline:singletonGroup' warning. The resulting features are
% still valid; the warning is informational. If many groups are singletons,
% consider setting NormalizationVar = '' (no group normalization) or merging
% small recordings before assembling the FeatureStore.

opts_norm = struct();
opts_norm.NormalizationVar = 'ChipID';
opts_norm.Algorithm        = 'rf';

result_norm = exp.classify('Unit', 'Mutation', opts_norm);

%% 12b  ComBat batch correction
%
% For multi-chip experiments with systematic technical offsets that go beyond
% mean/variance shifts, ComBat (Johnson et al. 2007) fits per-batch additive
% and multiplicative effects via empirical Bayes shrinkage, then removes them
% while preserving biological signal.
%
% Set NormalizationPipeline = 'combat' and NormalizationVar to the batch
% metadata field (e.g., ChipID). The classification target is automatically
% used as the biological covariate in ComBat's design matrix so that class
% separation is protected during batch correction.
%
% When to prefer ComBat over per-group z-score:
%   - Many chips/batches with different recording conditions (media lots,
%     electrode coating, temperature drift across days)
%   - The feature distributions look batch-stratified in PCA/UMAP before ML
%
% Note: ComBat already centers residuals per batch, so group_zscore before
% ComBat is unnecessary. The pipeline adds global z-score, clipping, and
% max-abs scaling after correction.

opts_combat = struct();
opts_combat.NormalizationVar      = 'ChipID';    % batch variable
opts_combat.NormalizationPipeline = 'combat';    % use combatThenGlobal pipeline
opts_combat.Algorithm             = 'rf';

result_combat = exp.classify('Unit', 'Mutation', opts_combat);
summary_combat = ClassificationResult.summarizeFolds(result_combat);
fprintf('ComBat-corrected accuracy: %.2f\n', summary_combat.mean_accuracy);

%% 13  Feature group selection for classification
%
% Ablation check: restrict FeatureGroups to a subset (here Activity +
% Waveform, dropping e.g. ACG/graph features) to see which categories of
% features actually carry the predictive signal, rather than treating the
% full feature set as a black box.

opts_fg = struct();
opts_fg.FeatureGroups = ["ActivityFeatures", "WaveformFeatures"];
opts_fg.Algorithm     = 'rf';

result_fg = exp.classify('Unit', 'Mutation', opts_fg);
summary_fg = ClassificationResult.summarizeFolds(result_fg);
fprintf('Activity+Waveform accuracy: %.2f\n', summary_fg.mean_accuracy);

%% 14  Confusion matrix visualization
%
% A single accuracy number hides *which* classes get confused for which.
% Pooling predictions across all folds (rather than plotting per-fold)
% gives a more stable picture of error structure with limited data.

% Aggregate predictions across all folds for a single confusion chart.
T_unit = ClassificationResult.results2table(result_unit);
if ~isempty(T_unit)
    figure;
    confusionchart(T_unit.Y_test, T_unit.Y_pred);
    title('Unit classification — confusion matrix (all folds)');
end

%% 15  Feature importance across folds
%
% Aggregates OOB permutation importance across all CV folds (RF only).
% Provides mean ± std per feature — more reliable than single-fold importance.

T_imp = ClassificationResult.summarizeImportance(result_unit);
if ~isempty(T_imp)
    disp(T_imp(1:min(10, height(T_imp)), :));   % top 10 features
end

%% 16  Minimum firing rate filter
%
% Exclude near-silent units (< 0.1 Hz) before classification.
% These units have unreliable ACG, waveform, and ISI features.

opts_filt = struct('Algorithm', 'rf', 'MinFiringRate', 0.1);
result_filt = exp.classify('Unit', 'Mutation', opts_filt);
summary_filt = ClassificationResult.summarizeFolds(result_filt);
fprintf('Accuracy (FR≥0.1 Hz): %.2f ± %.2f\n', summary_filt.mean_accuracy, summary_filt.std_accuracy);

%% 17  Permutation test — significance of classification accuracy
%
% Builds a null distribution by shuffling labels 100 times.
% p-value = fraction of null runs with accuracy >= observed.
% Recommended for confirming that accuracy exceeds chance, especially
% with small sample sizes or mild batch effects.

opts_perm = struct('Algorithm', 'rf', 'KFold', 5);
opts_perm.CVGroups = exp.FeatureStore.UnitTable.RecordingID;
[p_val, null_accs, obs_acc] = Classifier.permutationTest( ...
    exp.FeatureStore.unitMatrix('all'), exp.FeatureStore.UnitTable.Mutation, ...
    opts_perm, 100);
fprintf('Observed accuracy: %.3f,  permutation p = %.4f\n', obs_acc, p_val);
