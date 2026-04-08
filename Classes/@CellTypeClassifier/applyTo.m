function [target_labels, target_ctc] = applyTo(ctc, target_rg, opts)
% APPLYTO  Apply a trained CellTypeClassifier to a different RecordingGroup.
%
% Takes the training labels and stored normalization statistics from this
% CTC and classifies all units in the target RecordingGroup, producing a
% new CellTypeClassifier for the target.
%
% Key design guarantees:
%   - Normalization reuses ctc.NormStats exactly (same mu, sigma, scale,
%     nan_cols as used in generateTrainLabels and classifyUnits) so the
%     target data lands in the same feature space as the training data.
%   - Feature group weights (ACGWeight, WaveformWeight) are applied
%     identically to both train and target via applyFeatureWeights.
%   - Within-group z-score is applied independently to source train units
%     and target units using their respective NormalizationVar groups,
%     before the shared global normalization is applied.
%
% The trained CTC must have TrainLabels and NormStats set (i.e., you must
% have run generateTrainLabels() before calling applyTo()).
%
% USAGE:
%   ctc_source = CellTypeClassifier(rg_source, params);
%   ctc_source.identifyResponsiveUnits();
%   ctc_source.generateTrainLabels();
%   [labels, ctc_target] = ctc_source.applyTo(rg_target);
%
% INPUTS:
%   ctc       - trained CellTypeClassifier (source), NormStats must be set
%   target_rg - RecordingGroup to classify
%
% NAME-VALUE:
%   NormalizationVar - metadata field for within-group z-score on target
%                      (default: ctc.Parameters.UMAP.NormalizationVar)
%
% OUTPUTS:
%   target_labels - (1 x N_target) double, 1=excitatory, 2=inhibitory
%   target_ctc    - CellTypeClassifier for the target RG with UnitLabels,
%                   Reduction, HarmonizedWaveforms/ACGs populated

arguments
    ctc        CellTypeClassifier
    target_rg  RecordingGroup
    opts.NormalizationVar string = ctc.Parameters.UMAP.NormalizationVar
end

assert(~isempty(ctc.TrainLabels), ...
    'Source CTC has no TrainLabels — run generateTrainLabels() first.');
assert(~isempty(ctc.NormStats), ...
    'Source CTC has no NormStats — run generateTrainLabels() first.');
assert(~isempty(ctc.UnitList), ...
    'Source CTC has no UnitList — run identifyResponsiveUnits() first.');

labels  = ctc.TrainLabels;
p_umap  = ctc.Parameters.UMAP;
ns      = ctc.NormStats;

% ── Create target CTC (inherits all parameters including Harmonization) ───────
target_ctc = CellTypeClassifier(target_rg, ctc.Parameters);

% ── Build target unit list from baseline recordings ───────────────────────────
% For cultures with a single recording, sorted_recs(1) is the only recording.
% This mirrors identifyResponsiveUnits: unit list comes from the first
% (baseline) recording of each culture.
target_units = Unit.empty;
for c = 1:numel(target_rg.Cultures)
    culture = target_rg.Cultures(c);
    if isscalar(culture.Recordings)
        % Single-recording culture (e.g. in vivo) — take directly
        target_units = [target_units, culture.Recordings(1).Units]; %#ok<AGROW>
    else
        % Multi-recording culture — sort by concentration and take baseline
        concs_c = arrayfun(@(r) r.Metadata.Concentration, culture.Recordings);
        [~, si] = sort(concs_c);
        target_units = [target_units, culture.Recordings(si(1)).Units]; %#ok<AGROW>
    end
end
assert(~isempty(target_units), 'Target RecordingGroup has no units.');
target_ctc.UnitList = target_units;

fprintf('Applying trained classifier to %d target units\n', numel(target_units));

% ── Extract and normalize training features ───────────────────────────────────
% Use cached extraction from source CTC — avoids redundant recomputation.
% Apply the same normalization pipeline as generateTrainLabels/classifyUnits.
[wf_train, acg_train, sr_train]          = ctc.getOrExtract(ctc.UnitList);
[X_train_raw, feat_names, ~, ~, fg_train] = buildFeatureMatrix(ctc, wf_train, acg_train, sr_train);
X_train_raw(isnan(X_train_raw)) = 0;

% Step 1: within-group z-score on source units (matching generateTrainLabels)
[iG_train, G_train] = ctc.RecordingGroup.combineMetadataIndices( ...
    ctc.UnitList, p_umap.NormalizationVar);
for g = unique(iG_train)'
    mask = (iG_train == g);
    X_train_raw(mask, :) = normalize(X_train_raw(mask, :));
end

% Step 2: apply stored global z-score (fit during generateTrainLabels)
if ~isempty(ns.mu)
    X_train_all = normalize(X_train_raw, 'center', ns.mu, 'scale', ns.sigma);
else
    X_train_all = X_train_raw;
end

% Step 3: remove stored NaN columns
X_train_all(:, ns.nan_cols) = [];
feat_names_trimmed = feat_names(~ns.nan_cols);

% Step 4: apply stored scale
X_train_all = X_train_all ./ ns.scale;

% Step 5: apply feature group weights
X_train_all = CellTypeClassifier.applyFeatureWeights(X_train_all, ns.feature_groups, p_umap);

% Select training subset
X_train = X_train_all(labels.umap_train_idx, :);
y_train = labels.sorted_y_train;

% ── Extract and normalize target features ─────────────────────────────────────
[wf_target, acg_target, sr_target]                             = extractUnitWaveformsAndACGs(target_ctc, target_units);
[X_target_raw, feat_names_target, aligned_wf, norm_acgs, fg_target] = buildFeatureMatrix(target_ctc, wf_target, acg_target, sr_target);
X_target_raw(isnan(X_target_raw)) = 0;

% Validate feature count matches source
assert(numel(feat_names) == numel(feat_names_target), ...
    ['Feature count mismatch: source has %d features, target has %d. ' ...
     'Check Harmonization parameters are identical.'], ...
    numel(feat_names), numel(feat_names_target));

% Store harmonized data on target CTC for plotCellTypeFeatures
target_ctc.HarmonizedWaveforms = aligned_wf;
target_ctc.HarmonizedACGs      = norm_acgs;
target_ctc.HarmonizedSR        = ctc.Parameters.Harmonization.WaveformTargetSamplingRate;

% Step 1: within-group z-score on target units using target metadata
[iG_target, ~] = target_rg.combineMetadataIndices(target_units, opts.NormalizationVar);
for g = unique(iG_target)
    mask = (iG_target == g);
    X_target_raw(mask, :) = normalize(X_target_raw(mask, :));
end

% Step 2: apply stored global z-score (same stats as training data)
if ~isempty(ns.mu)
    X_target = normalize(X_target_raw, 'center', ns.mu, 'scale', ns.sigma);
else
    X_target = X_target_raw;
end

% Step 3: remove same NaN columns as training
X_target(:, ns.nan_cols) = [];

% Step 4: apply stored scale
X_target = X_target ./ ns.scale;

% Step 5: apply feature group weights (identical to training)
X_target = CellTypeClassifier.applyFeatureWeights(X_target, ns.feature_groups, p_umap);

% ── Supervised UMAP classification ───────────────────────────────────────────
[target_labels, ~, ~, target_reduction, train_reduction] = supervisedUMAP(ctc, ...
    X_train, y_train, feat_names_trimmed, X_target);

% ── Populate target CTC ───────────────────────────────────────────────────────
% Construct full label vector matching target UnitList
target_ctc.UnitLabels      = double(target_labels);
target_ctc.NormStats       = ns;   % propagate for downstream use

% Build TrainLabels struct for target CTC so plotUMAPSanityCheck works
% All target units are test units — none are training units
target_labels_struct.sorted_train_ids = labels.sorted_train_ids;
target_labels_struct.sorted_y_train   = labels.sorted_y_train;
target_labels_struct.umap_train_idx   = false(1, numel(target_units));
target_labels_struct.umap_test_idx    = true(1,  numel(target_units));
target_ctc.TrainLabels     = target_labels_struct;

target_ctc.Reduction.Train = train_reduction;
target_ctc.Reduction.Test  = target_reduction;

% ── Summary ───────────────────────────────────────────────────────────────────
n_exc = sum(target_labels == 1);
n_inh = sum(target_labels == 2);
fprintf('Target classification: %d excitatory, %d inhibitory (%.1f%% inhibitory)\n', ...
    n_exc, n_inh, 100 * n_inh / numel(target_labels));

% Per-recording inhibitory fraction
[rec_idx, ~] = target_rg.combineMetadataIndices(target_units, opts.NormalizationVar);
inh_count = histcounts(rec_idx(target_labels == 2), 1:max(rec_idx)+1);
all_count = histcounts(rec_idx,                     1:max(rec_idx)+1);
ie_per_rec = inh_count ./ max(all_count, 1);
fprintf('Inhibitory fraction per recording: ');
fprintf('%.2f  ', ie_per_rec);
fprintf('\n');
end