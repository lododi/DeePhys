function classifyUnits(ctc)
% CLASSIFYUNITS  Classify all units using supervised UMAP projection.
%
% Trains a supervised UMAP on the labelled training set, then projects all
% remaining (test) units into the same embedding space to predict cell type.
%
% Feature extraction uses ctc.Parameters.Harmonization so DeePhys ACGs are
% recomputed at the configured bin size and lag (default: 0.5 ms / 100 ms).
%
% Normalisation matches the main-branch pipeline:
%   1. Z-score within each NormalizationVar group (e.g. ChipID)
%   2. Global z-score: fit on train, apply train stats to test
%   3. Remove NaN columns, scale by max(abs(train))
%
% Requires ctc.TrainLabels to be set (run generateTrainLabels first).
% Sets: ctc.UnitLabels  (1 = excitatory, 2 = inhibitory, NaN = unclassified)

arguments
    ctc CellTypeClassifier
end

assert(~isempty(ctc.TrainLabels), ...
    'Run generateTrainLabels() before classifyUnits()');

labels = ctc.TrainLabels;
p_umap = ctc.Parameters.UMAP;
rg     = ctc.RecordingGroup;

% ── Extract features via harmonized path ─────────────────────────────────────
[wf, acg, sr]                        = ctc.getOrExtract(ctc.UnitList);
[X_raw, feat_names, aligned_wf, norm_acgs, ~] = buildFeatureMatrix(ctc, wf, acg, sr);

% Store harmonized data for downstream inspection / plotting
ctc.HarmonizedWaveforms = aligned_wf;
ctc.HarmonizedACGs      = norm_acgs;
ctc.HarmonizedSR        = ctc.Parameters.Harmonization.WaveformTargetSamplingRate;

% ── Step 1: Chip-level normalisation (matching generateTrainLabels) ───────────
X_raw(isnan(X_raw)) = 0;
[iG, G] = rg.combineMetadataIndices(ctc.UnitList, p_umap.NormalizationVar);
g_idx = unique(iG);
for g = 1:length(g_idx)
    mask = (iG == g_idx(g));
    X_raw(mask, :) = normalize(X_raw(mask, :));
end

% ── Steps 2-4: Apply stored normalization from generateTrainLabels ────────────
% Reusing ctc.NormStats ensures train and test data are on identical scales
% to the unsupervised embedding learned in generateTrainLabels.
% This prevents the normalization mismatch that causes ring artifacts.
assert(~isempty(ctc.NormStats), ...
    'ctc.NormStats is empty — run generateTrainLabels() before classifyUnits().');
ns = ctc.NormStats;

% Apply stored global z-score (fit on subset in generateTrainLabels)
if ~isempty(ns.mu)
    X_all = normalize(X_raw, 'center', ns.mu, 'scale', ns.sigma);
else
    X_all = X_raw;
end

% Remove the same NaN columns identified during label generation
X_all(:, ns.nan_cols) = [];

% Apply stored scale vector
X_all = X_all ./ ns.scale;

% Apply stored feature group weights
X_all = CellTypeClassifier.applyFeatureWeights(X_all, ns.feature_groups, p_umap);

% Split into train and test using stored label indices
train_idx = logical(labels.umap_train_idx);
test_idx  = logical(labels.umap_test_idx);
X_train   = X_all(train_idx, :);
X_test    = X_all(test_idx,  :);

% ── Supervised UMAP classification ───────────────────────────────────────────
[Y_pred, ~, ~, test_reduction, train_reduction] = supervisedUMAP(ctc, ...
    X_train, labels.sorted_y_train, feat_names, X_test);

ctc.Reduction.Train = train_reduction;
ctc.Reduction.Test  = test_reduction;

full_labels = nan(1, length(ctc.UnitList));
full_labels(labels.sorted_train_ids) = labels.sorted_y_train;
full_labels(labels.umap_test_idx)    = Y_pred;
ctc.UnitLabels = full_labels;

n_exc = sum(full_labels == 1, 'omitnan');
n_inh = sum(full_labels == 2, 'omitnan');
fprintf('Classified %i excitatory, %i inhibitory units (%i unclassified)\n', ...
    n_exc, n_inh, sum(isnan(full_labels)));
end
