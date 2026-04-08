function generateTrainLabels(ctc)
% GENERATETRAINLABELS  Build training labels using UMAP embedding and isolation forest.
%
% Projects units into a low-dimensional UMAP space, then uses an isolation
% forest to remove outliers from the responsive-unit candidates. Selects
% counterexamples (non-responsive units far from the responsive cluster centroid)
% to form a balanced training set (class 1 = excitatory, class 2 = inhibitory).
%
% Key design decisions:
%   - Counterexample selection is performed in NORMALIZED FEATURE SPACE
%     (not in UMAP space) so that training labels are independent of UMAP
%     hyperparameters. This is the prerequisite for stable co-optimization.
%   - Isolation forest is run on PCA-reduced feature space (15 components)
%     rather than 2D UMAP coordinates, which is more statistically robust.
%   - Normalization statistics (mu, sigma, scale, nan_cols) are stored in
%     ctc.NormStats and reused by classifyUnits() to ensure train and test
%     data are on identical scales.
%
% Normalisation pipeline:
%   1. Z-score within each NormalizationVar group (e.g. ChipID) — all units
%   2. Global z-score on the UMAP subset — fit here, stored for reuse
%   3. Remove NaN columns, scale by max(abs) — stored for reuse
%
% Requires ctc.ResponsiveUnitIdx to be set (run identifyResponsiveUnits first).
% Sets: ctc.TrainLabels, ctc.UMAP, ctc.NormStats

arguments
    ctc CellTypeClassifier
end

assert(~isempty(ctc.ResponsiveUnitIdx), ...
    'Run identifyResponsiveUnits() before generateTrainLabels()');

p_umap  = ctc.Parameters.UMAP;
p_outlr = ctc.Parameters.OutlierDetection;
rg      = ctc.RecordingGroup;

% ── Seed RNG for reproducibility ─────────────────────────────────────────────
rng(ctc.Parameters.RNGSeed, 'twister');

% ── Extract features via harmonized path ─────────────────────────────────────
[wf, acg, sr] = ctc.getOrExtract(ctc.UnitList);
[X_raw, ~, ~, ~, feature_groups] = buildFeatureMatrix(ctc, wf, acg, sr);

% ── Step 1: Chip-level normalisation (all units) ──────────────────────────────
% Z-score within each NormalizationVar group (e.g. ChipID).
% Applied to ALL units before any subsetting.
X_raw(isnan(X_raw)) = 0;
[iG, G] = rg.combineMetadataIndices(ctc.UnitList, p_umap.NormalizationVar);
g_idx   = unique(iG);
for g = 1:length(g_idx)
    mask = (iG == g_idx(g));
    X_raw(mask, :) = normalize(X_raw(mask, :));
end

% ── Determine training subset ─────────────────────────────────────────────────
if ~isempty(p_umap.TrainingCultureIdx)
    subset_mask = buildCultureUnitMask(ctc, p_umap.TrainingCultureIdx);
else
    subset_mask = true(1, numel(ctc.UnitList));
end

X_subset = X_raw(subset_mask, :);

% ── Step 2: Global z-score on subset — fit and store ─────────────────────────
% These statistics are stored in ctc.NormStats and reused by classifyUnits()
% to ensure train and test data are normalized identically.
if length(G) > 1
    [X_subset, mu_global, sigma_global] = normalize(X_subset);
else
    mu_global    = zeros(1, size(X_subset, 2));
    sigma_global = ones(1,  size(X_subset, 2));
end

% ── Step 3: Remove NaN columns and scale — store for reuse ───────────────────
nan_cols = any(isnan(X_subset), 1);
X_subset(:, nan_cols) = [];
scale = max(abs(X_subset), [], 1);
scale(scale == 0) = 1;
X_subset = X_subset ./ scale;

% ── Step 4: Apply feature group weights ──────────────────────────────────────
% Adjust n_acg and n_wf for removed NaN columns
acg_cols_kept = sum(~nan_cols(1:feature_groups.n_acg));
wf_cols_kept  = sum(~nan_cols(feature_groups.n_acg+1:end));
feature_groups_trimmed.n_acg = acg_cols_kept;
feature_groups_trimmed.n_wf  = wf_cols_kept;

X_subset = ctc.applyFeatureWeights(X_subset, feature_groups_trimmed, p_umap);

% Store all normalization stats for classifyUnits()
ctc.NormStats.mu             = mu_global;
ctc.NormStats.sigma          = sigma_global;
ctc.NormStats.scale          = scale;
ctc.NormStats.nan_cols       = nan_cols;
ctc.NormStats.feature_groups = feature_groups_trimmed;  % after NaN col removal

% ── UMAP embedding (unsupervised, on subset only) ─────────────────────────────
template_file = fullfile(p_umap.TemplateDir, 'ctc_umap_template.mat');
[reduction, umap_model, ~, ~] = run_umap(X_subset, ...
    'n_components',       p_umap.NDims, ...
    'n_neighbors',        p_umap.NNeighbors, ...
    'min_dist',           p_umap.MinDist, ...
    'spread',             p_umap.Spread, ...
    'sgd_tasks',          20, ...
    'method',             'Java', ...
    'verbose',            'none', ...
    'save_template_file', template_file);
ctc.UMAP = umap_model;
ctc.Reduction.Unsupervised = reduction;

% ── Map responsive units to subset-local indices ──────────────────────────────
subset_responsive    = ctc.ResponsiveUnitIdx(subset_mask);
responsive_in_subset = find(subset_responsive);

% ── Isolation forest in PCA-reduced feature space ────────────────────────────
% Running iforest on PCA components rather than 2D UMAP coordinates is more
% statistically robust and decouples outlier detection from UMAP geometry.
candidate_features = X_subset(subset_responsive, :);

% Remove near-zero variance columns before PCA to avoid singular matrix
col_var = var(candidate_features, 0, 1);
candidate_features_pca = candidate_features(:, col_var > 1e-10);

n_pca_components = min([15, size(candidate_features_pca, 1) - 1, size(candidate_features_pca, 2)]);
[~, pca_scores]  = pca(candidate_features_pca, 'NumComponents', n_pca_components);

[~, tf_forest, ~] = iforest(pca_scores, ...
    'ContaminationFraction',     p_outlr.ContaminationFraction, ...
    'NumObservationsPerLearner', p_outlr.NObsPerLearner);

clean_candidate_features = candidate_features(~tf_forest, :);
if isempty(clean_candidate_features)
    warning('CellTypeClassifier:generateTrainLabels', ...
        ['All %i responsive candidates flagged as outliers by iforest ' ...
         '— falling back to all candidates.'], size(candidate_features, 1));
    clean_candidate_features = candidate_features;
    tf_forest = false(size(tf_forest));
end

% ── Geometric consistency filter ─────────────────────────────────────────────
% Remove inhibitory candidates that are geometrically closer to the
% excitatory population centroid than to the inhibitory centroid.
% These are likely false positives from the bootstrap filter — units
% that showed a firing rate increase for network reasons rather than
% direct DREADD activation, but have excitatory-like waveform/ACG features.
non_responsive_mask  = ~ctc.ResponsiveUnitIdx(subset_mask);
exc_centroid         = mean(X_subset(non_responsive_mask, :), 1);
inh_centroid         = mean(clean_candidate_features, 1);

dist_to_inh = pdist2(inh_centroid, clean_candidate_features, 'correlation');
dist_to_exc = pdist2(exc_centroid, clean_candidate_features, 'correlation');

% Keep only candidates geometrically closer to inhibitory centroid
geom_consistent              = dist_to_inh < dist_to_exc;
n_removed_geom               = sum(~geom_consistent);
clean_candidate_features     = clean_candidate_features(geom_consistent, :);

% Update tf_forest to reflect geometric removal for downstream index mapping
non_outlier_idx              = find(~tf_forest);
additionally_removed         = non_outlier_idx(~geom_consistent);
tf_forest(additionally_removed) = true;

if isempty(clean_candidate_features)
    warning('CellTypeClassifier:generateTrainLabels', ...
        'Geometric consistency filter removed all candidates — reverting to pre-filter set.');
    tf_forest(:) = false;
    clean_candidate_features = candidate_features;
    n_removed_geom = 0;
end

fprintf('Geometric consistency filter: removed %d / %d inhibitory candidates\n', ...
        n_removed_geom, sum(~tf_forest) + n_removed_geom);

% ── Counterexample selection in feature space ─────────────────────────────────
% Distances computed using correlation metric (matches UMAP metric) directly
% in normalized feature space. Labels are therefore independent of UMAP
% hyperparameters, which is the prerequisite for stable co-optimization.
centroid = mean(clean_candidate_features, 1);

% Distance from centroid to ALL subset units (for counterexample pool)
dists_all = pdist2(centroid, X_subset, 'correlation');

% Distance from centroid to clean inhibitory candidates (for threshold)
dists_candidates = pdist2(centroid, clean_candidate_features, 'correlation');
dist_threshold   = prctile(dists_candidates, p_outlr.DistancePercentile);

% Sample counterexamples from units far from the inhibitory centroid
far_idx_local      = find(dists_all > dist_threshold);
n_clean_candidates = sum(~tf_forest);
n_counterexamples  = round(p_outlr.CounterexampleRatio * n_clean_candidates);
counterexample_idx_local = randsample( ...
    far_idx_local, min(n_counterexamples, length(far_idx_local)), false);

% ── Map local indices back to global UnitList indices ─────────────────────────
subset_global          = find(subset_mask);
clean_responsive_local = responsive_in_subset(~tf_forest(:)');

in_train_id        = subset_global(clean_responsive_local);
ex_train_id_global = subset_global(counterexample_idx_local);
ex_train_id        = ex_train_id_global(~ismember(ex_train_id_global, in_train_id));

train_ids = [ex_train_id, in_train_id];
y_train   = [ones(1, length(ex_train_id)), 2*ones(1, length(in_train_id))];
[sorted_train_ids, sort_idx] = sort(train_ids, 'ascend');

labels.sorted_train_ids = sorted_train_ids;
labels.sorted_y_train   = y_train(sort_idx);
labels.umap_train_idx   = false(1, length(ctc.UnitList));
labels.umap_train_idx(sorted_train_ids) = true;
labels.umap_test_idx    = ~labels.umap_train_idx;

% Store diagnostic info for plotUMAPSanityCheck
outlier_local = responsive_in_subset(tf_forest(:)');
labels.outlier_global_idx              = subset_global(outlier_local);
labels.excitatory_candidate_global_idx = subset_global(counterexample_idx_local);

ctc.TrainLabels = labels;
fprintf('Training set: %i excitatory, %i inhibitory candidates\n', ...
    sum(y_train == 1), sum(y_train == 2));
end

function mask = buildCultureUnitMask(ctc, culture_indices)
% Build a logical mask over ctc.UnitList selecting units from specified cultures.
    rg   = ctc.RecordingGroup;
    mask = false(1, numel(ctc.UnitList));
    offset = 0;
    for c = 1:length(rg.Cultures)
        n_units_c = numel(rg.Cultures(c).Units);
        if ismember(c, culture_indices)
            mask(offset+1 : offset+n_units_c) = true;
        end
        offset = offset + n_units_c;
    end
end