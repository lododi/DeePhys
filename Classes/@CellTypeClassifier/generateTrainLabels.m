function generateTrainLabels(ctc, opts)
% GENERATETRAINLABELS  Build training labels via outlier detection and feature-space operations.
%
% Pipeline:
%   1. Global normalization + optional feature selection  (computeGlobalNormalization)
%   2. Unsupervised UMAP on all unique units              (runUnsupervisedUMAP)
%   3. Outlier detection and CE selection                 (selectExplicitCE |
%                                                          detectCommunityOutliers |
%                                                          selectIforestCE | selectNoCE)
%   4. Assemble training labels                           (assembleTrainLabels)
%
% Two label-source strategies:
%   Drug response  — ground-truth label 1 = inhibitory ground truth; CEs inferred.
%   Metadata       — ctc.GroundTruthLabel2Idx provides explicit excitatory labels.
%
% Requires ctc.GroundTruthLabel1Idx (run identifyGroundTruthUnits first).
% Sets: ctc.TrainLabels, ctc.NormalizationParams, ctc.Reduction.Unsupervised

arguments
    ctc CellTypeClassifier
    opts.ReuseExistingUMAP (1,1) logical = true
end

assert(~isempty(ctc.GroundTruthLabel1Idx), ...
    'Run identifyGroundTruthUnits() before generateTrainLabels()');

p_umap  = ctc.Parameters.UMAP;
p_outlr = ctc.Parameters.OutlierDetection;
p_train = ctc.Parameters.TrainLabels;

rng(ctc.Parameters.RNGSeed, 'twister');
ctc.clearCache();
ctc.buildNormalizedFeatures();
nf = ctc.NormalizedFeatures;

% -- Stage 1: global normalization + feature selection ------------------------
[X_all, norm_params] = computeGlobalNormalization(nf, p_umap, ctc.Parameters.Harmonization);
ctc.NormalizationParams = norm_params;

% -- Stage 2: unsupervised UMAP -----------------------------------------------
[reduction, n_neighbors] = ctc.runUnsupervisedUMAP(X_all, p_umap, opts.ReuseExistingUMAP);

% -- Map ground-truth label-1 units to training-subset local indices ----------
unit_ids_table    = string(ctc.FeatureStore.UnitTable.UnitID);
unique_ud_ids     = string({nf.unique_ud.UnitID});
resp_uids         = unique(unit_ids_table(ctc.GroundTruthLabel1Idx));
resp_unique       = ismember(unique_ud_ids, resp_uids);
subset_responsive = resp_unique(nf.subset_mask);
responsive_local  = find(subset_responsive);
subset_global_idx = find(nf.subset_mask);

if isempty(responsive_local)
    n_resp_total = sum(ctc.GroundTruthLabel1Idx);
    error('CellTypeClassifier:noGroundTruthLabel1Units', ...
        ['No ground-truth label-1 units found in the training culture subset.\n' ...
         '  Total label-1 units across all cultures: %d\n' ...
         '  Training subset size: %d unique units\n' ...
         'Possible causes:\n' ...
         '  1. Parameters.UMAP.TrainingCultureIdx excludes cultures with label-1 units.\n' ...
         '  2. Parameters.UMAP.GroupingValues does not match the baseline dose in the data.\n' ...
         '  3. identifyGroundTruthUnits() was not run, or produced 0 label-1 units.\n' ...
         'Diagnostic: check sum(ctc.GroundTruthLabel1Idx) and ctc.Parameters.UMAP.TrainingCultureIdx.'], ...
        n_resp_total, sum(nf.subset_mask));
end

X_feat          = X_all(nf.subset_mask, :);
has_explicit_ce = ~isempty(ctc.GroundTruthLabel2Idx) && any(ctc.GroundTruthLabel2Idx);

% -- Stage 3: outlier detection and CE selection ------------------------------
if has_explicit_ce
    ce_uid_table = unique(unit_ids_table(ctc.GroundTruthLabel2Idx));
    ce_subset    = ismember(unique_ud_ids(nf.subset_mask), ce_uid_table);
    ce_local_all = find(ce_subset);
    det = selectExplicitCE(X_feat, responsive_local, subset_responsive, ...
        subset_global_idx, ce_local_all, reduction, p_outlr);
else
    odm = lower(string(p_outlr.Method));

    comm_diagnostics = [];
    if odm == "community"
        det = detectCommunityOutliers(X_feat, responsive_local, subset_responsive, ...
            subset_global_idx, reduction, size(X_all, 1), n_neighbors, ...
            ctc.Parameters.Community, ctc.UMAP, ctc.Parameters.RNGSeed);
        if ~det.success
            comm_diagnostics = det;   % preserve Louvain results even on fallback
            odm = "iforest";
        end
    end

    if odm == "iforest"
        det = selectIforestCE(X_feat, responsive_local, subset_responsive, ...
            subset_global_idx, reduction, p_outlr);
        if ~isempty(comm_diagnostics)
            det.community_ids   = comm_diagnostics.community_ids;
            det.inh_comm_ids    = comm_diagnostics.inh_comm_ids;
            det.community_qvals = comm_diagnostics.community_qvals;
            det.Q_modularity    = comm_diagnostics.Q_modularity;
        end
    elseif odm == "none"
        det = selectNoCE(X_feat, responsive_local, subset_responsive, p_outlr);
    end
end

det.has_explicit_ce = has_explicit_ce;

% -- Stage 4: assemble training labels ----------------------------------------
ctc.TrainLabels = assembleTrainLabels(det, subset_global_idx, nf, p_train, ...
    unit_ids_table, ctc.GroundTruthLabel1Strength);

resp_label    = p_train.GroundTruthLabel1;
counter_label = 3 - resp_label;
fprintf('Training set: %i %s, %i %s candidates\n', ...
    sum(ctc.TrainLabels.sorted_y_train == counter_label), localLabelName(counter_label), ...
    sum(ctc.TrainLabels.sorted_y_train == resp_label),    localLabelName(resp_label));

if ctc.Parameters.Diagnostics.Enable
    ctc.diagnosticTrainLabels();
end
end

function name = localLabelName(label)
    if label == 1; name = 'excitatory';
    elseif label == 2; name = 'inhibitory';
    else; name = sprintf('class_%d', label);
    end
end
