function best_params = optimizeCTCParameters(ctc, rg)
% OPTIMIZECTCPARAMETERS  Bayesian optimization of CellTypeClassifier UMAP parameters.
%
% Optimizes four parameters jointly using Bayesian optimization with a
% silhouette-based objective function evaluated on the supervised test embedding.
%
% Objective: weighted combination of inhibitory and excitatory silhouette scores
% on the out-of-sample test embedding (0.7 * inh_sil + 0.3 * exc_sil).
% The test embedding is used rather than train to get an honest out-of-sample
% measure unaffected by the supervised signal.
%
% Parameters optimized:
%   SupervisedNNeighbors — controls local vs global structure in supervised UMAP
%   TargetWeight         — balance between class separation and manifold geometry
%   ACGWeight            — relative contribution of ACG feature group
%   WaveformWeight       — relative contribution of waveform feature group
%
% INPUTS:
%   ctc - CellTypeClassifier after identifyResponsiveUnits() has been run
%         (ResponsiveUnitIdx and UnitList must be set)
%   rg  - RecordingGroup (used for per-chip consistency logging)
%
% OUTPUT:
%   best_params - optimizableVariable result struct with best parameter values
%
% USAGE:
%   best_params = optimizeCTCParameters(ctc, rg);
%   % Apply best parameters:
%   ctc.Parameters.UMAP.SupervisedNNeighbors = best_params.SupervisedNNeighbors;
%   ctc.Parameters.UMAP.TargetWeight         = best_params.TargetWeight;
%   ctc.Parameters.UMAP.ACGWeight            = best_params.ACGWeight;
%   ctc.Parameters.UMAP.WaveformWeight       = best_params.WaveformWeight;
%   ctc.generateTrainLabels();
%   ctc.classifyUnits();

arguments
    ctc CellTypeClassifier
    rg  RecordingGroup
end

assert(~isempty(ctc.ResponsiveUnitIdx), ...
    'Run identifyResponsiveUnits() before optimizeCTCParameters().');

% ── Search space ──────────────────────────────────────────────────────────────
% Bounds chosen based on:
%   SupervisedNNeighbors: 5–30 (ceiling = minority_class_size/5, enforced in objective)
%   TargetWeight:         0.05–0.4 (above 0.4 risks disconnected manifolds)
%   ACGWeight:            0.5–3.0 (allows ACG to dominate up to 3x or be halved)
%   WaveformWeight:       0.5–3.0 (symmetric range for waveform)
vars = [
    optimizableVariable('SupervisedNNeighbors', [5, 30], ...
        'Type', 'integer')
    optimizableVariable('TargetWeight',         [0.05, 0.40], ...
        'Type', 'real')
    optimizableVariable('ACGWeight',            [0.5, 3.0], ...
        'Type', 'real')
    optimizableVariable('WaveformWeight',       [0.5, 3.0], ...
        'Type', 'real')
];

% ── Objective function ────────────────────────────────────────────────────────
objective = @(params) ctcObjective(params, ctc, rg);

% ── Run Bayesian optimization ─────────────────────────────────────────────────
results = bayesopt(objective, vars, ...
    'MaxObjectiveEvaluations', 40, ...
    'IsObjectiveDeterministic', false, ...  % UMAP has stochastic elements
    'ExplorationRatio',         0.5, ...    % balanced explore/exploit
    'UseParallel',              false, ...  % UMAP is not thread-safe
    'Verbose',                  1, ...      % print progress
    'PlotFcn', {@plotObjectiveModel, @plotMinObjective});  % live plots

best_params = bestPoint(results);

fprintf('\n=== Optimization complete ===\n');
fprintf('Best SupervisedNNeighbors : %d\n',   best_params.SupervisedNNeighbors);
fprintf('Best TargetWeight         : %.3f\n', best_params.TargetWeight);
fprintf('Best ACGWeight            : %.3f\n', best_params.ACGWeight);
fprintf('Best WaveformWeight       : %.3f\n', best_params.WaveformWeight);
fprintf('Best objective value      : %.4f\n', results.MinObjective);
end


function score = ctcObjective(params, ctc, rg)
% CTCOBJECTIVE  Evaluate one parameter configuration.
%
% Runs generateTrainLabels + classifyUnits with the given parameters,
% then computes a weighted silhouette score on the test embedding.
% Returns Inf for degenerate solutions.

    % ── Apply parameters ──────────────────────────────────────────────────────
    ctc.Parameters.UMAP.SupervisedNNeighbors = params.SupervisedNNeighbors;
    ctc.Parameters.UMAP.TargetWeight         = params.TargetWeight;
    ctc.Parameters.UMAP.ACGWeight            = params.ACGWeight;
    ctc.Parameters.UMAP.WaveformWeight       = params.WaveformWeight;

    % ── Run pipeline ──────────────────────────────────────────────────────────
    try
        ctc.generateTrainLabels();

        % Enforce minority class ceiling on SupervisedNNeighbors
        n_inh = sum(ctc.TrainLabels.sorted_y_train == 2);
        safe_nn = floor(n_inh / 5);
        if params.SupervisedNNeighbors > safe_nn
            fprintf('  SupervisedNNeighbors=%d exceeds ceiling %d — invalid\n', ...
                params.SupervisedNNeighbors, safe_nn);
            score = Inf;
            return
        end

        ctc.classifyUnits();

    catch e
        warning('ctcObjective: pipeline failed — %s', e.message);
        score = Inf;
        return
    end

    % ── Check for degenerate solutions ────────────────────────────────────────
    n_inh_classified = sum(ctc.UnitLabels == 2, 'omitnan');
    n_exc_classified = sum(ctc.UnitLabels == 1, 'omitnan');

    if n_inh_classified < 10 || n_exc_classified < 10
        fprintf('  Degenerate: %d inh, %d exc — invalid\n', ...
            n_inh_classified, n_exc_classified);
        score = Inf;
        return
    end

    % ── Compute silhouette on test embedding ──────────────────────────────────
    % Test embedding is the honest out-of-sample measure — unaffected by
    % the supervised signal that shaped the training embedding directly.
    test_mask   = logical(ctc.TrainLabels.umap_test_idx);
    test_labels = ctc.UnitLabels(test_mask);
    test_embed  = ctc.Reduction.Test;
    valid       = ~isnan(test_labels);

    if sum(test_labels(valid) == 2) < 5 || sum(test_labels(valid) == 1) < 5
        fprintf('  Too few test points in one class — invalid\n');
        score = Inf;
        return
    end

    s = silhouette(test_embed(valid, :), test_labels(valid)');

    inh_sil = mean(s(test_labels(valid) == 2));
    exc_sil = mean(s(test_labels(valid) == 1));

    % Weighted combination — inhibitory weighted more heavily as the
    % minority class and harder classification target
    combined_sil = 0.7 * inh_sil + 0.3 * exc_sil;

    % Minimize negative silhouette
    score = -combined_sil;

    % ── Log this evaluation ───────────────────────────────────────────────────
    [chip_idx, ~] = rg.combineMetadataIndices(ctc.UnitList, "ChipID");
    inh_count  = histcounts(chip_idx(ctc.UnitLabels == 2), 1:max(chip_idx)+1);
    all_count  = histcounts(chip_idx, 1:max(chip_idx)+1);
    ie_per_chip = inh_count ./ all_count;

    fprintf(['  NNeigh=%2d  TW=%.2f  ACGW=%.2f  WFW=%.2f  ' ...
             'inh_sil=%.3f  exc_sil=%.3f  score=%.4f  ' ...
             'inh%%=[%s]\n'], ...
        params.SupervisedNNeighbors, params.TargetWeight, ...
        params.ACGWeight, params.WaveformWeight, ...
        inh_sil, exc_sil, score, ...
        num2str(round(ie_per_chip * 100), '%d '));
end