function diagnosticActivityConfound(ctc)
% DIAGNOSTICACTIVITYCONFOUND  Check whether curated activity/waveform-shape
%   features distinguish the classifier's predicted classes as strongly as
%   the raw ACG/waveform features the classifier was actually built from.
%
% Motivation: diagnosticClassification's "Mean ACG by class" panel can show a
% large baseline-height gap between predicted classes with almost no
% corresponding difference in mean waveform shape -- consistent with the
% classifier separating units by overall firing rate/regularity rather than
% true cell identity. Since identifyGroundTruthUnits' ground truth is itself
% firing-rate-based, that would be a real circularity, not just a cosmetic
% concern.
%
% This pulls FeatureCatalog's curated, literature-validated E/I descriptors
% (FiringRate, CV2InterSpikeInterval, RevisedLocalVariation, FanoFactor,
% T2Pdelay, HalfWidth, Asymmetry, RegularityFrequency, RegularityFit) straight
% from FeatureStore.UnitTable -- columns the classifier's own feature matrix
% never sees -- and compares them between predicted classes directly
% (Wilcoxon/t-test + Cohen's d via StatisticalTest.compareGroups).
%
% Effect sizes here comparable to or larger than the classifier's own top
% raw-feature Cohen's d (dashed reference line, tile 1) are evidence the
% classifier is riding a firing-rate proxy rather than waveform identity.
% Small, non-significant effect sizes here instead support the classifier
% capturing genuine, independent cell-type information.
%
% Requires ctc.UnitLabels (run classify() first).

arguments
    ctc CellTypeClassifier
end

C_INH = [0.8 0.1 0.1];
C_EXC = [0.1 0.3 0.8];

assert(~isempty(ctc.UnitLabels), 'Run classify() before diagnosticActivityConfound().');

unit_tbl = ctc.FeatureStore.UnitTable;
labels   = ctc.UnitLabels;
valid_lb = labels == 1 | labels == 2;

candidate_feats = [FeatureCatalog.features("ActivityFeatures", "core"), ...
                    FeatureCatalog.features("WaveformFeatures", "core"), ...
                    FeatureCatalog.features("RegularityFeatures", "core")];
present = ismember(candidate_feats, string(unit_tbl.Properties.VariableNames));
if ~any(present)
    warning('CellTypeClassifier:diagnosticActivityConfound', ...
        'None of the curated feature columns (%s) found in FeatureStore.UnitTable.', ...
        strjoin(candidate_feats, ', '));
    return
end
feats  = candidate_feats(present);
n_feat = numel(feats);

d_vals = nan(1, n_feat);
p_vals = nan(1, n_feat);
for fi = 1:n_feat
    vals = double(unit_tbl.(feats(fi)))';
    ok   = valid_lb & isfinite(vals);
    if sum(ok & labels == 1) < 2 || sum(ok & labels == 2) < 2
        continue
    end
    result = StatisticalTest.compareGroups(vals(ok)', labels(ok)');
    % compareGroups sorts groups ascending (exc=1 first), so effect_size is
    % (mu_exc - mu_inh)/pooled -- flip sign to match this codebase's
    % "inh - exc" Cohen's d convention used in diagnosticClassification.
    d_vals(fi) = -result.effect_size;
    p_vals(fi) = result.p;
end

% -- Reference: classifier's own top raw-feature |Cohen's d| ------------------
ref_d = NaN;
try
    np_h = ctc.NormalizationParams;
    nf_h = ctc.NormalizedFeatures;
    X = normalize(nf_h.X_pergroup, 'center', np_h.mu_global, 'scale', np_h.sigma_global);
    X(:, np_h.nan_cols) = [];
    X = X ./ np_h.scale;
    if isfield(np_h, 'feature_selection_mask') && ~all(np_h.feature_selection_mask)
        X = X(:, np_h.feature_selection_mask);
    end
    labels_unique = labels(nf_h.unique_to_rep);
    m1 = labels_unique == 1;
    m2 = labels_unique == 2;
    n1 = sum(m1);
    n2 = sum(m2);
    if n1 >= 2 && n2 >= 2
        mu1 = mean(X(m1, :), 1);
        mu2 = mean(X(m2, :), 1);
        s1  = std(X(m1, :), 0, 1);
        s2  = std(X(m2, :), 0, 1);
        pooled = sqrt(((n1-1).*s1.^2 + (n2-1).*s2.^2) / (n1+n2-2));
        pooled(pooled == 0) = 1;
        ref_d = max(abs((mu2 - mu1) ./ pooled));
    end
catch
    % leave ref_d as NaN if reconstruction fails -- reference line just omitted
end

fig = figure('Visible', 'on');
set(fig, 'Position', [100 100 1100 500]);
tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Activity/waveform confound check', 'FontWeight', 'bold');

% -- Tile 1: Cohen's d per curated feature, sorted by magnitude --------------
nexttile(1);
try
    valid_d   = ~isnan(d_vals);
    idx_valid = find(valid_d);
    [~, ord_rel] = sort(abs(d_vals(idx_valid)), 'descend');
    ord = idx_valid(ord_rel);

    hold on;
    for i = 1:numel(ord)
        fi  = ord(i);
        col = C_INH;
        if d_vals(fi) < 0; col = C_EXC; end
        barh(i, d_vals(fi), 'FaceColor', col, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
        text_offset = 0.02;
        if d_vals(fi) < 0; text_offset = -0.02; end
        text(d_vals(fi) + text_offset, i, sprintf('p=%.1e', p_vals(fi)), ...
            'FontSize', 7, 'VerticalAlignment', 'middle');
    end
    if isfinite(ref_d)
        xline(ref_d,  '--k', 'LineWidth', 1.5, 'DisplayName', 'Classifier top |d|');
        xline(-ref_d, '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    end
    hold off;
    yticks(1:numel(ord));
    yticklabels(feats(ord));
    xlabel("Cohen's d (inh - exc)");
    if isfinite(ref_d)
        title(sprintf('Curated feature effect sizes (classifier top |d|=%.2f)', ref_d));
    else
        title('Curated feature effect sizes (predicted classes)');
    end
    box off;
catch ME
    cla; axis off;
    text(0.5, 0.5, sprintf('Error: %s', ME.message), 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', [0.6 0 0]);
end

% -- Tile 2: distribution of the strongest curated feature by class ----------
nexttile(2);
try
    [~, best_i] = max(abs(d_vals));
    best_feat = feats(best_i);
    vals = double(unit_tbl.(best_feat))';
    ok   = valid_lb & isfinite(vals);

    v_exc = vals(ok & labels == 1);
    v_inh = vals(ok & labels == 2);

    edges = linspace(min(vals(ok)), max(vals(ok)), 30);
    hold on;
    histogram(v_exc, edges, 'FaceColor', C_EXC, 'FaceAlpha', 0.6, 'EdgeColor', 'none', ...
        'Normalization', 'probability', 'DisplayName', 'Excitatory (pred)');
    histogram(v_inh, edges, 'FaceColor', C_INH, 'FaceAlpha', 0.6, 'EdgeColor', 'none', ...
        'Normalization', 'probability', 'DisplayName', 'Inhibitory (pred)');
    hold off;
    xlabel(best_feat, 'Interpreter', 'none');
    ylabel('Fraction of units');
    title(sprintf('%s by predicted class (d=%.2f, p=%.1e)', ...
        best_feat, d_vals(best_i), p_vals(best_i)), 'Interpreter', 'none');
    legend('Location', 'best', 'Box', 'off', 'FontSize', 8);
    box off;
catch ME
    cla; axis off;
    text(0.5, 0.5, sprintf('Error: %s', ME.message), 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', [0.6 0 0]);
end

ctc.saveDiagnosticFigure(fig, 'activity_confound');
end
