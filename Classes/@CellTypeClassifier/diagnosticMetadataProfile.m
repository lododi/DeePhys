function diagnosticMetadataProfile(ctc, opts)
% DIAGNOSTICMETADATAPROFILE  Mean waveform + ACG panels grouped by a metadata field.
%
% One column per unique value of opts.Field (default "E_IRatio"); top row is
% mean waveform +/- SEM, bottom row is mean ACG +/- SEM. Independent of the
% classifier's predicted labels — groups units purely by a FeatureStore.UnitTable
% metadata column, so this checks whether raw waveform/ACG shape actually
% tracks a known ground-truth quantity (e.g. culture E/I mixing ratio) before
% any classification is applied.
%
%   ctc.diagnosticMetadataProfile();                    % groups by E_IRatio
%   ctc.diagnosticMetadataProfile('Field', 'Mutation');  % groups by any column

arguments
    ctc CellTypeClassifier
    opts.Field (1,1) string = "E_IRatio"
end

fs = ctc.FeatureStore;
assert(ismember(opts.Field, string(fs.UnitTable.Properties.VariableNames)), ...
    'CellTypeClassifier:diagnosticMetadataProfile', ...
    'Field "%s" not found in FeatureStore.UnitTable. Available: %s', ...
    opts.Field, strjoin(string(fs.UnitTable.Properties.VariableNames), ', '));

if isempty(ctc.HarmonizedWaveforms) || isempty(ctc.HarmonizedACGs)
    ctc.buildNormalizedFeatures();
end
wf  = ctc.HarmonizedWaveforms;   % (N_samples x N_units), UnitDataArray order
acg = ctc.HarmonizedACGs;        % (N_bins x N_units), UnitDataArray order
sr  = ctc.HarmonizedSR;
assert(~isempty(wf) && ~isempty(acg), ...
    'HarmonizedWaveforms/HarmonizedACGs not populated — run buildNormalizedFeatures() (or generateTrainLabels()) first.');

% -- Map FeatureStore-order metadata values onto UnitDataArray order ----------
% (mirrors the UnitID-matching pattern used in diagnosticGroundTruthUnits)
field_vals_table = string(fs.UnitTable.(opts.Field));
unit_ids_table   = string(fs.UnitTable.UnitID);
ud_ids_all       = string({ctc.UnitDataArray.UnitID});

[unique_uids, first_idx] = unique(unit_ids_table, 'stable');
uid_to_field = containers.Map(cellstr(unique_uids), cellstr(field_vals_table(first_idx)));

group_vals_ud = strings(1, numel(ud_ids_all));
for i = 1:numel(ud_ids_all)
    key = char(ud_ids_all(i));
    if isKey(uid_to_field, key)
        group_vals_ud(i) = string(uid_to_field(key));
    end
end

valid_group = group_vals_ud ~= "" & group_vals_ud ~= "NaN" & group_vals_ud ~= "<missing>";
unique_groups = unique(group_vals_ud(valid_group));
assert(~isempty(unique_groups), ...
    'No valid (non-empty) values found for field "%s".', opts.Field);

% Sort numerically by leading number (handles "0:100", "25:75", ... ratios)
% when every value parses; otherwise fall back to alphabetical (default 'unique' order).
sort_keys = nan(1, numel(unique_groups));
for gi = 1:numel(unique_groups)
    tok = regexp(unique_groups(gi), '^-?\d+\.?\d*', 'match', 'once');
    if ~isempty(tok)
        sort_keys(gi) = str2double(tok);
    end
end
if all(~isnan(sort_keys))
    [~, ord] = sort(sort_keys);
    unique_groups = unique_groups(ord);
end

n_groups = numel(unique_groups);

% -- Time/lag axes (mirrors diagnosticGroundTruthUnits) -----------------------
n_samp = size(wf, 1);
ph_h   = ctc.Parameters.Harmonization;
if isfield(ph_h, 'WaveformPreTrough') && isfield(ph_h, 'WaveformPostTrough')
    t_ms = linspace(-ph_h.WaveformPreTrough, ph_h.WaveformPostTrough, n_samp);
else
    t_ms = (0:(n_samp-1)) / sr * 1000;
end
n_bins     = size(acg, 1);
acg_bin_ms = ph_h.ACGBinSize * 1000;
lag_ms     = ((0:(n_bins-1)) - floor(n_bins/2)) * acg_bin_ms;

fig = figure('Visible', 'on');
set(fig, 'Position', [100 100 min(220 * n_groups + 150, 1800), 650]);
tl = tiledlayout(2, n_groups, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Waveform / ACG profile by %s', opts.Field), ...
    'FontWeight', 'bold', 'Interpreter', 'none');

cmap = lines(n_groups);
wf_ylims  = [Inf, -Inf];
acg_ylims = [Inf, -Inf];

% -- Row 1: mean waveform per group -------------------------------------------
for gi = 1:n_groups
    mask      = group_vals_ud == unique_groups(gi);
    n_units_g = sum(mask);

    nexttile(gi);
    try
        if n_units_g == 0
            error('No units');
        end
        mu  = mean(wf(:, mask), 2);
        sem = std(wf(:, mask), 0, 2) / sqrt(n_units_g);
        shadedErrorBar_local(t_ms, mu', sem', cmap(gi, :), 0.3);
        plot(t_ms, mu, '-', 'Color', cmap(gi, :), 'LineWidth', 1.5);
        wf_ylims = [min(wf_ylims(1), min(mu - sem)), max(wf_ylims(2), max(mu + sem))];
        xlabel('Time (ms)');
        if gi == 1; ylabel('Waveform (norm.)'); end
        title(sprintf('%s (n=%d)', unique_groups(gi), n_units_g), 'Interpreter', 'none');
        box off;
    catch ME
        cla; axis off;
        text(0.5, 0.5, sprintf('%s\n(n=%d)', ME.message, n_units_g), ...
            'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 8);
    end
end

% -- Row 2: mean ACG per group -------------------------------------------------
for gi = 1:n_groups
    mask      = group_vals_ud == unique_groups(gi);
    n_units_g = sum(mask);

    nexttile(n_groups + gi);
    try
        if n_units_g == 0
            error('No units');
        end
        mu  = mean(acg(:, mask), 2);
        sem = std(acg(:, mask), 0, 2) / sqrt(n_units_g);
        shadedErrorBar_local(lag_ms, mu', sem', cmap(gi, :), 0.3);
        plot(lag_ms, mu, '-', 'Color', cmap(gi, :), 'LineWidth', 1.5);
        acg_ylims = [min(acg_ylims(1), min(mu - sem)), max(acg_ylims(2), max(mu + sem))];
        xlabel('Lag (ms)');
        if gi == 1; ylabel('ACG (norm.)'); end
        box off;
    catch ME
        cla; axis off;
        text(0.5, 0.5, sprintf('%s\n(n=%d)', ME.message, n_units_g), ...
            'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 8);
    end
end

% Harmonize y-limits across columns within each row for visual comparability
if all(isfinite(wf_ylims)) && wf_ylims(1) < wf_ylims(2)
    for gi = 1:n_groups; nexttile(gi); ylim(wf_ylims); end
end
if all(isfinite(acg_ylims)) && acg_ylims(1) < acg_ylims(2)
    for gi = 1:n_groups; nexttile(n_groups + gi); ylim(acg_ylims); end
end

ctc.saveDiagnosticFigure(fig, sprintf('%s_profile', lower(opts.Field)));

end

%% ── Local helper: shaded error bar without toolbox dependency ────────────────

function shadedErrorBar_local(x, y, err, color, alpha_val)
% Draw a filled patch for mean +/- err. x, y, err are 1xN row vectors.
    x      = x(:)';
    y      = y(:)';
    err    = err(:)';
    x_patch = [x, fliplr(x)];
    y_patch = [y + err, fliplr(y - err)];
    patch(x_patch, y_patch, color, ...
        'FaceAlpha', alpha_val, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    hold on;
end
