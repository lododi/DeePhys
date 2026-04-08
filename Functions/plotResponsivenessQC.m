function plotResponsivenessQC(ctc, opts)
% PLOTRESPONSIVENESSQC  Diagnostic figure for identifyResponsiveUnits results.
%
% Produces a single figure with four panels:
%   1. Filter funnel       — units passing each stage, broken down by culture
%   2. Spearman rho dist   — histogram of rho values, threshold marked
%   3. Dose-response curves — mean FR curves for passing vs failing units
%   4. DR heatmap          — all units sorted by rho, cultures indicated
%
% INPUTS:
%   ctc  - CellTypeClassifier after identifyResponsiveUnits()
%   opts - (optional) name-value pairs:
%       MaxCurvesPerClass  (default 50)  max DR curves to plot per class
%       FigureTitle        (default '')  suptitle string
%
% REQUIRES: ctc.ResponsivenessDetail to be populated.

arguments
    ctc  CellTypeClassifier
    opts.MaxCurvesPerClass (1,1) double = 50
    opts.FigureTitle       (1,:) char   = ''
end

assert(~isempty(ctc.ResponsivenessDetail), ...
    'Run identifyResponsiveUnits() before plotResponsivenessQC().');

d  = ctc.ResponsivenessDetail;
p  = ctc.Parameters.Bootstrap;

% ── Concatenate across cultures ───────────────────────────────────────────────
n_cultures  = numel(d.GlobalUnitIdx);
dose_labels = arrayfun(@(x) sprintf('%.2g', x), p.DoseValues, ...
    'UniformOutput', false);
dose_labels{1} = 'BL';     % baseline label
n_epochs    = numel(p.DoseValues);

% Per-unit arrays (all cultures concatenated)
all_rho      = [];
all_fc       = [];
all_mono     = [];
all_fcpass   = [];
all_bspass   = [];
all_passall  = [];
all_curves   = [];   % (N_total x N_epochs)
all_culture  = [];   % culture index per unit

for c = 1:n_cultures
    if isempty(d.GlobalUnitIdx{c}); continue; end
    n_c = numel(d.PassAll{c});
    all_rho     = [all_rho,    d.SpearmanRho{c}];                       
    all_fc      = [all_fc,     d.FoldChange{c}];                        
    all_mono   = logical([all_mono,   logical(d.PassMonotonicity{c})]); 
    all_fcpass = logical([all_fcpass, logical(d.PassFoldChange{c})]);   
    all_bspass = logical([all_bspass, logical(d.PassBootstrap{c})]);    
    all_passall= logical([all_passall,logical(d.PassAll{c})]);          
    all_curves  = [all_curves; d.DoseResponseCurve{c}];                 
    all_culture = [all_culture, c * ones(1, n_c)];                      
end

n_total   = numel(all_rho);
n_passing = sum(all_passall);

% Normalise DR curves by each unit's max FR for display
% (keeps silent units from dominating the colour scale)
curve_max = max(all_curves, [], 2);
curve_max(curve_max == 0) = 1;
norm_curves = all_curves ./ curve_max;

% Culture colours — one per culture, used across all panels
cmap_cultures = lines(n_cultures);

% ── Figure layout ─────────────────────────────────────────────────────────────
fig = figure('Color', 'w', 'Position', [50 50 1400 900]);
if ~isempty(opts.FigureTitle)
    sgtitle(opts.FigureTitle, 'FontSize', 14, 'FontWeight', 'bold');
end

% Layout: 2 rows x 3 cols, panels span as follows:
%   [1,1]-[1,2]  filter funnel (wide)
%   [1,3]        rho histogram
%   [2,1]-[2,2]  DR heatmap (wide)
%   [2,3]        DR curves
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'tight');

% ── Panel 1: Filter funnel ────────────────────────────────────────────────────
ax1 = nexttile(tl, 1, [1 2]);
hold(ax1, 'on');

stage_labels = {'All units', 'Monotone', '+ Fold change', '+ Bootstrap'};
n_stages     = numel(stage_labels);

% Count per culture per stage
counts = zeros(n_cultures, n_stages);
for c = 1:n_cultures
    if isempty(d.GlobalUnitIdx{c}); continue; end
    mask = all_culture == c;
    counts(c, 1) = sum(mask);
    counts(c, 2) = sum(all_mono(mask));
    counts(c, 3) = sum(all_mono(mask) & all_fcpass(mask));
    counts(c, 4) = sum(all_passall(mask));
end

zero_cultures = find(sum(counts, 2) == 0);
if ~isempty(zero_cultures)
    fprintf('Cultures with no units in funnel: %s\n', num2str(zero_cultures'));
end

% Stacked bar — each culture is a colour
b = bar(ax1, counts', 'stacked');  % transpose: stages on x-axis, cultures stacked
for c = 1:n_cultures
    b(c).FaceColor = cmap_cultures(c, :);
    b(c).EdgeColor = 'none';
    b(c).DisplayName = sprintf('Culture %d', c);
end
stage_totals = sum(counts, 1);  % sum across cultures per stage
for s = 1:n_stages
    text(ax1, s, stage_totals(s) + 0.5, num2str(stage_totals(s)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 9, 'FontWeight', 'bold');
end

% Overlay total count labels on top of each bar group
stage_totals = sum(counts, 1);
for s = 1:n_stages
    text(ax1, s, stage_totals(s) + 0.5, num2str(stage_totals(s)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 9, 'FontWeight', 'bold');
end

set(ax1, 'XTick', 1:n_stages, 'XTickLabel', stage_labels, ...
    'XTickLabelRotation', 15, 'Box', 'off');
ylabel(ax1, 'Unit count');
title(ax1, 'Filter funnel');
legend(ax1, 'Location', 'northeast', 'Box', 'off');

% Add percentage labels on the last bar
pct_passing = 100 * stage_totals(end) / max(stage_totals(1), 1);
text(ax1, n_stages, stage_totals(end) / 2, ...
    sprintf('%.1f%% pass', pct_passing), ...
    'HorizontalAlignment', 'center', 'Color', 'w', ...
    'FontSize', 9, 'FontWeight', 'bold');

% ── Panel 2: Spearman rho distribution ───────────────────────────────────────
ax2 = nexttile(tl, 3);
hold(ax2, 'on');

rho_edges  = linspace(-1, 1, 41);
rho_fail   = all_rho(~all_passall);
rho_pass   = all_rho(all_passall);

histogram(ax2, rho_fail, rho_edges, ...
    'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none', ...
    'DisplayName', sprintf('Non-responsive (n=%d)', numel(rho_fail)));
histogram(ax2, rho_pass, rho_edges, ...
    'FaceColor', [0.85 0.33 0.10], 'EdgeColor', 'none', ...
    'DisplayName', sprintf('Responsive (n=%d)', n_passing));

% Threshold line
xline(ax2, p.MonotonicityThreshold, '--k', ...
    'LineWidth', 1.5, 'Label', ...
    sprintf('\\rho = %.2f', p.MonotonicityThreshold), ...
    'LabelVerticalAlignment', 'bottom');

xlabel(ax2, 'Spearman \rho (dose rank vs mean FR)');
ylabel(ax2, 'Unit count');
title(ax2, 'Monotonicity distribution');
legend(ax2, 'Location', 'northwest', 'Box', 'off');
xlim(ax2, [-1 1]);
box(ax2, 'off');

% ── Panel 3: DR heatmap ───────────────────────────────────────────────────────
ax3 = nexttile(tl, 4, [1 2]);

% Sort all units by Spearman rho (descending)
[sorted_rho, sort_idx] = sort(all_rho, 'descend');
sorted_curves  = norm_curves(sort_idx, :);
sorted_culture = all_culture(sort_idx);
sorted_pass    = all_passall(sort_idx);

imagesc(ax3, sorted_curves);
colormap(ax3, 'hot');
cb = colorbar(ax3);
cb.Label.String = 'Normalised FR';

% Mark the pass/fail boundary with a horizontal line
n_pass_sorted = sum(sorted_pass);   % passing units are at the top (highest rho)
if n_pass_sorted > 0 && n_pass_sorted < n_total
    yline(ax3, n_pass_sorted + 0.5, 'c-', 'LineWidth', 2, ...
        'Label', 'Pass threshold', ...
        'LabelHorizontalAlignment', 'left', ...
        'LabelVerticalAlignment', 'bottom');
end

% Culture indicator strip on the right y-axis
% Draw coloured markers at the right edge of the heatmap
ax3_pos = ax3.Position;
ax_strip = axes(fig, 'Position', ...
    [ax3_pos(1) + ax3_pos(3) + 0.005, ax3_pos(2), 0.012, ax3_pos(4)]);
culture_img = reshape(sorted_culture, [], 1);
imagesc(ax_strip, culture_img);
colormap(ax_strip, cmap_cultures(1:n_cultures, :));
clim(ax_strip, [1, n_cultures]);
axis(ax_strip, 'off');
title(ax_strip, 'Cul.', 'FontSize', 7);

set(ax3, 'XTick', 1:n_epochs, 'XTickLabel', dose_labels, ...
    'YDir', 'normal');
xlabel(ax3, 'Dose epoch');
ylabel(ax3, sprintf('Units (sorted by \\rho, n=%d total)', n_total));
title(ax3, 'Dose-response heatmap (normalised FR)');

% Overlay rho values as a right-side axis tick reference
% — add a text annotation showing rho range at top and bottom
text(ax3, n_epochs + 0.1, 1, ...
    sprintf('\\rho=%.2f', sorted_rho(1)), ...
    'FontSize', 7, 'Color', [0 0.6 0], 'Clipping', 'off');
text(ax3, n_epochs + 0.1, n_total, ...
    sprintf('\\rho=%.2f', sorted_rho(end)), ...
    'FontSize', 7, 'Color', [0.6 0 0], 'Clipping', 'off');

% ── Panel 4: DR curves (passing vs failing) ───────────────────────────────────
ax4 = nexttile(tl, 6);
hold(ax4, 'on');

x_doses = 1:n_epochs;

% Subsample for legibility
pass_idx = find(all_passall);
fail_idx = find(~all_passall);

rng(ctc.Parameters.RNGSeed, 'twister');   % reproducible subsample
n_pass_plot = min(opts.MaxCurvesPerClass, numel(pass_idx));
n_fail_plot = min(opts.MaxCurvesPerClass, numel(fail_idx));
pass_sample = pass_idx(randperm(numel(pass_idx), n_pass_plot));
fail_sample = fail_idx(randperm(numel(fail_idx), n_fail_plot));

% Individual traces (thin, transparent)
for u = fail_sample
    plot(ax4, x_doses, norm_curves(u, :), ...
        'Color', [0.6 0.6 0.9 0.15], 'LineWidth', 0.5);
end
for u = pass_sample
    plot(ax4, x_doses, norm_curves(u, :), ...
        'Color', [0.85 0.33 0.10 0.2], 'LineWidth', 0.5);
end

% Mean traces (thick)
if ~isempty(fail_sample)
    mean_fail = mean(norm_curves(fail_idx, :), 1);
    plot(ax4, x_doses, mean_fail, ...
        'Color', [0.3 0.3 0.7], 'LineWidth', 2.5, ...
        'DisplayName', sprintf('Non-responsive mean (n=%d)', numel(fail_idx)));
end
if ~isempty(pass_sample)
    mean_pass = mean(norm_curves(pass_idx, :), 1);
    plot(ax4, x_doses, mean_pass, ...
        'Color', [0.7 0.1 0.0], 'LineWidth', 2.5, ...
        'DisplayName', sprintf('Responsive mean (n=%d)', numel(pass_idx)));
end

set(ax4, 'XTick', x_doses, 'XTickLabel', dose_labels, ...
    'XTickLabelRotation', 15, 'Box', 'off');
xlabel(ax4, 'Dose epoch');
ylabel(ax4, 'Normalised FR');
title(ax4, sprintf('DR curves (showing %d/%d per class)', ...
    opts.MaxCurvesPerClass, max(numel(pass_idx), numel(fail_idx))));
legend(ax4, 'Location', 'northwest', 'Box', 'off');
ylim(ax4, [0, 1.05]);