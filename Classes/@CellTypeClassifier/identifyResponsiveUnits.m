function identifyResponsiveUnits(ctc, metadata_filter)
% IDENTIFYRESPONSIVEUNITS  Identify inhibitory candidates from a cumulative
% dose-response recording using three complementary criteria:
%
%   1. Monotonicity (primary): Spearman correlation between dose rank and
%      mean firing rate per epoch must exceed MonotonicityThreshold.
%      Exploits the full dose-response structure — the strongest biological
%      constraint available in this experimental design.
%
%   2. Effect size: fold change from baseline to maximum dose must exceed
%      MinFoldChange. Ensures biological relevance independent of
%      statistical significance.
%
%   3. Bootstrap confirmation (permutation test on extreme comparison):
%      Retained from original pipeline but with relaxed alpha
%      (ConfirmationAlpha), since monotonicity + effect size are already
%      doing the primary filtering.
%
% A unit must pass ALL THREE criteria to be flagged as a responsive
% (inhibitory) candidate.
%
% Additionally computes and stores per-unit dose-response curves
% (mean FR per epoch) in ctc.ResponsivenessDetail for downstream use
% (e.g. EC50 estimation, subtype analysis).
%
% INPUTS:
%   ctc             - CellTypeClassifier object
%   metadata_filter - (optional) {field, value} cell pair to restrict which
%                     cultures contribute responsive units.
%                     Default: all cultures.
%
% Sets:
%   ctc.UnitList              - (1 x N) Unit array across all cultures
%   ctc.ResponsiveUnitIdx     - (1 x N) logical, true = inhibitory candidate
%   ctc.ResponsivenessDetail  - struct with per-unit filter outcomes and
%                               dose-response curves

arguments
    ctc             CellTypeClassifier
    metadata_filter cell = {}
end

p  = ctc.Parameters.Bootstrap;
rg = ctc.RecordingGroup;

% ── Validate dose-response parameters ────────────────────────────────────────
n_epochs = size(p.DoseWindows, 1);
assert(numel(p.DoseValues) == n_epochs, ...
    'Bootstrap.DoseValues must have one entry per row of Bootstrap.DoseWindows.');
assert(p.DoseWindows(1, 1) == 0, ...
    'First DoseWindows row must start at t=0 (baseline window).');
assert(size(p.DoseWindows, 1) == numel(p.DoseValues), ...
    'DoseWindows and DoseValues must have the same number of entries.');

% Derive pre/post cutouts directly from DoseWindows
pre_cutout  = p.DoseWindows(1,   :);
post_cutout = p.DoseWindows(end, :);

% Ordinal rank for Spearman — rank-based so log transform not needed
dose_rank = (1:n_epochs)';

% ── Preallocate outputs ───────────────────────────────────────────────────────
responsive_idx = [];
unit_list      = Unit.empty;

detail = struct( ...
    'GlobalUnitIdx',    {{}}, ...
    'PassMonotonicity', {{}}, ...
    'PassFoldChange',   {{}}, ...
    'PassBootstrap',    {{}}, ...
    'PassAll',          {{}}, ...
    'SpearmanRho',      {{}}, ...
    'FoldChange',       {{}}, ...
    'DoseResponseCurve',{{}}, ...
    'EC50',             {{}});

for c = 1:length(rg.Cultures)
    culture        = rg.Cultures(c);
    n_units_before = numel(unit_list);

    % ── Apply metadata filter ─────────────────────────────────────────────────
    if ~isempty(metadata_filter) && isfield(culture.Metadata, metadata_filter{1})
        culture_matches = any([culture.Metadata.(metadata_filter{1})] == metadata_filter{2});
    else
        culture_matches = true;
    end

    if ~culture_matches
        % Still accumulate units so global indexing stays consistent
        % Use baseline recording units as the reference unit list
        concs_nc = arrayfun(@(r) r.Metadata.Concentration, culture.Recordings);
        [~, sort_idx_nc] = sort(concs_nc);
        unit_list = [unit_list, culture.Recordings(sort_idx_nc(1)).Units]; %#ok<AGROW>
        continue
    end

    % ── Sort recordings by concentration (= dose order) ──────────────────────
    concs = arrayfun(@(r) r.Metadata.Concentration, culture.Recordings);
    [sorted_concs, sort_idx] = sort(concs);
    sorted_recs = culture.Recordings(sort_idx);

    % Validate recording concentrations match DoseValues
    assert(numel(sorted_recs) == n_epochs, ...
        'Culture %d has %d recordings but DoseWindows has %d epochs.', ...
        c, numel(sorted_recs), n_epochs);
    assert(all(abs(sorted_concs(:) - p.DoseValues(:)) < 1e-9), ...
        'Culture %d concentrations [%s] do not match DoseValues [%s].', ...
        c, num2str(sorted_concs), num2str(p.DoseValues));

    % ── Build unit list from baseline recording ───────────────────────────────
    % Baseline recording is the reference — all epochs matched to these units
    baseline_units = sorted_recs(1).Units;
    n_units_c      = numel(baseline_units);
    baseline_ids   = arrayfun(@(u) char(u.StableID), baseline_units, ...
        'UniformOutput', false);

    % ── Compute per-unit FR for each dose epoch ───────────────────────────────
    % Each recording has LOCAL spike times (0 to ~epoch_dur seconds).
    % Match units across recordings by StableID.
    epoch_dur = diff(p.DoseWindows, 1, 2);   % (N_epochs x 1) duration in seconds
    dr_curves = zeros(n_units_c, n_epochs);

    for e = 1:n_epochs
        rec     = sorted_recs(e);
        rec_ids = arrayfun(@(u) char(u.StableID), rec.Units, ...
            'UniformOutput', false);
        rec_dur = epoch_dur(e);

        for u = 1:n_units_c
            match_idx = find(strcmp(rec_ids, baseline_ids{u}), 1);
            if isempty(match_idx)
                % Unit not detected in this epoch — FR stays 0
                continue
            end
            st = double(rec.Units(match_idx).SpikeTimes);
            % Spike times are local to each recording so count all spikes
            dr_curves(u, e) = numel(st) / rec_dur;
        end
    end

    % ── Stage 1: Monotonicity ─────────────────────────────────────────────────
    rho = zeros(n_units_c, 1);
    for u = 1:n_units_c
        rho(u) = corr(dose_rank, dr_curves(u, :)', 'Type', 'Spearman');
    end
    mono_pass = rho >= p.MonotonicityThreshold;

    % ── Stage 2: Effect size ──────────────────────────────────────────────────
    baseline_fr = dr_curves(:, 1);
    maxdose_fr  = dr_curves(:, end);

    fold_change = nan(n_units_c, 1);
    fc_pass     = false(n_units_c, 1);

    active               = baseline_fr > 0;
    fold_change(active)  = maxdose_fr(active) ./ baseline_fr(active);
    fc_pass(active)      = fold_change(active) >= p.MinFoldChange;

    % Silent-to-active: valid response if spikes appear at max dose
    silent              = ~active;
    fold_change(silent) = inf;
    fc_pass(silent)     = maxdose_fr(silent) > 0;

    % ── Stage 3: Bootstrap confirmation ──────────────────────────────────────
    % Only run on units passing stages 1 & 2 to save compute time.
    % Uses relaxed ConfirmationAlpha since monotonicity + effect size
    % are the primary filters.
    pre_post_pass  = false(n_units_c, 1);
    candidate_mask = mono_pass & fc_pass;

    if any(candidate_mask)
        response = culture.bootstrapResponse( ...
            pre_cutout, post_cutout, p.BinSize, p.NIter, p.ConfirmationAlpha);
        bootstrap_increase = false(n_units_c, 1);
        bootstrap_increase(response.increase) = true;
        pre_post_pass(candidate_mask) = bootstrap_increase(candidate_mask);
    end

    % ── Combine all three criteria ────────────────────────────────────────────
    passes_all     = mono_pass & fc_pass & pre_post_pass;
    local_resp_idx = find(passes_all)';
    responsive_idx = [responsive_idx, local_resp_idx + n_units_before]; %#ok<AGROW>

    % ── EC50 estimation ───────────────────────────────────────────────────────
    ec50          = nan(n_units_c, 1);
    nonzero_mask  = p.DoseValues > 0;
    nonzero_doses = p.DoseValues(nonzero_mask);

    for u = find(passes_all)'
        fr_nonzero = dr_curves(u, nonzero_mask);
        fr_max     = max(fr_nonzero);
        if fr_max == 0; continue; end
        fr_norm = fr_nonzero / fr_max;

        try
            hill = @(b, x) x.^b(2) ./ (b(1).^b(2) + x.^b(2));
            b0   = [1, 1];
            opts = statset('MaxIter', 200, 'Display', 'off');
            b    = nlinfit(nonzero_doses, fr_norm, hill, b0, opts);
            if b(1) > 0 && b(1) < max(nonzero_doses) * 10
                ec50(u) = b(1);
            end
        catch
            % nlinfit failed — leave as NaN
        end
    end

    % ── Store per-culture diagnostics ─────────────────────────────────────────
    global_idx = n_units_before + (1:n_units_c);
    detail.GlobalUnitIdx{c}     = global_idx;
    detail.PassMonotonicity{c}  = mono_pass';
    detail.PassFoldChange{c}    = fc_pass';
    detail.PassBootstrap{c}     = pre_post_pass';
    detail.PassAll{c}           = passes_all';
    detail.SpearmanRho{c}       = rho';
    detail.FoldChange{c}        = fold_change';
    detail.DoseResponseCurve{c} = dr_curves;
    detail.EC50{c}              = ec50';

    unit_list = [unit_list, baseline_units]; %#ok<AGROW>

    fprintf(['Culture %2d: %3d monotone | %3d + fold change ' ...
        '| %3d + bootstrap  →  %3d candidates\n'], ...
        c, sum(mono_pass), sum(mono_pass & fc_pass), ...
        sum(passes_all), sum(passes_all));
end

% ── Finalise ──────────────────────────────────────────────────────────────────
ctc.UnitList                          = unit_list;
ctc.ResponsiveUnitIdx                 = false(1, numel(unit_list));
ctc.ResponsiveUnitIdx(responsive_idx) = true;
ctc.ResponsivenessDetail              = detail;

fprintf('\nTotal inhibitory candidates: %i / %i units (%.1f%%)\n', ...
    sum(ctc.ResponsiveUnitIdx), numel(ctc.UnitList), ...
    100 * mean(ctc.ResponsiveUnitIdx));
end