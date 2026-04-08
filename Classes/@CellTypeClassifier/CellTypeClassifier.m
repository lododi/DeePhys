classdef CellTypeClassifier < handle
    % CELLTYPECLASSIFIER  Supervised cell-type classification pipeline for MEA perturbation experiments.
    %
    % Identifies inhibitory (interneuron) vs excitatory neurons in MEA cultures
    % using a bootstrap firing-rate response test followed by supervised UMAP classification.
    %
    % USAGE:
    %   ctc = CellTypeClassifier(rg, params);
    %   ctc.identifyResponsiveUnits();    % bootstrap pre/post firing rate comparison
    %   ctc.generateTrainLabels();        % UMAP embedding + isolation forest label generation
    %   ctc.classifyUnits();              % supervised UMAP projection and classification
    %   labels = ctc.UnitLabels;          % 1 = excitatory, 2 = inhibitory, NaN = unclassified

    properties
        RecordingGroup          % RecordingGroup with .Cultures already set
        Parameters              % struct — merged from returnDefaultParams + user overrides
        UnitList                % (1 x N) Unit array — all units across all cultures
        ResponsiveUnitIdx       % (1 x N) logical — units identified as responsive (inhibitory candidates) by bootstrap
        ResponsivenessDetail    % struct with per-unit filter outcomes and dose-response curves
        TrainLabels             % struct: .sorted_train_ids, .sorted_y_train, .umap_train_idx, .umap_test_idx
        UnitLabels              % (1 x N) double — 1=excitatory, 2=inhibitory, NaN=unclassified
        NormStats               % struct: .mu .sigma .scale .nan_cols — normalization handoff
        UMAP                    % trained UMAP model returned by run_umap
        Reduction = struct('Unsupervised', [], 'Train', [], 'Test', [], 'External', [])
        HarmonizedWaveforms     % (N_samples × N_units) processed waveforms (upsampled, aligned, trimmed)
        HarmonizedACGs          % (N_bins × N_units)   ACGs at Harmonization params
        HarmonizedSR            % scalar — waveform sampling rate after harmonization
        CachedExtraction        % struct caching extractUnitWaveformsAndACGs output (see getOrExtract)
        % struct storing UMAP embeddings:
        %   .Unsupervised — (N x D) from generateTrainLabels (all units, unsupervised)
        %   .Train        — (N_train x D) from classifyUnits (supervised training embedding)
        %   .Test         — (N_test  x D) from classifyUnits (supervised test projection)
        %   .External     — (N_ext  x D) from classifyExternalUnits
    end

    methods
        function ctc = CellTypeClassifier(rg, parameters)
            arguments
                rg RecordingGroup
                parameters struct = struct()
            end
            ctc.RecordingGroup = rg;
            ctc.Parameters = parseStructParameters(ctc.returnDefaultParams(), parameters);
            % UnitList is populated by identifyResponsiveUnits(), which also sets
            % the cumulative culture offsets needed for ResponsiveUnitIdx.
        end


        function [wf, acg, sr] = getOrExtract(ctc, unit_list)
            % getOrExtract  Return cached extraction or compute and cache.
            %   Avoids redundant FullACG recomputation across generateTrainLabels,
            %   classifyUnits, etc. Cache is invalidated when unit count or
            %   harmonization parameters change.
            ph = ctc.Parameters.Harmonization;
            cache = ctc.CachedExtraction;
            if ~isempty(cache) ...
                    && cache.N == numel(unit_list) ...
                    && abs(cache.ACGBinSize - ph.ACGBinSize) < 1e-12 ...
                    && abs(cache.ACGLag - ph.ACGLag) < 1e-12 ...
                    && cache.WaveformTargetSR == ph.WaveformTargetSamplingRate ...
                    && cache.ACGSource == ph.ACGSource
                wf  = cache.wf;
                acg = cache.acg;
                sr  = cache.sr;
                fprintf('Using cached waveform/ACG extraction (%d units)\n', cache.N);
                return
            end
            [wf, acg, sr] = extractUnitWaveformsAndACGs(ctc, unit_list);
            ctc.CachedExtraction = struct( ...
                'wf', wf, 'acg', acg, 'sr', sr, ...
                'N', numel(unit_list), ...
                'ACGBinSize', ph.ACGBinSize, ...
                'ACGLag', ph.ACGLag, ...
                'WaveformTargetSR', ph.WaveformTargetSamplingRate, ...
                'ACGSource', ph.ACGSource);
        end
    end

    methods (Static)
        function defaultParams = returnDefaultParams()
            defaultParams.Bootstrap.BinSize            = 20;
            defaultParams.Bootstrap.NIter              = 1000;
            defaultParams.Bootstrap.DoseWindows = ...
                [0 1200; 1200 2400; 2400 3600; 3600 4800; 4800 6000; 6000 7200];
            defaultParams.Bootstrap.DoseValues  = [0, 0.01, 0.1, 1, 10, 100];  % nM
            defaultParams.Bootstrap.MonotonicityThreshold = 0.8;   % Spearman rho
            defaultParams.Bootstrap.MinFoldChange         = 1.5;   % 100nM vs baseline
            defaultParams.Bootstrap.ConfirmationAlpha     = 1e-3;  % relaxed for confirmation

            defaultParams.UMAP.NDims                = 10;       % unsupervised UMAP output dimensions
            defaultParams.UMAP.NNeighbors           = 50;       % unsupervised UMAP n_neighbors (local vs global structure balance)
            defaultParams.UMAP.MinDist              = 0.1;      % minimum distance between embedded points (lower = tighter clusters)
            defaultParams.UMAP.Spread               = 1.0;      % scale of embedded points, used with MinDist
            defaultParams.UMAP.SupervisedNDims      = 2;        % supervised UMAP output dimensions
            defaultParams.UMAP.SupervisedNNeighbors = 100;      % supervised UMAP n_neighbors — must be < minority_class_size/5
            defaultParams.UMAP.UnitFeatures         = ["FullACG","ReferenceWaveform"]; % feature groups to extract per unit
            defaultParams.UMAP.NormalizationVar     = "ChipID"; % metadata field for within-group z-score normalisation
            defaultParams.UMAP.GroupingVar          = "Concentration"; % metadata field selecting which recordings enter UMAP
            defaultParams.UMAP.GroupingValues       = 0;        % value(s) of GroupingVar to include (0 = baseline only)
            defaultParams.UMAP.TrainingCultureIdx   = [];       % culture indices for label generation (empty = all cultures)
            defaultParams.UMAP.TemplateDir          = tempdir;  % directory for UMAP template files used in test projection
            defaultParams.UMAP.TargetWeight         = 0.5;      % supervised label influence (0-1); values >0.5 risk ring artifacts
            defaultParams.UMAP.Metric               = 'euclidean'; % UMAP distance metric — must match pdist2 metric in generateTrainLabels
            defaultParams.UMAP.MatchSupervisors     = '3';      % label propagation mode for test projection (3 = most aggressive)
            defaultParams.UMAP.ACGWeight            = 1.0;      % relative contribution of ACG feature group (dimension-normalised)
            defaultParams.UMAP.WaveformWeight       = 1.0;      % relative contribution of waveform feature group (dimension-normalised)

            defaultParams.OutlierDetection.ContaminationFraction = 0.5;
            defaultParams.OutlierDetection.NObsPerLearner        = 50;
            defaultParams.OutlierDetection.DistancePercentile    = 80;
            defaultParams.OutlierDetection.CounterexampleRatio   = 2;   % excitatory:inhibitory sampling ratio

            defaultParams.RNGSeed = 42;  % seed for reproducible label generation & classification

            % Harmonization parameters — shared target spec for DeePhys and external data.
            % buildFeatureMatrix and extractUnitWaveformsAndACGs both read these to
            % ensure DeePhys and external units are projected into the same feature space.
            % Change these to harmonize to a different ACG resolution or waveform rate.
            defaultParams.Harmonization.WaveformTargetSamplingRate = 120000;  % Hz — output rate after interpolation
            defaultParams.Harmonization.WaveformPreTrough          = 1.0;   % ms before trough to include
            defaultParams.Harmonization.WaveformPostTrough         = 1.0;   % ms after  trough to include
            defaultParams.Harmonization.WaveformEdgeMode           = "trim"; % "zero" | "edge" | "trim"
            % zero — zero-pad regions outside available input
            % edge — hold the boundary value (nearest-neighbor extrapolation)
            % trim — shrink output window to the shortest available span across all units
            defaultParams.Harmonization.ACGBinSize                 = 0.0005; % s — ACG bin width (0.5 ms)
            defaultParams.Harmonization.ACGLag                     = 0.1;    % s — one-sided ACG lag (100 ms → 401 bins)
            defaultParams.Harmonization.ACGSource                  = "FullACG"; % "FullACG" or "ACG" — which DeePhys ACG property to use
        end

        function X = applyFeatureWeights(X, feature_groups, p_umap)
            % APPLYFEATUREWEIGHTS  Scale ACG and waveform feature groups by explicit weights.
            %
            % Applies dimension-normalized weights so each group's total L2 contribution
            % is proportional to its weight parameter, independent of group size.
            % With ACGWeight=WaveformWeight=1 this equalizes group contributions.
            %
            % INPUTS:
            %   X             - (N x F) normalized feature matrix
            %   feature_groups - struct with .n_acg and .n_wf
            %   p_umap        - Parameters.UMAP struct with .ACGWeight and .WaveformWeight
            %
            % OUTPUT:
            %   X             - (N x F) reweighted feature matrix

            n_acg = feature_groups.n_acg;
            n_wf  = feature_groups.n_wf;

            % sqrt because UMAP distances are quadratic in feature values
            acg_scale = sqrt(p_umap.ACGWeight / n_acg);
            wf_scale  = sqrt(p_umap.WaveformWeight / n_wf);

            X(:, 1:n_acg)     = X(:, 1:n_acg)     * acg_scale;
            X(:, n_acg+1:end) = X(:, n_acg+1:end) * wf_scale;
        end
    end
end
