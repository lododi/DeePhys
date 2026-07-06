classdef ParentSpikeLoader
% PARENTSPIKELOADER  Shared session cache for parent-recording ACG computation.
%
% Static utility — no instantiation. For each unique combination of parent
% Kilosort directory + ACG parameters, computes a batch CCG once for ALL
% templates in that parent recording, extracts the auto-correlogram diagonal,
% and caches the result in memory.
%
% When multiple child RecordingProcessors share the same parent, this ensures
% the CCG is computed only once and each child simply looks up its own
% TemplateIDs from the returned map.
%
% USAGE:
%   params = struct('BinSize', 0.0005, 'Lag', 0.1);
%   acg_map = ParentSpikeLoader.loadFullACGs('/path/to/parent_ks', 30000, params);
%   acg = acg_map('42');   % ACG for TemplateID 42 (key is string)
%
%   ParentSpikeLoader.clearCache();  % free memory between experiments

    methods (Static)

        function acg_map = loadFullACGs(parent_path, sampling_rate, acg_params)
        % LOADFULLACGS  Return (or compute and cache) full-recording ACGs.
        %
        % INPUTS:
        %   parent_path   - absolute path to parent Kilosort directory
        %   sampling_rate - sampling rate in Hz (from child SpikeData)
        %   acg_params    - struct with fields BinSize (s) and Lag (s)
        %
        % OUTPUT:
        %   acg_map - containers.Map from TemplateID string → (n_bins×1) double ACG
            arguments
                parent_path   (1,1) string
                sampling_rate (1,1) double
                acg_params    (1,1) struct
            end

            cache = ParentSpikeLoader.getCache();
            key   = paramHash(struct( ...
                'ParentPath', char(parent_path), ...
                'BinSize',    acg_params.BinSize, ...
                'Lag',        acg_params.Lag));

            if cache.isKey(key)
                acg_map = cache(key);
                return
            end

            % Load full parent spike train
            npy_times = fullfile(parent_path, 'spike_times.npy');
            npy_units = fullfile(parent_path, 'spike_templates.npy');
            if ~isfile(npy_times) || ~isfile(npy_units)
                error('ParentSpikeLoader:missingFiles', ...
                    'spike_times.npy or spike_templates.npy not found in: %s', parent_path);
            end

            raw_times  = double(readNPY(npy_times));
            raw_units  = double(readNPY(npy_units));
            spike_times = raw_times / sampling_rate;
            spike_units = raw_units + 1;  % 0-indexed → 1-indexed

            % Sort by time (required by CCG)
            [spike_times, sort_idx] = sort(spike_times);
            spike_units = spike_units(sort_idx);

            % Compact local IDs (findgroups preserves ascending template order)
            [unit_ids_local, group_tids] = findgroups(spike_units);

            n_bins = round(2 * acg_params.Lag / acg_params.BinSize) + 1;
            fprintf('ParentSpikeLoader: computing ACG for %d templates in %s ...\n', ...
                numel(group_tids), parent_path);

            % Compute each template's autocorrelogram individually (its own
            % spike train against itself) rather than one all-pairs CCG call.
            % A single call requesting all templates at once returns a full
            % nBins x nGroups x nGroups cross-correlogram tensor -- only the
            % diagonal (autocorrelograms) is ever used below, but the tensor
            % itself scales as O(nGroups^2) and can demand tens of GB for a
            % parent recording with hundreds-to-thousands of templates,
            % crashing MATLAB outright. Per-template calls are O(nGroups).
            acg_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
            for t = 1:numel(group_tids)
                unit_times = spike_times(unit_ids_local == t);
                if numel(unit_times) <= 1
                    raw_acg = zeros(n_bins, 1);
                else
                    ccg_t = CCG(unit_times, ones(size(unit_times)), ...
                        'binSize', acg_params.BinSize, ...
                        'duration', 2 * acg_params.Lag, ...
                        'Fs', 1 / sampling_rate);
                    raw_acg = double(ccg_t(:, 1, 1));
                    if numel(raw_acg) ~= n_bins
                        warning('ParentSpikeLoader:binMismatch', ...
                            'Template %d: expected %d bins, got %d — ACG zeroed.', ...
                            group_tids(t), n_bins, numel(raw_acg));
                        raw_acg = zeros(n_bins, 1);
                    end
                end
                mx = max(raw_acg);
                if mx > 0
                    raw_acg = raw_acg / mx;
                end
                acg_map(num2str(group_tids(t))) = raw_acg;
            end

            cache(key) = acg_map;  %#ok<NASGU>
            fprintf('ParentSpikeLoader: done (%d templates cached).\n', numel(group_tids));
        end

        function clearCache()
        % CLEARCACHE  Remove all cached ACG maps to free memory.
            cache = ParentSpikeLoader.getCache();
            remove(cache, keys(cache));
        end

    end

    methods (Static, Access = private)

        function cache = getCache()
        % Persistent session-scoped cache (containers.Map of ACG maps).
            persistent c
            if isempty(c)
                c = containers.Map('KeyType', 'char', 'ValueType', 'any');
            end
            cache = c;
        end

    end
end
