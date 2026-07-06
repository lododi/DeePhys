function recTable = discoverRecordings(base_dir, excludeFolders)
% DISCOVERRECORDINGS Find every recording under base_dir and parse its
%   identity (ChipID, Well, RecordingDate) from the folder path.
%
%   Concatenated recordings: a sorter_output folder containing segment_0,
%   segment_1, ... subfolders is treated as a PARENT. Each segment_N
%   subfolder becomes its own job, with ParentPath pointing at the
%   parent sorter_output (used for computeParentFeatures / shared ACG
%   cache) and SegmentIndex0 giving the 0-indexed segment number, matching
%   the "segment_N" folder naming and the Excel sheet's SegmentIndex
%   convention.
%
%   Non-segmented recordings (a sorter_output with spike_times.npy
%   directly inside, no segment_* subfolders) become a single job with
%   IsSegment = false and SegmentIndex0 = NaN.
%
%   Only sorter_output folders that sit under a 'Network' folder are
%   parsed by parsePathMetadata (ChipID/Well/Date extraction assumes the
%   Network/well<NNN> layout). Any sorter_output folder whose path
%   contains one of excludeFolders as a path component (default:
%   {'AxonTracking'}) is skipped, since those use a different recording
%   layout not currently supported by this pipeline. A summary of how
%   many were skipped, and for what reason, is printed.
%
%   Returns a table with columns:
%       ks_path, save_dir, ChipID, Well, ChipWellKey, RecordingDate,
%       IsSegment, ParentPath, SegmentIndex0

    if nargin < 2 || isempty(excludeFolders)
        excludeFolders = {'AxonTracking'};
    end

    ks_dirs = findSorterOutputDirs(base_dir);
    rows = {};
    n_excluded = 0;
    n_no_network = 0;

    for k = 1:numel(ks_dirs)
        d = ks_dirs{k};
        parts = strsplit(d, filesep);

        if any(ismember(excludeFolders, parts))
            n_excluded = n_excluded + 1;
            continue
        end
        if ~any(strcmpi(parts, 'Network'))
            n_no_network = n_no_network + 1;
            continue
        end
        seg_dirs = dir(fullfile(d, 'segment_*'));
        seg_dirs = seg_dirs([seg_dirs.isdir]);

        if ~isempty(seg_dirs)
            meta = parsePathMetadata(d);
            for s = 1:numel(seg_dirs)
                child_path = fullfile(seg_dirs(s).folder, seg_dirs(s).name);
                if ~isfile(fullfile(child_path, 'spike_times.npy'))
                    continue
                end
                seg_match = regexp(seg_dirs(s).name, '^segment_(\d+)$', 'tokens', 'once');
                if isempty(seg_match)
                    continue
                end
                seg_idx = str2double(seg_match{1});
                rows(end+1,:) = {child_path, child_path, meta.ChipID, meta.Well, ...
                                  meta.ChipWellKey, meta.RecordingDate, true, ...
                                  string(d), seg_idx}; %#ok<AGROW>
            end
        elseif isfile(fullfile(d, 'spike_times.npy'))
            meta = parsePathMetadata(d);
            rows(end+1,:) = {d, d, meta.ChipID, meta.Well, meta.ChipWellKey, ...
                              meta.RecordingDate, false, "", NaN}; %#ok<AGROW>
        end
    end

    if isempty(rows)
        recTable = table('Size', [0 9], ...
            'VariableTypes', {'string','string','string','string','string','string','logical','string','double'}, ...
            'VariableNames', {'ks_path','save_dir','ChipID','Well','ChipWellKey','RecordingDate', ...
                               'IsSegment','ParentPath','SegmentIndex0'});
    else
        recTable = cell2table(rows, 'VariableNames', ...
            {'ks_path','save_dir','ChipID','Well','ChipWellKey','RecordingDate', ...
             'IsSegment','ParentPath','SegmentIndex0'});
        recTable.ks_path  = string(recTable.ks_path);
        recTable.save_dir = string(recTable.save_dir);
    end

    if n_excluded > 0
        fprintf('discoverRecordings: skipped %d sorter_output folder(s) under excluded folders (%s).\n', ...
            n_excluded, strjoin(excludeFolders, ', '));
    end
    if n_no_network > 0
        fprintf('discoverRecordings: skipped %d sorter_output folder(s) with no "Network" ancestor.\n', ...
            n_no_network);
    end
end