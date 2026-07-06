%% Batch Legacy Conversion — using RecordingProcessor.convertMany
%
% Discovers legacy MEArecording.mat files under a root directory and
% converts them to RecordingProcessor format in parallel.
%
% Fill in root_path, path_logic, and save_dir before running.

%% 1  Discover legacy recordings

root_path  = "/path/to/your/data";
path_logic = {'C*', '*', 'w*', 'sorter_output', 'segment_*'};
save_dir   = "/path/to/save/converted";

sorting_paths = generate_sorting_path_list(root_path, path_logic);
fprintf("Found %d sorting paths\n", numel(sorting_paths));

%% 2  Build list of MEArecording.mat paths

mat_paths = fullfile(string(sorting_paths), "MEArecording.mat");
exists    = isfile(mat_paths);
mat_paths = mat_paths(exists);
fprintf("%d MEArecording.mat files found\n", numel(mat_paths));

%% 3  Batch convert using convertMany

proc_paths = RecordingProcessor.convertMany(mat_paths, save_dir);

%% 4  Report results

n_ok     = sum(proc_paths ~= "");
n_failed = sum(proc_paths == "");
fprintf('\nConversion complete: %d succeeded, %d failed.\n', n_ok, n_failed);

if n_failed > 0
    failed_idx = find(proc_paths == "");
    fprintf('\nFailed files:\n');
    for i = 1:numel(failed_idx)
        fprintf('  %s\n', mat_paths(failed_idx(i)));
    end
end

%% 5  Optional: assemble FeatureStore from converted processors

% good_paths = proc_paths(proc_paths ~= "");
% fs = buildFeatureStoreInChunks(good_paths, save_dir, 20);
% fs.save(fullfile(save_dir, 'FeatureStore.mat'));
