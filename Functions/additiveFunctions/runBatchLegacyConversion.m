%%
root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Maxtwo/phornauer"; %Root path
path_logic = {'C*','*Week_1','w*','sorter_output','segment_*'}; %Variable parts
ei_path_list = generate_sorting_path_list(root_path, path_logic);
fprintf("Generated %i sorting paths\n",length(ei_path_list))

%%
root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Maxtwo/lododi2/Neuropixel_Dataset/251121"; %Root path
path_logic = {'T*','w*','*'}; %Variable parts
ei_path_list = generate_sorting_path_list(root_path, path_logic);
fprintf("Generated %i sorting paths\n",length(ei_path_list))

%%
for i = 1:length(ei_path_list)
    for k = 0:6
        save_dir = fullfile(ei_path_list(i),"segment_"+k,'test_proc');
        if ~exist(save_dir, 'dir')
            mkdir(save_dir)
        end
        proc = RecordingProcessor.fromLegacyMat(string(fullfile(ei_path_list(i),"segment_"+k,"MEArecording.mat")));
        fprintf('Migrated %d units from %s\n', numel(proc.Units), old_mat);
        fprintf('Status: QC=%s UnitFeatures=%s NetworkFeatures=%s\n', ...
            proc.Status.QC, proc.Status.UnitFeatures, proc.Status.NetworkFeatures);
        proc_file = fullfile(save_dir, 'RecordingProcessor.mat');
        proc.save(proc_file);
    end
end

%%
for i = 1:length(ei_path_list)
    save_dir = fullfile(ei_path_list(i),'test_proc');
    if ~exist(save_dir, 'dir')
            mkdir(save_dir)
    end
    proc = RecordingProcessor.fromLegacyMat(string(fullfile(ei_path_list(i),"MEArecording.mat")));
    fprintf('Migrated %d units from %s\n', numel(proc.Units), string(fullfile(ei_path_list(i),"MEArecording.mat")));
    fprintf('Status: QC=%s UnitFeatures=%s NetworkFeatures=%s\n', ...
        proc.Status.QC, proc.Status.UnitFeatures, proc.Status.NetworkFeatures);
    proc_file = fullfile(save_dir, 'RecordingProcessor.mat');
    proc.save(proc_file);
end
