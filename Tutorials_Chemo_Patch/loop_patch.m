
% Specify the DeePhys root path
addpath(genpath("/home/lododi/git/DeePhys"))

%% Generate sorting_path_list
filepath = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/tgaenswein/PathsSpikeSorted_fullSortings.mat";
sortingPaths = load(filepath);
sortingPaths.sortingPaths(8) = "/links/groups/hierlemann/Projects/EI_class/AnalyzedData/EIP-35-TG0084-B-1/datastore/.cache/EIP-35-TG0084-B-1/Cell_B_positive/KS-2.5-260218-1253";
sortingPathList = sortingPaths.sortingPaths;

% Set parameter values

emptyRec = MEArecording();
params = emptyRec.returnDefaultParams();
params.QC.Amplitude             = [];
params.QC.FiringRate            = [0.01 100];
params.QC.Axon                  = 0.8;
params.QC.Noise                 = 8;
params.QC.N_Units               = 2;
params.Analyses.SingleCell      = 1;
params.Analyses.Regularity      = 1;
params.Analyses.Catch22         = 1;
params.Analyses.Bursts          = 1; %For this tutorial we dont perform burst detection, as the cultures did not display network-wide synchronization
params.Analyses.Connectivity    = ["CCG"]; % Takes a LONG time if you have a lot of units
params.Outlier.Method           = []; %No outlier removal
params.Save.Flag                = 1; %Save individual MEArecordings to prevent data loss if the execution is interrupted
params.Save.Overwrite           = false;

% Set up metadata
metadata = struct();
metadata.RecordingDate = "231205";
metadata.PlateID = "2977";
metadata.WellID = compose("%02d", 1);
metadata.LookupPath = "/links/groups/hierlemann/Projects/ADHD/Recordings/metadata/Deephys_Patch.xlsx";
metadata.min_N_units = params.QC.N_Units; %Minimum number of spike-sorted units to perform feature extraction
parallel = false; %Perform feature extraction in parallel, HIGHLY RECOMMENDED IF YOUR SERVER/MACHINE HAS ENOUGH RAM

% Run full loop

% This function performs the actual feature extraction and returns a list
% of sortings for which it did not work. You can then run those again
% without parallelization to obtain proper error messages/use breakpoints
% to find out where things went wrong.

failed_sortings = generate_MEArecordings_from_sorting_list(sortingPathList, metadata, params, parallel);

%%
burst_peak_times = cell(1,length(sortingPaths.sortingPaths));
burst_peak_times_old = cell(1,length(sortingPaths.sortingPaths));
for iPath = 1:length(sortingPaths.sortingPaths)
    currentPath = sortingPaths.sortingPaths(iPath);
    try
        patch_mearec = load(fullfile(currentPath,"MEArecording.mat")).obj;
    catch
        continue
    end

    burst_peaks = 1:length(patch_mearec.Bursts.T_start);
    for iBurst = 1:length(patch_mearec.Bursts.T_start)
        t_start = patch_mearec.Bursts.T_start(iBurst);
        t_end = patch_mearec.Bursts.T_end(iBurst);
        burst_spikes = patch_mearec.Spikes.Times(t_start<patch_mearec.Spikes.Times&t_end>patch_mearec.Spikes.Times);
        sec_cutout = [t_start t_end];%[4800,6000];%[3600,4800];[0,1200];[2400,3600];
        bin_size = 0.01;
        binned_mat = generate_unit_spk_mat(patch_mearec,bin_size, sec_cutout);
        response.unchanged = 1:length(patch_mearec.Units);
        response.increase = [];
        response.decrease = [];
        activity = create_activity_cutout(binned_mat, bin_size, sec_cutout, response.increase, [response.unchanged; response.decrease]);
        activity_peaks = find(activity.total == max(activity.total));
        activity_peak = activity_peaks(1)/100;
        burst_peaks(iBurst) = t_start+activity_peak;
    end
    burst_peak_times{iPath} = burst_peaks;
    % old method
    sec_cutout = [0 patch_mearec.RecordingInfo.Duration];%[4800,6000];%[3600,4800];[0,1200];[2400,3600];
    bin_size = 0.01;
    binned_mat = generate_unit_spk_mat(patch_mearec, bin_size, sec_cutout);
    activity = create_activity_cutout(binned_mat, bin_size, sec_cutout, response.increase, []);
    response.increase = double((1:length(patch_mearec.Units))');
    threshold = 0.1; smooth_window = 9;
    z_binary = plot_EI_network_activity(activity,threshold, smooth_window);xlim([0 patch_mearec.RecordingInfo.Duration])
    z = plot_EI_network_activity(activity);xlim([0 patch_mearec.RecordingInfo.Duration])
    only_burst_activity = activity.total .* (z_binary'>=0);%(z'==2);

    peak_cutout = 150; pp = 0.1; ph = 0.1;
    [~,locs] = findpeaks(only_burst_activity,'MinPeakHeight',ph,'MinPeakProminence',pp,'MinPeakDistance',6);
    burst_peak_times_old{iPath} = locs/100;
end
