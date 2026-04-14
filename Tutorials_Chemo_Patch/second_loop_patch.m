burst_locations = cell(1,7);
tot_patch_rg = RecordingGroup(patch_array);
for sel_culture = 1:length(tot_patch_rg.Recordings)
    n_units = arrayfun(@(x) length(x.Units),patch_array);
    cumsum_units = [0, cumsum(n_units)] + 1;
    rg_params.Selection.Inclusion = {}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
    rg_params.Selection.Exclusion = {}; %Cell array of cell arrays with fieldname + value
    patch_rg = RecordingGroup(patch_array(sel_culture), rg_params);
    
    %% Generate actual spk mat
    sec_cutout = [0,patch_rg.Recordings.RecordingInfo.Duration];
    bin_size = 0.01;
    binned_mat = generate_unit_spk_mat(patch_rg,bin_size,sec_cutout);
    activity = create_activity_cutout(binned_mat, bin_size, sec_cutout, iN_idx(cumsum_units(sel_culture):cumsum_units(sel_culture+1) -1), E_idx(cumsum_units(sel_culture):cumsum_units(sel_culture+1) -1));
    
    %% Set cutoff at first minimum of IE ratio histogram
    figure('Color','w');histogram(activity.ratio(activity.ratio>0),0:0.1:max(activity.ratio));
    
    %%
    threshold = 2; smooth_window = 9;
    z = plot_EI_network_activity(activity); %continuous
    z = plot_EI_network_activity(activity,threshold, smooth_window);%xlim([0 100])
    
    %%
    peak_cutout = 100; pp = 0.1; ph = 0.1; %Peak prominence and peak height, respectively
    [burst_mats,burst_locations{sel_culture}] = generate_EI_burst_cutouts(activity, z, peak_cutout, ph, pp);
end