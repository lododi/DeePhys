% This tutorial shows how to perform the feature extraction from a list of
% sortings that are stored in a phy-compatible format
% (https://github.com/cortex-lab/phy).


% This tutorial assumes that you are working with the dataset provided at
% Zenodo (...). You will need to change paths etc. accordingly to work with
% your own data.

% This step will take some time, so if you want to play around with the
% data immediately you can make use of the existing files and skip to the
% next tutorials.

% I would recommend to copy these scripts into a separate folder, where you
% can change them without being overwritten if you pull a newer version of
% DeePhys.

% Specify the DeePhys root path
addpath(genpath("/home/lododi/git/DeePhys"))

% Or if your working directory is 1_FeatureExtraction we can also use a
% relative path
addpath(genpath("../../Functions"))
addpath(genpath("../../Classes"))
addpath(genpath("../../Toolboxes"))
%% %% Generate sorting_path_list 
% Next we need to specify the list of sortings for which we want to perform
% the feature extraction. Ideally, you have stored your data in a way that
% can make use of the 'generate_sorting_path_list' function. 

% The idea of this function is to find all folders that contain spike sorting
% output from a specified 'root_path'.
root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Maxtwo/phornauer/Torsten_2/241016";  %"/net/bs-filesvr02/export/group/hierlemann/recordings/Maxtwo/phornauer/Torsten_2/241016"; % INSERT YOUR OWN PATH HERE

% The path_logic is then used to selectively include certain subfolders.
% Each element in the cell array hereby represents one level in the folder
% structure (i.e., one subfolder). 
% Often, you will have similar folders with only slight differences
% (e.g. recording dates or chip ids), which can be found in one
% go by using the wildcard character *. 
path_logic = {'T0*','Network','w*','sorter*'}; 

sorting_path_list = generate_sorting_path_list(root_path, path_logic);
fprintf("Generated %i sorting paths\n",length(sorting_path_list))

%% Generate sorting_path_list 
% Next we need to specify the list of sortings for which we want to perform
% the feature extraction. Ideally, you have stored your data in a way that
% can make use of the 'generate_sorting_path_list' function. 

% The idea of this function is to find all folders that contain spike sorting
% output from a specified 'root_path'.
root_path = "/links/groups/hierlemann/Projects/ADHD/AnalyzedData/datastore/.cache/";  %"/net/bs-filesvr02/export/group/hierlemann/recordings/Maxtwo/phornauer/Torsten_2/241016"; % INSERT YOUR OWN PATH HERE

% The path_logic is then used to selectively include certain subfolders.
% Each element in the cell array hereby represents one level in the folder
% structure (i.e., one subfolder). 
% Often, you will have similar folders with only slight differences
% (e.g. recording dates or chip ids), which can be found in one
% go by using the wildcard character *. 
path_logic = {'*','data','w*','KS-2-*','spike_times.npy'}; 

sorting_path_list = generate_sorting_path_list(root_path, path_logic);
fprintf("Generated %i sorting paths\n",length(sorting_path_list))

% If your data is stored differently, you can also manually provide the
% sorting_path_list as a list of strings.

%% Extract cache path and uniqueIDs
% for i = 1:length(path_list)
%     current_path = path_list{i};
% 
%     % Split the full path into its individual components
%     path_parts = strsplit(current_path, filesep);
% 
%     % Filter out any empty strings that might result from splitting
%     path_parts = path_parts(~cellfun('isempty', path_parts));
% 
%     % Check if the path has at least 4 components
%     if length(path_parts) >= 4
%         % Get the last four components from the end of the list
%         last_four_components = path_parts(end-3:end);
% 
%         % Join the last four components back into a path string
%         trimmed_paths{i} = fullfile(last_four_components{:});
%     else
%         % If the path is too short, keep the original path or handle as needed
%         trimmed_paths{i} = current_path;
%     end
% end

%% --- Step 1: Extract the fourth-shallowest folder from each path ---
% for i = 1:length(sorting_path_list)
%     % Parse the full path and split it into its component parts
%     [folder_path, ~, ~] = fileparts(sorting_path_list{i});
% 
%     % Use filesep to split the folder path, which works on both Windows and Unix
%     path_parts = strsplit(folder_path, filesep);
% 
%     % Filter out any empty strings that might result from trailing separators
%     path_parts = path_parts(~cellfun('isempty', path_parts));
% 
%     % Check if the path has at least four folders
%     if length(path_parts) >= 3
%         % Get the fourth folder from the end of the list
%         IDs_unique{i} = path_parts{end-2};
%     else
%         % Use a placeholder for paths that are too short
%         IDs_unique{i} = 'TooShort'; 
%     end
% end
% 
% % --- Step 2: Count the occurrences of each unique folder ---
% % Find unique folder names and their indices
% [unique_folders, ~, folder_indices] = unique(IDs_unique);
% 
% % Use accumarray to count how many times each unique folder appears
% folder_counts = accumarray(folder_indices, 1);
% 
% % --- Step 3: Check for duplicates (counts > 1) ---
% % Find the indices where the count is greater than 1
% duplicate_indices = find(folder_counts > 1);
% 
% % Check if the list of duplicate indices is not empty
% if ~isempty(duplicate_indices)
%     % Get the names of the folders that have duplicates
%     duplicate_folders = unique_folders(duplicate_indices);
%     fprintf('The following fourth-shallowest folders have more than one path:\n');
%     disp(duplicate_folders);
% else
%     fprintf('No fourth-shallowest folders have more than one path in the list.\n');
% end
%% Set parameter values
% Next, we need to specify the relevant parameters. The default values
% should be a pretty good starting point, but you can adjust them here. To
% get a description of all parameter values look at the
% 'returnDefaultParams' function in the 'MEArecording' class.

emptyRec = MEArecording();
params = emptyRec.returnDefaultParams();

% Since we extract a variety of different features, we have a lot of
% parameters. Several of those are, however, the defaults used by the
% authors and should not need changing. Feel free to tinker with those as
% well, also if they are not explicitly listed here.

params.QC.Amplitude             = [];
params.QC.FiringRate            = [0.01 100];
params.QC.Axon                  = 0.8;
params.QC.Noise                 = 8;
params.QC.N_Units               = 10;
params.Analyses.SingleCell      = 1;
params.Analyses.Regularity      = 1;
params.Analyses.Catch22         = 1;
params.Analyses.Bursts          = 0; %For this tutorial we dont perform burst detection, as the cultures did not display network-wide synchronization
% As the automated burst detection tries to find optimal parameters for the
% burst detection even if none are present, performing burst analyses
% without bursts being present MIGHT result in errors, so be careful!
params.Analyses.Connectivity    = ["STTC","CCG"]; % Takes a LONG time if you have a lot of units
params.Outlier.Method           = []; %No outlier removal
params.Save.Flag                = 1; %Save individual MEArecordings to prevent data loss if the execution is interrupted
params.Save.Overwrite           = true; 

%% Set up metadata
metadata = struct();
[cacheIDs,uniqueIDs] = generate_IDs_from_sorting_list(sorting_path_list);
% This file provides the metadata necessary to identify/pool data later
% during analysis. The files in the 'metadata' folder can serve as a
% template, where you can add any other metadata information as a new
% column. 

% DONT CHANGE THE PlatingDate and Chip_IDs column headers when creating
% your own metadata lookup file! All other names can be chosen freely.
%metadata.LookupPath = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Maxtwo/phornauer/Torsten_2/metadata/Torsten_2.xlsx"; % INSERT YOUR OWN PATH HERE


[iuniqueID, ~, idx] = unique(uniqueIDs);
counts = histcounts(idx, 'BinMethod', 'integers');
repeatedValues = iuniqueID(counts > 1);
duplicateIndices = ismember(uniqueIDs, repeatedValues);
correspondingValues = cacheIDs(duplicateIndices);
if isempty(duplicateIndices)
    disp("No Duplicates!")
else
    disp(correspondingValues');
    disp(find(duplicateIndices)');
end

% Possible filtering out/selection of duplicates

uniqueIDs = uniqueIDs(2:end);

metadata.UniqueID = uniqueIDs;

metadata.Variables = ["RecordingDate","PlateID","WellID","PlatingDate","DIV","CellLine","Condition","ChipID"];

% project needs to be active

% DB = recordingsDB();
% 
% DB.filter('Main','UniqueID',fourth_shallowest_folders);
% metadata.RecordingDate = string(DB.Selected.Main.Date,'yyMMdd');
% metadata.PlateID = DB.Selected.Main.PlateID;
% metadata.WellID = compose("%02d", DB.Selected.Main.WellID);
% metadata.PlatingDate = string(DB.Selected.Main.PlatingDate,'yyMMdd');
% metadata.DIV = DB.Selected.Main.DIV;
% metadata.Patterning = DB.Selected.Main.CellLine;
% metadata.Density = DB.Selected.Main.Condition;
% metadata.Chip_IDs = DB.Selected.Main.ChipID;

metadata.LookupPath = "/links/groups/hierlemann/Projects/ADHD/Recordings/metadata/Deephys_45678_0904.xlsx"; % INSERT YOUR OWN PATH HERE

% The 'path_part_idx' variable exists to extract important metadata from the sorting path.
% Similarly to the previous function, this assumes that the sortings were
% saved in a systematic fashion and that the file paths contain the
% indicated information. The values should indicate the inidices after
% splitting the sorting_path at every '/'. 

% To give you an example:
%path_parts = strsplit('/mp/nfs4krb5/bs-nas02/bsse_group_hierlemann/Projects/ADHD/AnalyzedData/datastore/.cache/70000121/data/well_1/KS-2-250820-1149-70000121','/')
%path_parts = strsplit("/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Maxtwo/phornauer/Torsten_2/241023/T002523/Network/well000/sorter_output",'/')
% Now we only have to select the correct indices, in this example this would
% be 6 for the recording date, 7 for the plate ID and 9 for the well ID.

%% 
% If this approach is not applicable for your data, you can also specify these
% information separately for each recording by writing your own function or by doing it after the feature extraction.
% However, I would strongly recommend to store your data in a way that
% allows for this approach to work (also for a general ease of use).

%metadata.PathPartIndex = [11,12,14];%insert indices you chose above [recording_date, plate_id, well_id]

metadata.min_N_units = params.QC.N_Units; %Minimum number of spike-sorted units to perform feature extraction
parallel = true; %Perform feature extraction in parallel, HIGHLY RECOMMENDED IF YOUR SERVER/MACHINE HAS ENOUGH RAM

%% Run full loop

% This function performs the actual feature extraction and returns a list
% of sortings for which it did not work. You can then run those again
% without parallelization to obtain proper error messages/use breakpoints
% to find out where things went wrong.

failed_sortings = generate_MEArecordings_from_sorting_list(sorting_path_list, metadata, params, parallel);

%% (Hopefully) Optional: Find source of failed feature extractions
check_failed_analysis(sorting_path_list, metadata, params);

% This is just a convenience functions that iterates through all failed
% sortings to narrow down the source of the error.
% Play around with inserting breakpoints in the MEArecording code and
% hopefully you'll figure out the reason.