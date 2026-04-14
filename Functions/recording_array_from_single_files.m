function rec_array = recording_array_from_single_files(sorting_path_list, rm_con, parallel)
arguments
    sorting_path_list string
    rm_con logical = true %If true, removes connectivity data to save memory
    parallel logical = true
end
rec_array = cell(1,length(sorting_path_list));

if parallel
    parfor iPath = 1:length(sorting_path_list)
        mearec = fullfile(sorting_path_list(iPath), "MEArecording.mat");
        if exist(mearec,'file')
            opened = load(mearec,"obj");
            %parseMetadata(opened, metadata)
            if rm_con
                opened.obj.Connectivity = [];
            end
            rec_array{iPath} = opened.obj;
        end
    end
else
    for iPath = 1:length(sorting_path_list)
        mearec = fullfile(sorting_path_list(iPath), "MEArecording.mat");
        if exist(mearec,'file')
            opened = load(mearec,"obj");
            %parseMetadata(opened.obj, opened.obj.Metadata)
            if rm_con
                opened.obj.Connectivity = [];
            end
            rec_array{iPath} = opened.obj;
        end
    end
end


rec_array = [rec_array{:}];
empty_idx = arrayfun(@(x) isempty(x.Units),rec_array); %Recordings with no units
rec_array(empty_idx) = []; %Remove them

fprintf("%i/%i objects retrieved\n", length(rec_array),length(sorting_path_list))