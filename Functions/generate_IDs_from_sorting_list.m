function [cacheIDs,uniqueIDs] = generate_IDs_from_sorting_list(sorting_path_list) %% new LD
arguments
    sorting_path_list string
end
[cacheIDs,uniqueIDs] = deal(cell(1,length(sorting_path_list)));
for i = 1:length(sorting_path_list)
    path_parts = strsplit(sorting_path_list{i}, filesep);
    path_parts = path_parts(~cellfun('isempty', path_parts));
    last_four_components = path_parts(end-3:end);
    cacheIDs{i} = fullfile(last_four_components{:});
    uniqueIDs{i} = path_parts{end-3};
end
end