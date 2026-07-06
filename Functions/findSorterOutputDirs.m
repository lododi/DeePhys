function dirs = findSorterOutputDirs(base_dir)
% FINDSORTEROUTPUTDIRS Recursively find every folder named 'sorter_output'
%   under base_dir. Returns a cell array of full paths.

    dirs = {};
    listing = dir(base_dir);
    listing = listing([listing.isdir] & ~ismember({listing.name}, {'.', '..'}));

    for i = 1:numel(listing)
        full = fullfile(listing(i).folder, listing(i).name);
        if strcmp(listing(i).name, 'sorter_output')
            dirs{end+1} = full; %#ok<AGROW>
        else
            dirs = [dirs, findSorterOutputDirs(full)]; %#ok<AGROW>
        end
    end
end