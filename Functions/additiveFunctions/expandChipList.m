function longMeta = expandChipList(xlmeta)
% EXPANDCHIPLIST Explode each Excel row's comma-separated ChipList cell
%   (e.g. "T002523_01,T002523_02, T002523_03") into one row per
%   ChipWellKey token, duplicating all other metadata columns onto
%   each exploded row.
%
%   Requires xlmeta to have a 'ChipList' column (string). Returns a
%   table with the same columns as xlmeta minus 'ChipList', plus
%   'ChipWellKey' (string).

    metaCols = setdiff(xlmeta.Properties.VariableNames, {'ChipList'}, 'stable');
    outRows = {};
    keys = {};

    for r = 1:height(xlmeta)
        raw = xlmeta.ChipList(r);
        tokens = strtrim(strsplit(raw, ','));
        tokens = tokens(strlength(string(tokens)) > 0);
        for t = 1:numel(tokens)
            keys{end+1,1} = string(tokens{t}); %#ok<AGROW>
            outRows(end+1,:) = table2cell(xlmeta(r, metaCols)); %#ok<AGROW>
        end
    end

    longMeta = cell2table(outRows, 'VariableNames', metaCols);
    longMeta.ChipWellKey = string(keys);
end