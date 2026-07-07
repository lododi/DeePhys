function allSD = generate_spikedata_from_sorting_path_list(PATH_LIST, SHEET_PATH)
%GENERATE_SPIKEDATA_FROM_SORTING_PATH_LIST  Batch-load SpikeData objects from a list of qc_output paths.
%
%   allSD = generate_spikedata_from_sorting_path_list(PATH_LIST, SHEET_PATH)
%   builds a SpikeData object for each path in PATH_LIST, pulling metadata
%   from the Excel sheet at SHEET_PATH.
%
%   Each path is expected to follow the layout:
%       .../<RecordingDate>/<Plate>/Network/well###/sorter_output/qc_output
%   where <RecordingDate> is a yymmdd string, <Plate> is the chip-ID prefix
%   (e.g. 'T002523', 'M04248'), and well### is zero-indexed
%   (well000 -> chip suffix _01, well001 -> _02, ...).
%
%   The function reconstructs the chip ID as '<Plate>_<NN>' and matches it
%   against the comma-separated Chip_IDs column of the sheet, using the
%   RecordingDate column as an additional key. Paths with no matching row
%   are skipped; a one-line summary is printed at the end.
%
%   INPUTS
%       PATH_LIST  - cell array or string array of qc_output folder paths.
%       SHEET_PATH - path to the metadata .xlsx file. Expected columns:
%                    PlatingDate, RecordingDate (yymmdd strings),
%                    DIV, Drug, Patterning, CellLine, E_IRatio, Chip_IDs.
%
%   OUTPUT
%       allSD - 1xN cell array of SpikeData objects.

    PATH_LIST = cellstr(PATH_LIST);   % accept string array or cell

    opts = detectImportOptions(SHEET_PATH);
    for col = {'PlatingDate','RecordingDate','Chip_IDs'}
        opts = setvartype(opts, col{1}, 'string');
    end
    T = readtable(SHEET_PATH, opts);

    allSD     = {};
    nOk       = 0;
    nNoMatch  = 0;

    for i = 1:numel(PATH_LIST)
        ks_path = PATH_LIST{i};

        % Parse: .../<recDate>/<plate>/Network/well###/sorter_output[/qc_output|/segment_N]
        tok = regexp(ks_path, ...
            '/(\d{6})/([^/]+)/Network/well(\d{3})/sorter_output(?:/(?:qc_output|segment_\d+))?/?$', ...
            'tokens', 'once');

        if isempty(tok)
            nNoMatch = nNoMatch + 1;
            continue
        end

        recDateRaw = tok{1};                              % '241128'
        plate      = tok{2};                              % 'M04248'
        wellIdx    = str2double(tok{3});                  % 0
        chipSuffix = sprintf('%02d', wellIdx + 1);        % '01'
        chipID     = sprintf('%s_%s', plate, chipSuffix); % 'M04248_01'

        % Find the sheet row: matching RecordingDate AND chipID in Chip_IDs
        dateMask = T.RecordingDate == recDateRaw;
        rowIdx   = [];
        for r = find(dateMask)'
            chipList = strtrim(split(T.Chip_IDs(r), ','));
            if any(chipList == chipID)
                rowIdx = r;
                break
            end
        end

        if isempty(rowIdx)
            nNoMatch = nNoMatch + 1;
            continue
        end

        platDateRaw = char(T.PlatingDate(rowIdx));
        recDateISO  = char(yymmddToISO(recDateRaw));
        platDateISO = char(yymmddToISO(platDateRaw));

        metadata = struct();
        metadata.ChipID        = chipID;
        metadata.PlatingDate   = platDateISO;
        metadata.RecordingDate = recDateISO;
        metadata.DIV           = T.DIV(rowIdx);
        metadata.Mutation      = char(T.Drug(rowIdx));
        metadata.Concentration = 0;
        metadata.EI_Ratio      = char(T.E_IRatio(rowIdx));
        metadata.CellLine      = char(T.CellLine(rowIdx));
        metadata.Patterning    = char(T.Patterning(rowIdx));

        % Detect segment paths and set parent path.
        % Priority: sibling qc_output folder (preferred) > parent dir itself.
        [parent_dir, leaf] = fileparts(ks_path);
        qc_candidate = fullfile(parent_dir, 'qc_output');
        if startsWith(leaf, 'segment_') && isfolder(qc_candidate) && ...
                isfile(fullfile(qc_candidate, 'spike_times.npy'))
            metadata.ParentInputPath = qc_candidate;
        elseif startsWith(leaf, 'segment_') && isfile(fullfile(parent_dir, 'spike_times.npy'))
            metadata.ParentInputPath = parent_dir;
        end

        allSD{end+1} = SpikeData.fromKilosort(ks_path, metadata);
        nOk = nOk + 1;
    end

    fprintf('Loaded %d recordings, %d paths had no matching sheet row\n', ...
            nOk, nNoMatch);
end