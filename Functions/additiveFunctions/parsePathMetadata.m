function meta = parsePathMetadata(ks_path)
% PARSEPATHMETADATA Extract ChipID, Well, and RecordingDate from a
%   sorter_output path of the form:
%       .../<YYMMDD>/<ChipID>/Network/well<NNN>/sorter_output[/segment_<S>]
%   using the 'Network' folder as an anchor (robust to whatever sits
%   above the date folder).
%
%   Well numbering on disk is 0-indexed (well000, well001, ...); the
%   returned ChipWellKey uses 1-indexed numbering to match the
%   convention used in the Excel metadata sheet (ChipID_01, ChipID_02, ...).
%
%   Returns a struct with fields:
%       Well          - string, e.g. "well008"
%       WellNum0      - double, 0-indexed well number (as on disk)
%       WellNum1      - double, 1-indexed well number (as in the sheet)
%       ChipID        - string, e.g. "T002523"
%       RecordingDate - string, ISO 'YYYY-MM-DD'
%       ChipWellKey   - string, e.g. "T002523_09"

    parts = strsplit(ks_path, filesep);
    net_idx = find(strcmpi(parts, 'Network'), 1, 'last');
    if isempty(net_idx)
        error('parsePathMetadata:noAnchor', 'No "Network" folder found in: %s', ks_path);
    end

    well_tok = parts{net_idx + 1};
    chip_tok = parts{net_idx - 1};
    date_tok = parts{net_idx - 2};

    well_match = regexp(well_tok, '^well(\d+)$', 'tokens', 'once');
    if isempty(well_match)
        error('parsePathMetadata:badWell', 'Unexpected well folder name: %s', well_tok);
    end

    meta.WellNum0      = str2double(well_match{1});
    meta.WellNum1      = meta.WellNum0 + 1;
    meta.Well          = string(well_tok);
    meta.ChipID        = string(chip_tok);
    meta.RecordingDate = yymmddToISO(date_tok);
    meta.ChipWellKey   = sprintf('%s_%02d', chip_tok, meta.WellNum1);
end