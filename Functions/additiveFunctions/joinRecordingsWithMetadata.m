function jobs = joinRecordingsWithMetadata(recTable, xlmeta)
% JOINRECORDINGSWITHMETADATA Match each discovered recording (from
%   discoverRecordings) to its Excel metadata row.
%
%   Match key is (ChipWellKey, RecordingDate) together, since the same
%   chip/well typically recurs across multiple recording timepoints
%   (different DIV) with one Excel row per timepoint.
%
%   Within a matched (ChipWellKey, RecordingDate) group:
%     1. If the recording is a segment AND a metadata row exists with
%        SegmentsDiffer == true and SegmentIndex matching this
%        recording's segment index -> use that row (segment-specific).
%     2. Otherwise, use a row with SegmentsDiffer == false/blank
%        (applies uniformly to the whole chip/well, all segments alike).
%     3. Otherwise -> unmatched (Matched = false), metadata columns set
%        to <missing>.
%
%   xlmeta must already have ChipList exploded via expandChipList is NOT
%   required here -- this function calls expandChipList internally.
%   xlmeta must contain: ChipList, RecordingDate (ISO string), and
%   optionally SegmentsDiffer (logical) / SegmentIndex (double, 0-indexed).
%
%   Errors if the sheet contains duplicate (ChipWellKey, RecordingDate,
%   SegmentsDiffer, SegmentIndex) combinations, since that would make
%   the match ambiguous.
%
%   Returns recTable with all xlmeta metadata columns appended, plus a
%   logical 'Matched' column.

    if ~ismember('SegmentsDiffer', xlmeta.Properties.VariableNames)
        xlmeta.SegmentsDiffer = false(height(xlmeta), 1);
    end
    if ~ismember('SegmentIndex', xlmeta.Properties.VariableNames)
        xlmeta.SegmentIndex = nan(height(xlmeta), 1);
    end
    xlmeta.SegmentsDiffer = logical(xlmeta.SegmentsDiffer);

    longMeta = expandChipList(xlmeta);

    % Duplicate-key check: fail loudly rather than silently picking one
    keyStr = longMeta.ChipWellKey + "|" + longMeta.RecordingDate + "|" + ...
             string(longMeta.SegmentsDiffer) + "|" + string(longMeta.SegmentIndex);
    [~, ia] = unique(keyStr);
    if numel(ia) < height(longMeta)
        dupIdx = setdiff(1:height(longMeta), ia);
        error('joinRecordingsWithMetadata:duplicateKey', ...
            'Duplicate metadata rows for: %s', ...
            strjoin(unique(longMeta.ChipWellKey(dupIdx)), ', '));
    end

    % Exclude match-key/internal columns, AND any column name already
    % present in recTable (e.g. RecordingDate, which is parsed from the
    % path and used as part of the match key -- keeping it here too would
    % collide on horzcat below since the value is identical by construction).
    reservedCols = [{'ChipWellKey','SegmentsDiffer','SegmentIndex'}, ...
                     recTable.Properties.VariableNames];
    metaCols = setdiff(longMeta.Properties.VariableNames, reservedCols, 'stable');

    n = height(recTable);
    matched = false(n, 1);
    extra = table();

    for i = 1:n
        key  = recTable.ChipWellKey(i);
        date = recTable.RecordingDate(i);
        candRows = longMeta(longMeta.ChipWellKey == key & longMeta.RecordingDate == date, :);
        src = table();

        if recTable.IsSegment(i)
            segSpecific = candRows(candRows.SegmentsDiffer & ...
                candRows.SegmentIndex == recTable.SegmentIndex0(i), :);
            if height(segSpecific) >= 1
                src = segSpecific(1, metaCols);
                matched(i) = true;
            end
        end

        if ~matched(i)
            wholeChip = candRows(~candRows.SegmentsDiffer, :);
            if height(wholeChip) >= 1
                src = wholeChip(1, metaCols);
                matched(i) = true;
            end
        end

        if ~matched(i)
            src = longMeta(1, metaCols);
            src{1,:} = missing;
        end

        extra = [extra; src]; %#ok<AGROW>
    end

    jobs = [recTable, extra];
    jobs.Matched = matched;
end