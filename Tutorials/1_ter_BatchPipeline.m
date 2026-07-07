%% Tutorial 1 ter — Excel-Driven Batch Pipeline
%
% A more explicit alternative to 1_bis_DataProcessing.m's bulk pipeline.
% 1_bis relies on generate_spikedata_from_sorting_path_list to match each
% discovered path against the Excel sheet internally. This tutorial
% separates that into inspectable steps — useful when the automatic
% matching in 1_bis isn't enough, e.g. when:
%
%   - Segments of a concatenated recording need different metadata rows
%     (SegmentsDiffer / SegmentIndex columns in the sheet), not just one
%     row applied uniformly to the whole chip/well.
%   - You want a per-recording report of which paths matched, which
%     didn't, and why, before committing to a full (re)processing run.
%   - You want to skip recordings that were already processed in a
%     previous run rather than reprocessing everything from scratch.
%
% Pipeline: discover recordings on disk -> parse identity from each path
% -> join against Excel metadata (by ChipWellKey + RecordingDate) -> run
% RecordingProcessor on unmatched/new recordings only -> assemble a
% FeatureStore from the result via FeatureStore.fromProcessorPaths, which
% reads each recording's lightweight feature sidecar rather than the full
% processor -- safe to call on the whole path list at once even for
% hundreds of recordings (see RecordingProcessor.save/loadForFeatureStore).
%
% Prerequisites: DeePhys on the MATLAB path (run startup.m from repo root).
%
% Fill in the placeholder paths in Section 0 before running.

%% 0  Setup — fill in your paths

base_dir  = '/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Maxtwo/phornauer/EI_iNeurons';
xlsx_path = '/links/groups/hierlemann/Projects/lododi/EI_AdvancedSciences_Philipp/Recordings/DeePhys_Original_EI_241119';
out_dir   = '/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Maxtwo/phornauer/EI_iNeurons';

% Set true to reprocess every matched recording even if a
% RecordingProcessor.mat already exists in its save_dir. Leave false to
% skip already-processed recordings and only run new/failed ones.
force_reprocess = true;

% FeatureStore.fromProcessorPaths (used below) reads each recording's
% lightweight feature sidecar rather than the full processor, so no
% memory-bounded chunking is needed here even for hundreds of recordings.

%% 1  Read Excel metadata

xlmeta = readtable(xlsx_path, 'TextType', 'string', 'Sheet', 1, ...
    'VariableNamingRule', 'preserve');

% Rename to valid, consistent internal names.
% Adjust the left-hand list if your sheet's header text differs.
xlmeta = renamevars(xlmeta, ...
    ["PlatingDate","RecordingDate","DIV","Drug","Patterning","CellLine","E:I Ratio","Chip_IDs"], ...
    ["PlatingDate","RecordingDate","DIV","Drug","Patterning","CellLine","EI_Ratio","ChipList"]);

% Convert 6-digit dates to ISO strings for reliable matching against
% path-parsed dates
xlmeta.PlatingDate   = arrayfun(@yymmddToISO, xlmeta.PlatingDate);
xlmeta.RecordingDate = arrayfun(@yymmddToISO, xlmeta.RecordingDate);

xlmeta.ChipList = strtrim(xlmeta.ChipList);

%% 2  Discover recordings on disk
%
% discoverRecordings walks base_dir for sorter_output folders, splits
% concatenated recordings into one job per segment_N subfolder (each
% pointing back at its parent via ParentPath), and parses ChipID/Well/
% RecordingDate from each path via parsePathMetadata.

recTable = discoverRecordings(base_dir);
fprintf('Discovered %d recordings under %s\n', height(recTable), base_dir);
disp(recTable(:, ["ks_path","ChipID","Well","RecordingDate","IsSegment","ParentPath","SegmentIndex0"]));

%% 3  Join with Excel metadata
%
% joinRecordingsWithMetadata matches each discovered recording to its
% Excel row by (ChipWellKey, RecordingDate), preferring a segment-specific
% row (SegmentsDiffer + matching SegmentIndex) over a whole-chip row.
% Errors loudly on ambiguous duplicate rows rather than picking one silently.

jobs = joinRecordingsWithMetadata(recTable, xlmeta);

fprintf('%d/%d recordings matched to metadata.\n', sum(jobs.Matched), height(jobs));
if any(~jobs.Matched)
    warning('Unmatched recordings (will be skipped):');
    disp(jobs(~jobs.Matched, ["ks_path","ChipWellKey","RecordingDate"]));
end
jobs = jobs(jobs.Matched, :);

%% 4  Detect already-processed recordings

jobs.proc_file = fullfile(jobs.save_dir, 'RecordingProcessor.mat');
jobs.AlreadyProcessed = isfile(jobs.proc_file) & ~force_reprocess;

fprintf('%d/%d recordings already processed (found RecordingProcessor.mat), %d to (re)process.\n', ...
    sum(jobs.AlreadyProcessed), height(jobs), sum(~jobs.AlreadyProcessed));
if any(jobs.AlreadyProcessed)
    disp(jobs(jobs.AlreadyProcessed, ["ks_path","ChipWellKey","RecordingDate","proc_file"]));
end

todo_idx = find(~jobs.AlreadyProcessed);
n = numel(todo_idx);

proc_paths = strings(height(jobs), 1);
status     = strings(height(jobs), 1);

% Pre-fill already-processed entries so they flow straight into Section 6
proc_paths(jobs.AlreadyProcessed) = jobs.proc_file(jobs.AlreadyProcessed);
status(jobs.AlreadyProcessed) = "skipped (already processed)";

%% 5  Run batch (parfor) — only recordings not already processed

metaFieldMap = ["ChipWellKey","RecordingDate","PlatingDate","DIV","Drug", ...
    "Patterning","CellLine","EI_Ratio"];

% Slice out only the rows needing processing so parfor doesn't redo work
% or need to reason about the skip condition inside the loop
todoJobs = jobs(todo_idx, :);
todo_proc_paths = strings(n, 1);
todo_status     = strings(n, 1);

parfor i = 1:n
    try
        metadata = struct();
        for f = metaFieldMap
            if ismember(f, todoJobs.Properties.VariableNames)
                val = todoJobs.(f)(i);
                if iscell(val), val = val{1}; end     % unwrap cell -> underlying value
                if isstring(val), val = char(val); end % normalize string -> char
                metadata.(char(f)) = val;
            end
        end
        if ismember("ChipWellKey", todoJobs.Properties.VariableNames)
            val = todoJobs.ChipWellKey(i);
            if iscell(val), val = val{1}; end
            if isstring(val), val = char(val); end
            metadata.ChipID = val;
        end
        if todoJobs.IsSegment(i)
            metadata.ParentInputPath = char(todoJobs.ParentPath(i));
        end

        sd = SpikeData.fromKilosort(char(todoJobs.ks_path(i)), metadata);

        proc = RecordingProcessor(sd);
        proc.Parameters.QC.Amplitude = [0, 1000];
        proc.runAll();   % QC + unit features + parent features + network
        %   + connectivity + cell-type + spatial

        proc_file = fullfile(char(todoJobs.save_dir(i)), 'RecordingProcessor.mat');
        proc.save(proc_file);

        todo_proc_paths(i) = proc_file;
        todo_status(i) = "ok";
    catch ME
        todo_proc_paths(i) = "";
        todo_status(i) = "failed: " + string(ME.message);
    end
end

proc_paths(todo_idx) = todo_proc_paths;
status(todo_idx) = todo_status;

jobs.Status = status;
jobs.proc_file = proc_paths;   % keep proc_file authoritative post-run
writetable(removevars(jobs, intersect(jobs.Properties.VariableNames, {'AlreadyProcessed'})), ...
    fullfile(out_dir, 'batch_run_log.csv'));
disp(jobs(:, ["ks_path","ChipWellKey","RecordingDate","Status"]));

%% 6  Batch-load and assemble FeatureStore

good  = jobs.Status == "ok" | jobs.Status == "skipped (already processed)";
fprintf('Assembling FeatureStore from %d/%d recordings...\n', sum(good), height(jobs));

fs = FeatureStore.fromProcessorPaths(proc_paths(good));
fs.save(fullfile(out_dir, 'FeatureStore_batch.mat'));

fprintf('FeatureStore: %d units, %d recordings.\n', height(fs.UnitTable), height(fs.RecordingTable));
