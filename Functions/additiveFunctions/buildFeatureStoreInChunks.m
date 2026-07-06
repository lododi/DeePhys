function fs = buildFeatureStoreInChunks(proc_paths, out_dir, chunk_size)
% BUILDFEATURESTOREINCHUNKS Assemble a FeatureStore from many saved
%   RecordingProcessor .mat files without ever holding more than
%   chunk_size full RecordingProcessor objects in memory at once.
%
%   Loading everything via RecordingProcessor.loadMany + a single
%   FeatureStore.fromProcessors call keeps every processor's full data
%   (spike times, waveforms, CCG matrices, spatial data, etc.) resident
%   in RAM simultaneously. For a few hundred recordings this can exceed
%   available memory and cause a hard crash. This function instead:
%
%     1. Loads recordings in chunks of chunk_size via
%        RecordingProcessor.loadMany, builds a FeatureStore for that
%        chunk, saves it to disk, then clears the chunk's processors and
%        FeatureStore from memory before moving to the next chunk.
%     2. Reloads all the (much lighter, table-only) chunk FeatureStores
%        and vertically concatenates their UnitTable, RecordingTable,
%        and MetadataTable into one combined FeatureStore.
%
%   Inputs:
%       proc_paths  - string array of RecordingProcessor.mat paths
%       out_dir     - directory to write intermediate chunk FeatureStore
%                     files (FeatureStore_chunk_NNN.mat)
%       chunk_size  - number of RecordingProcessors to hold in memory at
%                     once (default 20; lower this if crashes persist,
%                     raise it if memory allows for faster runs)
%
%   Returns the combined FeatureStore. Intermediate chunk files are left
%   in out_dir (useful for resuming / inspection) and are NOT deleted.
%
%   NOTE: this function constructs the final FeatureStore by setting
%   UnitTable/RecordingTable/MetadataTable directly on a default-
%   constructed FeatureStore(). If your FeatureStore class does not
%   support a parameterless constructor or does not expose these as
%   settable properties, the final combination step will error -- in
%   that case use the returned combinedTables struct saved alongside
%   the chunk files (FeatureStore_tables_combined.mat) and adapt the
%   reconstruction to your class's actual API.

    if nargin < 3 || isempty(chunk_size)
        chunk_size = 20;
    end

    proc_paths = proc_paths(:);
    n = numel(proc_paths);
    n_chunks = ceil(n / chunk_size);
    chunk_files = strings(n_chunks, 1);

    fprintf('buildFeatureStoreInChunks: %d recordings in %d chunks of up to %d.\n', ...
        n, n_chunks, chunk_size);

    %% Phase 1: process each chunk, save, then release memory
    for c = 1:n_chunks
        idx0 = (c-1)*chunk_size + 1;
        idx1 = min(c*chunk_size, n);
        chunk_paths = proc_paths(idx0:idx1);

        fprintf('[chunk %d/%d] loading %d recordings...\n', c, n_chunks, numel(chunk_paths));

        procs_chunk = RecordingProcessor.loadMany(chunk_paths);
        fs_chunk = FeatureStore.fromProcessors(procs_chunk);

        chunk_file = fullfile(out_dir, sprintf('FeatureStore_chunk_%03d.mat', c));
        fs_chunk.save(chunk_file);
        chunk_files(c) = chunk_file;

        fprintf('[chunk %d/%d] saved %s (%d units, %d recordings)\n', ...
            c, n_chunks, chunk_file, height(fs_chunk.UnitTable), height(fs_chunk.RecordingTable));

        clear procs_chunk fs_chunk
    end

    %% Phase 2: reload chunk FeatureStores (lightweight) and combine
    fprintf('Combining %d chunk FeatureStores...\n', n_chunks);

    unitTables = cell(n_chunks, 1);
    recTables  = cell(n_chunks, 1);
    metaTables = cell(n_chunks, 1);

    for c = 1:n_chunks
        fs_c = FeatureStore.load(chunk_files(c));
        unitTables{c} = fs_c.UnitTable;
        recTables{c}  = fs_c.RecordingTable;
        metaTables{c} = fs_c.MetadataTable;
        clear fs_c
    end

    UnitTable      = vertcat(unitTables{:});
    RecordingTable = vertcat(recTables{:});
    MetadataTable  = vertcat(metaTables{:});

    % Always save the combined tables directly, independent of whether
    % FeatureStore reconstruction below succeeds -- this is your
    % guaranteed fallback if the class API doesn't match the assumption.
    combinedTablesFile = fullfile(out_dir, 'FeatureStore_tables_combined.mat');
    save(combinedTablesFile, 'UnitTable', 'RecordingTable', 'MetadataTable');
    fprintf('Combined tables saved to %s\n', combinedTablesFile);

    try
        fs = FeatureStore();
        fs.UnitTable = UnitTable;
        fs.RecordingTable = RecordingTable;
        fs.MetadataTable = MetadataTable;
        fprintf('Combined FeatureStore: %d units, %d recordings.\n', ...
            height(fs.UnitTable), height(fs.RecordingTable));
    catch ME
        warning(['Could not reconstruct a combined FeatureStore object directly ' ...
                 '(%s). Combined tables are available in %s -- load them and ' ...
                 'adapt reconstruction to your FeatureStore class API.'], ...
                 ME.message, combinedTablesFile);
        fs = struct('UnitTable', UnitTable, 'RecordingTable', RecordingTable, ...
                     'MetadataTable', MetadataTable);
    end
end