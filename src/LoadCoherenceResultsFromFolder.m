function [ResultsArray, SubjectIDs, LoadedFiles, SkippedFiles] = LoadCoherenceResultsFromFolder(output_dir, options)
    % ===============================================================================================================
    % LOAD PER-SUBJECT COHERENCE RESULTS SAVED BY script_runCoherenceAnalysis.m, READY FOR PLOTTING
    % ===============================================================================================================
    % script_runCoherenceAnalysis.m (and its ancestors, script_runTestWelchCoherence.m /
    % script_TestWelchCoherence_fewSubjects.m) save each subject's results to its OWN .mat file in
    % output_dir, rather than collecting them into one struct array in the workspace -- see that
    % script's header for the rationale (a crash or a later subject's error never costs you the
    % subjects already done). This function is the reload step: point it at that same output_dir
    % and it hands back an array of just the inner results structs (Coherence_Analysis_ShuffleTestResults
    % by default), in the shape PlotSubjectsAndGroupTopography.m (and other per-subject plotters/
    % analyses in this repo) expect -- so you don't have to re-derive this loop from scratch every
    % time you come back to plot a previous run.
    %
    % Typical use, right after re-opening MATLAB with no workspace state:
    %     ResultsArray = LoadCoherenceResultsFromFolder(output_dir);
    %     PlotSubjectsAndGroupTopography(ResultsArray, 'True_MSC', OutputDir=fullfile(output_dir, 'plots'), ...
    %         SubjectIDField='Subj_id', ColorbarLabel='True MSC');
    % (Only 'True_MSC' is guaranteed to exist -- it's present whether or not the shuffle test ran.
    % 'PValues'/'ZScores'/etc. are only present on files saved with NumShuffles >= 1 -- see
    % SpeechEEGCoherence_ShuffleTest.m and open-items.md item 18.)
    %
    % ARGUMENTS:
    %   output_dir: folder to search, e.g. the SAME output_dir used when script_runCoherenceAnalysis.m
    %       produced the files (printed at the end of that script's run, and recorded in
    %       coherence_control_analysis_parameters.txt inside that folder).
    %   options.FilePattern (default 'coherence_control_analysis_*.mat'): glob pattern (relative to
    %       output_dir) matching the per-subject files to load. Change this if you're loading from
    %       a differently-prefixed batch (e.g. a customized output_file_prefix) in the same folder.
    %   options.ResultsVarName (default 'Coherence_Analysis_ShuffleTestResults'): name of the
    %       variable inside each .mat file that holds the actual results struct (matches the
    %       variable name script_runCoherenceAnalysis.m calls `save` with).
    %   options.SkipMissingVar (default true): if a matched file doesn't contain
    %       options.ResultsVarName, skip it (with a warning, and the file path recorded in
    %       SkippedFiles) rather than erroring the whole load. Pass false to error on the first
    %       such file instead.
    %
    % RETURNS:
    %   ResultsArray: struct array (or, if the subjects' results structs don't all share identical
    %       fields, a cell array -- same struct-array-with-cell-array-fallback pattern used by
    %       ExtractStructArrayField.m and PlotSubjectsAndGroupTopography.m itself) of the inner
    %       results structs, one per successfully-loaded file, in the order the files were found
    %       (dir() -- typically alphabetical, i.e. subject-tag order for the S01/S02/... naming
    %       scheme).
    %   SubjectIDs: cell array of char subject tags parsed from each loaded filename (e.g. 'S04'
    %       from 'coherence_control_analysis_S04.mat'), same order as ResultsArray -- handy as
    %       options.SubjectIDs for PlotSubjectsAndGroupTopography.m if you'd rather not rely on
    %       options.SubjectIDField (e.g. SubjectIDField='Subj_id') instead.
    %   LoadedFiles, SkippedFiles: cell arrays of full file paths, for provenance/debugging.
    % ===============================================================================================================
    arguments
        output_dir (1,:) char
        options.FilePattern (1,:) char = 'coherence_control_analysis_*.mat'
        options.ResultsVarName (1,:) char = 'Coherence_Analysis_ShuffleTestResults'
        options.SkipMissingVar (1,1) logical = true
    end

    if ~isfolder(output_dir)
        error('LoadCoherenceResultsFromFolder:OutputDirNotFound', ...
            'output_dir does not exist: %s', output_dir);
    end

    file_listing = dir(fullfile(output_dir, options.FilePattern));
    file_listing = file_listing(~[file_listing.isdir]);   % just in case the pattern ever matches a folder

    if isempty(file_listing)
        error('LoadCoherenceResultsFromFolder:NoFilesFound', ...
            'No files matching ''%s'' found in %s.', options.FilePattern, output_dir);
    end

    ResultsCell = {};
    SubjectIDs = {};
    LoadedFiles = {};
    SkippedFiles = {};

    for k = 1:numel(file_listing)
        this_file = fullfile(output_dir, file_listing(k).name);
        S = load(this_file);

        if ~isfield(S, options.ResultsVarName)
            msg = sprintf('%s does not contain a variable named ''%s'' -- skipping.', ...
                this_file, options.ResultsVarName);
            if options.SkipMissingVar
                warning('LoadCoherenceResultsFromFolder:MissingVar', '%s', msg);
                SkippedFiles{end+1} = this_file; %#ok<AGROW>
                continue;
            else
                error('LoadCoherenceResultsFromFolder:MissingVar', '%s', msg);
            end
        end

        ResultsCell{end+1} = S.(options.ResultsVarName); %#ok<AGROW>
        LoadedFiles{end+1} = this_file; %#ok<AGROW>

        % Parse the subject tag back out of the filename, e.g. 'S04' from
        % 'coherence_control_analysis_S04.mat' -- strip the matched glob's fixed prefix/suffix by
        % reusing options.FilePattern's own literal (non-wildcard) parts rather than assuming the
        % specific 'coherence_control_analysis_' prefix, so this still works if FilePattern was
        % customized.
        [~, name_no_ext] = fileparts(file_listing(k).name);
        pattern_parts = strsplit(options.FilePattern, '*');
        tag = name_no_ext;
        if ~isempty(pattern_parts) && ~isempty(pattern_parts{1})
            tag = regexprep(tag, ['^' regexptranslate('escape', pattern_parts{1})], '');
        end
        if numel(pattern_parts) > 1 && ~isempty(pattern_parts{end})
            [~, pattern_suffix_noext] = fileparts(pattern_parts{end});
            if ~isempty(pattern_suffix_noext)
                tag = regexprep(tag, [regexptranslate('escape', pattern_suffix_noext) '$'], '');
            end
        end
        SubjectIDs{end+1} = tag; %#ok<AGROW>
    end

    if isempty(ResultsCell)
        error('LoadCoherenceResultsFromFolder:NothingLoaded', ...
            'Every matched file in %s was skipped (none contained a ''%s'' variable).', ...
            output_dir, options.ResultsVarName);
    end

    % Prefer a plain struct array (what most of this repo's plotters/analyses expect) when every
    % loaded results struct shares identical fields; fall back to a cell array otherwise -- same
    % fallback philosophy as ExtractStructArrayField.m and PlotSubjectsAndGroupTopography.m's own
    % input handling, since a struct array requires identical fields across elements.
    try
        ResultsArray = [ResultsCell{:}];
    catch
        warning('LoadCoherenceResultsFromFolder:FieldMismatch', ...
            ['Loaded results structs do not all share identical fields -- returning a cell array ' ...
             'instead of a struct array. PlotSubjectsAndGroupTopography.m and similar functions in ' ...
             'this repo accept either.']);
        ResultsArray = ResultsCell;
    end
end
