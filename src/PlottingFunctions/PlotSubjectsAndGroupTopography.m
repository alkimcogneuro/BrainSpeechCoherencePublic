function [fig_handle, PlotResults] = PlotSubjectsAndGroupTopography(ResultsArray, FieldName, options)
    % =========================================================================================
    % Plot Every Subject's Per-Channel Data as a Topography, Plus the Group-Average Topography
    % =========================================================================================
    % Takes an array of your full results structures (e.g. Coherence_Analysis_ShuffleTestResults,
    % one element per subject) UNCHANGED -- this function never strips, copies-out-only, or
    % otherwise throws away any of the fields on those structs. You tell it which field to plot
    % via FieldName (a string, e.g. 'True_MSC' or 'ZScores'); it reaches into that one field per
    % subject via dynamic field access, draws one topoplot per subject, averages across subjects
    % at each channel, and draws one more topoplot of that group-level average. Because the field
    % to plot is a parameter rather than something hardcoded, the SAME function covers MSC,
    % p-values, z-scores, RankZ, or any other per-subject/per-channel quantity that lives as a
    % field on your results structs -- call it again with a different FieldName to look at
    % something else from the exact same ResultsArray.
    %
    % This produces TWO things:
    %   1) One combined "overview" figure: every subject panel plus the group panel, tiled
    %      together in a single figure sharing one color scale, so every subject (and the group
    %      average) can be compared directly at a glance.
    %   2) Individual PNG files, saved to options.OutputDir: one per subject, plus one for the
    %      group average -- for looking more closely at any single result than the (necessarily
    %      small) tiles in the combined figure allow.
    %
    % Arguments:
    %   ResultsArray: cell array OR struct array of your results structures, one per subject
    %       (e.g. Coherence_Analysis_ShuffleTestResults, or a CoherenceResults_bandavgmsc array
    %       from BandAverageMSC.m -- anything with per-channel numeric fields and a .Chanlocs).
    %       Cell array or struct array are both accepted (mirrors GroupLevel_MSC_Analysis.m's
    %       handling of subject arrays), since a struct array requires every element to share
    %       identical fields, which may not always hold across subjects. Every field on every
    %       element -- band limits, subject ID, timestamps, whatever else you're tracking -- is
    %       left completely intact; this function only ever READS FieldName (and, optionally,
    %       .Chanlocs / a subject-ID field) off each element, nothing is copied out and discarded.
    %   FieldName: char/string naming the field to plot, e.g. 'True_MSC', 'PValues', 'ZScores'.
    %       Must be a numeric, per-channel field present on every element of ResultsArray
    %       ([1 x Num_channels] or [Num_channels x 1]).
    %
    % Optional name-value arguments:
    %   options.OutputDir     (REQUIRED): folder to save each subject's PNG, plus the group
    %                                        average's PNG, into. Created automatically if it
    %                                        doesn't already exist.
    %   options.Chanlocs      (default []): EEGLAB channel locations structure. If left empty,
    %                                        this function looks for a .Chanlocs field on the
    %                                        FIRST subject's results struct and uses that. An
    %                                        error is raised if neither is available.
    %   options.SubjectIDs    (default {}): cell array of char/string subject identifiers, one
    %                                        per element of ResultsArray, used as each subject
    %                                        panel's title AND (sanitized) filename. Takes
    %                                        priority over options.SubjectIDField if both are
    %                                        given. Defaults to {'Subject01', ...} if neither is
    %                                        supplied.
    %   options.SubjectIDField (default ''): name of a field on each element of ResultsArray to
    %                                        use as that subject's ID instead of auto-generating
    %                                        one (e.g. 'Subject_id', if that's what your results
    %                                        structs call it -- check the exact spelling in your
    %                                        data). Values are converted to char if numeric or
    %                                        string. Ignored if options.SubjectIDs is supplied.
    %   options.ColorbarLabel (default ''): colorbar y-axis label, shared by every panel. If left
    %                                        empty, defaults to FieldName itself (e.g. 'True_MSC'),
    %                                        so you get a reasonable label for free; pass your own
    %                                        (e.g. 'MSC', '-log10(p)') for something nicer.
    %   options.GroupLabel    (default 'Group Average'): title stem for the group panel (the
    %                                        actual title also appends '(N=<Num_subjects>)') and
    %                                        the base of the group PNG's filename.
    %   options.FigureTitle   (default ''): overall title for the COMBINED figure, shown above
    %                                        all tiles. Skipped if empty.
    %   options.Clim          (default []): [1 x 2] color limits shared by every panel -- the
    %                                        combined figure's tiles AND every saved PNG (so a
    %                                        PNG you open later still reads on the same color
    %                                        scale as the overview figure it came from). If left
    %                                        empty, limits are computed automatically from the
    %                                        combined range of every subject's FieldName values
    %                                        (see compute_shared_clim below).
    %   options.Colormap      (default []): [N x 3] colormap matrix, shared by every panel and
    %                                        every saved PNG. Defaults to parula(256) if left
    %                                        empty (and options.HighlightThreshold is also
    %                                        empty -- see below). Takes priority over
    %                                        options.HighlightThreshold if both are supplied.
    %   options.HighlightThreshold (default []): scalar. If supplied (and options.Colormap is
    %                                        NOT also supplied), the colormap is built
    %                                        automatically via BuildThresholdHighlightColormap.m
    %                                        instead of defaulting to parula(256): values below
    %                                        this threshold stay in a muted blue, values above it
    %                                        ramp from orange to bright yellow, so channels
    %                                        exceeding the threshold visually pop out instead of
    %                                        blending into a smooth, evenly-spread color scale.
    %                                        Built using the SAME Clim this call resolves to
    %                                        (see options.Clim above), so the threshold's visual
    %                                        position is always correct relative to what's
    %                                        actually plotted -- e.g. for ZScores, pass
    %                                        HighlightThreshold=2.0 to highlight channels with
    %                                        Z > ~2. For full control over the highlight/muted
    %                                        colors (not just the threshold), call
    %                                        BuildThresholdHighlightColormap.m yourself and pass
    %                                        its result via options.Colormap instead.
    %   options.NumCols       (default []): number of columns in the COMBINED figure's tile grid.
    %                                        If left empty, a roughly-square grid is chosen
    %                                        automatically based on the total number of panels
    %                                        (Num_subjects + 1).
    %   options.HighlightGroupPanel (default true): if true, draws a colored border around the
    %                                        group-average tile IN THE COMBINED FIGURE so it's
    %                                        visually distinguishable from the individual subject
    %                                        panels at a glance. (The standalone group PNG has no
    %                                        other panels to distinguish itself from, so this
    %                                        doesn't apply there.)
    %   options.Visible       (default true): if false, the COMBINED figure is created invisibly.
    %   options.SavePath      (default ''): if non-empty, the COMBINED figure is ALSO saved to
    %                                        this path (extension determines format, e.g. '.png',
    %                                        '.fig', '.pdf'), in addition to the per-subject/group
    %                                        PNGs that always get saved to options.OutputDir.
    %   options.IndividualVisible          (default false): if false, the per-subject/group PNG
    %                                        figures are created invisibly (avoids flooding your
    %                                        screen with Num_subjects+1 windows while saving).
    %   options.IndividualCloseAfterSaving (default true): if true, each per-subject/group PNG
    %                                        figure is closed immediately after it's saved.
    %
    % Returns:
    %   fig_handle:  handle to the COMBINED overview figure (all subjects + group, tiled).
    %   PlotResults: a structure containing:
    %     .FieldName           the FieldName you passed in, echoed back for provenance.
    %     .Data_matrix         [Num_subjects x Num_channels] each subject's ResultsArray(i).(FieldName)
    %                          (or ResultsArray{i}.(FieldName)), stacked into one matrix.
    %     .GroupMean           [1 x Num_channels] Data_matrix averaged across subjects, per
    %                          channel (simple mean, omitnan).
    %     .SubjectIDs          {1 x Num_subjects} cell array of subject IDs actually used
    %                          (options.SubjectIDs, options.SubjectIDField-derived, or the
    %                          auto-generated defaults, in that priority order).
    %     .Clim                [1 x 2] color limits actually used for every panel/PNG.
    %     .Colormap            [N x 3] colormap actually used for every panel/PNG.
    %     .HighlightThreshold  the options.HighlightThreshold you passed in, echoed back for
    %                          provenance (empty if you didn't use it).
    %     .Chanlocs            the channel locations structure actually used.
    %     .SubjectPNGFiles     {1 x Num_subjects} cell array of full paths to each subject's
    %                          saved PNG (in the same order as ResultsArray).
    %     .GroupPNGFile        full path to the saved group-average PNG.
    %     .SubjectFigureHandles [Num_subjects x 1] figure handles for the per-subject PNG
    %                          figures (invalid/closed if options.IndividualCloseAfterSaving).
    %     .GroupFigureHandle   figure handle for the group PNG figure (invalid/closed if
    %                          options.IndividualCloseAfterSaving).
    %   NOTE: PlotResults does NOT include a copy of ResultsArray itself -- you already have it,
    %   and it's untouched by this call (see "Arguments" above), so there's nothing new to return.
    %
    % -----------------------------------------------------------------------------------------
    % Design notes
    % -----------------------------------------------------------------------------------------
    % 1) Why a field-name selector instead of a plain matrix: an earlier version of this function
    %    (briefly, same day) took a bare [Num_subjects x Num_channels] matrix built by the caller,
    %    which meant the caller had to pre-extract whatever field they wanted and everything else
    %    on the results struct (band limits, subject ID, run history/provenance -- see the many
    %    other fields on Coherence_Analysis_ShuffleTestResults) got left behind outside this
    %    function's view entirely. Taking ResultsArray + FieldName instead means you keep passing
    %    your full, richly-annotated results structs around everywhere else in your code, and
    %    only tell THIS function, at call time, which one field of them to draw -- no separate
    %    extraction step, and nothing about your results structs' other fields is ever at risk of
    %    being dropped or forgotten. This is also still the single generic function that replaces
    %    the near-duplicate family in open-items.md item 7 (PlotShuffleTestSignificance_Topography
    %    etc., each hardcoded to a different field) -- now via FieldName instead of a pre-built
    %    matrix.
    %
    % 2) The group panel is always a plain per-channel mean (omitnan) of the FieldName values
    %    across subjects. If you want the group panel to show something else -- e.g. a
    %    t-statistic against zero (like GroupLevel_ShuffleTest_Analysis.m produces) rather than a
    %    mean of per-subject values -- compute that vector yourself and plot it as its own
    %    single-panel topography (e.g. via call_topoplot_custom.m directly); this function's group
    %    panel is deliberately just "the mean of what's in that field," not a general statistics
    %    engine.
    %
    % 3) call_topoplot_custom, not call_topoplot_eeglab: this function draws every panel (tiled
    %    AND standalone PNGs) with call_topoplot_custom.m rather than call_topoplot_eeglab.m.
    %    call_topoplot_custom.m is a drop-in replacement (same first 7 positional args, same
    %    "draws into the current/given axes and returns the axes handle" contract), documented
    %    as requiring no changes anywhere else -- so the tiledlayout/axes_in_tile_slot machinery
    %    below is UNCHANGED from what call_topoplot_eeglab would have needed. One consequence,
    %    per call_topoplot_custom's own header: because it scopes colormap(ax,...)/clim(ax,...)
    %    to an EXPLICIT axes handle every time (rather than relying on EEGLAB topoplot()'s
    %    current-axes/current-figure colormap state, which does not cleanly scope to one axes),
    %    the "final pass: reassert every tile's colormap/clim" loop below is now REDUNDANT for
    %    correctness -- each call already sets its own tile correctly. It's left in place anyway
    %    (harmless no-op) rather than removed, per call_topoplot_custom's own migration guidance
    %    ("safe to leave in place during migration and delete later").
    %
    % 4) One shared color scale, everywhere: the combined figure's tiles AND every standalone
    %    per-subject/group PNG all use the SAME options.Clim / options.Colormap (auto-computed
    %    once, up front, from every subject's FieldName values, unless supplied explicitly). This
    %    is deliberate: since the stated purpose of the individual PNGs is to look more closely at
    %    a panel you already saw in the combined overview, having the PNG land on a DIFFERENT
    %    (e.g. auto-rescaled-per-subject) color scale would make the zoomed-in view visually
    %    inconsistent with the overview it came from, and would make the saved PNGs uncomparable
    %    to each other too.
    %
    % 5) EEGLAB topoplot + tiledlayout + colormap (background, largely superseded by note 3):
    %    the axes_in_tile_slot() workaround exists because EEGLAB's topoplot() internally
    %    resizes/squares its own axes, which a TiledChartLayout tile rejects. call_topoplot_custom
    %    does not do this (it only ever calls axis(ax,'equal')/axis(ax,'off'), never repositions
    %    the axes), so this workaround is likely no longer strictly necessary either -- but it's
    %    kept as-is here too, for the same "minimal-diff, verify side-by-side first" reasoning as
    %    note 3.
    %
    % 6) Non-finite values: if a subject's FieldName value has non-finite (NaN/Inf) entries at
    %    some channel, that channel is simply that subject's own missing value; the group average
    %    is computed with 'omitnan' per channel, so different channels may effectively be
    %    averaged over different numbers of subjects. A warning is issued if this occurs.
    %
    % Migration note (2026-09-18): this function has changed signature TWICE today in the same
    % session -- first from (CoherenceResults_bandavgmsc_array) [reading .msc_band internally] to
    % (Data_matrix, chanlocs) [a plain matrix, no field-reading at all], and now to
    % (ResultsArray, FieldName, options) [full structs preserved, field read internally via
    % FieldName]. The old PlotSubjectsAndGroupMSC_Topography.m file still exists (not deleted,
    % header-noted as deprecated) for any call sites not yet updated. The one known caller
    % (dev/TestScripts/script_AnalyzeCoherence_Control30s.m) was updated again to match this final
    % signature -- see open-items.md item 13.
    % =========================================================================================
    arguments
        ResultsArray
        FieldName (1,:) char
        options.OutputDir (1,:) char
        options.Chanlocs = []
        options.SubjectIDs cell = {}
        options.SubjectIDField (1,:) char = ''
        options.ColorbarLabel (1,:) char = ''
        options.GroupLabel (1,:) char = 'Group Average'
        options.FigureTitle (1,:) char = ''
        options.Clim (:,:) double = []
        options.Colormap (:,:) double = []
        options.HighlightThreshold (:,:) double = []
        options.NumCols (:,:) double = []
        options.HighlightGroupPanel (1,1) logical = true
        options.Visible (1,1) logical = true
        options.SavePath (1,:) char = ''
        options.IndividualVisible (1,1) logical = false
        options.IndividualCloseAfterSaving (1,1) logical = true
    end

    % ---- Normalize subject count and per-subject access ------------------------------------
    Num_subjects = numel(ResultsArray);
    if Num_subjects < 1
        error('PlotSubjectsAndGroupTopography:NoSubjects', ...
            'ResultsArray has zero elements (no subjects).');
    end
    get_subj = @(arr, idx) get_subject_struct(arr, idx);

    % ---- Resolve channel locations -----------------------------------------------------------
    if isempty(options.Chanlocs)
        subj1 = get_subj(ResultsArray, 1);
        if isfield(subj1, 'Chanlocs')
            chanlocs = subj1.Chanlocs;
        else
            error('PlotSubjectsAndGroupTopography:NoChanlocs', ...
                ['options.Chanlocs was not supplied, and the first subject''s results struct ' ...
                 'has no .Chanlocs field. Either pass options.Chanlocs explicitly, or use ' ...
                 'results structs that have .Chanlocs populated.']);
        end
    else
        chanlocs = options.Chanlocs;
    end
    Num_channels = numel(chanlocs);

    % ---- Validate each subject and assemble the data matrix from FieldName ------------------
    Data_matrix = nan(Num_subjects, Num_channels);
    for subj_idx = 1:Num_subjects
        subj_struct = get_subj(ResultsArray, subj_idx);
        if ~isfield(subj_struct, FieldName)
            error('PlotSubjectsAndGroupTopography:MissingField', ...
                'Element %d of ResultsArray has no ''%s'' field.', subj_idx, FieldName);
        end
        subj_values = subj_struct.(FieldName);
        subj_values = subj_values(:)';  % row vector, [1 x Num_channels]
        if numel(subj_values) ~= Num_channels
            error('PlotSubjectsAndGroupTopography:ChanlocsMismatch', ...
                ['Element %d of ResultsArray has %d values in .%s, but chanlocs has %d ' ...
                 'elements. These must match.'], subj_idx, numel(subj_values), FieldName, ...
                Num_channels);
        end
        Data_matrix(subj_idx, :) = subj_values;
    end

    % ---- Resolve subject IDs -----------------------------------------------------------------
    if ~isempty(options.SubjectIDs)
        if numel(options.SubjectIDs) ~= Num_subjects
            error('PlotSubjectsAndGroupTopography:SubjectIDMismatch', ...
                'options.SubjectIDs has %d entries, but there are %d subjects. These must match.', ...
                numel(options.SubjectIDs), Num_subjects);
        end
        SubjectIDs = options.SubjectIDs;
    elseif ~isempty(options.SubjectIDField)
        SubjectIDs = cell(1, Num_subjects);
        for subj_idx = 1:Num_subjects
            subj_struct = get_subj(ResultsArray, subj_idx);
            if ~isfield(subj_struct, options.SubjectIDField)
                error('PlotSubjectsAndGroupTopography:MissingSubjectIDField', ...
                    ['Element %d of ResultsArray has no ''%s'' field (options.SubjectIDField). ' ...
                     'Check the exact field name on your results structs, or omit ' ...
                     'options.SubjectIDField to auto-generate subject IDs instead.'], ...
                    subj_idx, options.SubjectIDField);
            end
            SubjectIDs{subj_idx} = to_char_id(subj_struct.(options.SubjectIDField));
        end
    else
        pad_width = max(2, floor(log10(Num_subjects)) + 1);
        SubjectIDs = arrayfun(@(idx) sprintf('Subject%0*d', pad_width, idx), ...
            1:Num_subjects, 'UniformOutput', false);
    end

    % ---- Resolve colorbar label: default to FieldName itself if not supplied ----------------
    if isempty(options.ColorbarLabel)
        ColorbarLabel = FieldName;
    else
        ColorbarLabel = options.ColorbarLabel;
    end

    % ---- Warn (don't silently average over) any non-finite per-subject values ---------------
    nonfinite_mask = ~isfinite(Data_matrix);
    if any(nonfinite_mask(:))
        warning('PlotSubjectsAndGroupTopography:NonFiniteInputs', ...
            ['%d entries across %d subjects'' .%s values are non-finite. These will be ' ...
             'excluded via omitnan when computing the group average, which means different ' ...
             'channels may effectively be averaged over different numbers of subjects.'], ...
            sum(nonfinite_mask(:)), sum(any(nonfinite_mask, 2)), FieldName);
    end

    % ---- Compute the group-level average, per channel ----------------------------------------
    GroupMean = mean(Data_matrix, 1, 'omitnan');  % [1 x Num_channels]

    % ---- Determine the shared color scale (used for EVERY panel and EVERY saved PNG) --------
    if isempty(options.Clim)
        Clim = compute_shared_clim(Data_matrix);
    else
        if numel(options.Clim) ~= 2
            error('PlotSubjectsAndGroupTopography:BadClim', ...
                'options.Clim must be a 2-element vector [lo, hi].');
        end
        Clim = options.Clim;
    end

    % ---- Determine the shared colormap ---------------------------------------------------------
    % Priority: explicit options.Colormap, then options.HighlightThreshold (built via
    % BuildThresholdHighlightColormap.m using the Clim just resolved above), then parula(256).
    if ~isempty(options.Colormap)
        if size(options.Colormap, 2) ~= 3
            error('PlotSubjectsAndGroupTopography:BadColormap', ...
                'options.Colormap must be an [N x 3] RGB colormap matrix.');
        end
        cmap = options.Colormap;
        if ~isempty(options.HighlightThreshold)
            warning('PlotSubjectsAndGroupTopography:ColormapOverridesHighlightThreshold', ...
                ['Both options.Colormap and options.HighlightThreshold were supplied -- ' ...
                 'options.Colormap takes priority, so options.HighlightThreshold is ignored. ' ...
                 'Omit options.Colormap to have HighlightThreshold build the colormap instead.']);
        end
    elseif ~isempty(options.HighlightThreshold)
        if ~isscalar(options.HighlightThreshold)
            error('PlotSubjectsAndGroupTopography:BadHighlightThreshold', ...
                'options.HighlightThreshold must be a scalar.');
        end
        cmap = BuildThresholdHighlightColormap(Clim, options.HighlightThreshold);
    else
        cmap = parula(256);
    end

    % ---- Build the combined overview figure ----------------------------------------------------
    Num_panels = Num_subjects + 1;  % one per subject, plus the group average
    if isempty(options.NumCols)
        Num_cols = ceil(sqrt(Num_panels));
    else
        Num_cols = options.NumCols;
    end
    Num_rows = ceil(Num_panels / Num_cols);

    fig_handle = render_subjects_and_group_figure(Data_matrix, GroupMean, chanlocs, ...
        SubjectIDs, options.GroupLabel, ColorbarLabel, options.FigureTitle, ...
        Clim, cmap, Num_rows, Num_cols, options.HighlightGroupPanel, options.Visible);

    % ---- Save the combined figure, if requested ------------------------------------------------
    if ~isempty(options.SavePath)
        [~, ~, ext] = fileparts(options.SavePath);
        if strcmpi(ext, '.fig')
            savefig(fig_handle, options.SavePath);
        else
            exportgraphics(fig_handle, options.SavePath);
        end
    end

    % ---- Save one standalone PNG per subject, plus one for the group average -----------------
    if ~isfolder(options.OutputDir)
        mkdir(options.OutputDir);
    end

    SubjectPNGFiles      = cell(Num_subjects, 1);
    SubjectFigureHandles = gobjects(Num_subjects, 1);
    for subj_idx = 1:Num_subjects
        safe_subj_id           = sanitize_for_filename(SubjectIDs{subj_idx});
        subj_png_path          = fullfile(options.OutputDir, [safe_subj_id, '_topography.png']);
        [subj_fig, subj_saved] = render_and_save_single_topography(Data_matrix(subj_idx, :), ...
            chanlocs, SubjectIDs{subj_idx}, ColorbarLabel, Clim, cmap, subj_png_path, ...
            options.IndividualVisible);
        SubjectPNGFiles{subj_idx}    = subj_saved;
        SubjectFigureHandles(subj_idx) = subj_fig;
        if options.IndividualCloseAfterSaving
            close(subj_fig);
        end
    end

    group_title        = sprintf('%s (N=%d)', options.GroupLabel, Num_subjects);
    safe_group_label    = sanitize_for_filename(options.GroupLabel);
    group_png_path      = fullfile(options.OutputDir, [safe_group_label, '_topography.png']);
    [group_fig, GroupPNGFile] = render_and_save_single_topography(GroupMean, chanlocs, ...
        group_title, ColorbarLabel, Clim, cmap, group_png_path, options.IndividualVisible);
    GroupFigureHandle = group_fig;
    if options.IndividualCloseAfterSaving
        close(group_fig);
    end

    % ---- Assemble results structure -------------------------------------------------------------
    PlotResults = struct();
    PlotResults.FieldName            = FieldName;
    PlotResults.Data_matrix          = Data_matrix;
    PlotResults.GroupMean            = GroupMean;
    PlotResults.SubjectIDs           = SubjectIDs;
    PlotResults.Clim                 = Clim;
    PlotResults.Colormap             = cmap;
    PlotResults.HighlightThreshold   = options.HighlightThreshold;
    PlotResults.Chanlocs             = chanlocs;
    PlotResults.SubjectPNGFiles      = SubjectPNGFiles;
    PlotResults.GroupPNGFile         = GroupPNGFile;
    PlotResults.SubjectFigureHandles = SubjectFigureHandles;
    PlotResults.GroupFigureHandle    = GroupFigureHandle;
end

% =============================================================================================
% Local helper functions
% =============================================================================================

function s = get_subject_struct(arr, idx)
    % Allow ResultsArray to be passed in as either a cell array of structs or a struct array
    % (mirrors GroupLevel_MSC_Analysis.m's get_subject_struct).
    if iscell(arr)
        s = arr{idx};
    else
        s = arr(idx);
    end
end

function id_str = to_char_id(id_val)
    % Convert a subject-ID field's value to a plain char row, regardless of whether it was
    % stored as char, string, or numeric on the results struct.
    if ischar(id_val)
        id_str = id_val;
    elseif isstring(id_val)
        id_str = char(id_val);
    elseif isnumeric(id_val)
        id_str = num2str(id_val);
    else
        error('PlotSubjectsAndGroupTopography:BadSubjectIDValue', ...
            'options.SubjectIDField value must be char, string, or numeric (got %s).', ...
            class(id_val));
    end
end

function Clim = compute_shared_clim(Data_matrix)
    % Default clim logic when the caller doesn't supply options.Clim explicitly: computed from
    % the combined range of every subject's data (the group average, being an average, will
    % typically fall well within this range).
    finite_vals = Data_matrix(isfinite(Data_matrix));
    if isempty(finite_vals)
        % Degenerate case: every value is non-finite. Fall back to [0, 1] rather than erroring
        % outright.
        Clim = [0, 1];
        return;
    end
    clim_min = min(finite_vals);
    clim_max = max(finite_vals) * 1.01;  % small headroom, for better contrast
    if clim_max <= clim_min
        % Degenerate case: all values identical (or clim_min <= 0 so the 1.01 headroom trick
        % above didn't move clim_max). Avoid a zero-width color axis, which would error in
        % topoplot/clim.
        clim_max = clim_min + eps(clim_min) * 100 + 1e-6;
    end
    Clim = [clim_min, clim_max];
end

function safe_str = sanitize_for_filename(str)
    % Replace anything that isn't a letter, digit, hyphen, or underscore with '_' (handles
    % spaces, slashes, etc. in free-text subject IDs or GroupLabel), matching the convention
    % already used in PlotIndividualSubjectsMSC_Topography.m.
    safe_str = regexprep(str, '[^a-zA-Z0-9_-]', '_');
end

function [fig_handle, saved_path] = render_and_save_single_topography(data_vector, chanlocs, ...
        plot_title, colorbar_label, Clim, cmap, save_path, fig_visible)
    % Draw ONE topoplot into its own standalone figure and save it as a PNG. Used for both the
    % per-subject PNGs and the group-average PNG -- unlike the combined figure's tiles, a
    % standalone figure has only one axes, so none of the tiledlayout/axes_in_tile_slot machinery
    % is needed here regardless of which topoplot renderer is used.
    if fig_visible
        visible_str = 'on';
    else
        visible_str = 'off';
    end
    fig_handle = figure('Color', 'white', 'Position', [100 100 500 500], 'Visible', visible_str);
    ax = axes('Parent', fig_handle);
    call_topoplot_custom(data_vector, chanlocs, plot_title, colorbar_label, Clim(1), Clim(2), ...
        cmap, 'Axes', ax);
    exportgraphics(fig_handle, save_path);
    saved_path = save_path;
end

function fig_handle = render_subjects_and_group_figure(Data_matrix, GroupMean, chanlocs, ...
        SubjectIDs, group_label, colorbar_label, fig_title, Clim, cmap, Num_rows, Num_cols, ...
        highlight_group_panel, fig_visible)
    % Build one figure containing one topoplot tile per subject, plus one final tile for the
    % group average, all sharing one color scale. See design notes 3) and 5) above for why the
    % tiledlayout/topoplot workarounds (axes_in_tile_slot, final colormap-reassertion pass) are
    % kept even though call_topoplot_custom likely no longer strictly needs them.
    Num_subjects = size(Data_matrix, 1);
    Num_panels   = Num_subjects + 1;

    if fig_visible
        visible_str = 'on';
    else
        visible_str = 'off';
    end
    fig_width  = min(1800, max(500, 380 * Num_cols));
    fig_height = min(1400, max(450, 380 * Num_rows));
    fig_handle = figure('Color', 'white', 'Position', [80 80 fig_width fig_height], ...
        'Visible', visible_str);
    t = tiledlayout(fig_handle, Num_rows, Num_cols, 'TileSpacing', 'compact', 'Padding', 'compact');

    tile_axes = gobjects(Num_panels, 1);

    % ---- One tile per subject -----------------------------------------------------------------
    for subj_idx = 1:Num_subjects
        axes_in_tile_slot(fig_handle, t, subj_idx); % swapped for a standalone axes below
        tile_axes(subj_idx) = call_topoplot_custom(Data_matrix(subj_idx, :), chanlocs, ...
            SubjectIDs{subj_idx}, colorbar_label, Clim(1), Clim(2), cmap);
    end

    % ---- Final tile: the group average ---------------------------------------------------------
    group_tile_index = Num_panels;
    group_tile_position = axes_in_tile_slot(fig_handle, t, group_tile_index, true); %#ok<NASGU>
    group_title = sprintf('%s (N=%d)', group_label, Num_subjects);
    tile_axes(group_tile_index) = call_topoplot_custom(GroupMean, chanlocs, ...
        group_title, colorbar_label, Clim(1), Clim(2), cmap);

    if highlight_group_panel
        % Draw a colored border around the group tile so it stands out from the individual
        % subject panels at a glance. Uses a figure-level annotation (rather than styling the
        % topoplot axes itself), since the topoplot renderer manages its own axes appearance
        % internally.
        annotation(fig_handle, 'rectangle', group_tile_position, ...
            'Color', [0.75 0.1 0.1], 'LineWidth', 3);
    end

    % ---- Final pass: reassert every tile's colormap/clim. Redundant-but-harmless with
    % call_topoplot_custom (see design note 3 above) -- kept for minimal-diff continuity.
    for panel_idx = 1:Num_panels
        colormap(tile_axes(panel_idx), cmap);
        clim(tile_axes(panel_idx), Clim);
    end

    if ~isempty(fig_title)
        title(t, fig_title, 'FontSize', 16, 'FontWeight', 'bold');
    end
end

function varargout = axes_in_tile_slot(fig_handle, t, tile_index, return_position)
    % Reserve a SPECIFIC slot (by index) in tiledlayout t, then hand back a plain, standalone
    % axes occupying that exact same screen position -- instead of the tile axes itself. See
    % PlotThreePanelTopography.m's axes_in_tile_slot for the full explanation of why this is
    % needed for EEGLAB's topoplot() (its internal Position changes conflict with
    % TiledChartLayout-managed axes) and why tile_index must be passed explicitly rather than
    % relying on nexttile(t)'s auto-advance. call_topoplot_custom likely doesn't need this (see
    % design note 5 above), but it's kept for minimal-diff continuity.
    %
    % If return_position is true, this also returns the tile's [left bottom width height]
    % position (figure-normalized) as a second output, e.g. for drawing a figure-level
    % annotation at the same spot.
    if nargin < 4
        return_position = false;
    end
    nexttile(t, tile_index);
    tile_ax = gca;
    drawnow;   % force the layout to finalize this tile's Position before we read it
    tile_position = get(tile_ax, 'Position');
    delete(tile_ax);
    ax = axes('Parent', fig_handle, 'Position', tile_position);
    axes(ax);  % make it the current axes, since the topoplot renderer draws via gca

    varargout{1} = ax;
    if return_position
        varargout{1} = tile_position;
    end
end
