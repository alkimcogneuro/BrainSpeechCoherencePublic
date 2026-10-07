function [fig_handle, PlotResults] = PlotSubjectsAndGroupMSC_Topography(CoherenceResults_bandavgmsc_array, options)
    % =========================================================================================
    % DEPRECATED (2026-09-18): superseded by PlotSubjectsAndGroupTopography.m, which takes your
    % results struct array UNCHANGED (e.g. Coherence_Analysis_ShuffleTestResults, or the
    % CoherenceResults_bandavgmsc array this file expects) plus a FieldName string telling it
    % which field to plot ('msc_band', 'True_MSC', 'ZScores', etc.), instead of hardcoding
    % .msc_band internally. Every other field on your results structs (band limits, subject ID,
    % provenance) stays intact and untouched -- nothing is extracted into a separate matrix. This
    % file is left in place, unmodified below, only for any call sites not yet migrated -- prefer
    % the new function for all new code. See PlotSubjectsAndGroupTopography.m's header for the
    % migration note and open-items.md item 7 for why this consolidation happened.
    % =========================================================================================
    % =========================================================================================
    % Plot Every Subject's Band-Averaged MSC Topography, Plus the Group-Average Topography
    % =========================================================================================
    % Takes an array of single-subject CoherenceResults_bandavgmsc structures (the output of
    % BandAverageMSC.m, one element per subject), draws one topoplot per subject of that
    % subject's .msc_band field, averages .msc_band across all subjects at each channel, and
    % draws one more topoplot of that group-level average.
    %
    % This produces TWO things:
    %   1) One combined "overview" figure: every subject panel plus the group panel, tiled
    %      together in a single figure sharing one color scale, so every subject (and the group
    %      average) can be compared directly at a glance.
    %   2) Individual PNG files, saved to options.OutputDir: one per subject, plus one for the
    %      group average -- for looking more closely at any single result than the (necessarily
    %      small) tiles in the combined figure allow.
    %
    % This is the single-condition counterpart to PlotGroupMSC_Topography.m /
    % PlotIndividualSubjectsMSC_Topography.m (which both take two conditions and a difference).
    % Here there is only one quantity per subject (.msc_band from BandAverageMSC), so there is no
    % difference panel -- just N subject panels plus one group panel.
    %
    % Arguments:
    %   CoherenceResults_bandavgmsc_array: cell array (or struct array) of
    %       CoherenceResults_bandavgmsc structures, one per subject, as produced by
    %       BandAverageMSC.m. Each element must have a .msc_band field
    %       ([Num_channels x 1] or [1 x Num_channels]). Cell array or struct array are both
    %       accepted (mirrors GroupLevel_MSC_Analysis.m's handling of subject arrays), since a
    %       struct array requires every element to share identical fields, which may not always
    %       hold across subjects.
    %
    % Optional name-value arguments:
    %   options.OutputDir     (REQUIRED): folder to save each subject's PNG, plus the group
    %                                        average's PNG, into. Created automatically if it
    %                                        doesn't already exist.
    %   options.Chanlocs      (default []): EEGLAB channel locations structure. If left empty,
    %                                        this function looks for a .Chanlocs field on the
    %                                        FIRST subject's CoherenceResults_bandavgmsc struct
    %                                        (populated upstream by
    %                                        Apply2Dataset_CrossSpectralDensity / AnalyzeSpeechEEGCoherence_Dataset)
    %                                        and uses that. An error is raised if neither is
    %                                        available.
    %   options.SubjectIDs    (default {}): cell array of char/string subject identifiers, one
    %                                        per element of CoherenceResults_bandavgmsc_array,
    %                                        used as each subject panel's title AND (sanitized)
    %                                        filename. Defaults to {'Subject01', ...} if not
    %                                        supplied.
    %   options.ColorbarLabel (default 'MSC'): colorbar y-axis label, shared by every panel.
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
    %                                        combined range of every subject's .msc_band values
    %                                        (see compute_shared_clim below).
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
    %     .MSC_matrix          [Num_subjects x Num_channels] each subject's .msc_band, stacked
    %                          into one matrix (row i = CoherenceResults_bandavgmsc_array{i} or
    %                          (i), depending on which container type was passed in).
    %     .GroupMean_msc_band  [1 x Num_channels] .msc_band averaged across subjects, per channel.
    %     .SubjectIDs          {1 x Num_subjects} cell array of subject IDs actually used
    %                          (either what you passed in via options.SubjectIDs, or the
    %                          auto-generated defaults).
    %     .Clim                [1 x 2] color limits actually used for every panel/PNG.
    %     .Chanlocs            the channel locations structure actually used.
    %     .SubjectPNGFiles     {1 x Num_subjects} cell array of full paths to each subject's
    %                          saved PNG (in the same order as CoherenceResults_bandavgmsc_array).
    %     .GroupPNGFile        full path to the saved group-average PNG.
    %     .SubjectFigureHandles [Num_subjects x 1] figure handles for the per-subject PNG
    %                          figures (invalid/closed if options.IndividualCloseAfterSaving).
    %     .GroupFigureHandle   figure handle for the group PNG figure (invalid/closed if
    %                          options.IndividualCloseAfterSaving).
    %
    % -----------------------------------------------------------------------------------------
    % Design notes
    % -----------------------------------------------------------------------------------------
    % 1) call_topoplot_custom, not call_topoplot_eeglab: this function draws every panel (tiled
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
    % 2) One shared color scale, everywhere: the combined figure's tiles AND every standalone
    %    per-subject/group PNG all use the SAME options.Clim (auto-computed once, up front, from
    %    every subject's data). This is deliberate: since the stated purpose of the individual
    %    PNGs is to look more closely at a panel you already saw in the combined overview, having
    %    the PNG land on a DIFFERENT (e.g. auto-rescaled-per-subject) color scale would make the
    %    zoomed-in view visually inconsistent with the overview it came from, and would make the
    %    saved PNGs uncomparable to each other too.
    %
    % 3) EEGLAB topoplot + tiledlayout + colormap (background, largely superseded by note 1):
    %    the axes_in_tile_slot() workaround exists because EEGLAB's topoplot() internally
    %    resizes/squares its own axes, which a TiledChartLayout tile rejects. call_topoplot_custom
    %    does not do this (it only ever calls axis(ax,'equal')/axis(ax,'off'), never repositions
    %    the axes), so this workaround is likely no longer strictly necessary either -- but it's
    %    kept as-is here too, for the same "minimal-diff, verify side-by-side first" reasoning as
    %    note 1.
    %
    % 4) Non-finite values: if a subject's .msc_band has non-finite (NaN/Inf) entries at some
    %    channel (see BandAverageMSC's own NonFiniteResult warning for why this can happen), that
    %    channel is simply that subject's own missing value; the group average is computed with
    %    'omitnan' per channel, so different channels may effectively be averaged over different
    %    numbers of subjects. A warning is issued if this occurs, same as PlotGroupMSC_Topography.
    % =========================================================================================
    arguments
        CoherenceResults_bandavgmsc_array
        options.OutputDir (1,:) char
        options.Chanlocs = []
        options.SubjectIDs cell = {}
        options.ColorbarLabel (1,:) char = 'MSC'
        options.GroupLabel (1,:) char = 'Group Average'
        options.FigureTitle (1,:) char = ''
        options.Clim (:,:) double = []
        options.NumCols (:,:) double = []
        options.HighlightGroupPanel (1,1) logical = true
        options.Visible (1,1) logical = true
        options.SavePath (1,:) char = ''
        options.IndividualVisible (1,1) logical = false
        options.IndividualCloseAfterSaving (1,1) logical = true
    end

    % ---- Normalize subject count and per-subject access ------------------------------------
    Num_subjects = numel(CoherenceResults_bandavgmsc_array);
    if Num_subjects < 1
        error('PlotSubjectsAndGroupMSC_Topography:NoSubjects', ...
            'CoherenceResults_bandavgmsc_array has zero elements (no subjects).');
    end
    get_subj = @(arr, idx) get_subject_struct(arr, idx);

    % ---- Resolve channel locations -----------------------------------------------------------
    if isempty(options.Chanlocs)
        subj1 = get_subj(CoherenceResults_bandavgmsc_array, 1);
        if isfield(subj1, 'Chanlocs')
            chanlocs = subj1.Chanlocs;
        else
            error('PlotSubjectsAndGroupMSC_Topography:NoChanlocs', ...
                ['options.Chanlocs was not supplied, and the first subject''s ' ...
                 'CoherenceResults_bandavgmsc struct has no .Chanlocs field. Either pass ' ...
                 'options.Chanlocs explicitly, or use CoherenceResults_bandavgmsc structs that ' ...
                 'were derived from a CoherenceResults struct with .Chanlocs populated.']);
        end
    else
        chanlocs = options.Chanlocs;
    end

    % ---- Validate each subject and assemble the MSC matrix ----------------------------------
    Num_channels = numel(chanlocs);
    MSC_matrix = nan(Num_subjects, Num_channels);
    for subj_idx = 1:Num_subjects
        subj_struct = get_subj(CoherenceResults_bandavgmsc_array, subj_idx);
        if ~isfield(subj_struct, 'msc_band')
            error('PlotSubjectsAndGroupMSC_Topography:MissingMscBand', ...
                ['Element %d of CoherenceResults_bandavgmsc_array has no .msc_band field. ' ...
                 'Was this produced by the current version of BandAverageMSC?'], subj_idx);
        end
        subj_msc_band = subj_struct.msc_band(:)';  % row vector, [1 x Num_channels]
        if numel(subj_msc_band) ~= Num_channels
            error('PlotSubjectsAndGroupMSC_Topography:ChanlocsMismatch', ...
                ['Element %d of CoherenceResults_bandavgmsc_array has %d channels in ' ...
                 '.msc_band, but chanlocs has %d elements. These must match.'], ...
                subj_idx, numel(subj_msc_band), Num_channels);
        end
        MSC_matrix(subj_idx, :) = subj_msc_band;
    end

    % ---- Resolve subject IDs -----------------------------------------------------------------
    if isempty(options.SubjectIDs)
        pad_width = max(2, floor(log10(Num_subjects)) + 1);
        SubjectIDs = arrayfun(@(idx) sprintf('Subject%0*d', pad_width, idx), ...
            1:Num_subjects, 'UniformOutput', false);
    elseif numel(options.SubjectIDs) ~= Num_subjects
        error('PlotSubjectsAndGroupMSC_Topography:SubjectIDMismatch', ...
            'options.SubjectIDs has %d entries, but there are %d subjects. These must match.', ...
            numel(options.SubjectIDs), Num_subjects);
    else
        SubjectIDs = options.SubjectIDs;
    end

    % ---- Warn (don't silently average over) any non-finite per-subject values ---------------
    nonfinite_mask = ~isfinite(MSC_matrix);
    if any(nonfinite_mask(:))
        warning('PlotSubjectsAndGroupMSC_Topography:NonFiniteInputs', ...
            ['%d entries across %d subjects'' .msc_band are non-finite. These will be ' ...
             'excluded via omitnan when computing the group average, which means different ' ...
             'channels may effectively be averaged over different numbers of subjects. ' ...
             'Consider investigating the upstream cause (see BandAverageMSC warnings).'], ...
            sum(nonfinite_mask(:)), sum(any(nonfinite_mask, 2)));
    end

    % ---- Compute the group-level average, per channel ----------------------------------------
    GroupMean_msc_band = mean(MSC_matrix, 1, 'omitnan');  % [1 x Num_channels]

    % ---- Determine the shared color scale (used for EVERY panel and EVERY saved PNG) --------
    if isempty(options.Clim)
        Clim = compute_shared_clim(MSC_matrix);
    else
        if numel(options.Clim) ~= 2
            error('PlotSubjectsAndGroupMSC_Topography:BadClim', ...
                'options.Clim must be a 2-element vector [lo, hi].');
        end
        Clim = options.Clim;
    end

    % ---- Build the combined overview figure ----------------------------------------------------
    Num_panels = Num_subjects + 1;  % one per subject, plus the group average
    if isempty(options.NumCols)
        Num_cols = ceil(sqrt(Num_panels));
    else
        Num_cols = options.NumCols;
    end
    Num_rows = ceil(Num_panels / Num_cols);

    fig_handle = render_subjects_and_group_figure(MSC_matrix, GroupMean_msc_band, chanlocs, ...
        SubjectIDs, options.GroupLabel, options.ColorbarLabel, options.FigureTitle, ...
        Clim, Num_rows, Num_cols, options.HighlightGroupPanel, options.Visible);

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
    cmap = parula(256);

    SubjectPNGFiles      = cell(Num_subjects, 1);
    SubjectFigureHandles = gobjects(Num_subjects, 1);
    for subj_idx = 1:Num_subjects
        safe_subj_id           = sanitize_for_filename(SubjectIDs{subj_idx});
        subj_png_path          = fullfile(options.OutputDir, [safe_subj_id, '_topography.png']);
        [subj_fig, subj_saved] = render_and_save_single_topography(MSC_matrix(subj_idx, :), ...
            chanlocs, SubjectIDs{subj_idx}, options.ColorbarLabel, Clim, cmap, subj_png_path, ...
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
    [group_fig, GroupPNGFile] = render_and_save_single_topography(GroupMean_msc_band, chanlocs, ...
        group_title, options.ColorbarLabel, Clim, cmap, group_png_path, options.IndividualVisible);
    GroupFigureHandle = group_fig;
    if options.IndividualCloseAfterSaving
        close(group_fig);
    end

    % ---- Assemble results structure -------------------------------------------------------------
    PlotResults = struct();
    PlotResults.MSC_matrix           = MSC_matrix;
    PlotResults.GroupMean_msc_band   = GroupMean_msc_band;
    PlotResults.SubjectIDs           = SubjectIDs;
    PlotResults.Clim                 = Clim;
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
    % Allow CoherenceResults_bandavgmsc_array to be passed in as either a cell array of structs
    % or a struct array (mirrors GroupLevel_MSC_Analysis.m's get_subject_struct).
    if iscell(arr)
        s = arr{idx};
    else
        s = arr(idx);
    end
end

function Clim = compute_shared_clim(MSC_matrix)
    % Default clim logic when the caller doesn't supply options.Clim explicitly: computed from
    % the combined range of every subject's .msc_band values (the group average, being an
    % average, will typically fall well within this range).
    finite_vals = MSC_matrix(isfinite(MSC_matrix));
    if isempty(finite_vals)
        % Degenerate case: every value is non-finite. Fall back to [0, 1], the theoretical range
        % of MSC, rather than erroring outright.
        Clim = [0, 1];
        return;
    end
    clim_min = min(finite_vals);
    clim_max = max(finite_vals) * 1.01;  % small headroom, for better contrast
    if clim_max <= clim_min
        % Degenerate case: all values identical. Avoid a zero-width color axis, which would
        % error in topoplot/clim.
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

function fig_handle = render_subjects_and_group_figure(MSC_matrix, GroupMean_msc_band, chanlocs, ...
        SubjectIDs, group_label, colorbar_label, fig_title, Clim, Num_rows, Num_cols, ...
        highlight_group_panel, fig_visible)
    % Build one figure containing one topoplot tile per subject, plus one final tile for the
    % group average, all sharing one color scale. See design notes 1) and 3) above for why the
    % tiledlayout/topoplot workarounds (axes_in_tile_slot, final colormap-reassertion pass) are
    % kept even though call_topoplot_custom likely no longer strictly needs them.
    Num_subjects = size(MSC_matrix, 1);
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

    cmap = parula(256);
    tile_axes = gobjects(Num_panels, 1);

    % ---- One tile per subject -----------------------------------------------------------------
    for subj_idx = 1:Num_subjects
        axes_in_tile_slot(fig_handle, t, subj_idx); % swapped for a standalone axes below
        tile_axes(subj_idx) = call_topoplot_custom(MSC_matrix(subj_idx, :), chanlocs, ...
            SubjectIDs{subj_idx}, colorbar_label, Clim(1), Clim(2), cmap);
    end

    % ---- Final tile: the group average ---------------------------------------------------------
    group_tile_index = Num_panels;
    group_tile_position = axes_in_tile_slot(fig_handle, t, group_tile_index, true); %#ok<NASGU>
    group_title = sprintf('%s (N=%d)', group_label, Num_subjects);
    tile_axes(group_tile_index) = call_topoplot_custom(GroupMean_msc_band, chanlocs, ...
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
    % call_topoplot_custom (see design note 1 above) -- kept for minimal-diff continuity.
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
    % design note 3 above), but it's kept for minimal-diff continuity.
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
