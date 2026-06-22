function [fig_handle, PlotResults] = PlotGroupMSC_Topography(MSC_condition1, MSC_condition2, chanlocs, options)
    % =========================================================================================
    % Plot Group-Averaged MSC Topographies: Condition 1, Condition 2, and Their Difference
    % =========================================================================================
    % Takes per-subject band-averaged MSC values (e.g., the output of BandAverageMSC, one row
    % per subject) for two conditions, averages each condition across subjects per channel, and
    % plots three topoplots side by side: the Condition 1 group average, the Condition 2 group
    % average, and their difference (Condition 1 - Condition 2). Uses the user-supplied
    % call_topoplot_eeglab wrapper around EEGLAB's topoplot for the actual plotting.
    %
    % Arguments:
    %   MSC_condition1: [Num_subjects x Num_channels] matrix. Each row is one subject's
    %                   BandAverageMSC output for condition 1. Subjects must be in the SAME
    %                   ORDER as MSC_condition2 (row i in both matrices = same subject), since
    %                   the difference plot is a paired, per-subject-then-averaged comparison.
    %   MSC_condition2: [Num_subjects x Num_channels] matrix, condition 2, same subject order
    %                   and same Num_channels as MSC_condition1.
    %   chanlocs:       EEGLAB channel locations structure, passed through to
    %                   call_topoplot_eeglab. Must have Num_channels elements.
    %
    % Optional name-value arguments:
    %   options.Condition1_Label (default 'Condition 1'): title label for the first topoplot.
    %   options.Condition2_Label (default 'Condition 2'): title label for the second topoplot.
    %   options.ColorbarLabel    (default 'MSC'):          colorbar y-axis label for the two
    %                                                       condition plots (the difference plot
    %                                                       is always labeled 'MSC difference').
    %   options.FigureTitle      (default ''):              overall figure title (sgtitle), shown
    %                                                       above all three topoplots. Skipped if empty.
    %   options.SavePath         (default ''):              if non-empty, the GROUP figure is saved
    %                                                       to this path (extension determines
    %                                                       format, e.g. '.png', '.fig', '.pdf').
    %                                                       If empty, the figure is only displayed.
    %
    %   --- Per-subject plotting (one 3-panel figure per subject, each saved as its own PNG) ---
    %   options.PlotIndividualSubjects (default false):    if true, also produce and save a
    %                                                       3-panel (condition1 / condition2 /
    %                                                       difference) figure for EACH subject
    %                                                       (each row of MSC_condition1/2), in
    %                                                       addition to the group-average figure.
    %   options.IndividualOutputDir (default ''):          folder to save per-subject PNGs into.
    %                                                       REQUIRED if PlotIndividualSubjects is
    %                                                       true (created automatically if it
    %                                                       doesn't already exist).
    %   options.SubjectIDs (default {}):                   cell array of char/string subject
    %                                                       identifiers, one per row of
    %                                                       MSC_condition1/2, used both in each
    %                                                       per-subject figure's title and in its
    %                                                       saved filename (sanitized -- spaces
    %                                                       and other filesystem-unfriendly
    %                                                       characters are replaced with '_').
    %                                                       Defaults to {'Subject01','Subject02',...}
    %                                                       if not supplied.
    %   options.IndividualFigureVisible (default false):   if false, per-subject figures are
    %                                                       created invisibly (faster, and avoids
    %                                                       flooding your screen with dozens of
    %                                                       figure windows during batch saving).
    %                                                       Set true if you want to watch them
    %                                                       appear as they're generated.
    %   options.CloseIndividualFigures (default true):     if true, each per-subject figure is
    %                                                       closed (after saving) before moving on
    %                                                       to the next subject, so MATLAB doesn't
    %                                                       accumulate one open figure per subject.
    %                                                       Set false if you want to keep them all
    %                                                       open/accessible after the call returns
    %                                                       (PlotResults.IndividualFigureHandles
    %                                                       will then contain valid handles).
    %
    % Returns:
    %   fig_handle:  handle to the created GROUP-AVERAGE figure (per-subject figures, if
    %                requested, are handled separately -- see PlotResults below).
    %   PlotResults: a structure containing:
    %     .GroupMean_condition1  [1 x Num_channels] group average, condition 1
    %     .GroupMean_condition2  [1 x Num_channels] group average, condition 2
    %     .GroupMean_diff        [1 x Num_channels] condition1 - condition2, computed from the
    %                            group averages (NOT averaged per-subject differences -- see
    %                            note below)
    %     .clim_conditions  [1 x 2] shared color limits used for the two GROUP condition topoplots
    %     .clim_diff        [1 x 2] symmetric color limits used for the GROUP difference topoplot
    %     .IndividualSavedFiles   {Num_subjects x 1} cell array of full paths to each saved
    %                             per-subject PNG (only populated if PlotIndividualSubjects=true)
    %     .IndividualFigureHandles [Num_subjects x 1] array of figure handles for each per-subject
    %                              figure (only populated if PlotIndividualSubjects=true; handles
    %                              will be invalid/closed if CloseIndividualFigures=true)
    %     .IndividualClims        {Num_subjects x 1} cell array, each cell a struct with fields
    %                             .clim_conditions and .clim_diff -- the per-subject color limits
    %                             actually used for that subject's figure (only populated if
    %                             PlotIndividualSubjects=true; see design note 4 below on why
    %                             these are computed PER SUBJECT rather than reusing the group's)
    %
    % -----------------------------------------------------------------------------------------
    % Design notes
    % -----------------------------------------------------------------------------------------
    % 1) Color scales: the two condition topoplots share a common color scale (clim_conditions),
    %    computed from the combined data range across BOTH conditions, so the two maps are
    %    visually comparable to each other. The difference topoplot uses its own, separate,
    %    zero-centered (symmetric) color scale, since differences can be positive or negative and
    %    "no difference" should map to the neutral center of a diverging colormap. MSC itself is
    %    naturally bounded in [0, 1], but the shared scale is still computed from the data
    %    (rather than hardcoded to [0,1]) for better contrast, consistent with how
    %    call_topoplot_eeglab's own commented-out defaults already do this.
    %
    % 2) Order of operations for the difference: this function averages each condition across
    %    subjects FIRST, then takes the difference of the two group averages
    %    (GroupMean_diff = GroupMean_condition1 - GroupMean_condition2). This is equivalent to
    %    averaging the per-subject paired differences across subjects (the two are
    %    mathematically identical: mean(A) - mean(B) == mean(A - B) for paired rows), so either
    %    framing is valid -- we compute it via the group means since that's what's plotted.
    %    Note this function does NOT run a significance test on the difference (that's what
    %    GroupLevel_MSC_Analysis is for); this function is for visualization of the same
    %    quantities GroupLevel_MSC_Analysis tests.
    %
    % 3) Diverging colormap: MATLAB does not ship a built-in red-blue diverging colormap, so
    %    a simple one is generated locally (red = condition1 > condition2, blue = condition1 <
    %    condition2, white = no difference) via the local helper function diverging_redblue_cmap.
    %
    % 4) Per-subject color limits: when PlotIndividualSubjects is true, each subject's figure
    %    gets its OWN color limits (computed from that subject's data only), rather than reusing
    %    the group-average clims. Single-subject MSC data is noisier and can have a very
    %    different range than the group average (e.g. one subject's max MSC might be far above
    %    or below the group mean's max), so forcing every subject onto the group's scale would
    %    wash out or saturate many individual plots. The tradeoff is that per-subject figures are
    %    NOT directly comparable to each other by color alone -- if you need that, pass a fixed
    %    .maplimits-style scale in yourself by post-processing PlotResults.IndividualClims, or ask
    %    for that as a future option.
    % =========================================================================================
    arguments
        MSC_condition1 (:,:) double
        MSC_condition2 (:,:) double
        chanlocs
        options.Condition1_Label (1,:) char = 'Condition 1'
        options.Condition2_Label (1,:) char = 'Condition 2'
        options.ColorbarLabel (1,:) char = 'MSC'
        options.FigureTitle (1,:) char = ''
        options.SavePath (1,:) char = ''
        options.PlotIndividualSubjects (1,1) logical = false
        options.IndividualOutputDir (1,:) char = ''
        options.SubjectIDs cell = {}
        options.IndividualFigureVisible (1,1) logical = false
        options.CloseIndividualFigures (1,1) logical = true
    end

    % ---- Validate input shapes -------------------------------------------------------------
    [Num_subjects_1, Num_channels_1] = size(MSC_condition1);
    [Num_subjects_2, Num_channels_2] = size(MSC_condition2);

    if Num_subjects_1 ~= Num_subjects_2
        error('PlotGroupMSC_Topography:MismatchedSubjectCount', ...
            ['MSC_condition1 has %d subjects (rows) but MSC_condition2 has %d subjects (rows). ' ...
             'These must match, and rows must be in the same subject order, since the ' ...
             'difference plot assumes row i in both matrices is the same subject.'], ...
            Num_subjects_1, Num_subjects_2);
    end
    if Num_channels_1 ~= Num_channels_2
        error('PlotGroupMSC_Topography:MismatchedChannelCount', ...
            ['MSC_condition1 has %d channels (columns) but MSC_condition2 has %d channels ' ...
             '(columns). These must match.'], Num_channels_1, Num_channels_2);
    end
    if Num_channels_1 ~= numel(chanlocs)
        error('PlotGroupMSC_Topography:ChanlocsMismatch', ...
            ['MSC_condition1/MSC_condition2 have %d channels, but chanlocs has %d elements. ' ...
             'These must match for topoplot to place values at the correct electrode ' ...
             'positions.'], Num_channels_1, numel(chanlocs));
    end
    if Num_subjects_1 < 1
        error('PlotGroupMSC_Topography:NoSubjects', 'MSC_condition1/MSC_condition2 have zero rows (no subjects).');
    end

    % ---- Validate per-subject plotting options ----------------------------------------------
    Num_subjects = Num_subjects_1;
    if options.PlotIndividualSubjects
        if isempty(options.IndividualOutputDir)
            error('PlotGroupMSC_Topography:MissingOutputDir', ...
                ['options.PlotIndividualSubjects is true, but options.IndividualOutputDir was ' ...
                 'not supplied. Specify a folder to save the per-subject PNGs into.']);
        end
        if ~isfolder(options.IndividualOutputDir)
            % Create the output folder if it doesn't exist yet, rather than erroring -- this is
            % the common case (e.g. a fresh "figures/subjects" folder for a new analysis run).
            mkdir(options.IndividualOutputDir);
        end
        if isempty(options.SubjectIDs)
            % Default subject IDs: Subject01, Subject02, ... (zero-padded to the width needed
            % for the largest subject number, so filenames sort correctly even for 100+ subjects).
            pad_width = max(2, floor(log10(Num_subjects)) + 1);
            options.SubjectIDs = arrayfun(@(idx) sprintf('Subject%0*d', pad_width, idx), ...
                1:Num_subjects, 'UniformOutput', false);
        elseif numel(options.SubjectIDs) ~= Num_subjects
            error('PlotGroupMSC_Topography:SubjectIDMismatch', ...
                ['options.SubjectIDs has %d entries, but MSC_condition1/MSC_condition2 have %d ' ...
                 'subjects (rows). These must match.'], numel(options.SubjectIDs), Num_subjects);
        end
    end

    % ---- Warn (don't silently average over) any non-finite per-subject values --------------
    nonfinite_1 = ~isfinite(MSC_condition1);
    nonfinite_2 = ~isfinite(MSC_condition2);
    if any(nonfinite_1(:)) || any(nonfinite_2(:))
        warning('PlotGroupMSC_Topography:NonFiniteInputs', ...
            ['%d entries in MSC_condition1 and %d entries in MSC_condition2 are non-finite. ' ...
             'These will be excluded via omitnan when computing group averages, which means ' ...
             'different channels may effectively be averaged over different numbers of ' ...
             'subjects. Consider investigating the upstream cause (see BandAverageMSC warnings).'], ...
            sum(nonfinite_1(:)), sum(nonfinite_2(:)));
    end

    % ---- Compute group averages, per channel (omitting any non-finite subject values) ------
    GroupMean_condition1 = mean(MSC_condition1, 1, 'omitnan');  % [1 x Num_channels]
    GroupMean_condition2 = mean(MSC_condition2, 1, 'omitnan');  % [1 x Num_channels]

    % Equivalent to mean(MSC_condition1 - MSC_condition2, 1, 'omitnan') only when no row has a
    % NaN in just one condition; computing it from the already-omitnan'd group means here is the
    % simplest choice for plotting purposes (see design note 2 above for the algebraic identity
    % in the all-finite case).
    GroupMean_diff = GroupMean_condition1 - GroupMean_condition2;  % [1 x Num_channels]

    % ---- Determine color limits for the group figure ----------------------------------------
    [clim_conditions, clim_diff] = compute_topography_clims(GroupMean_condition1, ...
        GroupMean_condition2, GroupMean_diff);

    % ---- Build and save the group-average figure ---------------------------------------------
    fig_handle = render_three_panel_topography_figure(GroupMean_condition1, GroupMean_condition2, ...
        GroupMean_diff, chanlocs, options.Condition1_Label, options.Condition2_Label, ...
        options.ColorbarLabel, options.FigureTitle, clim_conditions, clim_diff, true);

    if ~isempty(options.SavePath)
        % exportgraphics handles most common formats (.png, .pdf, .jpg, .tif) based on the file
        % extension; for .fig specifically, use savefig instead, since exportgraphics does not
        % support MATLAB's native .fig format.
        [~, ~, ext] = fileparts(options.SavePath);
        if strcmpi(ext, '.fig')
            savefig(fig_handle, options.SavePath);
        else
            exportgraphics(fig_handle, options.SavePath);
        end
    end

    % ---- Per-subject figures, if requested -----------------------------------------------------
    % See design note 4 above: each subject's color limits are computed from that subject's own
    % data, not reused from the group, since a single subject's MSC range can differ substantially
    % from the group average's range.
    PlotResults = struct();
    if options.PlotIndividualSubjects
        IndividualSavedFiles    = cell(Num_subjects, 1);
        IndividualFigureHandles = gobjects(Num_subjects, 1);
        IndividualClims         = cell(Num_subjects, 1);

        for subj_idx = 1:Num_subjects
            subj_id    = options.SubjectIDs{subj_idx};
            subj_data1 = MSC_condition1(subj_idx, :);
            subj_data2 = MSC_condition2(subj_idx, :);
            subj_diff  = subj_data1 - subj_data2;

            [subj_clim_conditions, subj_clim_diff] = compute_topography_clims(subj_data1, ...
                subj_data2, subj_diff);

            subj_fig_title = subj_id;
            if ~isempty(options.FigureTitle)
                subj_fig_title = sprintf('%s -- %s', options.FigureTitle, subj_id);
            end

            subj_fig = render_three_panel_topography_figure(subj_data1, subj_data2, subj_diff, ...
                chanlocs, options.Condition1_Label, options.Condition2_Label, ...
                options.ColorbarLabel, subj_fig_title, subj_clim_conditions, subj_clim_diff, ...
                options.IndividualFigureVisible);

            % Sanitize the subject ID for use as a filename: replace anything that isn't a
            % letter, digit, hyphen, or underscore with '_' (handles spaces, slashes, etc. if a
            % caller supplies free-text SubjectIDs rather than the auto-generated 'SubjectNN').
            safe_subj_id = regexprep(subj_id, '[^a-zA-Z0-9_-]', '_');
            subj_filename = fullfile(options.IndividualOutputDir, [safe_subj_id, '_topography.png']);
            exportgraphics(subj_fig, subj_filename);

            IndividualSavedFiles{subj_idx}     = subj_filename;
            IndividualFigureHandles(subj_idx)  = subj_fig;
            IndividualClims{subj_idx}          = struct('clim_conditions', subj_clim_conditions, ...
                                                          'clim_diff', subj_clim_diff);

            if options.CloseIndividualFigures
                close(subj_fig);
            end
        end

        PlotResults.IndividualSavedFiles    = IndividualSavedFiles;
        PlotResults.IndividualFigureHandles = IndividualFigureHandles;
        PlotResults.IndividualClims         = IndividualClims;
    end

    % ---- Assemble results structure ----------------------------------------------------------
    PlotResults.GroupMean_condition1 = GroupMean_condition1;
    PlotResults.GroupMean_condition2 = GroupMean_condition2;
    PlotResults.GroupMean_diff       = GroupMean_diff;
    PlotResults.clim_conditions      = clim_conditions;
    PlotResults.clim_diff            = clim_diff;
end

% =============================================================================================
% Local helper functions
% =============================================================================================

function [clim_conditions, clim_diff] = compute_topography_clims(data1, data2, data_diff)
    % Shared clim logic for both the group-average figure and each per-subject figure (see
    % design notes 1 and 4 in the main function header for the rationale).
    %
    %   clim_conditions: shared color limits for the condition-1 and condition-2 topoplots,
    %                     computed from the combined range of BOTH inputs.
    %   clim_diff:        symmetric (zero-centered) color limits for the difference topoplot.
    combined_range = [data1, data2];
    clim_min_conditions = min(combined_range);
    clim_max_conditions = max(combined_range) * 1.1;  % 10% headroom, for better contrast
    if clim_max_conditions <= clim_min_conditions
        % Degenerate case: all values identical (e.g., all zero). Avoid a zero-width color axis,
        % which would error in topoplot/clim.
        clim_max_conditions = clim_min_conditions + eps(clim_min_conditions) * 100 + 1e-6;
    end
    clim_conditions = [clim_min_conditions, clim_max_conditions];

    max_abs_diff = max(abs(data_diff));
    if max_abs_diff == 0
        % Degenerate case: conditions are identical at every channel. Avoid a zero-width axis.
        max_abs_diff = 1e-6;
    end
    clim_diff = [-max_abs_diff, max_abs_diff];
end

function fig_handle = render_three_panel_topography_figure(data1, data2, data_diff, chanlocs, ...
        label1, label2, colorbar_label, fig_title, clim_conditions, clim_diff, fig_visible)
    % Build one figure containing three side-by-side topoplots (condition 1, condition 2, and
    % their difference) using the user-supplied call_topoplot_eeglab wrapper. Used for both the
    % group-average figure and each individual subject's figure.
    %
    % NOTE ON COLORMAPS: EEGLAB's topoplot() does not cleanly scope its colormap handling to the
    % current axes -- internally it reads/sets colormap state in a way that can affect later
    % calls (see https://github.com/sccn/eeglab, topoplot.m, and EEGLAB's icadefs.m default
    % colormap mechanism). In practice this means that if we set each tile's colormap
    % immediately after that tile's topoplot() call, a LATER tile's topoplot() call can still
    % retroactively clobber an EARLIER tile's colormap by the time the figure is actually
    % rendered. The reliable fix is to capture each tile's axes handle as we go, and reassert
    % colormap(ax, cmap) for every tile in one final pass AFTER all three topoplot() calls have
    % completed -- nothing runs after that final pass to overwrite it again.
    if fig_visible
        visible_str = 'on';
    else
        visible_str = 'off';
    end
    fig_handle = figure('Color', 'white', 'Position', [100 100 1500 500], 'Visible', visible_str);
    t = tiledlayout(fig_handle, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    cmap_conditions = parula(256);
    cmap_diff       = diverging_redblue_cmap(256);

    % NOTE ON THE "Unable to set 'Position' ... for objects in a TiledChartLayout" WARNING:
    % EEGLAB's topoplot() internally resizes/squares its own axes by calling something like
    % set(gca,'Position',...) (e.g. to make room for the head cartoon / keep the plot square).
    % An axes that lives directly inside a TiledChartLayout tile is managed by the layout and
    % will NOT accept an external Position change -- MATLAB just warns and ignores the request.
    % That means topoplot's own resizing silently fails to apply when called naively via
    % nexttile(t).
    %
    % Fix: use nexttile(t, idx) with an EXPLICIT tile index to ask the layout for that slot's
    % geometry, then replace that managed tile axes with an ordinary standalone axes placed at
    % the same Position. A standalone axes fully supports Position changes, so topoplot's
    % internal resizing works correctly and the warning goes away. The local helper
    % axes_in_tile_slot() below does this (see its comments for why the index must be explicit).
    ax1 = axes_in_tile_slot(fig_handle, t, 1); %#ok<NASGU> % overwritten by call_topoplot_eeglab below
    ax1 = call_topoplot_eeglab(data1, chanlocs, label1, ...
        colorbar_label, clim_conditions(1), clim_conditions(2), cmap_conditions);

    ax2 = axes_in_tile_slot(fig_handle, t, 2); %#ok<NASGU>
    ax2 = call_topoplot_eeglab(data2, chanlocs, label2, ...
        colorbar_label, clim_conditions(1), clim_conditions(2), cmap_conditions);

    ax3 = axes_in_tile_slot(fig_handle, t, 3); %#ok<NASGU>
    diff_title = sprintf('%s -minus- %s', label1, label2);
    ax3 = call_topoplot_eeglab(data_diff, chanlocs, diff_title, ...
        'MSC difference', clim_diff(1), clim_diff(2), cmap_diff);

    % Final pass: reassert each tile's colormap now that all three topoplot() calls are done,
    % so no subsequent topoplot() call can overwrite an earlier tile's colormap.
    colormap(ax1, cmap_conditions);
    colormap(ax2, cmap_conditions);
    colormap(ax3, cmap_diff);
    % Also reassert each axes' CLim, since some MATLAB/EEGLAB versions tie color-limit state to
    % colormap state internally and a late colormap(ax,...) call can occasionally reset CLim to
    % the axes' data range rather than the explicit limits we want.
    clim(ax1, clim_conditions);
    clim(ax2, clim_conditions);
    clim(ax3, clim_diff);

    if ~isempty(fig_title)
        % NOTE: sgtitle's target argument must be a Figure, Panel, or Tab -- it does NOT accept a
        % TiledChartLayout handle. The documented way to add an overall title to a tiledlayout
        % itself is title(t, ...), which is what we use here.
        title(t, fig_title, 'FontSize', 16, 'FontWeight', 'bold');
    end
end

function ax = axes_in_tile_slot(fig_handle, t, tile_index)
    % Reserve a SPECIFIC slot (by index) in tiledlayout t, then hand back a plain, standalone
    % axes occupying that exact same screen position -- instead of the tile axes itself.
    %
    % Why: a TiledChartLayout manages the Position of any axes placed directly in one of its
    % tiles, and will reject (with a warning) any external attempt to change that axes' Position.
    % EEGLAB's topoplot() does exactly that internally (e.g. to square up the plot / make room
    % for the head cartoon), which is the source of the
    %   "Unable to set 'Position', 'InnerPosition', 'OuterPosition', or 'PositionConstraint'
    %    for objects in a TiledChartLayout"
    % warning. A standalone axes (one that is NOT a tiledlayout child) has no such restriction,
    % so we ask the layout for a tile's geometry, then immediately swap in a normal axes at that
    % same spot.
    %
    % IMPORTANT: tile_index must be passed explicitly (nexttile(t, tile_index)), rather than
    % relying on nexttile(t)'s auto-advance. Auto-advance picks the next slot based on how many
    % tiles currently exist as children of the layout -- but since this function immediately
    % DELETES the tile axes after reading its Position (to replace it with a standalone axes),
    % the child count drops back down after every call. That fools auto-advance into handing out
    % the SAME slot (slot 1) every time, which is what caused all three topoplots to land on top
    % of each other. Requesting an explicit index sidesteps that entirely.
    nexttile(t, tile_index);
    tile_ax = gca;
    drawnow;   % force the layout to actually compute/finalize this tile's Position before we
               % read it -- tile geometry is not guaranteed to be set until a draw occurs
    tile_position = get(tile_ax, 'Position');   % figure-normalized [left bottom width height]
    delete(tile_ax);
    ax = axes('Parent', fig_handle, 'Position', tile_position);
    axes(ax);  % make it the current axes, since topoplot()/call_topoplot_eeglab() draw via gca
end

function cmap = diverging_redblue_cmap(n)
    % Generate a simple symmetric red-white-blue diverging colormap with n entries.
    % Blue (low/negative) -> white (zero/center) -> red (high/positive).
    % MATLAB does not ship a built-in diverging colormap in base MATLAB (parula, jet, hot, etc.
    % are all sequential), so this is constructed manually via simple linear interpolation
    % through three anchor colors.
    if mod(n, 2) == 0
        half = n/2;
        lower = [linspace(0.0, 1.0, half)', linspace(0.0, 1.0, half)', linspace(0.7, 1.0, half)'];  % blue -> white
        upper = [linspace(1.0, 0.8, half)', linspace(1.0, 0.0, half)', linspace(1.0, 0.0, half)'];  % white -> red
        cmap = [lower; upper];
    else
        half = (n-1)/2;
        lower = [linspace(0.0, 1.0, half)', linspace(0.0, 1.0, half)', linspace(0.7, 1.0, half)'];
        upper = [linspace(1.0, 0.8, half)', linspace(1.0, 0.0, half)', linspace(1.0, 0.0, half)'];
        cmap = [lower; [1 1 1]; upper];  % explicit pure-white center entry for odd n
    end
end
