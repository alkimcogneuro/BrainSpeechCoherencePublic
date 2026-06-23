function [fig_handle, clim_conditions, clim_diff] = PlotThreePanelTopography(data1, data2, chanlocs, options)
    % =========================================================================================
    % Plot Three Side-by-Side Topographies: Condition 1, Condition 2, and Their Difference
    % =========================================================================================
    % Core, reusable plotting primitive: given a single condition-1 vector and a single
    % condition-2 vector (one value per channel), plots three EEGLAB topoplots side by side --
    % condition 1, condition 2, and condition1-minus-condition2 -- in one figure, using the
    % user-supplied call_topoplot_eeglab wrapper around EEGLAB's topoplot.
    %
    % This function does NOT know or care whether data1/data2 came from a group average or a
    % single subject -- that distinction lives entirely in the caller. This is what makes it
    % possible to reuse the exact same plotting code for both:
    %   - PlotGroupMSC_Topography.m:        averages across subjects, then calls this function
    %                                        once on the group means.
    %   - a per-subject looping function:   calls this function once per subject, on each row
    %                                        of the subjects x channels data matrices.
    %
    % Arguments:
    %   data1:     [1 x Num_channels] (or [Num_channels x 1]) vector, condition 1 values.
    %   data2:     [1 x Num_channels] (or [Num_channels x 1]) vector, condition 2 values, same
    %              length and channel order as data1.
    %   chanlocs:  EEGLAB channel locations structure, passed through to call_topoplot_eeglab.
    %              Must have Num_channels elements.
    %
    % Optional name-value arguments:
    %   options.Label1         (default 'Condition 1'): title for the first topoplot.
    %   options.Label2         (default 'Condition 2'): title for the second topoplot.
    %   options.ColorbarLabel  (default 'MSC'):          colorbar y-axis label for the two
    %                                                     condition topoplots (the difference
    %                                                     topoplot is always labeled
    %                                                     'MSC difference').
    %   options.FigureTitle    (default ''):             overall figure title, shown above all
    %                                                     three topoplots. Skipped if empty.
    %   options.ClimConditions (default []):             [1 x 2] color limits to use for the two
    %                                                     condition topoplots. If left empty
    %                                                     (the default), limits are computed
    %                                                     automatically from data1/data2 (see
    %                                                     compute_topography_clims below). Pass
    %                                                     this explicitly if you want several
    %                                                     calls (e.g. one per subject) to all
    %                                                     share the same color scale.
    %   options.ClimDiff       (default []):             [1 x 2] color limits for the difference
    %                                                     topoplot. Same auto-vs-explicit behavior
    %                                                     as ClimConditions.
    %   options.Visible        (default true):           if false, the figure is created
    %                                                     invisibly (e.g. for batch saving many
    %                                                     per-subject figures without flooding
    %                                                     your screen with windows).
    %   options.SavePath       (default ''):              if non-empty, the figure is saved to
    %                                                     this path immediately (extension
    %                                                     determines format, e.g. '.png', '.fig',
    %                                                     '.pdf'). If empty, the figure is only
    %                                                     created/displayed, and it's the caller's
    %                                                     responsibility to save or close it.
    %
    % Returns:
    %   fig_handle:      handle to the created figure.
    %   clim_conditions: [1 x 2] color limits actually used for the two condition topoplots
    %                    (either what you passed in via options.ClimConditions, or what was
    %                    computed automatically).
    %   clim_diff:       [1 x 2] color limits actually used for the difference topoplot.
    %
    % -----------------------------------------------------------------------------------------
    % Design notes
    % -----------------------------------------------------------------------------------------
    % 1) Color scales: by default, the two condition topoplots share a common color scale
    %    (clim_conditions), computed from the combined data range across BOTH conditions, so the
    %    two maps are visually comparable to each other. The difference topoplot uses its own,
    %    separate, zero-centered (symmetric) color scale, since differences can be positive or
    %    negative and "no difference" should map to the neutral center of a diverging colormap.
    %
    % 2) Diverging colormap: MATLAB does not ship a built-in red-blue diverging colormap, so a
    %    simple one is generated locally (red = data1 > data2, blue = data1 < data2, white = no
    %    difference) via the local helper function diverging_redblue_cmap.
    %
    % 3) TiledChartLayout + topoplot: EEGLAB's topoplot() internally resizes/squares its own
    %    axes (e.g. to make room for the head cartoon), which conflicts with axes managed
    %    directly by a TiledChartLayout tile (MATLAB will warn
    %    "Unable to set 'Position' ... for objects in a TiledChartLayout" and ignore the resize).
    %    The local helper axes_in_tile_slot() works around this by reading each tile's intended
    %    position, then swapping in an ordinary standalone axes at that same spot, which fully
    %    supports Position changes.
    % =========================================================================================
    arguments
        data1 (1,:) double
        data2 (1,:) double
        chanlocs
        options.Label1 (1,:) char = 'Condition 1'
        options.Label2 (1,:) char = 'Condition 2'
        options.ColorbarLabel (1,:) char = 'MSC'
        options.FigureTitle (1,:) char = ''
        options.ClimConditions (:,:) double = []
        options.ClimDiff (:,:) double = []
        options.Visible (1,1) logical = true
        options.SavePath (1,:) char = ''
    end

    % ---- Validate input shapes -------------------------------------------------------------
    if numel(data1) ~= numel(data2)
        error('PlotThreePanelTopography:MismatchedLength', ...
            'data1 has %d elements but data2 has %d elements. These must match.', ...
            numel(data1), numel(data2));
    end
    if numel(data1) ~= numel(chanlocs)
        error('PlotThreePanelTopography:ChanlocsMismatch', ...
            ['data1/data2 have %d elements, but chanlocs has %d elements. These must match ' ...
             'for topoplot to place values at the correct electrode positions.'], ...
            numel(data1), numel(chanlocs));
    end
    if ~isempty(options.ClimConditions) && numel(options.ClimConditions) ~= 2
        error('PlotThreePanelTopography:BadClimConditions', ...
            'options.ClimConditions must be a 2-element vector [lo, hi].');
    end
    if ~isempty(options.ClimDiff) && numel(options.ClimDiff) ~= 2
        error('PlotThreePanelTopography:BadClimDiff', ...
            'options.ClimDiff must be a 2-element vector [lo, hi].');
    end

    data_diff = data1 - data2;

    % ---- Determine color limits (auto-compute any not explicitly supplied) -----------------
    [auto_clim_conditions, auto_clim_diff] = compute_topography_clims(data1, data2, data_diff);
    if isempty(options.ClimConditions)
        clim_conditions = auto_clim_conditions;
    else
        clim_conditions = options.ClimConditions;
    end
    if isempty(options.ClimDiff)
        clim_diff = auto_clim_diff;
    else
        clim_diff = options.ClimDiff;
    end

    % ---- Build the figure --------------------------------------------------------------------
    fig_handle = render_three_panel_topography_figure(data1, data2, data_diff, chanlocs, ...
        options.Label1, options.Label2, options.ColorbarLabel, options.FigureTitle, ...
        clim_conditions, clim_diff, options.Visible);

    % ---- Save the figure, if requested ------------------------------------------------------
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
end

% =============================================================================================
% Local helper functions
% =============================================================================================

function [clim_conditions, clim_diff] = compute_topography_clims(data1, data2, data_diff)
    % Default clim logic, used whenever the caller doesn't supply explicit ClimConditions /
    % ClimDiff (see design note 1 above for the rationale).
    %
    %   clim_conditions: shared color limits for the condition-1 and condition-2 topoplots,
    %                     computed from the combined range of BOTH inputs.
    %   clim_diff:        symmetric (zero-centered) color limits for the difference topoplot.
    combined_range = [data1, data2];
    clim_min_conditions = min(combined_range);
    clim_max_conditions = max(combined_range) * 1.01;  % 10%xx 5% headroom, for better contrast
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
    % their difference) using the user-supplied call_topoplot_eeglab wrapper.
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
    diff_title = format_diff_title(label1, label2);
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
    % the SAME slot (slot 1) every time, which would cause all three topoplots to land on top of
    % each other. Requesting an explicit index sidesteps that entirely.
    nexttile(t, tile_index);
    tile_ax = gca;
    drawnow;   % force the layout to actually compute/finalize this tile's Position before we
               % read it -- tile geometry is not guaranteed to be set until a draw occurs
    tile_position = get(tile_ax, 'Position');   % figure-normalized [left bottom width height]
    delete(tile_ax);
    ax = axes('Parent', fig_handle, 'Position', tile_position);
    axes(ax);  % make it the current axes, since topoplot()/call_topoplot_eeglab() draw via gca
end

function diff_title = format_diff_title(label1, label2)
    % Build the title for the difference panel, wrapping it onto multiple lines if it's long.
    %
    % Why: title() happily accepts a cell array of strings, rendering one line per cell -- but
    % if we always used a single long line like 'Eyes Open -minus- Eyes Closed', MATLAB expands
    % the title vertically as needed to avoid clipping it, and since this is the rightmost panel,
    % that vertical expansion crowds into (and can visually overlap) the colorbar legend sitting
    % just to its right. Splitting into a few shorter lines keeps the title's bounding box
    % narrower, which avoids that overlap.
    combined = sprintf('%s -minus- %s', label1, label2);

    % Threshold is a simple character-count heuristic, not a precise text-width measurement
    % (getting an exact rendered width would require a real text extent query, which depends on
    % font, DPI, and figure size at render time) -- but it works well in practice for the kinds
    % of condition labels this is meant for.
    LONG_TITLE_THRESHOLD = 25;

    if length(combined) <= LONG_TITLE_THRESHOLD
        diff_title = combined;
    else
        % Three short lines tend to stay narrower than two long ones for typical condition
        % labels, so split out '-minus-' onto its own middle line rather than appending it to
        % either label.
        diff_title = {label1, '-minus-', label2};
    end
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
