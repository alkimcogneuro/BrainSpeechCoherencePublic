function [ax] = call_topoplot_custom(coherence_vals, chanlocs, plot_title, y_label, clim_min, clim_max, cmap, options)
    % =========================================================================================
    % call_topoplot_custom.m -- Drop-In Replacement for call_topoplot_eeglab.m, No EEGLAB Needed
    % =========================================================================================
    % Same first 7 positional arguments, same "draws into the current axes and returns the axes
    % handle" contract as call_topoplot_eeglab.m -- every existing Plot*Topography.m file in this
    % repo should be able to switch by changing ONLY the function name at the call site. Nothing
    % else about how you build tiledlayouts, call axes_in_tile_slot, etc. needs to change.
    %
    % Renders the scalp map itself (interpolation, head outline, electrode markers, contours) --
    % no EEGLAB toolbox call anywhere in this file. Built on the same theta/radius -> x/y
    % projection and drawing conventions already used elsewhere in this repo (plot_topomap.m,
    % BuildChannelAdjacency.m, PlotChannelAdjacency.m), so all your custom scalp-map code now
    % shares one consistent coordinate system.
    %
    % Why this exists / what it fixes relative to call_topoplot_eeglab.m:
    %   1) THE MULTI-TILE COLORMAP BUG, structurally: call_topoplot_eeglab.m's own header
    %      documents that EEGLAB's topoplot() "does not cleanly scope its own colormap handling
    %      to a single axes," so a LATER topoplot() call (for a different tile) can retroactively
    %      change an EARLIER tile's rendered colormap -- which is why every Plot*Topography.m
    %      file in this repo has to duplicate a "final pass: reassert colormap(ax,cmap) and
    %      clim(ax,[...]) for every tile after all of them are drawn" workaround. This function
    %      never touches EEGLAB's rendering path at all: every colormap(...)/clim(...) call below
    %      is made with an EXPLICIT axes handle (colormap(ax, cmap), not colormap(cmap)), which is
    %      properly per-axes-scoped MATLAB behavior. That means the bug this function fixes can't
    %      happen here BY CONSTRUCTION -- not "usually doesn't happen," structurally can't. The
    %      existing "final pass" loops in your tiling code become REDUNDANT once you switch to
    %      this function (each call already sets its own tile's colormap/clim correctly, so a
    %      later re-assertion is a harmless no-op) -- safe to leave them in place during migration
    %      and delete them later once you've confirmed this function behaves as expected.
    %   2) Built-in cluster/significance overlays (options.HighlightChannels, options.ClusterChannels
    %      below) -- for marking individually-significant channels (e.g. from
    %      MaxStatisticPermutationTest_ShuffleTest.m) or outlining spatial clusters (e.g. from
    %      ClusterPermutationTest_ShuffleTest.m) directly on the topomap, without a separate
    %      post-processing step.
    %   3) No EEGLAB dependency: this only needs chanlocs(...).theta/.radius/.labels -- the same
    %      fields plot_topomap.m already relies on -- so plotting no longer requires the EEGLAB
    %      toolbox to be on the path at all.
    %
    % Arguments (first 7, positional, matching call_topoplot_eeglab.m exactly):
    %   coherence_vals: 1 x Num_channels (or Num_channels x 1) vector of values to plot.
    %   chanlocs:       chanlocs struct array (needs .theta, .radius, .labels per channel).
    %   plot_title:     axes title string.
    %   y_label:        colorbar label string.
    %   clim_min, clim_max: color axis limits.
    %   cmap:           Nx3 colormap matrix. Optional -- defaults to parula(256), same as
    %                   call_topoplot_eeglab.m's default, if omitted or passed as [].
    %
    % Optional name-value arguments (new -- these have no equivalent in call_topoplot_eeglab.m,
    % so omitting all of them reproduces that function's original visual behavior as closely as
    % a from-scratch renderer can):
    %   options.ShowLabels    (default true): draw channel-label text next to each electrode
    %                         (matches call_topoplot_eeglab.m's unconditional 'labelpoint' style).
    %   options.ShowContours  (default true): overlay black contour lines on the filled map
    %                         (matches EEGLAB topoplot()'s own 'both' default style).
    %   options.NumContours   (default 6): number of contour levels, if ShowContours is true.
    %   options.HeadRadius    (default 0.5): radius, in the same normalized units as chanlocs'
    %                         projected theta/radius coordinates, of the drawn head/nose/ears
    %                         CARTOON -- this is a fixed anatomical reference circle and does NOT
    %                         automatically grow to fit every electrode (see options.PlotRadius
    %                         and the design note below for why that's intentional).
    %   options.PlotRadius    (default NaN = auto): radius the INTERPOLATED COLOR MAP extends to
    %                         and is drawn within. Auto (NaN) computes this as the farthest valid
    %                         electrode's radius (plus a small margin), so every electrode in the
    %                         montage -- however peripheral -- ends up covered by colored/
    %                         interpolated data rather than floating in blank space beyond it.
    %                         Pass a number to override. See design note below: this is exactly
    %                         EEGLAB topoplot()'s own headrad-vs-plotrad distinction.
    %   options.GridResolution (default 300): interpolation grid is GridResolution x
    %                         GridResolution. Higher = smoother map, slower to draw.
    %   options.InterpMethod  (default 'v4'): passed to griddata(); 'v4' (biharmonic spline) is
    %                         the same choice plot_topomap.m already uses, and tends to look
    %                         smoother than 'linear' for sparse EEG electrode counts.
    %   options.HighlightChannels (default []): vector of channel INDICES to mark with an open
    %                         ring around their electrode dot (e.g. individually-significant
    %                         channels from MaxStatisticPermutationTest_ShuffleTest.m's
    %                         ChannelSignificant mask -- pass find(ms.ChannelSignificant)).
    %   options.HighlightColor (default [0 0 0], black): color of the HighlightChannels rings.
    %   options.ClusterChannels (default {}): cell array of channel-INDEX vectors, one per
    %                         cluster to outline (e.g. {ClusterResults.Clusters{ClusterResults.ClusterSignificant}}
    %                         to show only the significant clusters from
    %                         ClusterPermutationTest_ShuffleTest.m). Each cluster gets a halo
    %                         drawn around every one of its member electrodes (robust for a
    %                         cluster of any size, including a single channel), plus -- for
    %                         clusters of 3 or more channels -- a convex-hull outline around the
    %                         whole group, inflated slightly outward so it doesn't just trace
    %                         through the electrode dots themselves.
    %   options.ClusterColors (default []): Num_clusters x 3 color matrix, one row per entry in
    %                         ClusterChannels. Defaults to MATLAB's `lines` palette if not given.
    %   options.ClusterOutlineMargin (default 0.035): how far (in the same normalized units as
    %                         the head radius) the per-channel halos and the hull outline are
    %                         inflated outward from the electrode positions / cluster centroid.
    %   options.Axes          (default []): draw into this specific axes instead of gca. Same
    %                         effect as calling axes_in_tile_slot(...) before this function (the
    %                         pattern the rest of this repo already uses) -- provided as an
    %                         alternative for callers that would rather pass the axes explicitly.
    %
    % Returns:
    %   ax: the axes handle drawn into (options.Axes if given, otherwise gca) -- same contract as
    %       call_topoplot_eeglab.m, so ax_all(s) = call_topoplot_custom(...) inside a tiling loop
    %       works exactly like it did before.
    %
    % -----------------------------------------------------------------------------------------
    % Design notes
    % -----------------------------------------------------------------------------------------
    % 1) This is a genuinely different renderer from EEGLAB's topoplot() -- same overall visual
    %    language (filled color map, contour lines, head/nose/ears outline, labeled electrodes),
    %    but the specific interpolation/extrapolation math will not be pixel-identical to
    %    topoplot()'s. Compare a plot from this function against the equivalent
    %    call_topoplot_eeglab.m output side by side before fully switching over -- not run in an
    %    actual MATLAB session (none was reachable while writing this), so this has been
    %    syntax-checked only.
    % 2) HeadRadius vs. PlotRadius: EEG montages that extend beyond the classic 10-20 layout
    %    (e.g. 10-10/10-5 montages with sites like Fp1/Fp2, AF7/AF8, P9/P10, Iz) routinely include
    %    electrodes whose polar `radius` exceeds 0.5 -- they sit below the "equator" that the
    %    standard head-cartoon circle represents (near the inion, low frontal/temporal sites,
    %    etc.), not on top of the scalp as the 2-D head drawing depicts it. EEGLAB's topoplot()
    %    handles this with two SEPARATE radii: a fixed 'headrad' for the drawn head/nose/ears
    %    cartoon (an anatomical reference, not meant to grow), and a 'plotrad' that auto-expands
    %    to cover every electrode's actual position for the interpolated color map. This function
    %    does the same: options.HeadRadius stays fixed (the cartoon), while options.PlotRadius
    %    auto-expands (the color map + electrode markers/labels) so peripheral electrodes still
    %    have colored data under them and aren't just floating past the edge of a mostly-blank
    %    grid. Seeing electrode labels for sites like Fp1/Fp2, AF7/AF8, P9/P10, Iz sitting outside
    %    the drawn head circle is EXPECTED for an extended montage -- it reflects where those
    %    sites actually are relative to the standardized head outline, the same way EEGLAB's own
    %    topoplot() would show them. What earlier versions of this function got wrong was hard-
    %    coding the grid/mask to HeadRadius instead of a separate, auto-expanding PlotRadius,
    %    which left those channels outside the interpolated area entirely (blank space under
    %    their labels) rather than properly covered.
    % 3) Cluster outlines are drawn as: (a) a small halo circle around every member channel,
    %    always, regardless of cluster size, and (b) an inflated convex-hull polygon around the
    %    whole cluster, only when it has 3+ non-collinear channels (a hull isn't a meaningful
    %    shape for 1-2 points). This means a 2-channel bilateral cluster (e.g. one left-temporal,
    %    one right-temporal channel, if they end up as separate clusters) still shows clearly via
    %    its two halos, without a hull spuriously connecting unrelated territory between them.
    % =========================================================================================
    arguments
        coherence_vals (1,:) double
        chanlocs struct
        plot_title (1,:) char
        y_label (1,:) char
        clim_min (1,1) double
        clim_max (1,1) double
        cmap (:,3) double = []
        options.ShowLabels (1,1) logical = true
        options.ShowContours (1,1) logical = true
        options.NumContours (1,1) double {mustBeInteger, mustBePositive} = 6
        options.HeadRadius (1,1) double {mustBeReal, mustBePositive} = 0.5
        options.PlotRadius (1,1) double = NaN
        options.GridResolution (1,1) double {mustBeInteger, mustBePositive} = 300
        options.InterpMethod (1,:) char = 'v4'
        options.HighlightChannels (1,:) double = []
        options.HighlightColor (1,3) double = [0 0 0]
        options.ClusterChannels cell = {}
        options.ClusterColors (:,3) double = []
        options.ClusterOutlineMargin (1,1) double {mustBeReal, mustBePositive} = 0.035
        options.Axes = []
    end

    if isempty(cmap)
        cmap = parula(256);
    end

    if numel(coherence_vals) ~= length(chanlocs)
        error('call_topoplot_custom:SizeMismatch', ...
            'Length of coherence_vals (%d) does not match number of channels (%d).', ...
            numel(coherence_vals), length(chanlocs));
    end
    coherence_vals = coherence_vals(:)';
    Num_channels = numel(chanlocs);

    % ---- Project theta/radius to 2-D Cartesian (same convention as plot_topomap.m /
    %      BuildChannelAdjacency.m / PlotChannelAdjacency.m) --------------------------------------
    x = nan(1, Num_channels);
    y = nan(1, Num_channels);
    has_coords = false(1, Num_channels);
    for ch = 1:Num_channels
        th = chanlocs(ch).theta;
        rd = chanlocs(ch).radius;
        if ~isempty(th) && ~isempty(rd) && isfinite(th) && isfinite(rd)
            theta_rad = deg2rad(th);
            x(ch) = sin(theta_rad) * rd;
            y(ch) = cos(theta_rad) * rd;
            has_coords(ch) = true;
        end
    end

    if ~any(has_coords)
        error('call_topoplot_custom:NoCoordinates', ...
            'None of the %d channels have usable theta/radius coordinates.', Num_channels);
    end

    valid_idx = find(has_coords);
    xv = x(valid_idx);
    yv = y(valid_idx);
    valv = coherence_vals(valid_idx);

    % ---- Axes to draw into -----------------------------------------------------------------
    if isempty(options.Axes)
        ax = gca;
    else
        ax = options.Axes;
    end
    axes(ax); %#ok<LAXES>
    cla(ax);
    hold(ax, 'on');

    % ---- HeadRadius (fixed cartoon) vs. PlotRadius (auto-expanding color-map/grid extent) ----
    % See design note 2 above: PlotRadius must cover every electrode actually present, not just
    % HeadRadius, or peripheral 10-10/10-5 sites (Fp1/Fp2, AF7/AF8, P9/P10, Iz, etc.) end up with
    % their labels floating past the edge of the interpolated data.
    headrad = options.HeadRadius;
    channel_radii = sqrt(xv.^2 + yv.^2);
    max_channel_radius = max(channel_radii, [], 'omitnan');
    if isnan(options.PlotRadius)
        plotrad = max(headrad, max_channel_radius * 1.02);
    else
        plotrad = options.PlotRadius;
        if plotrad < max_channel_radius
            warning('call_topoplot_custom:PlotRadiusTooSmall', ...
                ['options.PlotRadius (%.3f) is smaller than the farthest electrode''s radius ' ...
                 '(%.3f) -- that electrode''s marker/label will be drawn outside the ' ...
                 'interpolated color map. Pass a larger PlotRadius, or leave it unset (NaN) to ' ...
                 'auto-fit every electrode.'], plotrad, max_channel_radius);
        end
    end

    % ---- Interpolation grid, sized to PlotRadius and masked to a circle of that radius --------
    lin = linspace(-plotrad * 1.05, plotrad * 1.05, options.GridResolution);
    [Xi, Yi] = meshgrid(lin, lin);
    Zi = griddata(xv, yv, valv, Xi, Yi, options.InterpMethod); %#ok<GRIDD>
    mask = sqrt(Xi.^2 + Yi.^2) > plotrad;
    Zi(mask) = NaN;

    % ---- Filled color surface -----------------------------------------------------------------
    h = pcolor(ax, Xi, Yi, Zi);
    set(h, 'EdgeColor', 'none', 'FaceColor', 'interp');
    shading(ax, 'interp');

    % ---- Contour lines --------------------------------------------------------------------
    if options.ShowContours
        contour(ax, Xi, Yi, Zi, options.NumContours, 'k', 'LineWidth', 0.5);
    end

    % ---- Head, nose, ears (same drawing convention as plot_topomap.m / PlotChannelAdjacency.m) -
    % Drawn at the fixed anatomical HeadRadius, deliberately NOT PlotRadius -- see design note 2.
    theta_circle = linspace(0, 2*pi, 360);
    plot(ax, headrad * cos(theta_circle), headrad * sin(theta_circle), 'k', 'LineWidth', 2);
    nose_x = headrad * [-0.05  0.00  0.05];
    nose_y = headrad * [ 1.00  1.12  1.00];
    plot(ax, nose_x, nose_y, 'k', 'LineWidth', 2);
    ear_w = 0.03; ear_h = 0.10;
    rectangle(ax, 'Position', [-headrad - ear_w, -ear_h/2, ear_w, ear_h], 'Curvature', 0.5, ...
        'EdgeColor', 'k', 'LineWidth', 2);
    rectangle(ax, 'Position', [ headrad + 0.015,  -ear_h/2, ear_w, ear_h], 'Curvature', 0.5, ...
        'EdgeColor', 'k', 'LineWidth', 2);

    % ---- Electrode markers + labels ----------------------------------------------------------
    plot(ax, xv, yv, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
    if options.ShowLabels
        for i = 1:numel(valid_idx)
            text(ax, xv(i), yv(i) + 0.02, chanlocs(valid_idx(i)).labels, 'FontSize', 6, ...
                'HorizontalAlignment', 'center', 'Color', 'k');
        end
    end

    % ---- Highlighted individual channels (e.g. individually-significant channels) -------------
    if ~isempty(options.HighlightChannels)
        hc = options.HighlightChannels;
        hc = hc(ismember(hc, valid_idx));
        if ~isempty(hc)
            plot(ax, x(hc), y(hc), 'o', 'MarkerSize', 10, 'LineWidth', 1.8, ...
                'MarkerEdgeColor', options.HighlightColor, 'MarkerFaceColor', 'none');
        end
    end

    % ---- Cluster outlines -----------------------------------------------------------------
    if ~isempty(options.ClusterChannels)
        Num_clusters = numel(options.ClusterChannels);
        if isempty(options.ClusterColors)
            cluster_colors = lines(max(Num_clusters, 1));
        else
            cluster_colors = options.ClusterColors;
            if size(cluster_colors, 1) < Num_clusters
                error('call_topoplot_custom:NotEnoughClusterColors', ...
                    'options.ClusterColors has %d rows but there are %d clusters.', ...
                    size(cluster_colors, 1), Num_clusters);
            end
        end
        for k = 1:Num_clusters
            draw_cluster_outline(ax, x, y, options.ClusterChannels{k}, cluster_colors(k, :), ...
                options.ClusterOutlineMargin);
        end
    end

    % ---- Colorbar, colormap, clim, title -- ALL explicitly scoped to `ax` -- this explicit
    %      per-axes scoping is what makes the multi-tile colormap bug structurally impossible
    %      here (see header). --------------------------------------------------------------
    colormap(ax, cmap);
    clim(ax, [clim_min, clim_max]);
    cb = colorbar(ax);
    ylabel(cb, y_label, 'FontSize', 12);
    title(ax, plot_title, 'FontSize', 14, 'FontWeight', 'bold');

    axis(ax, 'equal');
    axis(ax, 'off');
    hold(ax, 'off');
end

% =============================================================================================
% Local helper functions
% =============================================================================================

function draw_cluster_outline(ax, x, y, channel_indices, color, margin)
    % Draws a halo circle around every member channel of a cluster (robust for any cluster size,
    % including 1), plus an inflated convex-hull polygon around the whole group when it has 3 or
    % more non-collinear members. x, y are 1 x Num_channels coordinate vectors (NaN for channels
    % without usable coordinates); channel_indices indexes into them.
    channel_indices = channel_indices(:)';
    channel_indices = channel_indices(isfinite(x(channel_indices)) & isfinite(y(channel_indices)));
    if isempty(channel_indices)
        return;
    end
    cx = x(channel_indices);
    cy = y(channel_indices);

    theta_circle = linspace(0, 2*pi, 40);
    for i = 1:numel(channel_indices)
        plot(ax, cx(i) + margin * cos(theta_circle), cy(i) + margin * sin(theta_circle), ...
            '-', 'Color', color, 'LineWidth', 2);
    end

    if numel(channel_indices) >= 3
        try
            hull_idx = convhull(cx, cy);
            centroid_x = mean(cx);
            centroid_y = mean(cy);
            hull_x = cx(hull_idx);
            hull_y = cy(hull_idx);
            vec_x = hull_x - centroid_x;
            vec_y = hull_y - centroid_y;
            vec_len = sqrt(vec_x.^2 + vec_y.^2);
            vec_len(vec_len == 0) = 1;  % guard divide-by-zero for a hull point at the centroid
            hull_x = hull_x + margin * (vec_x ./ vec_len);
            hull_y = hull_y + margin * (vec_y ./ vec_len);
            plot(ax, hull_x, hull_y, '-', 'Color', color, 'LineWidth', 2.5);
        catch
            % Degenerate (e.g. collinear) point set -- convhull can error; the per-channel halos
            % above already convey cluster membership, so just skip the extra outline.
        end
    end
end
