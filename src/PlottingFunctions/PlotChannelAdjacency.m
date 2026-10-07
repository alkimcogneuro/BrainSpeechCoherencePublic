function fig_handle = PlotChannelAdjacency(chanlocs, Adjacency, options)
    % =========================================================================================
    % Sanity-Check Plot for a Channel Adjacency Structure (from BuildChannelAdjacency.m)
    % =========================================================================================
    % Draws the electrode montage (same theta/radius -> x/y projection convention as
    % PlottingFunctions/plot_topomap.m) with a line connecting every pair of channels that
    % BuildChannelAdjacency.m considered neighbors. ALWAYS look at this before trusting a
    % cluster test built on top of that adjacency structure -- a bad adjacency structure
    % (electrodes connected across the head, or a montage fragmented into disconnected islands)
    % will silently produce a bad cluster test, and this is the fastest way to catch that by eye.
    %
    % Arguments:
    %   chanlocs:   same chanlocs struct array passed to BuildChannelAdjacency.m.
    %   Adjacency:  the [Num_channels x Num_channels] logical matrix BuildChannelAdjacency.m
    %               returned for this chanlocs.
    %
    % Optional name-value arguments:
    %   options.ShowLabels (default true): draw channel labels next to each electrode.
    %   options.Title      (default 'Channel Adjacency'): figure title.
    %
    % Returns:
    %   fig_handle: handle to the created figure. Isolated channels (no neighbors) are circled
    %   in red -- expected for a handful of peripheral/reference channels, but if most of the
    %   montage comes back isolated, or the montage splits into more than one disconnected
    %   piece, options.Method/DistanceThreshold in BuildChannelAdjacency.m need adjusting.
    % =========================================================================================
    arguments
        chanlocs struct
        Adjacency logical
        options.ShowLabels (1,1) logical = true
        options.Title (1,:) char = 'Channel Adjacency'
    end

    Num_channels = numel(chanlocs);
    if ~isequal(size(Adjacency), [Num_channels, Num_channels])
        error('PlotChannelAdjacency:SizeMismatch', ...
            'Adjacency is %dx%d but chanlocs has %d channels.', ...
            size(Adjacency,1), size(Adjacency,2), Num_channels);
    end

    x = nan(Num_channels, 1);
    y = nan(Num_channels, 1);
    for ch = 1:Num_channels
        th = chanlocs(ch).theta;
        rd = chanlocs(ch).radius;
        if ~isempty(th) && ~isempty(rd) && isfinite(th) && isfinite(rd)
            theta_rad = deg2rad(th);
            x(ch) = sin(theta_rad) * rd;
            y(ch) = cos(theta_rad) * rd;
        end
    end

    fig_handle = figure('Color', 'white', 'Position', [100 100 600 600]);
    hold on;

    % Head outline (same drawing convention as plot_topomap.m)
    headrad = 0.5;
    theta_circle = linspace(0, 2*pi, 360);
    plot(headrad * cos(theta_circle), headrad * sin(theta_circle), 'k', 'LineWidth', 2);
    nose_x = headrad * [-0.05  0.00  0.05];
    nose_y = headrad * [ 1.00  1.12  1.00];
    plot(nose_x, nose_y, 'k', 'LineWidth', 2);

    % Adjacency edges
    [ii, jj] = find(triu(Adjacency, 1));
    for e = 1:numel(ii)
        plot([x(ii(e)), x(jj(e))], [y(ii(e)), y(jj(e))], '-', 'Color', [0.3 0.6 0.9], 'LineWidth', 1);
    end

    % Electrodes
    plot(x, y, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
    if options.ShowLabels
        for ch = 1:Num_channels
            if isfinite(x(ch))
                text(x(ch), y(ch) + 0.02, chanlocs(ch).labels, 'FontSize', 7, ...
                    'HorizontalAlignment', 'center', 'Color', [0.6 0 0]);
            end
        end
    end

    isolated = ~isfinite(x) | (sum(Adjacency, 2) == 0);
    if any(isolated)
        plot(x(isolated), y(isolated), 'ro', 'MarkerSize', 10, 'LineWidth', 1.5);
    end

    axis equal off;
    title(sprintf('%s  (%d channels, %d edges, %d isolated)', options.Title, Num_channels, ...
        numel(ii), sum(isolated)), 'FontSize', 13);
    hold off;
end
