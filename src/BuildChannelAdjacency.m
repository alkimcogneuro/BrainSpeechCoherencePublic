function [Adjacency, NeighborList, Coords] = BuildChannelAdjacency(chanlocs, options)
    % =========================================================================================
    % Build a Channel Adjacency ("Neighbor") Structure From EEGLAB chanlocs
    % =========================================================================================
    % Determines which channels count as spatial neighbors of which others, for use by
    % ClusterPermutationTest_ShuffleTest.m (and any other spatially-contiguous-cluster analysis).
    % This is the single-montage analog of what FieldTrip's ft_prepare_neighbours does for
    % MEG/EEG cluster-based permutation tests (Maris & Oostenveld, 2007).
    %
    % Electrode coordinates are taken from chanlocs.theta / chanlocs.radius (EEGLAB's standard
    % polar projection), converted to 2-D Cartesian with the SAME convention already used
    % elsewhere in this repo (see PlottingFunctions/plot_topomap.m):
    %     x = sin(theta_rad) .* radius;   y = cos(theta_rad) .* radius;
    % Using the same convention as plot_topomap.m means the adjacency structure lines up with
    % what you actually SEE when you look at a topoplot of these channels.
    %
    % Arguments:
    %   chanlocs: EEGLAB chanlocs struct array (same one passed to ComputeRankBasedZScore /
    %             call_topoplot_eeglab / etc. throughout this repo). Must have .theta, .radius,
    %             .labels for every channel that should participate in clustering. A channel
    %             missing .theta/.radius is treated as spatially isolated (no neighbors, so it
    %             can only ever form a cluster of size 1) rather than causing an error, since
    %             this can legitimately happen for a non-scalp reference/EOG channel included in
    %             the data.
    %
    % Optional name-value arguments:
    %   options.Method (default 'triangulation'): 'triangulation' or 'distance'.
    %       'triangulation': Delaunay-triangulate the 2-D electrode layout, then DROP edges
    %           longer than options.MaxEdgeLengthFactor times the triangulation's median edge
    %           length. Plain Delaunay triangulation alone tends to also connect far-apart
    %           electrodes around the outer boundary of the montage (an artifact of the
    %           triangulation, not a real spatial relationship) -- pruning those long edges is a
    %           standard fix. Mostly parameter-free; a reasonable default for a roughly evenly
    %           spaced montage (10-20, 10-10, etc.).
    %       'distance': connect every pair of channels within options.DistanceThreshold of each
    %           other (in the same normalized units as the projected x/y, where the head circle
    %           sits at radius 1). Simpler and fully under your control, but the right threshold
    %           depends on your montage's electrode density.
    %   options.MaxEdgeLengthFactor (default 1.5): used only for Method='triangulation'.
    %   options.DistanceThreshold (default 0.35): used only for Method='distance'. For reference,
    %       adjacent electrodes in a standard 10-10 montage are roughly 0.15-0.25 apart in these
    %       normalized units near the vertex, a bit more near the periphery -- start around there
    %       and use PlotChannelAdjacency.m to check the result looks right before trusting it.
    %
    % Returns:
    %   Adjacency:    [Num_channels x Num_channels] logical, symmetric, zero diagonal.
    %   NeighborList: 1 x Num_channels cell array; NeighborList{ch} = row vector of channel
    %                 indices adjacent to channel ch. (Same information as Adjacency, in a
    %                 format that's sometimes more convenient.)
    %   Coords:       [Num_channels x 2] the projected (x, y) coordinates actually used, in case
    %                 you want to inspect or reuse them (e.g. for your own plotting). Channels
    %                 missing theta/radius come back as NaN rows.
    %
    % ALWAYS sanity-check the result with PlotChannelAdjacency.m before trusting it for a
    % cluster test -- this function was written and syntax-checked without a MATLAB install
    % available to actually run it, so a visual check against your real montage is important.
    % =========================================================================================
    arguments
        chanlocs struct
        options.Method (1,:) char {mustBeMember(options.Method, {'triangulation', 'distance'})} = 'triangulation'
        options.MaxEdgeLengthFactor (1,1) double {mustBeReal, mustBePositive} = 1.5
        options.DistanceThreshold (1,1) double {mustBeReal, mustBePositive} = 0.35
    end

    Num_channels = numel(chanlocs);
    if Num_channels < 2
        error('BuildChannelAdjacency:TooFewChannels', ...
            'Need at least 2 channels to build an adjacency structure (got %d).', Num_channels);
    end

    % ---- Project theta/radius to 2-D Cartesian (same convention as plot_topomap.m) -----------
    has_coords = false(Num_channels, 1);
    x = nan(Num_channels, 1);
    y = nan(Num_channels, 1);
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
        error('BuildChannelAdjacency:NoCoordinates', ...
            'None of the %d channels have usable theta/radius coordinates.', Num_channels);
    end
    if any(~has_coords)
        missing_labels = {chanlocs(~has_coords).labels};
        warning('BuildChannelAdjacency:MissingCoordinates', ...
            ['%d channel(s) are missing theta/radius and will be treated as spatially isolated ' ...
             '(no neighbors): %s'], sum(~has_coords), strjoin(missing_labels, ', '));
    end

    valid_idx = find(has_coords);
    Num_valid = numel(valid_idx);

    Adjacency = false(Num_channels, Num_channels);

    if Num_valid < 2
        % Nothing to connect -- every channel isolated. Skip triangulation/distance entirely.
        NeighborList = arrayfun(@(ch) zeros(1,0), 1:Num_channels, 'UniformOutput', false);
        Coords = [x, y];
        return;
    end

    xv = x(valid_idx);
    yv = y(valid_idx);

    switch options.Method
        case 'triangulation'
            % Guard against degenerate configurations (e.g. all points collinear), which makes
            % delaunay() error or return something unusable.
            try
                tri = delaunay(xv, yv);
            catch ME
                error('BuildChannelAdjacency:TriangulationFailed', ...
                    ['Delaunay triangulation of the electrode layout failed (%s). This can ' ...
                     'happen with a degenerate (e.g. collinear or near-duplicate) set of ' ...
                     'electrode coordinates. Try Method=''distance'' instead.'], ME.message);
            end
            % Collect every triangle edge exactly once, as pairs of indices INTO valid_idx.
            edge_pairs = [tri(:,[1 2]); tri(:,[2 3]); tri(:,[1 3])];
            edge_pairs = unique(sort(edge_pairs, 2), 'rows');
            edge_lengths = sqrt((xv(edge_pairs(:,1)) - xv(edge_pairs(:,2))).^2 + ...
                                 (yv(edge_pairs(:,1)) - yv(edge_pairs(:,2))).^2);
            keep = edge_lengths <= options.MaxEdgeLengthFactor * median(edge_lengths);
            edge_pairs = edge_pairs(keep, :);

        case 'distance'
            [ii, jj] = find(triu(true(Num_valid), 1));
            d = sqrt((xv(ii) - xv(jj)).^2 + (yv(ii) - yv(jj)).^2);
            keep = d <= options.DistanceThreshold;
            edge_pairs = [ii(keep), jj(keep)];
    end

    % Map edge_pairs (indices into xv/yv, i.e. into valid_idx) back to real channel indices, and
    % fill in the symmetric adjacency matrix.
    for e = 1:size(edge_pairs, 1)
        a = valid_idx(edge_pairs(e,1));
        b = valid_idx(edge_pairs(e,2));
        Adjacency(a, b) = true;
        Adjacency(b, a) = true;
    end

    NeighborList = cell(1, Num_channels);
    for ch = 1:Num_channels
        NeighborList{ch} = find(Adjacency(ch, :));
    end

    Coords = [x, y];
end
