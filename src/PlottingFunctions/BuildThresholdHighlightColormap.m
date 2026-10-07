function cmap = BuildThresholdHighlightColormap(Clim, Threshold, options)
    % =========================================================================================
    % Build an [N x 3] Colormap That Stays Muted Below a Threshold and Ramps Through a "Hot"
    % Multi-Color Scale Above It -- for Spotting "Where Are the Large Values" (and How Large) at
    % a Glance
    % =========================================================================================
    % A plain sequential colormap (e.g. parula(256), the default everywhere in this repo) spreads
    % its color change evenly across the WHOLE data range. For something like ZScores, where most
    % channels sit near 0 (the "nothing interesting here" bulk of the distribution) and you only
    % care about the handful of channels standing out above some threshold, that even spread means
    % the interesting values barely change color from the uninteresting ones -- which is exactly
    % the "lots of blue, can't tell where Z is above threshold" symptom that motivated this file.
    %
    % This function builds a colormap with the color change concentrated where it's useful: values
    % below Threshold stay in one muted color family (a small dark-navy-to-slate-blue gradient, so
    % you're not throwing away ALL of the below-threshold structure, just de-emphasizing it), and
    % values from Threshold up to Clim(2) ramp through a "hot" thermal scale -- deep red, through
    % orange, through yellow, up to near-white at the very top -- by default. Crossing Threshold is
    % where the visible "pop" happens; the multi-stop hot ramp above it then gives you visible
    % separation between "just interesting" and "extremely large," rather than every above-
    % threshold value reading as roughly the same shade (2026-09-18 update -- see open-items.md
    % item 15's follow-up: with only two above-threshold colors, everything interesting ended up
    % some shade of yellow, with no way to tell a large value from a very large one).
    %
    % This is a THRESHOLD/highlight colormap, not a diverging one -- it does not assume Threshold is
    % meaningfully the midpoint of anything, and values below Threshold are not treated as mirror
    % images of values above it (unlike, say, a symmetric red-blue diverging map centered at 0).
    % That matches how you'd read ZScores/RankZ for this purpose: you're hunting for "unusually
    % high," not for "unusually high vs. unusually low, expressed the same way."
    %
    % Use it together with an EXPLICIT options.Clim when calling PlotSubjectsAndGroupTopography.m
    % (or call_topoplot_custom.m directly) -- the colormap's row-to-color mapping only makes sense
    % relative to the SAME [lo hi] range it was built with, since a colormap has no idea what data
    % range it will be stretched across at draw time. If you let Clim auto-compute instead, the
    % threshold's visual position will silently drift depending on that call's actual data range.
    %
    % Arguments:
    %   Clim: [1 x 2] the [lo, hi] color axis limits you intend to plot with (e.g. what you'll pass
    %       as options.Clim to PlotSubjectsAndGroupTopography.m). Must have Clim(2) > Clim(1).
    %   Threshold: scalar value, in the same units as your data, where the color should start
    %       visibly changing toward the highlight colors (e.g. a one-tailed significance cutoff
    %       for a Z-score-style "starting to look interesting" boundary). Should fall strictly
    %       between Clim(1) and Clim(2) -- if it doesn't, a warning is issued (see below) because
    %       the highlight effect won't be visible.
    %
    % Optional name-value arguments:
    %   options.N               (default 256): number of rows in the returned colormap.
    %   options.BelowColor      (default [0.16 0.20 0.45], a muted slate blue): the color at
    %                           Threshold itself, for the "below" side -- values well below
    %                           Threshold fade toward options.BelowColorDim instead (see below),
    %                           so the below-threshold region isn't perfectly flat/uninformative.
    %   options.BelowColorDim   (default [0.05 0.07 0.20], a near-black navy): the color at Clim(1)
    %                           -- the darkest, most muted end of the below-threshold gradient.
    %   options.AboveColors     (default a 5-stop "hot" ramp -- deep red, red-orange, orange,
    %                           yellow, near-white; see the matrix below): an [M x 3] RGB matrix
    %                           of color waypoints for the ABOVE-threshold region, in order from
    %                           "just above Threshold" (row 1) to "at Clim(2)" (row M). Values in
    %                           between are interpolated piecewise-linearly across the M-1
    %                           segments, so more rows = more visually distinct "how large is
    %                           large" bands. Pass your own [M x 3] matrix (M >= 2) to use a
    %                           different above-threshold palette while keeping the same
    %                           threshold-highlight structure -- e.g. a blue-green-white "cool"
    %                           ramp, if red/orange/yellow isn't the right metaphor for your data.
    %
    % Returns:
    %   cmap: [options.N x 3] RGB colormap, suitable for options.Colormap on
    %       PlotSubjectsAndGroupTopography.m, or the cmap argument of call_topoplot_custom.m /
    %       call_topoplot_eeglab.m directly.
    %
    % Example (ZScores, highlighting above a one-tailed 10% significance cutoff, z ~= 1.2816):
    %   my_clim = [-3, 6];   % pick something that comfortably covers your actual ZScores range
    %   cmap = BuildThresholdHighlightColormap(my_clim, 1.2816);
    %   [fig_z, PlotResults_Z] = PlotSubjectsAndGroupTopography(CoherenceResultsArray, 'ZScores', ...
    %       OutputDir=topography_output_dir, SubjectIDField='Subj_id', ColorbarLabel='ZScore', ...
    %       Clim=my_clim, Colormap=cmap);
    %
    % This is generic -- not specific to ZScores/MSC/this repo's field names -- so the same idea
    % (and this same function) works for RankZ with a different Threshold, or for -log10(p) with a
    % threshold at -log10(0.05), etc.
    % =========================================================================================
    arguments
        Clim (1,2) double
        Threshold (1,1) double
        options.N (1,1) double {mustBeInteger, mustBePositive} = 256
        options.BelowColor (1,3) double = [0.16 0.20 0.45]
        options.BelowColorDim (1,3) double = [0.05 0.07 0.20]
        options.AboveColors (:,3) double = [ ...
            0.55 0.00 0.00; ...   % deep red    -- just above Threshold
            0.85 0.25 0.05; ...   % red-orange
            0.97 0.55 0.10; ...   % orange
            1.00 0.85 0.25; ...   % yellow
            1.00 1.00 0.92   ...  % near-white  -- at Clim(2), the most extreme values
            ]
    end

    if Clim(2) <= Clim(1)
        error('BuildThresholdHighlightColormap:BadClim', ...
            'Clim must satisfy Clim(2) > Clim(1) (got [%g, %g]).', Clim(1), Clim(2));
    end
    if size(options.AboveColors, 1) < 2
        error('BuildThresholdHighlightColormap:TooFewAboveColors', ...
            ['options.AboveColors must have at least 2 rows (got %d) -- one row for the color ' ...
             'just above Threshold, and one row for the color at Clim(2), at minimum.'], ...
            size(options.AboveColors, 1));
    end

    N = options.N;
    frac = (Threshold - Clim(1)) / (Clim(2) - Clim(1));  % where Threshold falls within [0, 1]
    frac_clamped = min(max(frac, 0), 1);

    if frac <= 0
        warning('BuildThresholdHighlightColormap:ThresholdAtOrBelowClimLo', ...
            ['Threshold (%g) is at or below Clim(1) (%g) -- the entire colormap will be the ' ...
             'above-threshold "hot" ramp, with no muted "below" region. Lower Threshold or ' ...
             'raise Clim(1) if that''s not what you want.'], Threshold, Clim(1));
    elseif frac >= 1
        warning('BuildThresholdHighlightColormap:ThresholdAtOrAboveClimHi', ...
            ['Threshold (%g) is at or above Clim(2) (%g) -- the entire colormap will be the ' ...
             'muted "below" color, so nothing will visibly pop. Raise Clim(2) above Threshold, ' ...
             'or lower Threshold, for the highlight effect to actually show up.'], Threshold, Clim(2));
    end

    N_below = round(frac_clamped * N);
    N_above = N - N_below;

    below_cmap = interp_multistop_colors([options.BelowColorDim; options.BelowColor], N_below);
    above_cmap = interp_multistop_colors(options.AboveColors, N_above);

    cmap = [below_cmap; above_cmap];
end

function out = interp_multistop_colors(colors, N)
    % Piecewise-linear interpolation of an [M x 3] RGB waypoint matrix into an [N x 3] gradient.
    % M == 1 (a single color) returns that color repeated N times. Evenly spaces the M waypoints
    % across the output regardless of N, so this works the same way whether N is small or large.
    if N <= 0
        out = zeros(0, 3);
        return;
    end
    M = size(colors, 1);
    if M == 1
        out = repmat(colors, N, 1);
        return;
    end
    query_points = linspace(0, 1, N)';
    waypoint_positions = linspace(0, 1, M);
    out = zeros(N, 3);
    for channel = 1:3
        out(:, channel) = interp1(waypoint_positions, colors(:, channel), query_points);
    end
end
