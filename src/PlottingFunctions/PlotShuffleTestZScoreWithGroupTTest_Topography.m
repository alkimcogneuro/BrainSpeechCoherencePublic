function [fig_handle, PlotResults, GroupResults] = PlotShuffleTestZScoreWithGroupTTest_Topography(ShuffleTestResultsArray, options)
    % =========================================================================================
    % Plot Per-Subject Rank-Based Z-Score Topographies, Plus a Group-Level T-Statistic Panel
    % =========================================================================================
    % This supersedes the -log10(p)-based approach in PlotShuffleTestSignificance_Topography.m /
    % PlotShuffleTestSignificanceWithGroupAvg_Topography.m (kept as separate files, not deleted,
    % so you can compare). Two changes, both driven by the same underlying issue: averaging a
    % p-value-derived (or raw-mean/std-derived) quantity across subjects does not have a clean
    % statistical interpretation, and doesn't account for how CONSISTENT the effect is across
    % subjects. Instead:
    %   1) Per-subject panels plot a RANK-BASED z-score (from ComputeRankBasedZScore.m), not
    %      -log10(p) and not ShuffleTestResults.ZScores. ShuffleTestResults.ZScores standardizes
    %      True_MSC using its null distribution's raw mean/std, which is distorted by the fact
    %      that MSC is non-negative (bounded at 0) and its null distribution isn't guaranteed to
    %      be symmetric -- so a raw Z of -3 and +3 are NOT equally-extreme evidence, even though
    %      they look that way on a symmetric colormap. The rank-based transform instead uses
    %      True_MSC's PERCENTILE RANK within its own null distribution, which is valid regardless
    %      of the null's shape -- see ComputeRankBasedZScore.m and GroupLevel_ShuffleTest_Analysis.m
    %      for the full rationale. This makes the "large negative value" concern that motivated
    %      this rewrite actually meaningful: positive and negative RankZ of equal magnitude now
    %      really are comparable in magnitude, so the diverging colormap below is honest to use.
    %   2) The final "group" panel is the group-level one-sample t-statistic (per channel, across
    %      subjects, testing whether RankZ is reliably greater than 0), computed by
    %      GroupLevel_ShuffleTest_Analysis.m -- not a plain average.
    %
    % Arguments:
    %   ShuffleTestResultsArray: a struct array, one element per subject, where each element is a
    %                            ShuffleTestResults struct as returned by
    %                            ShuffleTest_SpeechEEGCoherence.m. Build this from
    %                            script_AnalyzeCoherence_withPermTest.m's All_ShuffleTestResults via:
    %                                ShuffleTestResultsArray = [All_ShuffleTestResults.Coherence_Analysis_ShuffleTestResults];
    %
    % Optional name-value arguments:
    %   options.FigureTitle (default '', auto-generated if left empty)
    %   options.Clim         (default [], auto-computed as a SYMMETRIC (zero-centered) range
    %                        covering both the per-subject RankZ values and the group t-statistic,
    %                        with 1% headroom. A single shared, zero-centered scale is used for
    %                        ALL panels -- including the t-stat panel -- so that "no effect"
    %                        always sits at the same neutral color. Pass an explicit [lo hi] (or,
    %                        since it's forced symmetric internally, just the abs bound) to override.)
    %   options.Colormap     (default: local diverging red-white-blue, blue=negative,
    %                        red=positive, matching PlotThreePanelTopography.m's diverging_redblue_cmap)
    %   options.FDR_alpha    (default 0.05): passed through to GroupLevel_ShuffleTest_Analysis.
    %   options.Tail         (default 'right'): passed through to GroupLevel_ShuffleTest_Analysis
    %                        (see that file for why 'right' is the default).
    %   options.Visible      (default true)
    %   options.SavePath     (default ''): if non-empty, the combined figure is saved to this path.
    %
    % Returns:
    %   fig_handle:   handle to the created figure.
    %   PlotResults:  struct with fields .RankZ [Num_subjects x Num_channels], .SubjectIDs,
    %                 .GroupLabel, .Clim (the shared, symmetric limits actually used).
    %   GroupResults: the full struct returned by GroupLevel_ShuffleTest_Analysis.m (tstat,
    %                 p_ttest, FDR-corrected p-values, sig_ttest_fdr, etc.) -- returned here too
    %                 so you have the underlying numbers/significance, not just the picture.
    %
    % -----------------------------------------------------------------------------------------
    % Design notes
    % -----------------------------------------------------------------------------------------
    % 1) Why the t-stat panel isn't just "another subject": it's visually distinguished with an
    %    ALL-CAPS label ("GROUP T-STAT...") the same way PlotShuffleTestSignificanceWithGroupAvg_
    %    Topography.m distinguished its group-average panel, for the same reason -- nothing else
    %    about the tiling code knows or needs to know one panel is special, but a human glancing
    %    at the figure should be able to tell immediately.
    % 2) Tiling/colormap mechanics (axes_in_tile_slot, the "reassert colormap in one final pass"
    %    workaround) are duplicated locally from PlotThreePanelTopography.m, same as the other
    %    Plot*Topography.m files in this folder -- see that file's design notes for why.
    % =========================================================================================
    arguments
        ShuffleTestResultsArray struct
        options.FigureTitle (1,:) char = ''
        options.Clim (:,:) double = []
        options.Colormap (:,3) double = []
        options.FDR_alpha (1,1) double {mustBeReal, mustBePositive} = 0.05
        options.Tail (1,:) char {mustBeMember(options.Tail, {'right', 'left', 'both'})} = 'right'
        options.Visible (1,1) logical = true
        options.SavePath (1,:) char = ''
    end

    if isempty(options.Colormap)
        cmap = diverging_redblue_cmap(256);
    else
        cmap = options.Colormap;
    end

    % ---- Run the group-level analysis (validates ShuffleTestResultsArray, computes RankZ per
    % subject via ComputeRankBasedZScore.m, and runs the group t-test) --------------------------
    GroupResults = GroupLevel_ShuffleTest_Analysis(ShuffleTestResultsArray, ...
        'FDR_alpha', options.FDR_alpha, 'Tail', options.Tail);

    RankZ      = GroupResults.RankZ;        % [Num_subjects x Num_channels]
    SubjectIDs = GroupResults.SubjectIDs;
    tstat      = GroupResults.tstat;        % [1 x Num_channels]
    Num_subjects = GroupResults.Num_subjects;

    GroupLabel = sprintf('GROUP T-STAT (one-sample, N=%d)', Num_subjects);

    % Combine per-subject RankZ panels and the group t-stat panel into one data/label set.
    PanelData   = [RankZ; tstat];
    PanelLabels = [SubjectIDs; {GroupLabel}];

    % ---- Determine shared, symmetric (zero-centered) color limits ------------------------------
    if isempty(options.Clim)
        max_abs = max(abs(PanelData(:)), [], 'omitnan') * 1.01;
        if ~isfinite(max_abs) || max_abs <= 0
            max_abs = 1e-6;  % degenerate case: everything is exactly 0
        end
        Clim = [-max_abs, max_abs];
    else
        if numel(options.Clim) == 2
            Clim = [min(options.Clim), max(options.Clim)];
        elseif numel(options.Clim) == 1
            Clim = [-abs(options.Clim), abs(options.Clim)];
        else
            error('PlotShuffleTestZScoreWithGroupTTest_Topography:BadClim', ...
                'options.Clim must be a scalar (symmetric bound) or a 2-element vector [lo, hi].');
        end
    end

    % ---- Figure title -----------------------------------------------------------------------
    fig_title = options.FigureTitle;
    if isempty(fig_title)
        if isfield(GroupResults, 'band_low') && isfield(GroupResults, 'band_high')
            fig_title = sprintf('Speech-EEG Coherence Rank Z + Group T-Test, [%.1f %.1f] Hz (%s-tailed, N=%d)', ...
                GroupResults.band_low, GroupResults.band_high, options.Tail, Num_subjects);
        else
            fig_title = sprintf('Speech-EEG Coherence Rank Z + Group T-Test (%s-tailed, N=%d)', ...
                options.Tail, Num_subjects);
        end
    end

    % ---- Build the figure ---------------------------------------------------------------------
    fig_handle = render_subject_grid_topography(PanelData, PanelLabels, GroupResults.Chanlocs, ...
        'Rank Z / t', fig_title, Clim, cmap, options.Visible);

    if ~isempty(options.SavePath)
        [~, ~, ext] = fileparts(options.SavePath);
        if strcmpi(ext, '.fig')
            savefig(fig_handle, options.SavePath);
        else
            exportgraphics(fig_handle, options.SavePath);
        end
    end

    PlotResults = struct();
    PlotResults.RankZ = RankZ;
    PlotResults.SubjectIDs = SubjectIDs;
    PlotResults.GroupLabel = GroupLabel;
    PlotResults.Clim = Clim;
end

% =============================================================================================
% Local helper functions
% =============================================================================================
%
% All three below are duplicated (not shared/imported) from PlotThreePanelTopography.m and/or
% PlotShuffleTestSignificanceWithGroupAvg_Topography.m -- see those files, and
% PreprocessSpeechEEGForCoherence.m's header, for the rationale on why this codebase duplicates
% small, stable, already-validated plotting-primitive logic rather than sharing it across files.

function fig_handle = render_subject_grid_topography(data_matrix, panel_labels, chanlocs, ...
        colorbar_label, fig_title, clim_shared, cmap, fig_visible)
    Num_panels = size(data_matrix, 1);
    Num_cols = ceil(sqrt(Num_panels));
    Num_rows = ceil(Num_panels / Num_cols);

    if fig_visible
        visible_str = 'on';
    else
        visible_str = 'off';
    end
    fig_handle = figure('Color', 'white', 'Position', [100 100 350*Num_cols 350*Num_rows], ...
        'Visible', visible_str);
    t = tiledlayout(fig_handle, Num_rows, Num_cols, 'TileSpacing', 'compact', 'Padding', 'compact');

    ax_all = gobjects(Num_panels, 1);
    for s = 1:Num_panels
        ax_all(s) = axes_in_tile_slot(fig_handle, t, s); %#ok<NASGU> % overwritten below
        ax_all(s) = call_topoplot_eeglab(data_matrix(s, :), chanlocs, panel_labels{s}, ...
            colorbar_label, clim_shared(1), clim_shared(2), cmap);
    end

    for s = 1:Num_panels
        colormap(ax_all(s), cmap);
        clim(ax_all(s), clim_shared);
    end

    if ~isempty(fig_title)
        title(t, fig_title, 'FontSize', 16, 'FontWeight', 'bold');
    end
end

function ax = axes_in_tile_slot(fig_handle, t, tile_index)
    nexttile(t, tile_index);
    tile_ax = gca;
    drawnow;
    tile_position = get(tile_ax, 'Position');
    delete(tile_ax);
    ax = axes('Parent', fig_handle, 'Position', tile_position);
    axes(ax);
end

function cmap = diverging_redblue_cmap(n)
    % Duplicated from PlotThreePanelTopography.m's identical local helper. Blue (negative) ->
    % white (zero) -> red (positive).
    if mod(n, 2) == 0
        half = n/2;
        lower = [linspace(0.0, 1.0, half)', linspace(0.0, 1.0, half)', linspace(0.7, 1.0, half)'];
        upper = [linspace(1.0, 0.8, half)', linspace(1.0, 0.0, half)', linspace(1.0, 0.0, half)'];
        cmap = [lower; upper];
    else
        half = (n-1)/2;
        lower = [linspace(0.0, 1.0, half)', linspace(0.0, 1.0, half)', linspace(0.7, 1.0, half)'];
        upper = [linspace(1.0, 0.8, half)', linspace(1.0, 0.0, half)', linspace(1.0, 0.0, half)'];
        cmap = [lower; [1 1 1]; upper];
    end
end
