function [fig_handle, PlotResults] = PlotShuffleTestSignificance_Topography(ShuffleTestResultsArray, options)
    % =========================================================================================
    % Plot -log10(p) Topographies, One Panel Per Subject, From ShuffleTest_SpeechEEGCoherence Output
    % =========================================================================================
    % Takes an array of ShuffleTestResults structs -- one per subject, as produced by
    % ShuffleTest_SpeechEEGCoherence.m -- and plots one topoplot per subject of
    % -log10(PValues), a "hot = more significant" significance map: low p-values (strong
    % evidence that band-averaged MSC exceeds the shuffle-derived null) map to hot colors
    % (bright/yellow), high p-values (no evidence of above-chance coherence) map to cool/dark
    % colors.
    %
    % Why -log10(p) rather than p itself: p ranges from ~1/(NumShuffles+1) (most significant) up
    % to 1 (least significant), so smaller p = more significant, which is the OPPOSITE of what a
    % "hot = significant" colormap wants without a sign flip. -log10(p) inverts that: it's largest
    % when p is smallest, so mapping -log10(p) onto a sequential "hot" colormap (black/red at the
    % low end, yellow/white at the high end) puts the most significant channels at the hottest end
    % of the scale, with no need for a diverging/zero-centered colormap (this is a one-sided test,
    % see ShuffleTest_SpeechEEGCoherence.m's docstring -- there's no "significant in the other
    % direction" to represent).
    %
    % Arguments:
    %   ShuffleTestResultsArray: a struct array, one element per subject, where each element is a
    %                            ShuffleTestResults struct as returned by
    %                            ShuffleTest_SpeechEEGCoherence.m (must contain at minimum the
    %                            fields .PValues, .Chanlocs, .Subj_id, .band_low, .band_high,
    %                            .NumShuffles). Build this from
    %                            script_AnalyzeCoherence_withPermTest.m's All_ShuffleTestResults via:
    %                                ShuffleTestResultsArray = [All_ShuffleTestResults.Coherence_Analysis_ShuffleTestResults];
    %
    % Optional name-value arguments:
    %   options.FigureTitle (default '', auto-generated from the first subject's band_low/band_high
    %                        and NumShuffles if left empty)
    %   options.Clim         (default [], auto-computed as [0, max(-log10(p)) across ALL subjects
    %                        and channels, with 1% headroom] -- shared across every panel so
    %                        subjects are directly comparable to each other on one scale. Pass an
    %                        explicit [lo hi] to override.)
    %   options.Colormap     (default hot(256))
    %   options.Visible      (default true)
    %   options.SavePath     (default ''): if non-empty, the combined figure is saved to this path
    %                        (extension determines format, e.g. '.png', '.pdf', '.fig').
    %
    % Returns:
    %   fig_handle:  handle to the created figure.
    %   PlotResults: struct with fields:
    %     .NegLog10P    [Num_subjects x Num_channels] the plotted values, one row per subject, in
    %                   the same subject order as ShuffleTestResultsArray.
    %     .SubjectIDs   {Num_subjects x 1} cell array of Subj_id, same order.
    %     .Clim         [1 x 2] color limits actually used (shared across all panels).
    %
    % -----------------------------------------------------------------------------------------
    % Design notes
    % -----------------------------------------------------------------------------------------
    % 1) Tiling technique: EEGLAB's topoplot() resizes its own axes internally, which conflicts
    %    with axes managed by a TiledChartLayout tile. This function reuses the same
    %    axes_in_tile_slot() workaround, and the same "reassert colormap/clim in one final pass"
    %    fix for EEGLAB's cross-tile colormap clobbering, as PlotThreePanelTopography.m -- see
    %    that file's design notes for the full explanation of why. The two copies are intentionally
    %    kept in sync rather than sharing code, following the same rationale
    %    PreprocessSpeechEEGForCoherence.m gives for its own intentional duplication: these are
    %    small, stable, already-validated pieces of plotting-primitive logic, not pairing-dependent
    %    analysis code, so duplicating them here avoids coupling this file's release cadence to
    %    PlotThreePanelTopography.m's.
    % 2) One shared color scale, not one per subject: unlike
    %    PlotIndividualSubjectsMSC_Topography.m (which deliberately uses each subject's OWN color
    %    range, since raw MSC magnitude varies a lot subject to subject), here every subject's
    %    -log10(p) lives on the SAME meaningful scale (a p-value is a p-value, regardless of
    %    subject), so a shared Clim across all panels is the right default -- it lets you compare
    %    significance at a glance across subjects, which is presumably the point of plotting them
    %    together.
    % 3) No significance markers (yet): ShuffleTest_SpeechEEGCoherence.m's own docstring flags that
    %    testing every channel independently needs a multiple-comparisons correction (e.g.
    %    Benjamini-Hochberg FDR) before declaring anything "significant". This function only plots
    %    the continuous -log10(p) map; it does not threshold or mark specific channels as
    %    significant. That would be a natural follow-on (overlay markers at FDR-significant
    %    channels), but doing it correctly requires matching this file's marker coordinates to
    %    topoplot()'s internal electrode-placement geometry exactly, which needs to be checked
    %    visually against a real figure -- left for a follow-up once you can see the output.
    % =========================================================================================
    arguments
        ShuffleTestResultsArray struct
        options.FigureTitle (1,:) char = ''
        options.Clim (:,:) double = []
        options.Colormap (:,3) double = hot(256)
        options.Visible (1,1) logical = true
        options.SavePath (1,:) char = ''
    end

    % ---- Validate inputs --------------------------------------------------------------------
    Num_subjects = numel(ShuffleTestResultsArray);
    if Num_subjects < 1
        error('PlotShuffleTestSignificance_Topography:NoSubjects', ...
            'ShuffleTestResultsArray is empty -- nothing to plot.');
    end

    required_fields = {'PValues', 'Chanlocs', 'Subj_id'};
    for s = 1:Num_subjects
        for f = 1:numel(required_fields)
            if ~isfield(ShuffleTestResultsArray(s), required_fields{f})
                error('PlotShuffleTestSignificance_Topography:MissingField', ...
                    ['Element %d of ShuffleTestResultsArray is missing required field "%s". ' ...
                     'Was this produced by ShuffleTest_SpeechEEGCoherence.m?'], s, required_fields{f});
            end
        end
    end

    Num_channels = numel(ShuffleTestResultsArray(1).Chanlocs);
    if numel(ShuffleTestResultsArray(1).PValues) ~= Num_channels
        error('PlotShuffleTestSignificance_Topography:ChanlocsMismatch', ...
            'Subject 1 (%s) has %d PValues but %d Chanlocs entries. These must match.', ...
            ShuffleTestResultsArray(1).Subj_id, numel(ShuffleTestResultsArray(1).PValues), Num_channels);
    end
    for s = 2:Num_subjects
        if numel(ShuffleTestResultsArray(s).PValues) ~= Num_channels
            error('PlotShuffleTestSignificance_Topography:InconsistentChannelCount', ...
                ['Subject %d (%s) has %d channels, but subject 1 (%s) has %d. All subjects must ' ...
                 'have the same number of channels, in the same order, for a shared color scale ' ...
                 'and consistent panel layout to make sense.'], s, ShuffleTestResultsArray(s).Subj_id, ...
                numel(ShuffleTestResultsArray(s).PValues), ShuffleTestResultsArray(1).Subj_id, Num_channels);
        end
    end

    if any(cellfun(@(p) any(p <= 0 | p > 1), {ShuffleTestResultsArray.PValues}))
        error('PlotShuffleTestSignificance_Topography:InvalidPValues', ...
            ['Found PValues outside (0, 1]. ShuffleTest_SpeechEEGCoherence.m''s p-value estimator ' ...
             'should never return exactly 0 or a value above 1 -- check that ShuffleTestResultsArray ' ...
             'actually came from that function.']);
    end

    % ---- Compute -log10(p) for every subject --------------------------------------------------
    NegLog10P = nan(Num_subjects, Num_channels);
    SubjectIDs = cell(Num_subjects, 1);
    for s = 1:Num_subjects
        NegLog10P(s, :) = -log10(reshape(ShuffleTestResultsArray(s).PValues, 1, []));
        SubjectIDs{s} = ShuffleTestResultsArray(s).Subj_id;
    end

    % ---- Determine shared color limits (see Design note 2) ------------------------------------
    if isempty(options.Clim)
        clim_max = max(NegLog10P(:)) * 1.01;
        if clim_max <= 0
            % Degenerate case: every p-value was 1 everywhere (no evidence of coherence anywhere,
            % for any subject). Avoid a zero-width color axis.
            clim_max = 1e-6;
        end
        Clim = [0, clim_max];
    else
        if numel(options.Clim) ~= 2
            error('PlotShuffleTestSignificance_Topography:BadClim', ...
                'options.Clim must be a 2-element vector [lo, hi].');
        end
        Clim = options.Clim;
    end

    % ---- Figure title (see Arguments doc above) ------------------------------------------------
    fig_title = options.FigureTitle;
    if isempty(fig_title) && isfield(ShuffleTestResultsArray(1), 'band_low') && isfield(ShuffleTestResultsArray(1), 'band_high')
        fig_title = sprintf('Speech-EEG Coherence Significance: -log10(p), [%.1f %.1f] Hz (%d shuffles)', ...
            ShuffleTestResultsArray(1).band_low, ShuffleTestResultsArray(1).band_high, ...
            ShuffleTestResultsArray(1).NumShuffles);
    end

    % ---- Build the figure: one tiled topoplot per subject --------------------------------------
    fig_handle = render_subject_grid_topography(NegLog10P, SubjectIDs, ShuffleTestResultsArray(1).Chanlocs, ...
        '-log10(p)', fig_title, Clim, options.Colormap, options.Visible);

    if ~isempty(options.SavePath)
        [~, ~, ext] = fileparts(options.SavePath);
        if strcmpi(ext, '.fig')
            savefig(fig_handle, options.SavePath);
        else
            exportgraphics(fig_handle, options.SavePath);
        end
    end

    PlotResults = struct();
    PlotResults.NegLog10P = NegLog10P;
    PlotResults.SubjectIDs = SubjectIDs;
    PlotResults.Clim = Clim;
end

% =============================================================================================
% Local helper functions
% =============================================================================================

function fig_handle = render_subject_grid_topography(data_matrix, subject_ids, chanlocs, ...
        colorbar_label, fig_title, clim_shared, cmap, fig_visible)
    % Build one figure with an auto-sized grid of tiles (roughly square), one topoplot per row
    % of data_matrix, all sharing clim_shared and cmap. Mirrors PlotThreePanelTopography.m's
    % render_three_panel_topography_figure, generalized from a fixed 3 panels to N panels -- see
    % that function for the reasoning behind the tiling/colormap workarounds reused here.
    Num_subjects = size(data_matrix, 1);
    Num_cols = ceil(sqrt(Num_subjects));
    Num_rows = ceil(Num_subjects / Num_cols);

    if fig_visible
        visible_str = 'on';
    else
        visible_str = 'off';
    end
    fig_handle = figure('Color', 'white', 'Position', [100 100 350*Num_cols 350*Num_rows], ...
        'Visible', visible_str);
    t = tiledlayout(fig_handle, Num_rows, Num_cols, 'TileSpacing', 'compact', 'Padding', 'compact');

    ax_all = gobjects(Num_subjects, 1);
    for s = 1:Num_subjects
        ax_all(s) = axes_in_tile_slot(fig_handle, t, s); %#ok<NASGU> % overwritten below
        ax_all(s) = call_topoplot_eeglab(data_matrix(s, :), chanlocs, subject_ids{s}, ...
            colorbar_label, clim_shared(1), clim_shared(2), cmap);
    end

    % Final pass: reassert colormap/CLim for every tile now that all topoplot() calls are done
    % (see PlotThreePanelTopography.m's design notes for why this has to happen in a final pass,
    % not immediately after each individual topoplot() call).
    for s = 1:Num_subjects
        colormap(ax_all(s), cmap);
        clim(ax_all(s), clim_shared);
    end

    if ~isempty(fig_title)
        title(t, fig_title, 'FontSize', 16, 'FontWeight', 'bold');
    end
end

function ax = axes_in_tile_slot(fig_handle, t, tile_index)
    % Duplicated from PlotThreePanelTopography.m -- see that file's copy for the full explanation
    % of why an explicit tile_index and a standalone axes swap-in are both necessary (EEGLAB's
    % topoplot() resizes its own axes, which a TiledChartLayout-managed axes won't allow; and
    % nexttile(t)'s auto-advance gets fooled once we start deleting tile axes). Keep this in sync
    % with that copy if the underlying workaround ever needs to change.
    nexttile(t, tile_index);
    tile_ax = gca;
    drawnow;
    tile_position = get(tile_ax, 'Position');
    delete(tile_ax);
    ax = axes('Parent', fig_handle, 'Position', tile_position);
    axes(ax);
end
