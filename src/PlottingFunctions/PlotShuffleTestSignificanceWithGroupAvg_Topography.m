function [fig_handle, PlotResults] = PlotShuffleTestSignificanceWithGroupAvg_Topography(ShuffleTestResultsArray, options)
    % =========================================================================================
    % Plot -log10(p) Topographies (Per Subject, Plus a Group-Average Panel)
    % =========================================================================================
    % This is a variant of PlotShuffleTestSignificance_Topography.m that appends ONE extra panel
    % to the grid: the across-subject average of -log10(p), channel by channel. It's kept as a
    % SEPARATE file (rather than added to PlotShuffleTestSignificance_Topography.m in place) at
    % the user's request, so both versions are available to compare/keep. Everything else about
    % this function -- inputs, per-subject panels, color scale approach, tiling mechanics -- is
    % identical to that file; see it for the fuller design-notes discussion (in particular why
    % -log10(p) rather than raw p, and why a shared color scale across panels makes sense here).
    %
    % Takes an array of ShuffleTestResults structs -- one per subject, as produced by
    % ShuffleTest_SpeechEEGCoherence.m -- and plots one topoplot per subject of -log10(PValues),
    % PLUS one final topoplot of the across-subject average of -log10(PValues) per channel.
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
    %                        AND the group-average panel, with 1% headroom] -- shared across every
    %                        panel, including the group-average one, so all panels are directly
    %                        comparable on one scale. Pass an explicit [lo hi] to override.)
    %   options.Colormap     (default hot(256))
    %   options.GroupLabel   (default '', auto-generated as 'GROUP AVERAGE (N=<Num_subjects>)' if
    %                        left empty)
    %   options.Visible      (default true)
    %   options.SavePath     (default ''): if non-empty, the combined figure is saved to this path
    %                        (extension determines format, e.g. '.png', '.pdf', '.fig').
    %
    % Returns:
    %   fig_handle:  handle to the created figure.
    %   PlotResults: struct with fields:
    %     .NegLog10P            [Num_subjects x Num_channels] per-subject plotted values, one row
    %                           per subject, same subject order as ShuffleTestResultsArray.
    %     .SubjectIDs           {Num_subjects x 1} cell array of Subj_id, same order.
    %     .GroupAvg_NegLog10P   [1 x Num_channels] across-subject average of NegLog10P, per
    %                           channel -- this is what's plotted in the final panel.
    %     .GroupLabel           the label used for the group-average panel.
    %     .Clim                 [1 x 2] color limits actually used (shared across all panels,
    %                           including the group-average one).
    %
    % -----------------------------------------------------------------------------------------
    % Design note: averaging -log10(p) itself, not re-deriving a group-level p-value
    % -----------------------------------------------------------------------------------------
    % The group panel is the channel-wise MEAN of the already-computed per-subject -log10(p)
    % values (mean(NegLog10P, 1)) -- the same "average the already-plotted quantity across
    % subjects" approach PlotGroupMSC_Topography.m uses for MSC. This is a quick visual summary
    % of where -log10(p) tends to run high across your subjects, NOT a new group-level
    % significance test. In particular, it is NOT the same thing as, and should not be read as,
    % a proper group-level test of whether coherence exceeds chance across subjects (e.g. a
    % one-sample t-test of each subject's True_MSC or Z-score against zero/chance across the
    % group -- that's what GroupLevel_MSC_Analysis.m is for). If you want an actual group-level
    % statistical test rather than a descriptive average of per-subject p-values, that's a
    % different computation and a good candidate for a follow-up.
    % =========================================================================================
    arguments
        ShuffleTestResultsArray struct
        options.FigureTitle (1,:) char = ''
        options.Clim (:,:) double = []
        options.Colormap (:,3) double = hot(256)
        options.GroupLabel (1,:) char = ''
        options.Visible (1,1) logical = true
        options.SavePath (1,:) char = ''
    end

    % ---- Validate inputs --------------------------------------------------------------------
    Num_subjects = numel(ShuffleTestResultsArray);
    if Num_subjects < 1
        error('PlotShuffleTestSignificanceWithGroupAvg_Topography:NoSubjects', ...
            'ShuffleTestResultsArray is empty -- nothing to plot.');
    end

    required_fields = {'PValues', 'Chanlocs', 'Subj_id'};
    for s = 1:Num_subjects
        for f = 1:numel(required_fields)
            if ~isfield(ShuffleTestResultsArray(s), required_fields{f})
                error('PlotShuffleTestSignificanceWithGroupAvg_Topography:MissingField', ...
                    ['Element %d of ShuffleTestResultsArray is missing required field "%s". ' ...
                     'Was this produced by ShuffleTest_SpeechEEGCoherence.m?'], s, required_fields{f});
            end
        end
    end

    Num_channels = numel(ShuffleTestResultsArray(1).Chanlocs);
    if numel(ShuffleTestResultsArray(1).PValues) ~= Num_channels
        error('PlotShuffleTestSignificanceWithGroupAvg_Topography:ChanlocsMismatch', ...
            'Subject 1 (%s) has %d PValues but %d Chanlocs entries. These must match.', ...
            ShuffleTestResultsArray(1).Subj_id, numel(ShuffleTestResultsArray(1).PValues), Num_channels);
    end
    for s = 2:Num_subjects
        if numel(ShuffleTestResultsArray(s).PValues) ~= Num_channels
            error('PlotShuffleTestSignificanceWithGroupAvg_Topography:InconsistentChannelCount', ...
                ['Subject %d (%s) has %d channels, but subject 1 (%s) has %d. All subjects must ' ...
                 'have the same number of channels, in the same order, for a group average and a ' ...
                 'shared color scale to make sense.'], s, ShuffleTestResultsArray(s).Subj_id, ...
                numel(ShuffleTestResultsArray(s).PValues), ShuffleTestResultsArray(1).Subj_id, Num_channels);
        end
    end

    if any(cellfun(@(p) any(p <= 0 | p > 1), {ShuffleTestResultsArray.PValues}))
        error('PlotShuffleTestSignificanceWithGroupAvg_Topography:InvalidPValues', ...
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

    % ---- Group-average panel: channel-wise mean of the per-subject -log10(p) values -----------
    % (see the Design note above for what this is and, importantly, what it is NOT)
    GroupAvg_NegLog10P = mean(NegLog10P, 1, 'omitnan');  % [1 x Num_channels]

    GroupLabel = options.GroupLabel;
    if isempty(GroupLabel)
        GroupLabel = sprintf('GROUP AVERAGE (N=%d)', Num_subjects);
    end

    % Combine per-subject panels and the group-average panel into one data/label set, so the
    % shared tiling+colormap machinery below just sees "Num_subjects+1 panels" and doesn't need
    % to know one of them is special.
    PanelData   = [NegLog10P; GroupAvg_NegLog10P];
    PanelLabels = [SubjectIDs; {GroupLabel}];

    % ---- Determine shared color limits, computed across ALL panels including the group avg ----
    if isempty(options.Clim)
        clim_max = max(PanelData(:)) * 1.01;
        if clim_max <= 0
            % Degenerate case: every p-value was 1 everywhere, for every subject and therefore
            % the group average too. Avoid a zero-width color axis.
            clim_max = 1e-6;
        end
        Clim = [0, clim_max];
    else
        if numel(options.Clim) ~= 2
            error('PlotShuffleTestSignificanceWithGroupAvg_Topography:BadClim', ...
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

    % ---- Build the figure: one tiled topoplot per subject, plus the group-average panel --------
    fig_handle = render_subject_grid_topography(PanelData, PanelLabels, ShuffleTestResultsArray(1).Chanlocs, ...
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
    PlotResults.GroupAvg_NegLog10P = GroupAvg_NegLog10P;
    PlotResults.GroupLabel = GroupLabel;
    PlotResults.Clim = Clim;
end

% =============================================================================================
% Local helper functions
% =============================================================================================
%
% Both helpers below are duplicated (not shared/imported) from PlotShuffleTestSignificance_Topography.m
% / PlotThreePanelTopography.m, following the same "small, stable, already-validated
% plotting-primitive logic can be duplicated rather than coupled across files" rationale those
% files already document (see PreprocessSpeechEEGForCoherence.m for the original statement of
% that rationale). Keep in sync with those copies if the underlying EEGLAB workarounds ever need
% to change.

function fig_handle = render_subject_grid_topography(data_matrix, panel_labels, chanlocs, ...
        colorbar_label, fig_title, clim_shared, cmap, fig_visible)
    % Build one figure with an auto-sized grid of tiles (roughly square), one topoplot per row
    % of data_matrix (here: one per subject, plus one for the group-average row appended by the
    % caller), all sharing clim_shared and cmap.
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

    % Final pass: reassert colormap/CLim for every tile now that all topoplot() calls are done
    % (see PlotThreePanelTopography.m's design notes for why this has to happen in a final pass,
    % not immediately after each individual topoplot() call).
    for s = 1:Num_panels
        colormap(ax_all(s), cmap);
        clim(ax_all(s), clim_shared);
    end

    if ~isempty(fig_title)
        title(t, fig_title, 'FontSize', 16, 'FontWeight', 'bold');
    end
end

function ax = axes_in_tile_slot(fig_handle, t, tile_index)
    % See PlotThreePanelTopography.m's copy of this same helper for the full explanation of why
    % an explicit tile_index and a standalone axes swap-in are both necessary.
    nexttile(t, tile_index);
    tile_ax = gca;
    drawnow;
    tile_position = get(tile_ax, 'Position');
    delete(tile_ax);
    ax = axes('Parent', fig_handle, 'Position', tile_position);
    axes(ax);
end
