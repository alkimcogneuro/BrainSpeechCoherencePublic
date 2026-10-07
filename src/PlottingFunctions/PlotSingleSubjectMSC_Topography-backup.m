function [fig_handle, clim_used] = PlotSingleSubjectMSC_Topography(MSC_values, chanlocs, options)
    % =========================================================================================
    % Plot a Single Subject's Band-Averaged MSC as One Topographic Map
    % =========================================================================================
    % Takes ONE subject's band-averaged MSC values (one value per channel -- e.g. the output of
    % BandAverageMSC.m for a single subject, or the .True_MSC field returned by
    % ShuffleTest_SpeechEEGCoherence.m) and plots a single EEGLAB topoplot of those values.
    %
    % This is a separate, simpler primitive from PlotThreePanelTopography.m / PlotGroupMSC_
    % Topography.m / PlotIndividualSubjectsMSC_Topography.m, rather than a wrapper around them:
    % those three all exist to compare TWO conditions side by side (plus their difference), which
    % needs three panels and a diverging difference colormap. Plotting a single subject's single
    % MSC vector has no second condition and no difference to show, so it doesn't fit that
    % 3-panel shape -- forcing it through PlotThreePanelTopography with a dummy second condition
    % would draw (and throw away) two unwanted panels for no reason. This function calls the same
    % underlying call_topoplot_eeglab wrapper directly, on a single standalone axes, so it stays
    % visually consistent with the other topography plots in this codebase while only drawing
    % what's actually being shown.
    %
    % Arguments:
    %   MSC_values: a vector of band-averaged MSC values, one per channel. Accepts either
    %               orientation -- [Num_channels x 1] (the shape BandAverageMSC.m and
    %               ShuffleTest_SpeechEEGCoherence.m's .True_MSC both return) or
    %               [1 x Num_channels] -- and is reshaped internally as needed.
    %   chanlocs:   EEGLAB channel locations structure, passed through to call_topoplot_eeglab.
    %               Must have Num_channels elements, in the same channel order as MSC_values.
    %
    % Optional name-value arguments:
    %   options.Title         (default ''):    title shown above the topoplot. Skipped if empty.
    %   options.ColorbarLabel (default 'MSC'): colorbar y-axis label.
    %   options.Clim          (default []):    [1 x 2] color limits [lo hi]. If left empty, limits
    %                                          are computed automatically from MSC_values (data
    %                                          min to data max, with a small margin on the high
    %                                          end -- see compute below), matching the same
    %                                          auto-clim convention used for the condition panels
    %                                          in PlotThreePanelTopography.m. Since MSC is bounded
    %                                          in [0, 1], you may instead want a fixed scale (e.g.
    %                                          Clim=[0 1], or Clim=[0 X] to match a group plot's
    %                                          scale for direct visual comparison) -- pass it
    %                                          explicitly to override the data-driven default.
    %   options.Colormap      (default parula(256)): colormap for the topoplot. MSC is a single
    %                                          non-negative quantity (not a difference), so this
    %                                          uses the same sequential colormap as the condition
    %                                          panels elsewhere in this codebase, not the
    %                                          diverging red/blue map used only for differences.
    %   options.SavePath      (default ''):    if non-empty, the figure is saved to this path
    %                                          immediately (extension determines format, e.g.
    %                                          '.png', '.fig', '.pdf'). If empty, the figure is
    %                                          only created/displayed.
    %   options.Visible       (default true):  if false, the figure is created invisibly (e.g.
    %                                          for batch-saving many subjects without flooding
    %                                          your screen with windows).
    %
    % Returns:
    %   fig_handle: handle to the created figure.
    %   clim_used:  [1 x 2] color limits actually used (either what you passed in via
    %               options.Clim, or what was computed automatically).
    %
    % Example:
    %   CoherenceResults = AnalyzeSpeechEEGCoherence_Dataset(EEG_struct, Speech_RawData, 2, 35);
    %   MSC_band = BandAverageMSC(CoherenceResults, 3, 8);   % [Num_channels x 1], one subject
    %   PlotSingleSubjectMSC_Topography(MSC_band, EEG_struct.Chanlocs, ...
    %       Title=sprintf('%s: 3-8 Hz MSC', EEG_struct.Subj_id));
    %
    %   % Or directly from the shuffle test:
    %   ShuffleTestResults = ShuffleTest_SpeechEEGCoherence(EEG_struct, Speech_RawData, 3, 8);
    %   PlotSingleSubjectMSC_Topography(ShuffleTestResults.True_MSC, ShuffleTestResults.Chanlocs);
    % =========================================================================================
    arguments
        MSC_values (:,:) double
        chanlocs
        options.Title (1,:) char = ''
        options.ColorbarLabel (1,:) char = 'MSC'
        options.Clim (:,:) double = []
        options.Colormap (:,:) double = parula(256)
        options.SavePath (1,:) char = ''
        options.Visible (1,1) logical = true
    end

    % ---- Validate and normalize MSC_values --------------------------------------------------
    if ~isvector(MSC_values)
        error('PlotSingleSubjectMSC_Topography:NotAVector', ...
            'MSC_values must be a vector (one value per channel); got a %s array.', ...
            mat2str(size(MSC_values)));
    end
    MSC_values = MSC_values(:)';  % normalize to a row vector regardless of input orientation

    if numel(MSC_values) ~= numel(chanlocs)
        error('PlotSingleSubjectMSC_Topography:ChanlocsMismatch', ...
            ['MSC_values has %d channels but chanlocs has %d elements. These must match, in ' ...
             'the same channel order, for the topoplot to place values at the correct ' ...
             'electrode positions.'], numel(MSC_values), numel(chanlocs));
    end

    nonfinite = ~isfinite(MSC_values);
    if any(nonfinite)
        warning('PlotSingleSubjectMSC_Topography:NonFiniteInputs', ...
            ['%d of %d channels in MSC_values are non-finite (NaN/Inf). These channels may not ' ...
             'render correctly in the topoplot. Consider investigating the upstream cause ' ...
             '(see BandAverageMSC warnings).'], sum(nonfinite), numel(MSC_values));
    end

    if ~isempty(options.Clim) && numel(options.Clim) ~= 2
        error('PlotSingleSubjectMSC_Topography:BadClim', ...
            'options.Clim must be a 2-element vector [lo, hi].');
    end

    % ---- Determine color limits (auto-compute if not explicitly supplied) ------------------
    if isempty(options.Clim)
        clim_lo = min(MSC_values);
        clim_hi = max(MSC_values) * 1.01;  % small headroom, for better contrast
        if clim_hi <= clim_lo
            % Degenerate case: all channels identical (e.g., all zero). Avoid a zero-width color
            % axis, which would error in topoplot/clim.
            clim_hi = clim_lo + eps(clim_lo) * 100 + 1e-6;
        end
        clim_used = [clim_lo, clim_hi];
    else
        clim_used = options.Clim;
    end

    % ---- Build the figure ---------------------------------------------------------------------
    if options.Visible
        visible_str = 'on';
    else
        visible_str = 'off';
    end
    fig_handle = figure('Color', 'white', 'Position', [100 100 550 500], 'Visible', visible_str);
    ax = axes('Parent', fig_handle);
    axes(ax);  % make it the current axes, since topoplot()/call_topoplot_eeglab() draw via gca

    ax = call_topoplot_eeglab(MSC_values, chanlocs, options.Title, ...
        options.ColorbarLabel, clim_used(1), clim_used(2), options.Colormap);

    % Reassert colormap/clim on the returned axes, matching the defensive pattern used in
    % PlotThreePanelTopography.m -- some MATLAB/EEGLAB versions tie color-limit state to
    % colormap state internally, so a colormap(ax,...) call can occasionally reset CLim.
    colormap(ax, options.Colormap);
    clim(ax, clim_used);

    % ---- Save the figure, if requested ------------------------------------------------------
    if ~isempty(options.SavePath)
        [~, ~, ext] = fileparts(options.SavePath);
        if strcmpi(ext, '.fig')
            savefig(fig_handle, options.SavePath);
        else
            exportgraphics(fig_handle, options.SavePath);
        end
    end
end
