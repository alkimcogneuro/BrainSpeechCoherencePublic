function [fig_handle, PlotResults] = PlotGroupMSC_Topography(MSC_condition1, MSC_condition2, chanlocs, options)
    % =========================================================================================
    % Plot Group-Averaged MSC Topographies: Condition 1, Condition 2, and Their Difference
    % =========================================================================================
    % Takes per-subject band-averaged MSC values (e.g., the output of BandAverageMSC, one row
    % per subject) for two conditions, averages each condition across subjects per channel, and
    % plots three topoplots side by side: the Condition 1 group average, the Condition 2 group
    % average, and their difference (Condition 1 - Condition 2).
    %
    % This function is a thin wrapper: it computes the group averages and the difference, then
    % delegates the actual three-panel plotting to PlotThreePanelTopography.m. That function is
    % the reusable plotting primitive -- it doesn't know or care whether the data it's given is
    % a group average or a single subject. If you want to plot INDIVIDUAL subjects (one 3-panel
    % figure per subject, e.g. for batch-saving to PNGs), write a small loop that calls
    % PlotThreePanelTopography.m directly on each row of MSC_condition1/MSC_condition2, rather
    % than adding that logic here -- see PlotIndividualSubjectsMSC_Topography.m for a worked
    % example of exactly that pattern.
    %
    % Arguments:
    %   MSC_condition1: [Num_subjects x Num_channels] matrix. Each row is one subject's
    %                   BandAverageMSC output for condition 1. Subjects must be in the SAME
    %                   ORDER as MSC_condition2 (row i in both matrices = same subject), since
    %                   the difference plot is a paired, per-subject-then-averaged comparison.
    %   MSC_condition2: [Num_subjects x Num_channels] matrix, condition 2, same subject order
    %                   and same Num_channels as MSC_condition1.
    %   chanlocs:       EEGLAB channel locations structure, passed through to
    %                   PlotThreePanelTopography (and from there to call_topoplot_eeglab). Must
    %                   have Num_channels elements.
    %
    % Optional name-value arguments:
    %   options.Condition1_Label (default 'Condition 1'): title label for the first topoplot.
    %   options.Condition2_Label (default 'Condition 2'): title label for the second topoplot.
    %   options.ColorbarLabel    (default 'MSC'):          colorbar y-axis label for the two
    %                                                       condition plots (the difference plot
    %                                                       is always labeled 'MSC difference').
    %   options.FigureTitle      (default ''):              overall figure title, shown above
    %                                                       all three topoplots. Skipped if empty.
    %   options.SavePath         (default ''):              if non-empty, the figure is saved to
    %                                                       this path (extension determines
    %                                                       format, e.g. '.png', '.fig', '.pdf').
    %                                                       If empty, the figure is only displayed.
    %   options.Visible          (default true):            if false, the figure is created
    %                                                       invisibly.
    %
    % Returns:
    %   fig_handle:  handle to the created figure.
    %   PlotResults: a structure containing:
    %     .GroupMean_condition1  [1 x Num_channels] group average, condition 1
    %     .GroupMean_condition2  [1 x Num_channels] group average, condition 2
    %     .GroupMean_diff        [1 x Num_channels] condition1 - condition2, computed from the
    %                            group averages (NOT averaged per-subject differences -- see
    %                            note below)
    %     .clim_conditions  [1 x 2] shared color limits used for the two condition topoplots
    %     .clim_diff        [1 x 2] symmetric color limits used for the difference topoplot
    %
    % -----------------------------------------------------------------------------------------
    % Design notes
    % -----------------------------------------------------------------------------------------
    % 1) Color scales: see PlotThreePanelTopography.m's design notes -- the two condition
    %    topoplots share a common color scale computed from the combined data range across both
    %    conditions, and the difference topoplot uses its own zero-centered scale.
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
        options.Visible (1,1) logical = true
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

    % ---- Delegate the actual plotting to the modular three-panel plotting function ---------
    [fig_handle, clim_conditions, clim_diff] = PlotThreePanelTopography( ...
        GroupMean_condition1, GroupMean_condition2, chanlocs, ...
        'Label1', options.Condition1_Label, ...
        'Label2', options.Condition2_Label, ...
        'ColorbarLabel', options.ColorbarLabel, ...
        'FigureTitle', options.FigureTitle, ...
        'Visible', options.Visible, ...
        'SavePath', options.SavePath);

    % ---- Assemble results structure ----------------------------------------------------------
    PlotResults = struct();
    PlotResults.GroupMean_condition1 = GroupMean_condition1;
    PlotResults.GroupMean_condition2 = GroupMean_condition2;
    PlotResults.GroupMean_diff       = GroupMean_diff;
    PlotResults.clim_conditions      = clim_conditions;
    PlotResults.clim_diff            = clim_diff;
end
