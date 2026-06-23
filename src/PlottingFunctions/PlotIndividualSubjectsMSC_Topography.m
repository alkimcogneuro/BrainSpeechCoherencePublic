function PlotResults = PlotIndividualSubjectsMSC_Topography(MSC_condition1, MSC_condition2, chanlocs, options)
    % =========================================================================================
    % Plot Per-Subject MSC Topographies and Save Each as a PNG
    % =========================================================================================
    % Loops over each subject (each row of MSC_condition1/MSC_condition2) and produces one
    % 3-panel topography figure per subject (condition 1 / condition 2 / difference), saving
    % each as its own PNG. This is the per-subject counterpart to PlotGroupMSC_Topography.m --
    % both functions delegate the actual plotting to the same modular primitive,
    % PlotThreePanelTopography.m, so the two are guaranteed to render identically apart from
    % which data/title/color-limits get passed in.
    %
    % Arguments:
    %   MSC_condition1: [Num_subjects x Num_channels] matrix. Each row is one subject's
    %                   BandAverageMSC output for condition 1.
    %   MSC_condition2: [Num_subjects x Num_channels] matrix, condition 2, same subject order
    %                   and same Num_channels as MSC_condition1.
    %   chanlocs:       EEGLAB channel locations structure. Must have Num_channels elements.
    %
    % Optional name-value arguments:
    %   options.Condition1_Label (default 'Condition 1'): title label for the first topoplot.
    %   options.Condition2_Label (default 'Condition 2'): title label for the second topoplot.
    %   options.ColorbarLabel    (default 'MSC'):          colorbar y-axis label for the two
    %                                                       condition plots.
    %   options.OutputDir        (required):                folder to save each subject's PNG
    %                                                       into. Created automatically if it
    %                                                       doesn't already exist.
    %   options.SubjectIDs       (default {}):              cell array of char/string subject
    %                                                       identifiers, one per row of
    %                                                       MSC_condition1/2, used in each
    %                                                       figure's title and filename
    %                                                       (sanitized for filesystem safety).
    %                                                       Defaults to {'Subject01', ...} if
    %                                                       not supplied.
    %   options.FigureTitlePrefix (default ''):             if non-empty, prepended to each
    %                                                       subject's figure title as
    %                                                       '<prefix> -- <SubjectID>'.
    %   options.Visible          (default false):           if false, each figure is created
    %                                                       invisibly (faster, avoids flooding
    %                                                       your screen during batch saving).
    %   options.CloseAfterSaving (default true):            if true, each figure is closed
    %                                                       immediately after saving.
    %
    % Returns:
    %   PlotResults: a structure containing:
    %     .SavedFiles      {Num_subjects x 1} cell array of full paths to each saved PNG.
    %     .FigureHandles   [Num_subjects x 1] array of figure handles (invalid/closed if
    %                      CloseAfterSaving is true).
    %     .Clims           {Num_subjects x 1} cell array, each cell a struct with fields
    %                      .clim_conditions and .clim_diff -- the per-subject color limits
    %                      actually used for that subject's figure.
    %
    % -----------------------------------------------------------------------------------------
    % Design note: why per-subject color limits, not the group's
    % -----------------------------------------------------------------------------------------
    % Each subject's color limits are computed from that subject's OWN data (the default
    % behavior of PlotThreePanelTopography when ClimConditions/ClimDiff aren't supplied), rather
    % than reused from a group average. Single-subject MSC data is noisier and can have a very
    % different range than the group average, so forcing every subject onto a shared scale would
    % wash out or saturate many individual plots. If you DO want all subjects (and/or the group)
    % on one common scale for direct visual comparison, compute a shared [lo hi] pair yourself
    % (e.g. from PlotGroupMSC_Topography's PlotResults.clim_conditions, or from the min/max across
    % all of MSC_condition1/MSC_condition2) and pass it through to PlotThreePanelTopography via
    % its ClimConditions/ClimDiff options inside the loop below.
    % =========================================================================================
    arguments
        MSC_condition1 (:,:) double
        MSC_condition2 (:,:) double
        chanlocs
        options.Condition1_Label (1,:) char = 'Condition 1'
        options.Condition2_Label (1,:) char = 'Condition 2'
        options.ColorbarLabel (1,:) char = 'MSC'
        options.OutputDir (1,:) char
        options.SubjectIDs cell = {}
        options.FigureTitlePrefix (1,:) char = ''
        options.Visible (1,1) logical = false
        options.CloseAfterSaving (1,1) logical = true
    end

    % ---- Validate input shapes -------------------------------------------------------------
    [Num_subjects_1, Num_channels_1] = size(MSC_condition1);
    [Num_subjects_2, Num_channels_2] = size(MSC_condition2);

    if Num_subjects_1 ~= Num_subjects_2
        error('PlotIndividualSubjectsMSC_Topography:MismatchedSubjectCount', ...
            'MSC_condition1 has %d subjects (rows) but MSC_condition2 has %d. These must match.', ...
            Num_subjects_1, Num_subjects_2);
    end
    if Num_channels_1 ~= Num_channels_2
        error('PlotIndividualSubjectsMSC_Topography:MismatchedChannelCount', ...
            'MSC_condition1 has %d channels (columns) but MSC_condition2 has %d. These must match.', ...
            Num_channels_1, Num_channels_2);
    end
    if Num_channels_1 ~= numel(chanlocs)
        error('PlotIndividualSubjectsMSC_Topography:ChanlocsMismatch', ...
            ['MSC_condition1/MSC_condition2 have %d channels, but chanlocs has %d elements. ' ...
             'These must match.'], Num_channels_1, numel(chanlocs));
    end
    Num_subjects = Num_subjects_1;
    if Num_subjects < 1
        error('PlotIndividualSubjectsMSC_Topography:NoSubjects', ...
            'MSC_condition1/MSC_condition2 have zero rows (no subjects).');
    end

    if ~isfolder(options.OutputDir)
        mkdir(options.OutputDir);
    end

    if isempty(options.SubjectIDs)
        % Default subject IDs: Subject01, Subject02, ... (zero-padded so filenames sort
        % correctly even for 100+ subjects).
        pad_width = max(2, floor(log10(Num_subjects)) + 1);
        options.SubjectIDs = arrayfun(@(idx) sprintf('Subject%0*d', pad_width, idx), ...
            1:Num_subjects, 'UniformOutput', false);
    elseif numel(options.SubjectIDs) ~= Num_subjects
        error('PlotIndividualSubjectsMSC_Topography:SubjectIDMismatch', ...
            'options.SubjectIDs has %d entries, but there are %d subjects. These must match.', ...
            numel(options.SubjectIDs), Num_subjects);
    end

    % ---- Loop over subjects, plotting and saving each one -----------------------------------
    SavedFiles    = cell(Num_subjects, 1);
    FigureHandles = gobjects(Num_subjects, 1);
    Clims         = cell(Num_subjects, 1);

    for subj_idx = 1:Num_subjects
        subj_id    = options.SubjectIDs{subj_idx};
        subj_data1 = MSC_condition1(subj_idx, :);
        subj_data2 = MSC_condition2(subj_idx, :);

        subj_fig_title = subj_id;
        if ~isempty(options.FigureTitlePrefix)
            subj_fig_title = sprintf('%s -- %s', options.FigureTitlePrefix, subj_id);
        end

        % Sanitize the subject ID for use as a filename: replace anything that isn't a letter,
        % digit, hyphen, or underscore with '_' (handles spaces, slashes, etc. if SubjectIDs are
        % free-text rather than the auto-generated 'SubjectNN').
        safe_subj_id  = regexprep(subj_id, '[^a-zA-Z0-9_-]', '_');
        subj_filename = fullfile(options.OutputDir, [safe_subj_id, '_topography.png']);

        % This is the actual reuse: same plotting primitive as PlotGroupMSC_Topography.m calls,
        % just given one subject's row instead of the group mean. Per-subject color limits are
        % left to auto-compute (see the design note above for why).
        [subj_fig, subj_clim_conditions, subj_clim_diff] = PlotThreePanelTopography( ...
            subj_data1, subj_data2, chanlocs, ...
            'Label1', options.Condition1_Label, ...
            'Label2', options.Condition2_Label, ...
            'ColorbarLabel', options.ColorbarLabel, ...
            'FigureTitle', subj_fig_title, ...
            'Visible', options.Visible, ...
            'SavePath', subj_filename);

        SavedFiles{subj_idx}    = subj_filename;
        FigureHandles(subj_idx) = subj_fig;
        Clims{subj_idx}         = struct('clim_conditions', subj_clim_conditions, ...
                                          'clim_diff', subj_clim_diff);

        if options.CloseAfterSaving
            close(subj_fig);
        end
    end

    PlotResults = struct();
    PlotResults.SavedFiles    = SavedFiles;
    PlotResults.FigureHandles = FigureHandles;
    PlotResults.Clims         = Clims;
end
