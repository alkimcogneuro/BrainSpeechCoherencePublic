% filepath: /Users/aakim/MATLAB_CODE/BrainSpeechCoherence/dev/PlottingFunctions/PlotSingleSubjectMSC_Topography.m
function [fig_handle, clim_used] = PlotSingleSubjectMSC_Topography(MSC_values, MSC_values_control, chanlocs, options)
    % =========================================================================================
    % Plot a Single Subject's Band-Averaged MSC vs Control as Three Topographic Maps
    % =========================================================================================
    % Takes ONE subject's band-averaged MSC values and corresponding control MSC values
    % (with randomly shuffled EEG-speech pairings) and plots three EEGLAB topoplots:
    % 1. MSC values (true pairing)
    % 2. MSC control values (shuffled pairing)
    % 3. Difference (MSC - Control)
    %
    % Arguments:
    %   MSC_values:         vector of band-averaged MSC values, one per channel
    %   MSC_values_control: vector of control MSC values (shuffled), same shape as MSC_values
    %   chanlocs:           EEGLAB channel locations structure
    %   options:            name-value options structure
    %
    % Optional name-value arguments:
    %   options.Title           (default ''): overall figure title
    %   options.Condition1_Label (default 'True MSC'): label for MSC panel
    %   options.Condition2_Label (default 'Control MSC'): label for control panel
    %   options.ColorbarLabel   (default 'MSC'): colorbar label for condition panels
    %   options.Clim            (default []): [1 x 2] color limits for condition panels
    %   options.SavePath        (default ''): path to save figure
    %   options.Visible         (default true): figure visibility
    
    arguments
        MSC_values
        MSC_values_control
        chanlocs
        options.Title {mustBeTextScalar} = ''
        options.Condition1_Label {mustBeTextScalar} = 'True MSC'
        options.Condition2_Label {mustBeTextScalar} = 'Control (shuffled pairings) MSC'
        options.ColorbarLabel {mustBeTextScalar} = 'MSC'
        %%options.Clim (1,2) {mustBeNumeric} = []
        options.Clim (1,:) {mustBeNumeric} = [nan nan]
        options.SavePath {mustBeTextScalar} = ''
        options.Visible (1,1) logical = true
    end
    
    % Reshape to row vectors if needed
    MSC_values = reshape(MSC_values, 1, []);
    MSC_values_control = reshape(MSC_values_control, 1, []);
    
    % Compute shared color limits for the two condition panels
    if any(isnan(options.Clim))
        data_min = min([MSC_values, MSC_values_control]);
        data_max = max([MSC_values, MSC_values_control]);
        clim_conditions = [data_min, data_max * 1.05];  % 5% margin on high end
    else
        clim_conditions = options.Clim;
    end
    
    % Compute difference
    MSC_diff = MSC_values - MSC_values_control;
    
    % Call PlotThreePanelTopography to create the figure
    [fig_handle, PlotResults] = PlotThreePanelTopography(...
        MSC_values, MSC_values_control, chanlocs, ...
        Label1=options.Condition1_Label, ...
        Label2=options.Condition2_Label, ...
        ColorbarLabel=options.ColorbarLabel, ...
        FigureTitle=options.Title, ...
        SavePath=options.SavePath, ...
        Visible=options.Visible, ...
        ClimConditions=clim_conditions);
    
    clim_used = clim_conditions;
end