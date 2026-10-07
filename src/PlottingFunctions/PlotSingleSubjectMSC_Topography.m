% filepath: /Users/aakim/MATLAB_CODE/BrainSpeechCoherence/dev/PlottingFunctions/PlotSingleSubjectMSC_Topography.m
function [fig_handle, clim_used] = PlotSingleSubjectMSC_Topography(MSC_values, chanlocs, options)
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
    
    % Compute shared color limits for the two condition panels
    if any(isnan(options.Clim))
        data_min = min(MSC_values);
        data_max = max(MSC_values);
        clim_conditions = [data_min, data_max * 1.05];  % 5% margin on high end
    else
        clim_conditions = options.Clim;
    end
    
    plot_title = 'title for one subjects topomap';
    y_label = 'y label';
    clim_min = clim_conditions(1);
    clim_max = clim_conditions(2);
    cmap = parula(256);  % default colormap; can be changed to a diverging colormap if needed

    ax = call_topoplot_custom(MSC_values, chanlocs, plot_title, y_label, clim_min, clim_max, cmap);
    %% ax = call_topoplot_eeglab(MSC_values, chanlocs, plot_title, y_label, clim_min, clim_max, cmap);
    % return the figure handle for the figure
    fig_handle = gcf;
    
    clim_used = clim_conditions;
end