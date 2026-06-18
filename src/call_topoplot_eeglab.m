function [] = call_topoplot_eeglab(coherence_vals, chanlocs, plot_title, y_label)
    % topoplot_coherence.m
    % Plot magnitude squared coherence as a scalp topography using EEGLAB's topoplot.
    % Assumes:
    %   - EEG struct is already loaded in the workspace (with chanlocs)
    %   - coherence_vals is a 1 x 64 (or 64 x 1) vector of pre-computed MSC values
    
    % -------------------------------------------------------------------------
    % INPUTS — edit these
    % -------------------------------------------------------------------------
    %% plot_title     = 'Magnitude Squared Coherence';
    
    % Optional: set color axis limits (MSC ranges 0–1, but narrow if needed)
    clim_min = 0;
    %%clim_max = 1;
    %% note:  if we want to compare two topoplots, we should use the same color axis limits for both, to make them comparable.
    %%clim_max  = max(coherence_vals) * 1.1;  % set max to 10% above the max value in the data, for better contrast
    clim_max = 0.04;  % set max to 0.5, which is a reasonable upper bound for MSC values in EEG data, and allows for better contrast in the plot.
    % -------------------------------------------------------------------------
    % Validate inputs
    % -------------------------------------------------------------------------
    % if ~exist('EEG', 'var') || ~isfield(EEG, 'chanlocs')
    %     error('EEG struct with chanlocs not found in workspace.');
    % end
    
    if numel(coherence_vals) ~= length(chanlocs)
        error('Length of coherence_vals (%d) does not match number of channels (%d).', ...
        numel(coherence_vals), length(chanlocs));
    end
    
    % Ensure it's a row vector
    coherence_vals = coherence_vals(:)';
    
    % -------------------------------------------------------------------------
    % Plot
    % -------------------------------------------------------------------------
    figure('Color', 'white', 'Position', [100 100 600 500]);
    
    topoplot(coherence_vals, chanlocs, ...
    'maplimits',  [clim_min, clim_max], ...  % color axis limits
    'electrodes', 'labelpoint', ...           % show channel labels + dots
    'colormap',   parula(256), ...            % colormap (parula is more perceptually uniform than jet)
    'whitebk',    'on');                      % white background
    %%    'colormap',   jet(256), ...            % colormap (jet is traditional for EEG)
    % Colorbar and title
    cb = colorbar;
    ylabel(cb, y_label, 'FontSize', 12);
    clim([clim_min clim_max]);   % MATLAB R2022a+; use caxis() if on older MATLAB
    title(plot_title, 'FontSize', 14, 'FontWeight', 'bold');

end
