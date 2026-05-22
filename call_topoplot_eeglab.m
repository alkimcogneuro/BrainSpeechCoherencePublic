function [] = call_topoplot_eeglab(coherence_vals, chanlocs)
    %% topoplot_coherence.m
    % Plot magnitude squared coherence as a scalp topography using EEGLAB's topoplot.
    % Assumes:
    %   - EEG struct is already loaded in the workspace (with chanlocs)
    %   - coherence_vals is a 1 x 64 (or 64 x 1) vector of pre-computed MSC values
    
    % -------------------------------------------------------------------------
    % INPUTS — edit these
    % -------------------------------------------------------------------------
    plot_title     = 'Magnitude Squared Coherence';
    
    % Optional: set color axis limits (MSC ranges 0–1, but narrow if needed)
    clim_min = 0;
    %%clim_max = 1;
    clim_max  = max(coherence_vals) * 1.1;  % set max to 10% above the max value in the data, for better contrast
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
    'colormap',   jet(256), ...               % colormap (jet is traditional for EEG)
    'whitebk',    'on');                      % white background
    
    % Colorbar and title
    cb = colorbar;
    ylabel(cb, 'MSC', 'FontSize', 12);
    clim([clim_min clim_max]);   % MATLAB R2022a+; use caxis() if on older MATLAB
    title(plot_title, 'FontSize', 14, 'FontWeight', 'bold');

end
