function [ax] = call_topoplot_eeglab(coherence_vals, chanlocs, plot_title, y_label, clim_min, clim_max, cmap)
    % topoplot_coherence.m
    % Plot magnitude squared coherence as a scalp topography using EEGLAB's topoplot.
    % Assumes:
    %   - EEG struct is already loaded in the workspace (with chanlocs)
    %   - coherence_vals is a 1 x 64 (or 64 x 1) vector of pre-computed MSC values
    %
    % cmap (optional): an Nx3 colormap matrix (e.g., parula(256), or a diverging colormap for
    %                   difference plots). Defaults to parula(256) if not supplied, preserving
    %                   the original behavior of this function for any existing callers.
    %
    % Returns ax: the axes handle this function drew into. IMPORTANT -- if you are calling this
    %   function multiple times to draw several topoplots into the same figure (e.g., tiles of a
    %   tiledlayout) and want each tile to keep ITS OWN colormap, do not rely on this function's
    %   internal colormap(gca, cmap) call alone. EEGLAB's topoplot() does not cleanly scope its
    %   own colormap handling to a single axes, so a LATER call to this function (for a
    %   different tile) can retroactively change an EARLIER tile's rendered colormap. Instead,
    %   capture the ax output from each call, and once ALL topoplots have been drawn, run a
    %   final pass re-applying colormap(ax, cmap) and clim(ax, [clim_min clim_max]) for every
    %   tile. See PlotGroupMSC_Topography.m for a working example of this pattern.
    
    % -------------------------------------------------------------------------
    % INPUTS — edit these
    % -------------------------------------------------------------------------
    %% plot_title     = 'Magnitude Squared Coherence';
    
    % Optional: set color axis limits (MSC ranges 0–1, but narrow if needed)
    %clim_min = 0;
    %%clim_max = 0.05;
    %% note:  if we want to compare two topoplots, we should use the same color axis limits for both, to make them comparable.
    %clim_max  = max(coherence_vals) * 1.1;  % set max to 10% above the max value in the data, for better contrast
    %%% clim_max = 0.04;  % set max to 0.5, which is a reasonable upper bound for MSC values in EEG data, and allows for better contrast in the plot.
    % -------------------------------------------------------------------------
    % Validate inputs
    % -------------------------------------------------------------------------
    % if ~exist('EEG', 'var') || ~isfield(EEG, 'chanlocs')
    %     error('EEG struct with chanlocs not found in workspace.');
    % end

    if nargin < 7 || isempty(cmap)
        cmap = parula(256);  % preserve original default colormap if not specified
    end
    
    if numel(coherence_vals) ~= length(chanlocs)
        error('Length of coherence_vals (%d) does not match number of channels (%d).', ...
        numel(coherence_vals), length(chanlocs));
    end
    
    % Ensure it's a row vector
    coherence_vals = coherence_vals(:)';
    
    % -------------------------------------------------------------------------
    % Plot
    % -------------------------------------------------------------------------
     %% commenting out this call to figure(), so that the topoplot will be drawn in the current figure, rather than creating a new figure each time.
     %% moreover, it will be part of a tiled layout, if we are plotting multiple topoplots in a single figure.
     %% that tiledlayout is outside this function, so we don't want to create a new figure here.  
     %  figure('Color', 'white', 'Position', [100 100 600 500]);
    
    topoplot(coherence_vals, chanlocs, ...
    'maplimits',  [clim_min, clim_max], ...  % color axis limits
    'electrodes', 'labelpoint', ...           % show channel labels + dots
    'colormap',   cmap, ...                   % colormap (parula by default; pass a diverging colormap for difference plots)
    'whitebk',    'on');                      % white background
    % Colorbar and title
    ax = gca;  % capture the axes topoplot just drew into, for the caller to use in a final colormap pass if needed
    cb = colorbar;
    ylabel(cb, y_label, 'FontSize', 12);
    clim(ax, [clim_min clim_max]);   % MATLAB R2022a+; use caxis() if on older MATLAB
    colormap(ax, cmap);               % best-effort immediate attempt; NOT sufficient on its own if
                                       % more topoplot() calls will follow in the same figure -- see
                                       % the function-level comment above and PlotGroupMSC_Topography.m
    title(ax, plot_title, 'FontSize', 14, 'FontWeight', 'bold');

end
