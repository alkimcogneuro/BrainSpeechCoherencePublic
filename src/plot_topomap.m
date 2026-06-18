function plot_topomap(voltages, chanlocs, varargin)
    fprintf('running ak plot_topomap with %d channels\n', length(chanlocs));
    
% PLOT_TOPOMAP  Plot a topographic map of voltage values at scalp locations
%% --- Parse optional inputs ---
p = inputParser;
addParameter(p, 'colormap',     'jet',   @ischar);
addParameter(p, 'numcontours',  6,       @isnumeric);
addParameter(p, 'showlabels',   false,   @islogical);
addParameter(p, 'headrad',      0.5,     @isnumeric);
addParameter(p, 'interpmethod', 'v4',    @ischar);
parse(p, varargin{:});
opts = p.Results;

voltages = voltages(:);  % ensure column vector

%% --- Filter Missing Coordinates ---
% Identify which channels actually have spatial coordinates
has_coords = arrayfun(@(c) ~isempty(c.theta) && ~isempty(c.radius), chanlocs);

% Filter out the channels missing coordinates (Ensure indexing is column-wise)
valid_chanlocs = chanlocs(has_coords);
valid_voltages = voltages(has_coords(:)); 

%% --- Convert EEGLAB polar coords to 2D Cartesian ---
theta_rad = deg2rad([valid_chanlocs.theta]');
radius    = [valid_chanlocs.radius]';

% Standard EEGLAB convention
x =  sin(theta_rad) .* radius;
y =  cos(theta_rad) .* radius;

%% --- Create interpolation grid ---
headrad   = opts.headrad;
lin       = linspace(-headrad * 1.1, headrad * 1.1, 300);
[Xi, Yi]  = meshgrid(lin, lin);

% Interpolate scattered electrode values onto grid
Zi = griddata(x, y, valid_voltages, Xi, Yi, opts.interpmethod);

% Mask everything outside the head circle
mask         = sqrt(Xi.^2 + Yi.^2) > headrad;
Zi(mask)     = NaN;

%% --- Plot ---
figure('Color', 'w', 'Position', [100 100 500 500]);
hold on;

% Filled colour surface
h = pcolor(Xi, Yi, Zi);
set(h, 'EdgeColor', 'none', 'FaceColor', 'interp');
shading interp;
colormap(opts.colormap);
colorbar;

% Contour lines
contour(Xi, Yi, Zi, opts.numcontours, 'k', 'LineWidth', 0.8);

%% --- Draw head, nose, and ears ---
theta_circle = linspace(0, 2*pi, 360);
% Head
plot(headrad * cos(theta_circle), headrad * sin(theta_circle), 'k', 'LineWidth', 2);
% Nose
nose_x = headrad * [-0.05  0.00  0.05];
nose_y = headrad * [ 1.00  1.12  1.00];
plot(nose_x, nose_y, 'k', 'LineWidth', 2);
% Ears
ear_w = 0.03; ear_h = 0.10;
rectangle('Position', [-headrad - ear_w, -ear_h/2, ear_w, ear_h], 'Curvature', 0.5, 'EdgeColor', 'k', 'LineWidth', 2);
rectangle('Position', [ headrad + 0.015,  -ear_h/2, ear_w, ear_h], 'Curvature', 0.5, 'EdgeColor', 'k', 'LineWidth', 2);

%% --- Plot electrode locations ---
plot(x, y, 'k.', 'MarkerSize', 6);
if opts.showlabels
    labels = {valid_chanlocs.labels};
    for i = 1:numel(labels)
        text(x(i), y(i) + 0.02, labels{i}, 'FontSize', 6, 'HorizontalAlignment', 'center', 'Color', 'k');
    end
end

%% --- Formatting ---
axis equal off;
clim_val = max(abs(valid_voltages)); % Use max of valid voltages
clim([-clim_val  clim_val]);         % symmetric color scale around zero
title('Topographic Voltage Map', 'FontSize', 13);
hold off;
end