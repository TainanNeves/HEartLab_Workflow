function plot_potential_map(data, sample, varargin)
% PLOT_POTENTIAL_MAP - Visualize pre-interpolated electric potential data
%
% Usage:
%   plot_potential_map(data, sample)
%   plot_potential_map(data, sample, 'lim', [-100 100], 'title', 'My Plot', 'type', 'MEA')
%
% Inputs:
%   - data: 3D potential data matrix (grid x grid x samples)
%   - sample: Sample index/indices for plotting (single number or array like 1:10:100)
%
% Optional Parameters (name-value pairs):
%   - 'lim': Color axis limits (default: [], auto-scaling)
%   - 'title': Plot title string (default: 'Potential Map')
%   - 'type': Data type - 'MEA' for 11x11 or 'TANK' for 25x25 (default: auto-detected)
%   - 'colormap': Colormap to use (default: 'jet')
%   - 'colorbar': Show colorbar (default: true)
%   - 'grid': Show grid coordinates (default: false)

% Default parameters
defaultLim = [];
defaultTitle = 'Potential Map';
defaultType = '';
defaultColormap = 'jet';
defaultColorbar = true;
defaultGrid = false;

% Parse input arguments
p = inputParser;
addRequired(p, 'data', @(x) isnumeric(x) && ndims(x) == 3);
addRequired(p, 'sample', @(x) isnumeric(x) && all(x > 0));
addParameter(p, 'lim', defaultLim, @(x) isempty(x) || (isnumeric(x) && length(x) == 2));
addParameter(p, 'title', defaultTitle, @ischar);
addParameter(p, 'type', defaultType, @ischar);
addParameter(p, 'colormap', defaultColormap, @ischar);
addParameter(p, 'colorbar', defaultColorbar, @islogical);
addParameter(p, 'grid', defaultGrid, @islogical);
parse(p, data, sample, varargin{:});

% Extract parameters
lim = p.Results.lim;
title_str = p.Results.title;
data_type = p.Results.type;
cmap = p.Results.colormap;
show_colorbar = p.Results.colorbar;
show_grid = p.Results.grid;

% Auto-detect data type if not provided
if isempty(data_type)
    [grid_size, ~, ~] = size(data);
    if grid_size == 11
        data_type = 'MEA';
    elseif grid_size == 25
        data_type = 'TANK';
    else
        data_type = 'UNKNOWN';
        warning('Unknown grid size. Using default settings.');
    end
end

% Handle multiple samples
if length(sample) > 1
    plotMultipleSamples(data, sample, title_str, lim, data_type, cmap, show_colorbar, show_grid);
else
    plotSingleSample(data, sample, title_str, lim, data_type, cmap, show_colorbar, show_grid);
end
end

function plotSingleSample(data, sample, title_str, lim, data_type, cmap, show_colorbar, show_grid)
% PLOTSINGLESAMPLE - Plot a single sample

% Determine figure size based on data type
if strcmp(data_type, 'TANK')
    fig_size = [160 40 800 400]; % Wider figure for TANK
else
    fig_size = [40 40 500 400];
end

f1 = figure('color','white','Position', fig_size);

% Extract single frame
I = squeeze(data(:,:,sample));

% Create 2D surface plot
I = flipud(I); % Vertical Flip because surf() invert the y-axis
surf(I, 'EdgeColor', 'none', 'FaceColor', 'interp');
view(2); % 2D view from above

% Set rectangular aspect ratio for TANK (2:1 width:height)
if strcmp(data_type, 'TANK')
    aspect_ratio = [2 1 1]; % [x y z] aspect ratio
    pbaspect(aspect_ratio);
else
    axis equal; % Keep square for MEA
end

axis tight;

% Set title and color limits
title([title_str, ' | Sample: ', num2str(sample)], 'Interpreter', 'none');
if ~isempty(lim)
    caxis(lim);
end

% Set colormap
colormap(cmap);

% Add colorbar if requested
if show_colorbar
    hBar1 = colorbar('Location', 'eastoutside', 'Position', [0.92 0.1 0.02 0.8]);
    ylabel(hBar1, 'Potential [$\mu$V]', 'FontSize', 14, 'Interpreter', 'latex');
end

% Show grid coordinates if requested
if show_grid
    xlabel('X Grid');
    ylabel('Y Grid');
else
    set(gca, 'XTick', [], 'YTick', []);
end

set(gca,'fontsize', 14);
end

function plotMultipleSamples(data, samples, title_str, lim, data_type, cmap, show_colorbar, show_grid)
% PLOTMULTIPLESAMPLES - Plot multiple samples in subplots

num_samples = length(samples);

% Determine figure size based on data type
if strcmp(data_type, 'TANK')
    fig_size = [100 100 1400 900];
else
    fig_size = [100 100 1000 800];
end

% Calculate subplot arrangement (always 8 columns)
cols = 8;
rows = ceil(num_samples / cols);

% Create figure
f1 = figure('color','white','Position', fig_size);

for i = 1:num_samples
    subplot(rows, cols, i);
    
    % Extract data
    I = squeeze(data(:,:,samples(i)));
    
    % Create 2D surface plot
    I = flipud(I); % Vertical Flip because surf() invert the y-axis
    surf(I, 'EdgeColor', 'none', 'FaceColor', 'interp');
    view(2);
    
    % Set rectangular aspect ratio for TANK
    if strcmp(data_type, 'TANK')
        aspect_ratio = [2 1 1];
        pbaspect(aspect_ratio);
    else
        axis equal;
    end
    
    axis tight;
    
    if ~isempty(lim)
        caxis(lim);
    end
    
    title(['S: ', num2str(samples(i))], 'FontSize', 10);
    
    if ~show_grid
        set(gca, 'XTick', [], 'YTick', []);
    end
end

% Set colormap for all subplots
colormap(cmap);

% Add single colorbar for the entire figure if requested
if show_colorbar
    hBar1 = colorbar('Location', 'eastoutside', 'Position', [0.93 0.1 0.02 0.8]);
    ylabel(hBar1, 'Potential [$\mu$V]', 'FontSize', 14, 'Interpreter', 'latex');
end

% Main title
sgtitle([title_str, ' | Multiple Samples'], 'FontSize', 16, 'Interpreter', 'none');
end






