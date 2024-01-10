function plotTank(DF_TANK, Plane, tank_plane_indx, title_str, lim)
% PLOTTANK Plots the Tank using patch with dominant frequency information.
%
%   plotTank(DF_TANK, Plane, tank_plane_indx, title_str, lim) takes the following inputs:
%   - DF_TANK: Matrix representing the dominant frequency values for each vertex of the Tank.
%   - Plane: Structure containing Tank information with fields 'faces' and 'vertices'.
%   - tank_plane_indx: Index of electrodes to be highlighted on the Tank plot.
%   - title_str: String for the title of the plot.
%   - lim: Limits for the color axis.
%
%   Example:
%   DF_TANK = rand(1201, 1) * 10; % Example dominant frequency values
%   Plane = load('your_Tank_data.mat'); % Replace with actual Tank data
%   tank_plane_indx = 1:10; % Replace with actual electrode indices
%   title_str = 'Tank Plot';
%   lim = [0, 10]; % Replace with actual color axis limits
%   plotTank(DF_TANK, Plane, tank_plane_indx, title_str, lim);
%
%   Note: Make sure to replace the placeholder values with your actual data.

    % Set the sampling frequency
    Fs = 4000; % Hz

    % Create a new figure for the Tank plot
    f1 = figure('color','white','Position', [160 40 800 400]);

    % Plot the Tank using patch
    patch('faces', Plane.faces, 'vertices', Plane.vertices, ...
        'FaceVertexCData', DF_TANK(:), 'FaceColor', 'interp', ...
        'FaceAlpha', 1, 'EdgeAlpha', 0, 'FaceLighting','gouraud');

    % Set the title for the plot
    title([title_str, ' | Dominant Frequency']);

    % Adjust axis properties
    axis equal;
    axis off;
    caxis(lim);

    % Add colorbar outside the plot figure
    hBar1 = colorbar('Location', 'eastoutside', 'Position', [0.92 0.1 0.02 0.8]);
    c = jet(256);
    colormap(c);
    
    % Adjust the position of the colorbar to make room for the y-axis label
    pos = get(hBar1, 'Position');
    set(hBar1, 'Position', [pos(1)-0.02 pos(2) pos(3) pos(4)]);
    
    % Add y-axis label
    ylabel(hBar1, 'Frequency [Hz]', 'FontSize', 14);

    % Add a transparent overlay to highlight electrode positions
    hold on;
    set(gcf, 'color', 'white');
    patch('faces', Plane.faces, 'vertices', Plane.vertices, ...
        'FaceVertexCData', ones(1201,1), 'FaceColor', 'interp', ...
        'FaceAlpha', 0, 'EdgeAlpha', 0.05, 'FaceLighting','gouraud');

    % Highlight electrode positions with red circles
    axis equal;
    hold on;
    plot3(Plane.vertices(tank_plane_indx,1), Plane.vertices(tank_plane_indx,2), Plane.vertices(tank_plane_indx,3),'or', 'color', 'k', 'LineWidth', 2)

end