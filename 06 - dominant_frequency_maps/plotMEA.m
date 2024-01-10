function plotMEA(DF_MEA, MEA, MEA_plane_indx, title_str, lim)
% PLOTMEA Plots the MEA using patch with dominant frequency information.
%
%   plotMEA(DF_MEA, MEA, MEA_plane_indx, title_str, lim) takes the following inputs:
%   - DF_MEA: Matrix representing the dominant frequency values for each vertex of the MEA.
%   - MEA: Structure containing MEA information with fields 'faces' and 'vertices'.
%   - MEA_plane_indx: Index of electrodes to be highlighted on the MEA plot.
%   - title_str: String for the title of the plot.
%   - lim: Limits for the color axis.
%
%   Example:
%   DF_MEA = rand(221, 1) * 10; % Example dominant frequency values
%   MEA = load('your_MEA_data.mat'); % Replace with actual MEA data
%   MEA_plane_indx = 1:10; % Replace with actual electrode indices
%   title_str = 'MEA Plot';
%   lim = [0, 10]; % Replace with actual color axis limits
%   plotMEA(DF_MEA, MEA, MEA_plane_indx, title_str, lim);
%
%   Note: Make sure to replace the placeholder values with your actual data.

    % Set the sampling frequency
    Fs = 4000; % Hz

    % Create a new figure for the MEA plot
    f1 = figure('color','white','Position', [40 40 500 400]);

    % Plot the MEA using patch
    patch('faces', MEA.faces, 'vertices', MEA.vertices, ...
        'FaceVertexCData', DF_MEA(:), 'FaceColor', 'interp', ...
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
    patch('faces', MEA.faces, 'vertices', MEA.vertices, ...
        'FaceVertexCData', ones(221,1), 'FaceColor', 'interp', ...
        'FaceAlpha', 0, 'EdgeAlpha', 0.05, 'FaceLighting','gouraud');

    % Highlight electrode positions with red circles
    axis equal;
    hold on;
    plot3(MEA.vertices(MEA_plane_indx,1), MEA.vertices(MEA_plane_indx,2), MEA.vertices(MEA_plane_indx,3),'or', 'Color', 'k', 'LineWidth', 2);

    % Set font size for the axis
    set(gca,'fontsize', 14);

end
