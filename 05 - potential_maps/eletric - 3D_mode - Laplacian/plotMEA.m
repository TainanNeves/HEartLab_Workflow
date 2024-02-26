function plotMEA(V_MEA, sample, MEA, MEA_plane_indx, title_str, lim)
% PLOTMEA - Visualize electric potential data on a 3D mesh for a specific electrode configuration.
%
% Usage:
%   plotMEA(V_MEA, sample, MEA, MEA_plane_indx, title_str, lim)
%
% Inputs:
%   - V_MEA: Matrix containing electric potential data. Rows represent electrodes,
%            and columns represent samples.
%   - sample: Sample index for plotting.
%   - MEA: Structure containing mesh information (faces and vertices).
%   - MEA_plane_indx: Indices representing electrode positions on the mesh.
%   - title_str: String specifying the title of the plot.
%   - lim: Two-element vector specifying the color axis limits.
%
% Output:
%   The function generates a 3D plot visualizing electric potential data on the
%   specified 3D mesh for a specific electrode configuration.

% Create a new figure for the MEA plot
    f1 = figure('color','white','Position', [40 40 500 400]);

    % Plot the MEA using patch
    patch('faces', MEA.faces, 'vertices', MEA.vertices, ...
        'FaceVertexCData', V_MEA(:,sample), 'FaceColor', 'interp', ...
        'FaceAlpha', 1, 'EdgeAlpha', 0, 'FaceLighting','gouraud');

    % Set the title for the plot
    title([title_str, ' | Potential Map | Sample: ',num2str(sample)]);

    % Adjust axis properties
    axis equal;
    axis off;
    caxis(lim);

    % Add colorbar outside the plot figure
    hBar1 = colorbar('Location', 'eastoutside', 'Position', [0.92 0.1 0.02 0.8]);
    c = jet(256);
    colormap(c);
    
    % Add y-axis label
    ylabel(hBar1, 'Potential [$\mu$V]', 'FontSize', 14, 'Interpreter', 'latex');

    % Add a transparent overlay to highlight electrode positions
    hold on;
    set(gcf, 'color', 'white');
    patch('faces', MEA.faces, 'vertices', MEA.vertices, ...
        'FaceVertexCData', ones(221,1), 'FaceColor', 'interp', ...
        'FaceAlpha', 0, 'EdgeAlpha', 0.05, 'FaceLighting','gouraud');

    % Highlight electrode positions with black circles
    axis equal;
    hold on;
    plot3(MEA.vertices(MEA_plane_indx,1), MEA.vertices(MEA_plane_indx,2), MEA.vertices(MEA_plane_indx,3),'or', 'Color', 'k', 'LineWidth', 2);

    % Set font size for the axis
    set(gca,'fontsize', 14);

end