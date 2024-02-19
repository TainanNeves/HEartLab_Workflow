function plotTank(V_TANK, sample, Plane, tank_plane_indx, title_str, lim)
% PLOTTANK - Visualize electric potential data on a 3D tank-shaped mesh.
%
% Usage:
%   plotTank(V_TANK, sample, Plane, tank_plane_indx, title_str, lim)
%
% Inputs:
%   - V_TANK: Matrix containing electric potential data. Rows represent electrodes,
%             and columns represent samples.
%   - sample: Sample index for plotting.
%   - Plane: Structure containing mesh information (faces and vertices) for the tank.
%   - tank_plane_indx: Indices representing electrode positions on the tank mesh.
%   - title_str: String specifying the title of the plot.
%   - lim: Two-element vector specifying the color axis limits.
%
% Output:
%   The function generates a 3D plot visualizing electric potential data on the
%   specified tank-shaped mesh for a specific electrode configuration.

    % Create a new figure for the Tank plot
    f1 = figure('color','white','Position', [160 40 800 400]);

    % Plot the Tank using patch
    patch('faces', Plane.faces, 'vertices', Plane.vertices, ...
        'FaceVertexCData', V_TANK(:,sample), 'FaceColor', 'interp', ...
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
    patch('faces', Plane.faces, 'vertices', Plane.vertices, ...
        'FaceVertexCData', ones(1201,1), 'FaceColor', 'interp', ...
        'FaceAlpha', 0, 'EdgeAlpha', 0.05, 'FaceLighting','gouraud');

    % Highlight electrode positions with black circles
    axis equal;
    hold on;
    plot3(Plane.vertices(tank_plane_indx,1), Plane.vertices(tank_plane_indx,2), Plane.vertices(tank_plane_indx,3),'or', 'color', 'k', 'LineWidth', 2)

end