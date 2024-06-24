function plotTank_CL(CL_TANK, Plane, tank_plane_indx, title_str, lim)


    % Set the sampling frequency
    Fs = 4000; % Hz

    % Create a new figure for the Tank plot
    f1 = figure('color','white','Position', [160 40 800 400]);

    % Plot the Tank using patch
    patch('faces', Plane.faces, 'vertices', Plane.vertices, ...
        'FaceVertexCData', CL_TANK(:), 'FaceColor', 'interp', ...
        'FaceAlpha', 1, 'EdgeAlpha', 0, 'FaceLighting','gouraud');

    % Set the title for the plot
    title([title_str, ' | Cycle Length']);

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
    ylabel(hBar1, 'Cycle Length [ms]', 'FontSize', 14);

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