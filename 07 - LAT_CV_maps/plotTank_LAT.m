function plotTank_LAT(LAT_TANK, Plane, tank_plane_indx, title_str, lim, mode)
    % PLOTTANK_LAT plots the Local Activation Time (LAT) on a tank mesh.
    %
    % Syntax:
    %   plotTank_LAT(LAT_TANK, Plane, tank_plane_indx, title_str, lim, mode)
    %
    % Inputs:
    %   LAT_TANK:           Vector containing the LAT values for each vertex of the tank mesh.
    %   Plane:              Struct containing the vertices and faces information of the tank mesh.
    %   tank_plane_indx:    Indices of vertices corresponding to electrode positions on the tank mesh.
    %   title_str:          Title string for the plot.
    %   lim:                Color axis limits for the plot.
    %   mode:               Mode selection flag:
    %                           0: No adjustment.
    %                           1: Subtract the minimum LAT value from LAT_TANK (Local reference).
    %
    % Note:
    %   - The function assumes a certain format for the input 'Plane', where 'Plane.vertices'
    %     contains the coordinates of the vertices and 'Plane.faces' contains the connectivity
    %     information of the mesh faces.
    %   - 'LAT_TANK' and 'tank_plane_indx' should have matching dimensions, where each entry in
    %     'LAT_TANK' corresponds to the LAT value for the corresponding vertex in 'tank_plane_indx'.
    %   - 'lim' should be a 2-element vector specifying the color axis limits for the plot.
    %   - 'mode' determines whether to adjust LAT_TANK values. If mode is 1, the minimum LAT value
    %     will be subtracted from all LAT_TANK values.
    %
    % Example:
    %   % Define LAT_TANK, Plane, tank_plane_indx, title_str, lim, and mode
    %   plotTank_LAT(LAT_TANK, Plane, tank_plane_indx, 'Tank LAT Plot', [0 100], 1);
    %
    % See also:
    %   plot, patch, colorbar

    % Set the sampling frequency
    Fs = 4000; % Hz

    % Create a new figure for the Tank plot
    f1 = figure('color','white','Position', [160 40 800 400]);

    % Mode selection
    if mode == 1 % Local reference
        mini = min(LAT_TANK); % <-- This line references 'LAT_MEA', but it's not defined within this function
        LAT_TANK = LAT_TANK - mini;
    end

    % Plot the Tank using patch
    patch('faces', Plane.faces, 'vertices', Plane.vertices, ...
        'FaceVertexCData', LAT_TANK(:), 'FaceColor', 'interp', ...
        'FaceAlpha', 1, 'EdgeAlpha', 0, 'FaceLighting','gouraud');

    % Set the title for the plot
    title([title_str, ' | Local Activation Time']);

    % Adjust axis properties
    axis equal;
    axis off;
    caxis(lim);

    % Add colorbar outside the plot figure
    hBar1 = colorbar('Location', 'eastoutside', 'Position', [0.92 0.1 0.02 0.8]);
    c = parula(256);
    colormap(c);
    
    % Adjust the position of the colorbar to make room for the y-axis label
    pos = get(hBar1, 'Position');
    set(hBar1, 'Position', [pos(1)-0.02 pos(2) pos(3) pos(4)]);
    
    % Add y-axis label
    ylabel(hBar1, 'LAT [ms]', 'FontSize', 14);

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
