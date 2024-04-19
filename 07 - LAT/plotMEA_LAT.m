function plotMEA_LAT(LAT_MEA, MEA, MEA_plane_indx, title_str, lim, mode)
    % PLOTMEA_LAT plots the Local Activation Time (LAT) on a Microelectrode Array (MEA) mesh.
    %
    % Syntax:
    %   plotMEA_LAT(LAT_MEA, MEA, MEA_plane_indx, title_str, lim, mode)
    %
    % Inputs:
    %   LAT_MEA:            Vector containing the LAT values for each vertex of the MEA mesh.
    %   MEA:                Struct containing the vertices and faces information of the MEA mesh.
    %   MEA_plane_indx:     Indices of vertices corresponding to electrode positions on the MEA mesh.
    %   title_str:          Title string for the plot.
    %   lim:                Color axis limits for the plot.
    %   mode:               Mode selection flag:
    %                           0: No adjustment.
    %                           1: Subtract the minimum LAT value from LAT_MEA (Local reference).
    %
    % Note:
    %   - The function assumes a certain format for the input 'MEA', where 'MEA.vertices'
    %     contains the coordinates of the vertices and 'MEA.faces' contains the connectivity
    %     information of the mesh faces.
    %   - 'LAT_MEA' and 'MEA_plane_indx' should have matching dimensions, where each entry in
    %     'LAT_MEA' corresponds to the LAT value for the corresponding vertex in 'MEA_plane_indx'.
    %   - 'lim' should be a 2-element vector specifying the color axis limits for the plot.
    %   - 'mode' determines whether to adjust LAT_MEA values. If mode is 1, the minimum LAT value
    %     will be subtracted from all LAT_MEA values.
    %
    % Example:
    %   % Define LAT_MEA, MEA, MEA_plane_indx, title_str, lim, and mode
    %   plotMEA_LAT(LAT_MEA, MEA, MEA_plane_indx, 'MEA LAT Plot', [0 100], 1);
    %
    % See also:
    %   plot, patch, colorbar

    % Set the sampling frequency
    Fs = 4000; % Hz

    % Create a new figure for the MEA plot
    f1 = figure('color','white','Position', [40 40 500 400]);

    % Mode selection
    if mode == 1 % Local reference
        mini = min(LAT_MEA);
        LAT_MEA = LAT_MEA - mini;
    end
    
    % Plot the MEA using patch
    patch('faces', MEA.faces, 'vertices', MEA.vertices, ...
        'FaceVertexCData', LAT_MEA(:), 'FaceColor', 'interp', ...
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
