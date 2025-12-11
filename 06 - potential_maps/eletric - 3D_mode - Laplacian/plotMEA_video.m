function plotMEA_video(V_MEA, sample, MEA, MEA_plane_indx, lim)
    % Plot the MEA using patch
    patch('faces', MEA.faces, 'vertices', MEA.vertices, ...
        'FaceVertexCData', V_MEA(:,sample), 'FaceColor', 'interp', ...
        'FaceAlpha', 1, 'EdgeAlpha', 0, 'FaceLighting','gouraud');

    % Adjust axis properties
    axis equal;
    axis off;
    caxis(lim);

    % Add colorbar outside the plot figure
    hBar1 = colorbar();
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
    
    hold off
    hold off
end