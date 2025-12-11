function plotTank_video(V_TANK, sample, Plane, tank_plane_indx, lim)

% Plot the Tank using patch
patch('faces', Plane.faces, 'vertices', Plane.vertices, ...
    'FaceVertexCData', V_TANK(:,sample), 'FaceColor', 'interp', ...
    'FaceAlpha', 1, 'EdgeAlpha', 0, 'FaceLighting','gouraud');

% Adjust axis properties
axis equal;
axis off;
caxis(lim);

% Add colorbar
hBar1 = colorbar;
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

hold off;
hold off;

end