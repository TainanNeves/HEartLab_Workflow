function plotTank(V_TANK, Plane, tank_plane_indx, frame, lim)
% PLOTTANK - Visualizes Tank electrode data.
%
% Syntax:
%   plotTank(V_TANK, Plane, tank_plane_indx, frame, lim)
%
% Description:
%   This function plots the phase data for Tank electrodes at a specific frame index.
%
% Input:
%   - V_TANK: Interpolated and filtered phase data for the Tank electrodes.
%   - Plane: Tank geometry structure with 'faces' and 'vertices'.
%   - tank_plane_indx: Indices representing the positions of Tank electrodes.
%   - frame: Frame index for visualization.
%   - lim: Color axis limits for the phase map.
%
% Author:
%   Tainan Neves, HEartLab - UFABC
%
% Example:
%   plotTank(data, Plane, tank_plane_indx, 50, [-pi pi]);

% Set the sampling frequency
Fs = 4000; % Hz

% Create a new figure for the Tank plot
f1 = figure('color','white','Position', [160 40 800 600]);

% Plot the Tank using patch
patch('faces', Plane.faces, 'vertices', Plane.vertices, ...
    'FaceVertexCData', V_TANK(:, frame), 'FaceColor', 'interp', ...
    'FaceAlpha', 1, 'EdgeAlpha', 0, 'FaceLighting','gouraud');

% Set the title for the plot
title(['TANK | time: ', num2str(frame)]);

% Adjust axis properties
axis equal;
axis off;
caxis(lim);

% Add colorbar
load("PS_colormaps.mat");
hBar1 = colorbar('southoutside');
colormap(mycmap);
ylabel(hBar1, ['Radians [' char(960) ']'],'FontSize',14);
hBar1.Ticks = [-2, 0, 2];


% Add a transparent overlay to highlight electrode positions
hold on;
set(gcf, 'color', 'white');
patch('faces', Plane.faces, 'vertices', Plane.vertices, ...
    'FaceVertexCData', ones(1201,1), 'FaceColor', 'interp', ...
    'FaceAlpha', 0, 'EdgeAlpha', 0.05, 'FaceLighting','gouraud');

% Highlight electrode positions with red circles
axis equal;
hold on;
plot3(Plane.vertices(tank_plane_indx,1), Plane.vertices(tank_plane_indx,2), Plane.vertices(tank_plane_indx,3),'or', 'color', 'k','LineWidth',2);

% Set font size for the axis
set(gca,'fontsize', 14);

end
