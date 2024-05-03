function plotMEA(V_MEA, MEA, MEA_plane_indx, title_str, frame, lim)
% PLOTMEA - Visualizes Multi-Electrode Array (MEA) data.
%
% Syntax:
%   plotMEA(V_MEA, MEA, MEA_plane_indx, title_str, frame, lim)
%
% Description:
%   This function plots the phase data for a Multi-Electrode Array (MEA) at a
%   specific frame index.
%
% Input:
%   - V_MEA: Interpolated and filtered phase data for the MEA electrodes.
%   - MEA: MEA geometry structure with 'faces' and 'vertices'.
%   - MEA_plane_indx: Indices representing the positions of MEA electrodes.
%   - title_str: Title string for the plot.
%   - frame: Frame index for visualization.
%   - lim: Color axis limits for the phase map.
%
% Author:
%   Tainan Neves, HEartLab - UFABC
%
% Example:
%   plotMEA(data, MEA, MEA_plane_indx, 'MEA 1 (Right Atrium)', 50, [-pi pi]);

% Set the sampling frequency
Fs = 4000; % Hz

% Create a new figure for the MEA plot
f1 = figure('color','white','Position', [40 40 500 550]);

% Plot the MEA using patch
patch('faces', MEA.faces, 'vertices', MEA.vertices, ...
    'FaceVertexCData', V_MEA(:, frame), 'FaceColor', 'interp', ...
    'FaceAlpha', 1, 'EdgeAlpha', 0, 'FaceLighting','gouraud');

% Set the title for the plot
title([title_str, ' | time: ', num2str(frame)]);

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
