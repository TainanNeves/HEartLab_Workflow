%% Plot 3D geometry and meas electrodes positions
%
% Author: Angélica Quadros
% Affiliation: HeartLab, UFABC
% Year: 2024
%
% Description:
% This script plots the 3D heart geometry and electrode positions from the MEA
% on the surface of the heart geometry at specific time instants.
%

%% Electrodes position plot on 3D geometry

% Load heart geometry
heart_geo = projections.geometry_HR;

% Define time for plotting (seconds)
inst1 = 0.096;
fs = 4000; % Sampling frequency (Hz)

% Conversion to samples
inst1 = inst1 * fs + 1;

% Organizing signal
% sinal1 = x_hat(:, inst1);

% Create the figure
figure;
x = heart_geo.vertices(:,1);
y = heart_geo.vertices(:,2);
z = heart_geo.vertices(:,3);
faces = heart_geo.faces;

% Choose heart resolution
heart_resolution = 'MEAS_HR_IDX';
% MEAS_HR_IDX
% MEAS_IDX_20000
% MEAS_IDX_10000
% MEAS_IDX_2500
% MEAS_IDX_1200

idx1 = projections.(heart_resolution).MEA1;
idx2 = projections.(heart_resolution).MEA2;
idx3 = projections.(heart_resolution).MEA3;

% Plot heart surface
trisurf(faces, x, y, z, 'facecolor', 'interp', 'LineStyle', 'none');
hold on;
% Plot electrodes at their respective positions with red color
scatter3(x(idx1), y(idx1), z(idx1), 'ro', 'filled');
scatter3(x(idx2), y(idx2), z(idx2), 'ro', 'filled');
scatter3(x(idx3), y(idx3), z(idx3), 'ro', 'filled');
colormap('jet');

hold off;
grid off;
axis off;

xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
view(3);
axis equal;

%% Electrodes position and number plot in a specific instant
% This section plots the 3D geometry again, but this time displays the electrode
% positions and their indices at a specific time instant in the signal.

% Create a new figure
figure;
x = heart_geo.vertices(:,1);
y = heart_geo.vertices(:,2);
z = heart_geo.vertices(:,3);
faces = heart_geo.faces;

fs = 4000; % Sampling frequency (Hz)

% Define time for plotting (seconds)
inst1 = 0.5;
inst_samples = inst1 * fs;

% Choose heart resolution 
heart_resolution = 'MEAS_HR_IDX';
% MEAS_HR_IDX
% MEAS_IDX_20000
% MEAS_IDX_10000
% MEAS_IDX_2500
% MEAS_IDX_1200

idx1 = projections.(heart_resolution).MEA1;
idx2 = projections.(heart_resolution).MEA2;
idx3 = projections.(heart_resolution).MEA3;

% Plot the heart surface with signal values at the given instant
trisurf(faces, x, y, z, estimated_signal(:,inst_samples),'facecolor', 'interp', 'LineStyle', 'none');
title(['Instant ', num2str(inst1), 's Frame', num2str(inst_samples)]);
colormap('jet');

% Adding light
lighting gouraud;           % Smooth lighting
light('Position',[10 10 0]);  % Add a light source
camlight('headlight');  % Light from the camera direction
material shiny;             % Make the surface shiny

hold on;

% Plot the electrodes at their positions
scatter3(x(idx1), y(idx1), z(idx1), 'ro', 'filled');
scatter3(x(idx2), y(idx2), z(idx2), 'ro', 'filled');
scatter3(x(idx3), y(idx3), z(idx3), 'ro', 'filled');

% Calculate the centroid of the heart geometry
center_x = mean(x);
center_y = mean(y);
center_z = mean(z);

% Scale factor to position the text slightly outside the surface
scale_factor = 1.15;

% Annotate MEA1 electrodes with their index numbers
for i = 1:length(idx1)
    text(scale_factor * (x(idx1(i)) - center_x) + center_x, ...
         scale_factor * (y(idx1(i)) - center_y) + center_y, ...
         scale_factor * (z(idx1(i)) - center_z) + center_z, ...
         num2str(i), ...
        'Color', 'k', 'FontSize', 17, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
         'Margin', 1);
end

% Annotate MEA2 electrodes with their index numbers
for i = 1:length(idx2)
   text(scale_factor * (x(idx2(i)) - center_x) + center_x, ...
         scale_factor * (y(idx2(i)) - center_y) + center_y, ...
         scale_factor * (z(idx2(i)) - center_z) + center_z, ...
         num2str(i+16), ...
        'Color', 'k', 'FontSize', 17, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
         'Margin', 1);
end

% Annotate MEA3 electrodes with their index numbers
for i = 1:length(idx3)
    text(scale_factor * (x(idx3(i)) - center_x) + center_x, ...
         scale_factor * (y(idx3(i)) - center_y) + center_y, ...
         scale_factor * (z(idx3(i)) - center_z) + center_z, ...
         num2str(i+64), ...
        'Color', 'k', 'FontSize', 15, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
         'Margin', 1);
end

hold off;
grid off;
axis off;

xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
view(3);
axis equal;

%% Plot tank electrodes position
tank_electrodes = [129:174, 177:190];

x = tank_geo.vertices(:,1);
y = tank_geo.vertices(:,2);
z = tank_geo.vertices(:,3);
faces = tank_geo.faces;

trisurf(tank_geo.faces, tank_geo.vertices(:,1), tank_geo.vertices(:,2), tank_geo.vertices(:,3), ...
    'FaceAlpha', 0.2, 'FaceColor', 'g');
hold on;
scatter3(x(idx), y(idx), z(idx), 'ro', 'filled');

% Scale factor to position the text slightly outside the surface
scale_factor = 1.15;

% Annotate MEA1 electrodes with their index numbers
for i = 1:length(idx)
    text(scale_factor * (x(idx(i)) - center_x) + center_x, ...
         scale_factor * (y(idx(i)) - center_y) + center_y, ...
         scale_factor * (z(idx(i)) - center_z) + center_z, ...
         num2str(tank_electrodes(i)), ...
        'Color', 'k', 'FontSize', 17, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
         'Margin', 1);
end
