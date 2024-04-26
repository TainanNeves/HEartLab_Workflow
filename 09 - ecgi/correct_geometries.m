%% Data Loading

% Define file paths
heart_geo_file = 'C:\Users\HeartLAB\Documents\Documents\CinC 2024\ECGi\Dados\projected_signals_exp14.mat';
tank_geo_file = 'C:\Users\HeartLAB\Documents\Documents\CinC 2024\ECGi\Dados\LR_smoothed_tank.mat';

% Loading data
tank_data = load(tank_geo_file);
tank_geo = tank_data.(subsref(fieldnames(tank_data),substruct('{}',{1})));

heart_data = load(heart_geo_file);
heart_geo = heart_data.geometry_20000;

% Clear unnecessary variables
clear heart_geo_file tank_geo_file heart_data tank_data;

%% Correction of the faces vector (removing zero indices)
% Some heart geometries might contain zero indices. Check it and correct it
% if necessary.

for j=1:3
    for i=1:length(heart_geo.faces)
        heart_geo.faces(i,j) = heart_geo.faces(i,j) + 1;
    end
end

%% Shifting atrium geometry to the center of the tank
% Plot both tank and heart and verify if the heart is placed in the center
% of the tank before estimation. If not, use this code block for
% displacement.


for dim = 1:3
    if dim == 3
        dc = mean(heart_geo.vertices(:,dim)) - mean(tank_geo.vertices(:,dim));
    else
        dc = mean(heart_geo.vertices(:,dim));
    end

    % Shifting vertices
    for vertice_index = 1:length(heart_geo.vertices)
        vertices_new(vertice_index, dim) = heart_geo.vertices(vertice_index, dim) - dc;
    end
end

% Saving new geometry with displacement
heart_geo = struct('vertices', vertices_new, 'faces' , heart_geo.faces);

% Clear temporary variables
clear dim vertice_index dc vertices_new;

%% Plot
% Plot both tank and heart geometries

figure(1)
trisurf(heart_geo.faces, heart_geo.vertices(:,1), heart_geo.vertices(:,2), heart_geo.vertices(:,3), 'FaceAlpha',0.2,'FaceColor','y');
hold on;
trisurf(tank_geo.faces, tank_geo.vertices(:,1), tank_geo.vertices(:,2), tank_geo.vertices(:,3), 'FaceAlpha',0.2,'FaceColor','g');
