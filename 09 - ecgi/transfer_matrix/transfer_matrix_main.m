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
%% Define Conductivity

% Conductivity of blood, tissues, and air, respectively
conductivities = [0.3 0.2 0] * 0.01; % Old
%conductivities = [0 0.2 0] * 0.01; % New

%% Organize Geometries

% Store geometry data in a struct array
Model(1).vertices = heart_geo.vertices;
Model(1).faces = heart_geo.faces;
Model(2).vertices = tank_geo.vertices;
Model(2).faces = tank_geo.faces;

%% Create Matrix

disp('Calculating Transfer Matrix [Direct Problem]...')

% Save transfer matrix in a data matrix
[MTransfer] = BEM_TransfMat(Model, conductivities);

%% Saving

% Define file path for saving the transfer matrix
transfer_matrix_file = 'transfer_matrix.mat';

% Save the transfer matrix
save(transfer_matrix_file, 'MTransfer');

disp(['Transfer matrix saved to: ', transfer_matrix_file]);
