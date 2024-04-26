%% Data Loading

% ECG signal
signal_file = 'C:\Users\HeartLAB\Documents\Documents\Current codes\filtered_recordings\data_filtered_sync_E14_F4_R10.mat';

% Electrodes indices
electrodes_idx_file = 'C:\Users\HeartLAB\Documents\Documents\CinC 2024\ECGi\Data\electrodes_LR.mat';

% Heart geometry
heart_geo_file = 'C:\Users\HeartLAB\Documents\Documents\CinC 2024\ECGi\Data\heart_geometry_exp14.mat';

% Tank geometry
tank_geo_file = 'C:\Users\HeartLAB\Documents\Documents\CinC 2024\ECGi\Data\LR_smoothed_tank.mat';

% Transfer matrix
mtransfer_file = 'C:\Users\HeartLAB\Documents\Documents\CinC 2024\ECGi\Data\MTransfer_exp14_LR_20000.mat';

%% Loading Data

% Extracting tank signals
signal_data = load(signal_file);
signal_data = signal_data.(subsref(fieldnames(signal_data),substruct('{}',{1})));
signal = signal_data.EL([129:174,177:190],:); % keeping only 60 electrodes

% Electrodes
electrodes_data = load(electrodes_idx_file);
electrodes = electrodes_data.(subsref(fieldnames(electrodes_data),substruct('{}',{1})));

% Tank
tank_data = load(tank_geo_file);
tank_geo = tank_data.(subsref(fieldnames(tank_data),substruct('{}',{1})));

% Heart
heart_data = load(heart_geo_file);
heart_geo = heart_data.(subsref(fieldnames(heart_data),substruct('{}',{1})));

% Transfer Matrix
mtransfer_data = load(mtransfer_file);
MTransfer = mtransfer_data.(subsref(fieldnames(mtransfer_data),substruct('{}',{1})));

clear heart_geo_file tank_geo_file mtransfer_file signal_file electrodes_idx_file heart_data tank_data mtransfer_data electrodes_data;

%% Converting to mV (if needed)

% signal = 0.00019*signal;
% signal = signal/1000;

%% Defining Estimation Parameters

lambda = 10.^(-0.5:-0.5:-12.5); % arbitrary
order = 0; % 0, 1 or 2 
reg_param_method = 'i';
compute_params = 1;
model = 'SAF';
fs = 4000;
SNR = 100;

%% Interpolating ECG Signal

y = signal;

idx = electrodes';

[lap, edge] = mesh_laplacian(tank_geo.vertices, tank_geo.faces); 
[int, keepindex, repindex] = mesh_laplacian_interp_current(lap, idx);

% Keep measured data
interp_signal = int *  y;  
for i = 1:60
    interp_signal(idx(1, i), :) = y(i, :);  
end

%% Defining Transfer Matrix

A = MTransfer;

%% Estimation Calculation
Fs = 4000;

% Estimation time setting
init_time = 3.5;
final_time = 3.8;

% Conversion to samples
init_sample = init_time * Fs + 1;
final_sample = final_time * Fs;

% Estimation
[AA, L, LL] = precompute_matrices(A, order, heart_geo);

% Regularization
% Tikhonov
[x_hat, lambda_opt] = tikhonov(A, L, AA, LL, interp_signal(:,:), lambda, SNR, order, reg_param_method, compute_params);

% TSVD
%[x_hat, k] = tsvd(A, L, interp_signal(:,5*4000:6*4000), lambda, SNR, order, compute_params);

%% Plotting 3 Maps

% Define times for plotting (seconds)
inst1 = 1;
inst2 = 2;
inst3 = 3;

% Conversion to samples
inst1 = inst1 * fs + 1;
inst2 = inst2 * fs + 1;
inst3 = inst3 * fs + 1;

% Organizing geometry
faces = heart_geo.faces;
x = heart_geo.vertices(:,1);
y = heart_geo.vertices(:,2);
z = heart_geo.vertices(:,3);

% Organizing signal
sinal1 = x_hat(:, inst1);
sinal2 = x_hat(:, inst2);
sinal3 = x_hat(:, inst3);

% Plot
tlo = tiledlayout(1, 3);

f(1) = nexttile(tlo);
trisurf(faces, x, y, z, sinal1, 'facecolor', 'interp', 'LineStyle', 'none');
grid off; axis off;
title(['Instant ', num2str(inst1 / fs), 's']);

f(2) = nexttile(tlo);
trisurf(faces, x, y, z, sinal2, 'facecolor', 'interp', 'LineStyle', 'none');
grid off; axis off;
title(['Instant ', num2str(inst2 / fs), 's']);

f(3) = nexttile(tlo);
trisurf(faces, x, y, z, sinal3, 'facecolor', 'interp', 'LineStyle', 'none');
grid off; axis off;
title(['Instant ', num2str(inst3 / fs), 's']);

set(f, 'Colormap', jet);
c = colorbar(f(end));
c.Layout.Tile = 'east';
c.Label.String = 'Amplitude (\muV)';

%% Plotting 1 Map

% Define time for plotting (seconds)
inst1 = 1;

% Conversion to samples
inst1 = inst1 * fs + 1;

% Organizing geometry
faces = heart_geo.faces;
x = heart_geo.vertices(:,1);
y = heart_geo.vertices(:,2);
z = heart_geo.vertices(:,3);

% Organizing signal
sinal1 = x_hat(:, inst1);

% Plot
figure();
trisurf(faces, x, y, z, sinal1, 'facecolor', 'interp', 'LineStyle', 'none');
grid off; axis off;
title(['Instant ', num2str(inst1 / fs), 's Frame', num2str(inst1)]);
colormap('jet');

% Add colorbar with label
c = colorbar;
c.Label.String = 'Amplitude (\muV)';
