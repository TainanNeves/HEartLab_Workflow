% ECGi Reconstruction Script
%
% Author: Angélica Drielly Santos de Quadros
% Affiliation: Federal University of ABC (UFABC), Brazil
% Year: 2025
%
% This script performs electrocardiographic imaging (ECGi) reconstruction
% by estimating epicardial potentials from tank signals using
% inverse problem methods.
%
% Input files:
%   - Preprocessed tank signals (.mat)
%   - Tank electrode indices (.mat)
%   - Tank geometry (.mat)
%   - Heart geometry (.mat)
%   - Transfer matrix (.mat)
%
% Reconstruction methods available (uncomment to use):
%   - Tikhonov regularization (zero, first, or second order)
%   - TSVD (Truncated Singular Value Decomposition)
%   - DSVD (Damped Singular Value Decomposition)
%   - Bayes (requires epicardial reference data)
%
% Output:
%   - Estimated epicardial potentials saved as a .mat struct with fields:
%     Data, Time, Regularization, Order, Heart_geometry, Filtering, Sync
%
% Notes:
%   - Only one reconstruction method should be active at a time.
%   - Tank signals are expected at fs = 200 Hz (downsampled during preprocessing).
%   - The regularization parameter lambda is optimized internally.

%% Loading Files

% Define file paths
tank_signals_file = "../../01 - data/electric_data_Exx_Fxx_Rxx_filtered.mat";
heart_geo_file = "../../01 - data/heart_geometry.mat";
tank_geo_file = "../../01 - data/tank_geometry.mat";
mtransfer_file = "../../01 - data/projections/MTransfer.mat";

%% Reading Files

% Tank signals
tank_signals = load(tank_signals_file);

% Tank geometry
tank_data = load(tank_geo_file);
tank_geo = tank_data.(subsref(fieldnames(tank_data), substruct('{}', {1})));

% Heart geometry
heart_data = load(heart_geo_file);
heart_geo = heart_data.(subsref(fieldnames(heart_data), substruct('{}', {1})));

% Transfer Matrix
mtransfer_data = load(mtransfer_file);
MTransfer = mtransfer_data.(subsref(fieldnames(mtransfer_data), substruct('{}', {1})));

% Clear unnecessary variables
clear heart_geo_file tank_geo_file mtransfer_file tank_signals_file electrodes_idx_file heart_data tank_data mtransfer_data electrodes_data;

%% Convert to mV (if needed)
% The signals are exported in microvolts

% Uncomment if conversion is needed
% signal = 0.00019 * signal;
% signal = signal / 1000;

%% Define Transfer Matrix

A = MTransfer;

%% Define Estimation Parameters

lambda = logspace(-0.5, -12, 10); % it will be optimized inside of the regularization code
order = 0; % Order of regularization (0, 1, or 2)
reg_param_method = 'g'; % global (g) or by instant (i) calculation
compute_params = 1; % 1 if reg_params need to be calculated
model = 'SAF'; % Model type
fs = 200; % Sampling frequency of the tank signals
SNR = 100; % Signal-to-noise ratio; it will only be used in the l-curve; set it to any number different of 40.

%% Precompute matrices for the regularization method

[AA, L, LL] = precompute_matrices(A, order, heart_geo);

%% Estimation time
% Set estimation time window (in seconds)
est_start = 2;
est_end = 3;

% Convert time to samples
est_start_sample = est_start * fs + 1;
est_end_sample = est_end * fs;

%% Tikhonov

[x_hat, lambda_opt,x,z] = tikhonov(A, L, AA, LL, tank_signals(:, est_start_sample:est_end_sample), lambda, SNR, order, reg_param_method, compute_params);


%% TSVD

% [x_hat, k] = tsvd(A, L, tank_signals(:, est_start_sample:est_end_sample), lambda, SNR, order, compute_params);


%% DSVD

% [x_hat, lambda_opt] = dsvd (A, tank_signals(:, est_start_sample:est_end_sample), lambda, SNR, compute_params);

%% Bayes
% 
% ss_random = epicardial_data(:,1:150); %random samples outside of the analysis
% Nsamples = size(ss_random,2);
% 
% %covariance matrix
% Cn = estimate_Cn(tank_signals, 2*fs+1:6*fs); % select time windows outside of the analysis
% 
% x_hat = bayes (A, tank_signals(:, est_start_sample:est_end_sample), epicardial_data, Cn, Nsamples, ss_random);

%% Plot All Signals

figure();
plot(x_hat');
title('ECGi Signals', 'FontSize', 15);
ylabel('Amplitude (\muV)', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Samples', 'FontSize', 15);

%% Plot Vertices in Time

data_est = x_hat; 

% Define vertices to plot
v1 = 1;
v2 = 1000;
v3 = 2900; 
v4 = 7000;

% Define time window for plot (in seconds)
t_start = 0;
t_end = 1;

% Convert time to samples
t_start_sample = t_start * fs + 1;
t_end_sample = t_end * fs;

% Convert to time array
time = linspace(t_start, t_end, length(t_start_sample:t_end_sample));

figure();
subplot(4, 1, 1);
plot(time, data_est(v1, t_start_sample:t_end_sample), 'LineWidth', 1, 'Color', 'black');
title(['Vertex ', num2str(v1)], 'FontSize', 15);
ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex', 'FontSize', 14);

subplot(4, 1, 2);
plot(time, data_est(v2, t_start_sample:t_end_sample), 'LineWidth', 1, 'Color', 'black');
title(['Vertex ', num2str(v2)], 'FontSize', 15);
ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex', 'FontSize', 14);

subplot(4, 1, 3);
plot(time, data_est(v3, t_start_sample:t_end_sample), 'LineWidth', 1, 'Color', 'black');
title(['Vertex ', num2str(v3)], 'FontSize', 15);
ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex', 'FontSize', 14);

subplot(4, 1, 4);
plot(time, data_est(v4, t_start_sample:t_end_sample), 'LineWidth', 1, 'Color', 'black');
title(['Vertex ', num2str(v4)], 'FontSize', 15);
ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Time (s)', 'FontSize', 15);

%% Upsampling

estimated_down = x_hat;

% Downsampling factor
factor = 20;

% Estimated signals upsampling
estimated_up = resample(estimated_down', factor, 1)'; % Upsample estimated signals

% Update the sampling frequency
fs_up = fs * factor;

%% Normalizing signals

row_max = max(abs(estimated_up), [], 2);
estimated_signal = bsxfun(@rdivide, estimated_up, row_max);

%% Saving

% Define the regularization method used
reg_method = 'Tikhonov'; % set the regularization method used

% Filename to save
FileName = 'EXX_FXX_RXX';

% Get the current directory
current_dir = pwd;

% Determine the parent directory
[parent_dir, ~, ~] = fileparts(current_dir);

% Create a timestamp
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

% Define the path to save the file in the parent directory
output_file = fullfile(parent_dir, ['estimated_', FileName, '_', timestamp, '.mat']);

% Prepare the data structure
estimated = struct(); % Initialize structure
estimated.Data = x_hat; 
estimated.Time = [est_start, est_end];
estimated.Regularization = reg_method;
estimated.Order = order;
estimated.Heart_geometry = heart_geo;
estimated.Filtering = ''; %set applied filtering
estimated.Sync = 'Yes'; %set if it is synchronized with optical and electrical
estimated.Data_upsampled = estimated_up;
estimated.Data_norm = estimated_signal;

% Save the file to the specified path
save(output_file, 'estimated');

fprintf('File saved to: %s\n', output_file);
