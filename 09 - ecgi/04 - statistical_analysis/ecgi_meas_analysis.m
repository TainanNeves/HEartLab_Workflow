%% ECGi and MEAs comparison
%
% Author: Angélica Quadros
% Affiliation: HeartLab, UFABC
% Year: 2024
%
% This script performs a comprehensive analysis of ECGi and MEA signals to evaluate the quality of non-invasive reconstructions. 
% The workflow includes:
% - Loading measured and estimated signals.
% - Extracting MEA potentials and optionally normalizing them.
% - Computing key metrics such as RMSE, MSE, MAE, and correlation coefficients to assess reconstruction accuracy.
% - Analyzing correlation coefficients to identify electrodes with strong relationships between MEA and ECGi signals.
% - Displaying statistical summaries of metrics, including minimum and average values.
% - Generating boxplots to compare metrics across potentials and dominant frequency values.
% - Visualizing MEA signals alongside corresponding vertices on a 3D heart geometry.
%
%
%% Loading data

% Load the filtered signal data
signal_file = load("C:\Users\HeartLAB\Documents\Documents\Conferences\CinC 2024\Version2\Dados\data_filtered_sync_E20_F01_R01.mat");
% signal_file = load("C:\Users\HeartLAB\Documents\Documents\Codes\Current codes\extracted_filtered_recordings\data_filtered_sync_E14_F4_R10.mat"); % Synchronized data
% signal_file = load("C:\Users\HeartLAB\Documents\Documents\Conferences\CinC 2024\Version1\Recordings\data_filtered_sync_E18_F2_R2.mat");


% Extracting only the MEA signals from the loaded data
% meas_signal_raw = signal_file.D_SYNC.EL(1:80,:);
meas_signal = signal_file.D_SYNC.EL(1:80,:);

% Load the MEA indices (which map the electrodes to vertices in the heart geometry)
% projections = load("C:\Users\HeartLAB\Documents\Documents\Projections\EXP14\projected_signals_Rec_15_43_10_Bin=8\projected_signals_Rec_15_43_10_Bin=8.mat");
% projections = load("C:\Users\HeartLAB\Documents\Documents\Projections\projected_signals_exp18.mat");
projections = load("C:\Users\HeartLAB\Documents\Documents\Projections\projected_signals_exp20_2507\projected_signals_Rec_14_33_05_Bin=8_Cam1.mat");

% Extract the heart geometry
heart_geo = projections.geometry_HR;

% Load the estimated potentials
estimated_file_path = "C:\Users\HeartLAB\Documents\Documents\Codes\Current codes\ECGi\estimated_E20_F01_R01_20240809_145729.mat";
estimated_file = load(estimated_file_path);
estimated_signal = estimated_file.estimated.Data;
% load("C:\Users\HeartLAB\Documents\Documents\Codes\Current codes\ECGi\estimated_E14_F04_R10_20240819_154102.mat");
% estimated_file = load("C:\Users\HeartLAB\Documents\Documents\Codes\Current codes\ECGi\estimated_E18_F02_R02_20240814_123749.mat");
% estimated_signal_raw = estimated.Data;

% Flag for normalization
norm = 'No';

%% Normalizing data (optional)
% This section is useful if you want to standardize the data before performing further analyses or metrics calculations
%
% Normalizes the data along each row using the Z-score method.
% The robust option is used to make the normalization less sensitive to outliers in the data.
meas_signal_norm = normalize(meas_signal_raw,2,'zscore','robust');
estimated_signal_norm = normalize(estimated_signal_raw,2,'zscore','robust');
norm = 'Yes'; % Flag for normalization

%% Plot before and after normalization
% This section of the code plots the raw and normalized signals for both the MEA and estimated signals.

el = 20;  % Define the MEA electrode number
tin = 1.15; % Initial time (s)
tfin = 1.3; % Final time (s)
Fsampling = 4000;  % Sampling frequency in Hz
time = linspace(tin, tfin, length(tin*Fsampling:tfin*Fsampling));

% Create a new figure
figure();

% Plot the MEA electrode signal before and after normalization
subplot(2,1,1)
plot(time, meas_signal_raw(el,tin*Fsampling+init_est:tfin*Fsampling+init_est),'LineWidth', 1, 'Color','black');
hold on;
plot(time, meas_signal_norm(el,tin*Fsampling+init_est:tfin*Fsampling+init_est),'LineWidth', 1.5, 'Color','red');
title(['MEA - Electrode ', num2str(el)], 'FontSize',15);
ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex','FontSize',14);
legend('Raw','Normalized');

% Plot the estimated signal before and after normalization
subplot(2,1,2)
plot(time,estimated_signal_raw(id_est_mea,tin*Fsampling:tfin*Fsampling),'LineWidth', 1, 'Color','black');
hold on;
plot(time,estimated_signal_norm(id_est_mea,tin*Fsampling:tfin*Fsampling),'LineWidth', 1.5, 'Color','red');
title(['Estimated - Vertex ', num2str(id_est_mea)], 'FontSize',15);
ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex','FontSize',14);
legend('Raw','Normalized');

%% Calculating the statistics
% This section calculates various performance metrics to compare the measured
% (from MEAs) and estimated (ECGi) signals.
% These metrics include Root Mean Squared Error (RMSE), Mean Squared Error (MSE),
% Mean Absolute Error (MAE), Standard Deviation, Correlation, and Relative Error.
% The metrics are calculated within a defined time window and stored in the 'metrics' matrix.

% metric_type: metric to plot (1=RMSE, 2=MSE, 3=MAE, 4=Standard Deviation of MEAs,
% 5=Standard Deviation of ECGi, 6=Correlation, 7=Relative Error)

% meas_signal = meas_signal_norm;
% estimated_signal = estimated_signal_norm;

% meas_signal = meas_signal_raw;
% estimated_signal = estimated_signal_raw;
% estimated_signal = estimated.Data;

% Setting the frequency sampling
Fs = 4000;

% Initialize the metrics matrix
metrics = zeros(80, 7);

% Resolution of the heart geometry (number of vertices)
resolution = 20000;

% Define the estimation time window
init_time = 0.01; % seconds
final_time = 0.05; % seconds

init_sample = init_time * Fs + 1;
final_sample = final_time * Fs;

% init_est = estimated.Time(1)*Fs;
init_est = 0;

% Calculate metrics for each electrode
for electrode = 1:80
    if electrode < 33 || electrode > 64
        mea_vertex = get_electrode_position(electrode, projections, resolution);

        % Check the estimation time to ensure the same time window will be compared
        pot_estimated = estimated_signal(mea_vertex, init_sample:final_sample); % signal(vertex, time)
        pot_measured = meas_signal(electrode, init_sample+init_est:final_sample+init_est); % signal(electrode, time)

        % Calculate metrics using the new function
        metrics(electrode, :) = calculate_metrics(pot_measured, pot_estimated);
    end
end

%% Find the minimum and average values
% This section finds the minimum and average values of the metrics
% (RMSE, MSE, MAE, Relative Error) for each electrode and for the regions
% MEA1, MEA2, and MEA3. It displays the electrode with the smallest values
% for these metrics as well as the average values for the different regions.

% Find the electrode with the minimum RMSE across all regions
min_rmse = min(metrics([1:11, 14:32, 65:79], 1)); % (exclude bad electrodes if necessary)
electrode_min_rmse = find(metrics(:, 1) == min_rmse);

% Find the electrode with the minimum RMSE in MEA1 (RA)
min_rmse_mea1 = min(metrics([1:11, 14:16], 1)); % (exclude bad electrodes if necessary)
electrode_min_rmse_mea1 = find(metrics(:, 1) == min_rmse_mea1);

% Find the electrode with the minimum RMSE in MEA2 (V)
min_rmse_mea2 = min(metrics([17:32], 1)); % (exclude bad electrodes if necessary)
electrode_min_rmse_mea2 = find(metrics(:, 1) == min_rmse_mea2);

% Find the electrode with the minimum RMSE in MEA3 (LA)
min_rmse_mea3 = min(metrics([65:79], 1)); % (exclude bad electrodes if necessary)
electrode_min_rmse_mea3 = find(metrics(:, 1) == min_rmse_mea3);

% Find the electrode with the minimum MSE
min_mse = min(metrics([1:32, 65:end], 2));
electrode_min_mse = find(metrics(:, 2) == min_mse);

% Find the electrode with the minimum MAE
min_mae = min(metrics([1:32, 65:end], 3));
electrode_min_mae = find(metrics(:, 3) == min_mae);

% Find the electrode with the minimum relative error
min_rel_error = min(metrics([1:32, 65:79], 7));
electrode_min_rel_error = find(metrics(:, 7) == min_rel_error);

% Find the electrode with the minimum relative error in MEA1 (RA)
min_rel_error_mea1 = min(metrics([1:11, 14:16], 7));
electrode_min_rel_error_mea1 = find(metrics(:, 7) == min_rel_error_mea1);

% Find the electrode with the minimum relative error in MEA2 (V)
min_rel_error_mea2 = min(metrics([17:32], 7));
electrode_min_rel_error_mea2 = find(metrics(:, 7) == min_rel_error_mea2);

% Find the electrode with the minimum relative error in MEA3 (LA)
min_rel_error_mea3 = min(metrics([65:79], 7));
electrode_min_rel_error_mea3 = find(metrics(:, 7) == min_rel_error_mea3);

% Calculate average RMSE
mean_rmse = mean(metrics([1:11, 14:16, 17:32, 65:79], 1)); % MEA1, MEA2 and MEA3
mean_rmse_mea1 = mean(metrics([1:11, 14:16], 1)); % MEA1 (exclude bad electrodes if necessary)
mean_rmse_mea2 = mean(metrics([17:32], 1)); % MEA2
mean_rmse_mea3 = mean(metrics([65:79], 1)); % MEA3

% Calculate average mean relative error
mean_mre = mean(metrics([1:11, 14:16, 17:32, 65:79], 7)); % MEA1, MEA2 and MEA3
mean_mre_mea1 = mean(metrics([1:11, 14:16], 7)); % MEA1 (exclude bad electrodes if necessary)
mean_mre_mea2 = mean(metrics([17:32], 7)); % MEA2
mean_mre_mea3 = mean(metrics([65:79], 7)); % MEA3

% Display the results for RMSE, MRE, and their averages
% fprintf('Smallest MRE in MEA1: %.1f in electrode %d.\n', min_rel_error_mea1, electrode_min_rel_error_mea1);
% fprintf('Smallest MRE: %.1f in electrode %d.\n', min_rel_error_mea2, electrode_min_rel_error_mea2);
% fprintf('Smallest MRE in MEA3: %.1f in electrode %d.\n', min_rel_error_mea3, electrode_min_rel_error_mea3);
% fprintf('Average MRE: %.4f\n', mean_mre);
% fprintf('Average MRE for MEA1 (RA): %.4f\n', mean_mre_mea1);
% fprintf('Average MRE for MEA2 (RA): %.4f\n', mean_mre_mea2);
% fprintf('Average MRE for MEA3 (RA): %.4f\n', mean_mre_mea3);
% fprintf('Smallest RMSE: %.1f in electrode %d.\n', min_rmse, electrode_min_rmse);
fprintf('Smallest RMSE in MEA1: %.1f in electrode %d.\n', min_rmse_mea1, electrode_min_rmse_mea1);
fprintf('Smallest RMSE in MEA2: %.1f in electrode %d.\n', min_rmse_mea2, electrode_min_rmse_mea2);
fprintf('Smallest RMSE in MEA3: %.1f in electrode %d.\n', min_rmse_mea3, electrode_min_rmse_mea3);
%fprintf('Smallest MSE: %.1f in electrode %d.\n', min_mse, electrode_min_mse);
%fprintf('Smallest MAE: %.1f in electrode %d.\n', min_mae, electrode_min_mae);
% fprintf('Smallest mean relative error: %.1f in electrode %d.\n', min_rel_error, electrode_min_rel_error);
fprintf('Average RMSE: %.4f\n', mean_rmse);
fprintf('Average RMSE for MEA1 (RA): %.4f\n', mean_rmse_mea1);
fprintf('Average RMSE for MEA2 (V): %.4f\n', mean_rmse_mea2);
fprintf('Average RMSE for MEA3 (LA): %.4f\n -----\n', mean_rmse_mea3);

%% Saving

% Create a unique identifier using the init and final time
metrics_field_name = sprintf('Metrics_%s_%s', strrep(num2str(init_time), '.', '_'), strrep(num2str(final_time), '.', '_'));

% Check if Metrics field already exists, if not, create it
if ~isfield(estimated_file.estimated, 'Metrics')
    estimated_file.estimated.Metrics = struct();
end

% Append the new metrics to the existing Metrics field
estimated_file.estimated.Metrics.(metrics_field_name).Data = metrics;
estimated_file.estimated.Metrics.(metrics_field_name).Time = [init_time, final_time];
estimated_file.estimated.Metrics.(metrics_field_name).Columns = {'RMSE', 'MSE', 'MAE', 'Standard Deviation of MEAs', 'Standard Deviation of ECGi', 'Correlation', 'Relative Error'};
estimated_file.estimated.Metrics.(metrics_field_name).Normalization = norm;

% Save the updated estimated file
save(estimated_file_path, '-struct', 'estimated_file');

%% Correlation coefficient analysis

% Define the correlation threshold
threshold = 0.8; % Adjust as needed

% Identify electrodes with positive correlation above the threshold
positive_corr_el = find(metrics(:, 6) > threshold);
positive_corr_values = metrics(positive_corr_el, 6);

% Identify electrodes with negative correlation below the negative threshold
negative_corr_el = find(metrics(:, 6) < -threshold);
negative_corr_values = metrics(negative_corr_el, 6);

% Calculate the average correlation coefficient across all electrodes
average_correlation = mean(metrics(:, 6));

% Display results
fprintf('Number of electrodes with correlation above %.1f: %d\n', threshold, numel(positive_corr_el));
fprintf('Electrodes with correlation above %.1f: %s\n', threshold, num2str(positive_corr_el'));
fprintf('Number of electrodes with correlation below %.1f: %d\n', threshold, numel(negative_corr_el));
fprintf('Electrodes with correlation below %.1f: %s\n', threshold, num2str(negative_corr_el'));
fprintf('Average correlation coefficient: %.4f\n', average_correlation);

% Create a matrix containing electrode numbers and correlation values
% positive_corr_matrix = [positive_corr_el, positive_corr_values];
% negative_corr_matrix = [negative_corr_el, negative_corr_values];

% Save the matrices separately
% save('positive_corr_values.mat', 'positive_corr_matrix');
% save('negative_corr_values.mat', 'negative_corr_matrix');
% 
% Save all in one matrix 
% save('correlation_values.mat', 'positive_corr_matrix', 'negative_corr_matrix');

%% boxplot for metrics
% This section generates a boxplot to compare the selected metric across different conditions.

% metrics_data: cell array containing matrices of metrics, where each row corresponds to an electrode and each column corresponds to a metric. 
metrics_data = {metrics_1, metrics_2};

% metric_type: index of the metric to plot. The values range from 1 to 7, each corresponding to a different metric.
metric_type = 6; % metric

% rhythm: A cell array containing labels for the different conditions being compared.
% These labels could represent different types of normalization, regularization methods, or any other experimental characteristic.
rhythm = {'Standard norm','Robust norm'};

% call the plotMetricsBoxplot function
plotMetricsBoxplot(metrics_data, metric_type, rhythm);

%% boxplot for DF
% Generates a boxplot comparing DF values between MEAs and ECGi signals.

df_boxplot(df_values_meas_tank, df_values_ecgi, meas_signal, projections, 'HR');


%% Selecting only one MEA electrode and its corresponding vertex
% This section extracts the measured and estimated signals for a selected MEA electrode. 
% It retrieves the position of the chosen electrode in the 3D geometry and extracts 
% the corresponding vertex signal from the ECGi estimation, as well as the measured 
% signal from the corresponding MEA electrode.

% meas_signal = meas_signal_raw;
% estimated_signal = estimated_signal_raw;

% Choose one electrode number
electrode = 21;

% Get MEA electrode position in 3D geometry
[mea, id_mea, electrode_mea] = get_electrode_position(electrode, meas_signal, projections, 'HR');

% Extract the estimated signal for the corresponding vertex
id_est_mea = id_mea(electrode_mea); % MEA(electrode) vertex index
pot_estimated = estimated_signal(id_est_mea, :); % Estimated signal (vertex, time)

% Extract the MEA signal
pot_measured = mea(electrode_mea, :); % Measured signal (MEA electrode, time)

%% Comparison plot MEA vs ECGi
% This section plots the measured signals from one selected MEA electrode with the 
% estimated signal at the corresponding vertex of the heart geometry.

fs = 4000; % Sampling frequency (Hz)

% Define time window for plotting (in seconds)
t_start = 0;
t_end = 2;

% Convert time window to sample indices
t_start_sample = t_start * fs + 1;
t_end_sample = t_end * fs;

% Get the start time of the estimated signal
init_est = estimated_file.Time(1);

% Generate time array for plotting
time = linspace(t_start, t_end, length(t_start_sample:t_end_sample));

% Determine the MEA region based on electrode number
if ismember(electrode, 1:16)
    mea_region = 'Right Atrium (MEA1)';
elseif ismember(electrode, 17:32)
    mea_region = 'Ventricle (MEA2)';
elseif ismember(electrode, 65:80)
    mea_region = 'Left Atrium (MEA3)';
else
    mea_region = 'Unknown Region'; % In case of an invalid electrode number
end

% Create a new figure for plotting
figure();

% Plot the measured signal
subplot(2, 1, 1);
plot(time, pot_measured(t_start_sample+init_est:t_end_sample+init_est), 'LineWidth', 1, 'Color', 'black');
title(['MEA Electrode ', num2str(electrode), ' - ', mea_region], 'FontSize', 15);
ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex', 'FontSize', 14);

% Plot the estimated signal
subplot(2, 1, 2);
plot(time, pot_estimated(t_start_sample:t_end_sample), 'LineWidth', 1, 'Color', 'black');
title(['Vertex ', num2str(id_est_mea)], 'FontSize', 15);
ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex', 'FontSize', 14);

% Adjust layout for better visualization
xlabel('Time (s)', 'FontSize', 14);

%% Plotting all meas electrodes and their correspondent vertices
% This section generates a series of plots comparing the measured signal from multiple MEA electrodes 
% with the estimated signals at the corresponding vertices of the heart geometry.


% Define the electrodes for each MEA
mea1 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]; % RA
mea2 = [17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32]; % V
mea3 = [65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80]; % LA

% Define the MEA electrodes to plot (only one MEA per time)
electrodes_plot = [mea1];

fs = 4000; % Sampling frequency (Hz)

init_est = estimated_file.Time(1); % Starting time for the estimated signal

% Initialize matrices to store positions and signals
vertices_plot = zeros(length(electrodes_plot), 1);
pot_measured_all = cell(length(electrodes_plot), 1);
pot_estimated_all = cell(length(electrodes_plot), 1);

% Loop through each selected electrode
for i = 1:length(electrodes_plot)
    electrode = electrodes_plot(i);

    % Get the MEA electrode position and the corresponding vertex heart
    % geometry
    [mea, id_mea, electrode_mea] = get_electrode_position(electrode, meas_signal, projections, 'HR');
    
    % Store the vertex index and signal data for both measured and estimated signals
    vertices_plot(i) = id_mea(electrode_mea);
    pot_measured_all{i} = mea(electrode_mea, :);
    pot_estimated_all{i} = estimated_signal(vertices_plot(i), :);
end

% Create a figure
figure;

% Plot the MEA and estimated signals side by side for the selected electrodes
for i = 1:length(electrodes_plot)
    electrode = electrodes_plot(i);
    
    % Define time window for plotting (in seconds)
    t_start = 0;
    t_end = 0.8;
    t_start_sample = t_start * fs + 1;
    t_end_sample = t_end * fs;
    
    time = linspace(t_start, t_end, length(t_start_sample:t_end_sample));

    inst = 0.282;
    
    % Determine the MEA region based on electrode number
    if ismember(electrode, 1:16)
        mea_region = 'Right Atrium (MEA1)';
    elseif ismember(electrode, 17:32)
        mea_region = 'Ventricle (MEA2)';
    elseif ismember(electrode, 65:80)
        mea_region = 'Left Atrium (MEA3)';
    else
        mea_region = 'Unknown Region';
    end
    
    % Plot MEA signal
    subplot(8, 4, 2*i-1);
    plot(time, pot_measured_all{i}(t_start_sample+init_est:t_end_sample+init_est), 'LineWidth', 1, 'Color', 'black');
    title(['Electrode ', num2str(electrode)], 'FontSize', 8);
    ylabel('Amp ($\mu$V)', 'Interpreter', 'latex', 'FontSize', 8);
    xlabel('Time (s)', 'FontSize', 8);
%     xline(inst,'blue','LineWidth', 1.5);
    
    % Plot estimated signal
    subplot(8, 4, 2*i);
    plot(time, pot_estimated_all{i}(t_start_sample:t_end_sample), 'LineWidth', 1, 'Color', 'red');
    title(['Reconstructed  - Electrode ', num2str(electrode)], 'FontSize', 8);
%     title(['Reconstructed  - Electrode ', num2str(electrode),' - Vertex ', num2str(vertices_plot(i))], 'FontSize', 8);
    ylabel('Amp ($\mu$V)', 'Interpreter', 'latex', 'FontSize', 8);
    xlabel('Time (s)', 'FontSize', 8);
%     xline(inst,'blue','LineWidth', 1.5);
end
