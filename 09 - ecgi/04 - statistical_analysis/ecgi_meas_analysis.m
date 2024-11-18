%% ECGi and MEAs comparison

%% Loading data

% Load the signals
signal_file = load("C:\Users\HeartLAB\Documents\Documents\Conferences\CinC 2024\Version1\Recordings\data_filtered_sync_E14_F3_R4.mat");

% Getting only MEAs signals
meas_signal = signal_file.D_SYNC.EL(1:80,:);

% Load the MEA indices
file_id_meas = load("C:\Users\HeartLAB\Documents\Documents\Conferences\CinC 2024\Version1\ECGi\Dados\projected_signals_exp14.mat");

% Load the heart geometry
heart_geo = file_id_meas.geometry_HR;

% Load the estimated potentials
estimated_file_path = "C:\Users\HeartLAB\Documents\Documents\Codes\Current codes\ECGi\estimation\estimated_E14_F03_R04.mat";
estimated_file = estimated_file.(subsref(fieldnames(estimated_file), substruct('{}', {1})));
estimated_signal = estimated_file.Data;

norm = 'No';

clear estimated_file_path;
%% Normalizing data
% Only if you need to normalize the signals, don't run it as default

meas_signal_raw = signal_file.D_SYNC.EL(1:80,:);
estimated_signal_raw = estimated_file.Data;

meas_signal_norm = normalize(meas_signal_raw,2,'zscore','robust');
estimated_signal_norm = normalize(estimated_signal_raw,2,'zscore','robust');
norm = 'Yes';

% if you want to work with normalized signals from now on
% meas_signal = meas_signal_norm;
% estimated_signal = estimated_signal_norm;

%% Plot before and after normalization

el = 20;
tin = 1.15; % Initial time (s)
tfin = 1.3; % Final time (s)
Fsampling = 4000;
time = linspace(tin, tfin, length(tin*Fsampling:tfin*Fsampling));

figure();
subplot(2,1,1)
plot(time, meas_signal_raw(el,tin*Fsampling+init_est:tfin*Fsampling+init_est),'LineWidth', 1, 'Color','black');
hold on;
plot(time, meas_signal_norm(el,tin*Fsampling+init_est:tfin*Fsampling+init_est),'LineWidth', 1.5, 'Color','red');
title(['MEA - Electrode ', num2str(el)], 'FontSize',15);
ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex','FontSize',14);
legend('Raw','Normalized');

subplot(2,1,2)
plot(time,estimated_signal_raw(id_est_mea,tin*Fsampling:tfin*Fsampling),'LineWidth', 1, 'Color','black');
hold on;
plot(time,estimated_signal_norm(id_est_mea,tin*Fsampling:tfin*Fsampling),'LineWidth', 1.5, 'Color','red');
title(['Estimated - Vertex ', num2str(id_est_mea)], 'FontSize',15);
ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex','FontSize',14);
legend('Raw','Normalized');


%% Calculating the statistics

% Setting the frequency sampling
Fs = 4000;

% Initialize the metrics matrix
metrics = zeros(80, 7);

% Resolution of the heart geometry (number of vertices)
resolution = 'HR';

% Define the estimation time window
init_time = 0.1; % seconds
final_time = 0.3; % seconds

init_sample = init_time * Fs + 1;
final_sample = final_time * Fs;

init_est = estimated_file.Time(1)*Fs;

% Calculate metrics for each electrode
for electrode = 1:80
    if electrode < 33 || electrode > 64
        [mea, id_mea, electrode_mea] = get_electrode_position(electrode, meas_signal, file_id_meas, resolution);

        % Check the estimation time to ensure the same time window will be compared
        id_pot_est = id_mea(electrode_mea);
        pot_estimated = estimated_signal(id_pot_est, init_sample:final_sample); % signal(vertex, time)
        pot_measured = mea(electrode_mea, init_sample+init_est:final_sample+init_est); % mea(electrode, time)

        % Calculate metrics using the new function
        metrics(electrode, :) = calculate_metrics(pot_measured, pot_estimated);
    end
end

%% Find the minimum and average values

% Find the electrode with the minimum RMSE
min_rmse = min(metrics([1:11, 14:32, 65:79], 1)); % (exclude bad electrodes if necessary)
electrode_min_rmse = find(metrics(:, 1) == min_rmse);

% Find the electrode with the minimum RMSE in MEA1
min_rmse_mea1 = min(metrics([1:11, 14:16], 1)); % (exclude bad electrodes if necessary)
electrode_min_rmse_mea1 = find(metrics(:, 1) == min_rmse_mea1);

% Find the electrode with the minimum RMSE in MEA2
min_rmse_mea2 = min(metrics([17:32], 1)); % (exclude bad electrodes if necessary)
electrode_min_rmse_mea2 = find(metrics(:, 1) == min_rmse_mea2);

% Find the electrode with the minimum RMSE in MEA3
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

% Find the electrode with the minimum relative error in MEA1
min_rel_error_mea1 = min(metrics([1:11, 14:16], 7));
electrode_min_rel_error_mea1 = find(metrics(:, 7) == min_rel_error_mea1);

% Find the electrode with the minimum relative error in MEA2
min_rel_error_mea2 = min(metrics([17:32], 7));
electrode_min_rel_error_mea2 = find(metrics(:, 7) == min_rel_error_mea2);

% Find the electrode with the minimum relative error in MEA3
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

%% Correlation analysis

threshold = 0.8; % correlation threshold
count = 0;

% Matrix to store electrodes and correlation values within the threshold
high_corr = [];

for electrode = 1:80
    if metrics(electrode, 6) >= threshold || metrics(electrode, 6) <= -threshold
        count = count + 1;
        high_corr(count, 1) = electrode; % Store electrode number
        high_corr(count, 2) = metrics(electrode, 6); % Store corresponding correlation value
    end
end

fprintf('Number of electrodes with correlation higher %.3f and lower than -%.3f: %d\n', threshold, threshold, count);

% Print the electrodes and their values
% fprintf('Electrodes and correlation values:\n');
% for i = 1:count
%     fprintf('Electrode %d: %.4f\n', high_corr(i, 1), high_corr(i, 2));
% end

% Save the matrix containing electrodes and correlation values
% save('high_corr_values.mat', 'high_corr');

%% Call boxplot function
% metric_type: metric to plot (1=RMSE, 2=MSE, 3=MAE, 4=Standard Deviation of MEAs,
% 5=Standard Deviation of ECGi, 6=Correlation, 7=Relative Error)

metrics_data = {metrics_1; metrics_2};

metric_type = 3; % metric
rhythm_names = {'Sinus','Arrhythmia'};

plotMetricsBoxplot(metrics_data, metric_type, rhythm_names);

%% boxplot for DF

df_boxplot(df_values_meas_tank, df_values_ecgi, meas_signal, file_id_meas, resolution);


%% Selecting only one MEA electrode and its corresponding vertex

% Choose one electrode number
electrode = 10;

% Get MEA electrode position in 3D geometry
[mea, id_mea, electrode_mea] = get_electrode_position(electrode, meas_signal, file_id_meas, 'HR');

% Extract the estimated signal for the corresponding vertex
id_est_mea = id_mea(electrode_mea); % MEA(electrode) vertex index
pot_estimated = estimated_signal(id_est_mea, :); % Estimated signal (vertex, time)

% Extract the MEA signal
pot_measured = mea(electrode_mea, :); % Measured signal (MEA electrode, time)

%% Comparison plot MEA vs ECGi

fs = 4000; % Sampling frequency (Hz)

% Define time window for plotting (in seconds)
t_start = 1;
t_end = 2;

% Convert time window to sample indices
t_start_sample = t_start * fs + 1;
t_end_sample = t_end * fs;

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
plot(time, pot_measured(t_start_sample:t_end_sample), 'LineWidth', 1, 'Color', 'black');
title(['Electrode ', num2str(electrode), ' - ', mea_region], 'FontSize', 15);
ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex', 'FontSize', 14);

% Plot the estimated signal
subplot(2, 1, 2);
plot(time, pot_estimated(t_start_sample:t_end_sample), 'LineWidth', 1, 'Color', 'black');
title(['Correspondent signal - Vertex ', num2str(id_est_mea)], 'FontSize', 15);
ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex', 'FontSize', 14);

% Adjust layout for better visualization
xlabel('Time (s)', 'FontSize', 14);

%% Plotting all meas electrodes and their correspondent vertices
% only one mea per time

mea1 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16];
mea2 = [17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32];
mea3 = [65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80];

% Define the MEA electrodes
electrodes_plot = [mea1];

fs = 4000;

init_est = estimated_file.Time(1);

% Initialize matrices to store positions and signals
vertices_plot = zeros(length(electrodes_plot), 1);
pot_measured_all = cell(length(electrodes_plot), 1);
pot_estimated_all = cell(length(electrodes_plot), 1);

% Get MEA electrode positions and signals
for i = 1:length(electrodes_plot)
    electrode = electrodes_plot(i);
    [mea, id_mea, electrode_mea] = get_electrode_position(electrode, meas_signal, file_id_meas, 'HR');
    
    % Store corresponding vertex and signals
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




