%% ECGi and MEAs comparison

%% Loading data

% Load the signals
signal_file = load("C:\Users\HeartLAB\Documents\Documents\Conferences\CinC 2024\Version1\Recordings\data_filtered_sync_E14_F3_R4.mat");

% Getting only MEAs signals
meas_signal = signal_file.D_SYNC.EL(1:80,:);

% Load the MEA indices
file_id_meas = load("C:\Users\HeartLAB\Documents\Documents\Conferences\CinC 2024\Version1\ECGi\Dados\projected_signals_exp14.mat");

% Load the estimated potentials
estimated_file_path = "C:\Users\HeartLAB\Documents\Documents\Codes\Current codes\ECGi\estimation\estimated_E14_F03_R04.mat";
estimated_file = load(estimated_file_path);
estimated_signal = estimated_file.estimated.Data;

%% Calculating the statistics

% Setting the frequency sampling
Fs = 4000;

% Initialize the metrics matrix
metrics = zeros(80, 7);

% Resolution of the heart geometry (number of vertices)
resolution = 20000;

% Define the estimation time window
init_time = 0.1; % seconds
final_time = 0.3; % seconds

init_sample = init_time * Fs + 1;
final_sample = final_time * Fs;

% Calculate metrics for each electrode
for electrode = 1:80
    if electrode < 33 || electrode > 64
        [mea, id_mea, electrode_mea] = get_mea_electrode(electrode, meas_signal, file_id_meas, resolution);

        % Check the estimation time to ensure the same time window will be compared
        id_pot_est = id_mea(electrode_mea);
        pot_estimated = estimated_signal(id_pot_est, init_sample:final_sample); % signal(vertex, time)
        pot_measured = mea(electrode_mea, init_sample:final_sample); % mea(electrode, time)

        % Calculate metrics using the new function
        metrics(electrode, :) = calculate_metrics(pot_measured, pot_estimated);
    end
end

%% Find the minimum and average values

% Find the electrode with the minimum RMSE
min_rmse = min(metrics([1:11, 14:32, 65:79], 1)); % (exclude bad electrodes if necessary)
electrode_min_rmse = find(metrics(:, 1) == min_rmse);

% Find the electrode with the minimum MSE
min_mse = min(metrics([1:32, 65:end], 2));
electrode_min_mse = find(metrics(:, 2) == min_mse);

% Find the electrode with the minimum MAE
min_mae = min(metrics([1:32, 65:end], 3));
electrode_min_mae = find(metrics(:, 3) == min_mae);

% Find the electrode with the minimum relative error
min_rel_error = min(metrics([1:32, 65:end], 7));
electrode_min_rel_error = find(metrics(:, 7) == min_rel_error);

% Calculate average RMSE
mean_rmse = mean(metrics([1:11, 14:16, 17:32, 65:79], 1)); % MEA1, MEA2 and MEA3
mean_rmse_mea1 = mean(metrics([1:11, 14:16], 1)); % MEA1 (exclude bad electrodes if necessary)
mean_rmse_mea2 = mean(metrics([17:32], 1)); % MEA2
mean_rmse_mea3 = mean(metrics([65:79], 1)); % MEA3

fprintf('Smallest RMSE: %.1f in electrode %d.\n', min_rmse, electrode_min_rmse);
fprintf('Smallest MSE: %.1f in electrode %d.\n', min_mse, electrode_min_mse);
fprintf('Smallest MAE: %.1f in electrode %d.\n', min_mae, electrode_min_mae);
fprintf('Smallest relative error: %.1f in electrode %d.\n', min_rel_error, electrode_min_rel_error);
fprintf('Average RMSE: %.4f\n', mean_rmse);
fprintf('Average RMSE for MEA1 (RA): %.4f\n', mean_rmse_mea1);
fprintf('Average RMSE for MEA2 (V): %.4f\n', mean_rmse_mea2);
fprintf('Average RMSE for MEA3 (LA): %.4f\n', mean_rmse_mea3);

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

% Save the updated estimated file
save(estimated_file_path, '-struct', 'estimated_file');

%% Correlation analysis

threshold = 0.1; % correlation threshold
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

metrics_data = {metrics_04s; metrics_1s};

metric_type = 3; % metric
rhythm_names = {'Sinus','Arrhythmia'};

plotMetricsBoxplot(metrics_data, metric_type, rhythm_names);

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
title(['MEA Electrode ', num2str(electrode), ' - ', mea_region], 'FontSize', 15);
ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex', 'FontSize', 14);

% Plot the estimated signal
subplot(2, 1, 2);
plot(time, pot_estimated(t_start_sample:t_end_sample), 'LineWidth', 1, 'Color', 'black');
title(['Vertex ', num2str(id_est_mea)], 'FontSize', 15);
ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex', 'FontSize', 14);

% Adjust layout for better visualization
xlabel('Time (s)', 'FontSize', 14);