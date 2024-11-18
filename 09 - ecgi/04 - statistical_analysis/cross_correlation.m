%% Cross correlation analysis

%% Loading data

% Load the signals
signal_file = load("C:\Users\HeartLAB\Documents\Documents\Codes\Current codes\02 - extraction_filtering\electric\open_ephys_extraction_matlab_filtering\electric_data_E20_F01_R01_filtered.mat");

% Getting only MEAs signals
meas_signal_raw = signal_file.D_SYNC.EL.Data(1:80,:);

% Load the MEA indices
file_id_meas = load("C:\Users\HeartLAB\Documents\Documents\Projections\projected_signals_exp20_2507\projected_signals_Rec_14_33_05_Bin=8_Cam1.mat");

% Load the heart geometry
heart_geo = file_id_meas.geometry_HR;

% Load the estimated potentials
estimated_file_path = "C:\Users\HeartLAB\Documents\Documents\Codes\Current codes\ECGi\estimated_E20_F01_R01_20240731_185335.mat";
estimated_file = estimated_file.(subsref(fieldnames(estimated_file), substruct('{}', {1})));
estimated_signal = estimated_file.Data;

norm = 'No';

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

%% Choose to use normalized or raw signals

if strcmp(norm, 'Yes')
    meas_signal = meas_signal_norm;
    estimated_signal = estimated_signal_norm;
else
    meas_signal = meas_signal_raw;
    estimated_signal = estimated_signal_raw;
end

%% Cross-correlation calculation

% Setting the frequency sampling
Fs = 4000;

% Initialize the metrics matrix
metrics = zeros(80, 7);

% Resolution of the heart geometry (number of vertices)
resolution = 'HR';

% Define the estimation time window
init_time = 0; % seconds
final_time = 2; % seconds

init_sample = init_time * Fs + 1;
final_sample = final_time * Fs;

init_est = estimated_file.Time(1) * Fs;

% Calculate metrics for each electrode
for electrode = 1:80
    if electrode < 33 || electrode > 64
        [mea, id_mea, electrode_mea] = get_electrode_position(electrode, meas_signal, file_id_meas, resolution);

        % Check the estimation time to ensure the same time window will be compared
        id_pot_est = id_mea(electrode_mea);
        pot_estimated = estimated_signal(id_pot_est, init_sample:final_sample); % signal(vertex, time)
        pot_measured = mea(electrode_mea, init_sample + init_est:final_sample + init_est); % mea(electrode, time)

        % Calculate cross-correlation
        [C, lags] = xcorr(pot_measured, pot_estimated, 'coeff');
        [max_corr_coeff, max_index] = max(C);
        max_lag = lags(max_index);

        % Store the max correlation coefficient and lag in the metrics matrix
        metrics(electrode, 1) = max_corr_coeff;
        metrics(electrode, 2) = max_lag / Fs * 1000; % Convert lag to milliseconds
    end
end

%% Plotting all meas electrodes and their correspondent vertices along with cross correlation

% estimated_signal = x_hat;
% meas_signal = signal_file.D_SYNC.EL(1:80,:);
% estimated_signal = estimated_signal_raw;
% meas_signal = meas_signal_raw;

mea1 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16];
mea2 = [17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32];
mea3 = [65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80];

% Define the MEA electrodes
electrodes_plot = [mea3];

Fs = 4000;

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
% axis square;

% Plot the MEA and estimated signals side by side for the selected electrodes
for i = 1:length(electrodes_plot)
    electrode = electrodes_plot(i);
    
    % Define time window for plotting (in seconds)
    init_time = 1.75;
    final_time = 1.95;
    init_sample = init_time * Fs + 1;
    final_sample = final_time * Fs;
    init_est = estimated_file.Time(1)* Fs;
%     init_est = 2 * Fs;

    time = linspace(init_time, final_time, length(init_sample:final_sample));

    inst = 1.22;
    
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
%     axis square;
    % Plot MEA signal
    subplot(8, 4, 2*i-1);
    plot(time, pot_measured_all{i}(init_sample + init_est : final_sample + init_est), 'LineWidth', 1, 'Color', 'black');
    title(['Electrode ', num2str(electrode), ' - ', mea_region], 'FontSize', 8);
    ylabel('Amp ($\mu$V)', 'Interpreter', 'latex', 'FontSize', 8);
    xlabel('Time (s)', 'FontSize', 8);
%     xline(inst, 'blue', 'LineWidth', 1.5);
%     axis square;
    % Plot estimated signal
    subplot(8, 4, 2*i);
    plot(time, pot_estimated_all{i}(init_sample:final_sample), 'LineWidth', 1, 'Color', 'red');
    title(['Electrode ', num2str(electrode), ' - Vertex ', num2str(vertices_plot(i))], 'FontSize', 8);
    ylabel('Amp ($\mu$V)', 'Interpreter', 'latex', 'FontSize', 8);
    xlabel('Time (s)', 'FontSize', 8);
%     xline(inst, 'blue', 'LineWidth', 1.5);
%     axis square;
    % Calculate and plot cross-correlation for the selected time window
    [C, lags] = xcorr(pot_measured_all{i}(init_sample + init_est : final_sample + init_est), pot_estimated_all{i}(init_sample:final_sample), 'coeff');
    [max_corr_coeff, max_index] = max(C);
    max_lag = lags(max_index) / Fs * 1000; % Convert lag to milliseconds

    % Get subplot position
    subplot_pos = get(gca, 'Position');
    annotation_x = subplot_pos(1) + subplot_pos(3) * 0.81;
    annotation_y = subplot_pos(2) + subplot_pos(4) * 0.1;
    
    % Display max correlation coefficient and lag on the plot
    annotation('textbox', [annotation_x, annotation_y, 0.1, 0.1], 'String', ...
        {['Corr: ', num2str(max_corr_coeff, '%.2f')], ['Lag: ', num2str(max_lag, '%.2f'), ' ms']}, ...
        'FitBoxToText', 'on', 'FontSize', 8, 'Color','blue', 'BackgroundColor', 'none', 'EdgeColor', 'none');
end


%% Plotting all meas electrodes and their correspondent vertices after movemean along with cross correlation

% setting meas electrodes
mea1 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16];
mea2 = [17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32];
mea3 = [65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80];

% Define the MEA electrodes
electrodes_plot = [mea1];

window_size = 3; % Define the window size for smoothing
Fs = 4000;

% Initialize matrices to store positions and signals
vertices_plot = zeros(length(electrodes_plot), 1);
pot_measured_all = cell(length(electrodes_plot), 1);
pot_estimated_smoothed_all = cell(length(electrodes_plot), 1);
cross_corr_coeffs = zeros(length(electrodes_plot), 1);
cross_corr_lags = zeros(length(electrodes_plot), 1);

% Get MEA electrode positions and signals
for i = 1:length(electrodes_plot)
    electrode = electrodes_plot(i);
    [mea, id_mea, electrode_mea] = get_electrode_position(electrode, meas_signal, file_id_meas, 'HR');
    
    % Store corresponding vertex and signals
    vertices_plot(i) = id_mea(electrode_mea);
    pot_measured_all{i} = mea(electrode_mea, :);
    pot_estimated_all{i} = estimated_signal(vertices_plot(i), :);
    pot_estimated_smoothed_all{i} = movmean(estimated_signal(vertices_plot(i), :), window_size); % Apply smoothing
end

% Create a figure
figure;
axis square;

% Plot the MEA and estimated signals side by side for the selected electrodes
for i = 1:length(electrodes_plot)
    electrode = electrodes_plot(i);
    
    % Define time window for plotting (in seconds)
    init_time = 0.52;
    final_time = 0.72;
    init_sample = init_time * Fs + 1;
    final_sample = final_time * Fs;
    init_est = estimated_file.Time(1) * Fs;

    time = linspace(init_time, final_time, length(init_sample:final_sample));

    inst = 1.22;
    
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
    plot(time, pot_measured_all{i}(init_sample + init_est : final_sample + init_est), 'LineWidth', 1, 'Color', 'black');
    title(['Electrode ', num2str(electrode), ' - ', mea_region], 'FontSize', 8);
    ylabel('Amp ($\mu$V)', 'Interpreter', 'latex', 'FontSize', 8);
    xlabel('Time (s)', 'FontSize', 8);
    axis square;

    % Plot smoothed estimated signal
    subplot(8, 4, 2*i);
    plot(time, pot_estimated_smoothed_all{i}(init_sample:final_sample), 'LineWidth', 1, 'Color', 'red');
    title(['Electrode ', num2str(electrode), ' - Vertex ', num2str(vertices_plot(i))], 'FontSize', 8);
    ylabel('Amp ($\mu$V)', 'Interpreter', 'latex', 'FontSize', 8);
    xlabel('Time (s)', 'FontSize', 8);
    axis square;

    % Calculate and plot cross-correlation for the selected time window
    [C, lags] = xcorr(pot_measured_all{i}(init_sample + init_est : final_sample + init_est),pot_estimated_smoothed_all{i}(init_sample:final_sample), 'coeff');
    [max_corr_coeff, max_index] = max(C);
    max_lag = lags(max_index) / Fs * 1000; % Convert lag to milliseconds

    % Store the maximum correlation coefficient and lag
    cross_corr_coeffs(i) = max_corr_coeff;
    cross_corr_lags(i) = max_lag;

    % Calculate the percentage of electrodes with correlation coefficient < 0.7 and >= 0.7
    below_07 = sum(cross_corr_coeffs < 0.7);
    above_07 = sum(cross_corr_coeffs >= 0.7);
    total_electrodes = length(electrodes_plot);

    percentage_below_07 = (below_07 / total_electrodes) * 100;
    percentage_above_or_equal_07 = (above_07 / total_electrodes) * 100;

    % Adjusted position for the correlation annotation
    subplot_pos = get(gca, 'Position');
    annotation_x = subplot_pos(1) + subplot_pos(3) * 0.65;
    annotation_y = subplot_pos(2) + subplot_pos(4) * -0.45;
    
    % Display max correlation coefficient and lag on the plot
    annotation('textbox', [annotation_x, annotation_y, 0.1, 0.1], 'String', ...
        {['Corr: ', num2str(max_corr_coeff, '%.2f')], ['Lag: ', num2str(max_lag, '%.2f'), ' ms']}, ...
        'FitBoxToText', 'on', 'FontSize', 8, 'Color','blue', 'BackgroundColor', 'none', 'EdgeColor', 'none');
end

% Calculate mean and standard deviation of the correlation coefficients
mean_corr_coeff = mean(cross_corr_coeffs);
std_corr_coeff = std(cross_corr_coeffs);

% Display the results
disp('-----------');
disp(['Mean Correlation Coefficient: ', num2str(mean_corr_coeff)]);
disp(['Standard Deviation of Correlation Coefficient: ', num2str(std_corr_coeff)]);
fprintf('Percentage of electrodes with CC < 0.7: %.2f%%\n', percentage_below_07);
fprintf('Percentage of electrodes with CC >= 0.7: %.2f%%\n', percentage_above_or_equal_07);
disp('-----------');



%% Just a few electrodes

% setting meas electrodes
mea1 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16];
mea2 = [17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32];
mea3 = [65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80];

% Define the MEA electrodes
electrodes_plot = [72 77];

window_size = 3;

% Initialize matrices to store positions and signals
vertices_plot = zeros(length(electrodes_plot), 1);
pot_measured_all = cell(length(electrodes_plot), 1);
pot_estimated_smoothed_all = cell(length(electrodes_plot), 1); % For smoothed signals

% Get MEA electrode positions and signals
for i = 1:length(electrodes_plot)
    electrode = electrodes_plot(i);
    [mea, id_mea, electrode_mea] = get_electrode_position(electrode, meas_signal, file_id_meas, 'HR');
    
    % Store corresponding vertex and signals
    vertices_plot(i) = id_mea(electrode_mea);
    pot_measured_all{i} = mea(electrode_mea, :);
    pot_estimated_smoothed_all{i} = movmean(estimated_signal(vertices_plot(i), :), window_size); % Apply smoothing
end

% Create a tiled layout for the plots
num_rows = length(electrodes_plot); % Number of rows for electrodes
num_columns = 2; % Two columns: 1 for measured and 1 for smoothed estimated

% Create a figure and tiled layout
figure;
tiledlayout(num_rows, num_columns, 'Padding', 'compact', 'TileSpacing', 'compact');


% Plot the measured and smoothed estimated signals
for i = 1:length(electrodes_plot)
    electrode = electrodes_plot(i);

    % Define time window for plotting (in seconds)
    init_time = 1.75;
    final_time = 1.9;
    init_sample = init_time * Fs + 1;
    final_sample = final_time * Fs;
    init_est = estimated_file.Time(1) * Fs;

    time = (linspace(init_time, final_time, length(init_sample:final_sample)) - init_time) * 1000;

    % Plot measured signal in the first column
    ax1 = nexttile;
    
    plot(time, pot_measured_all{i}(init_sample + init_est : final_sample + init_est), 'LineWidth', 1.8, 'Color', 'black');
    title(['EGM - Electrode ', num2str(electrode)], 'FontSize', 12);
%     title('EGM','FontSize', 12);
    ylabel('$\mu$V', 'Interpreter', 'latex', 'FontSize', 14);
    xlabel('Time (ms)', 'FontSize', 8);
    grid off;
    axis square;
    axis tight; % Fit axis to the data
    ylim padded;
    set(gca,'FontSize', 14, 'FontWeight', 'bold');

    xlim([0, time(end)]);
 

    % Plot smoothed estimated signal in the second column
    ax2 = nexttile;
   
    plot(time, pot_estimated_smoothed_all{i}(init_sample:final_sample), 'LineWidth', 1.8, 'Color', 'red');
    title(['rEGM - Electrode ', num2str(electrode)], 'FontSize', 12);
%     title('rEGM', 'FontSize', 12);
    ylabel('$\mu$V', 'Interpreter', 'latex', 'FontSize', 14);
    xlabel('Time (ms)', 'FontSize', 8);
    grid off;
    axis square;
    axis tight; % Fit axis to the data
    ylim padded;
    set(gca,'FontSize', 14, 'FontWeight', 'bold');

    xlim([0, time(end)]);

    % Calculate and plot cross-correlation for the selected time window
    [C, lags] = xcorr(pot_measured_all{i}(init_sample + init_est : final_sample + init_est), pot_estimated_smoothed_all{i}(init_sample:final_sample), 'coeff');
    [max_corr_coeff, max_index] = max(C);
    max_lag = lags(max_index) / Fs * 1000; % Convert lag to milliseconds

    % Get subplot position
    subplot_pos = get(gca, 'Position');
    annotation_x = subplot_pos(1) + subplot_pos(3) * 1.1;
    annotation_y = subplot_pos(2) + subplot_pos(4) * 0.2;
    
    % Display max correlation coefficient and lag on the plot
%     annotation('textbox', [annotation_x, annotation_y, 0.1, 0.1], 'String', ...
%         {['Cross Correlation: ', num2str(max_corr_coeff, '%.2f')], ['Lag: ', num2str(max_lag, '%.2f'), ' ms']}, ...
%         'FitBoxToText', 'on', 'FontSize', 13, 'Color','blue', 'BackgroundColor', 'none', 'EdgeColor', 'none',FontWeight='bold');
end
