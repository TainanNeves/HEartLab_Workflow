%% Cross correlation analysis
%
% Author: Angélica Quadros
% Affiliation: HeartLab, UFABC
% Year: 2024
%
% This script performs a cross-correlation analysis between the measured
% signals from MEAs and the estimated potentials 
% obtained through ECGi. It loads the measured
% signal data, the heart geometry, and the estimated signal data. The 
% signals are optionally normalized, and the maximum cross-correlation 
% coefficient and corresponding lag are calculated for each electrode. 
% The results are visualized by plotting the measured and estimated signals 
% side by side, along with the cross-correlation values and lags.
%
% In addition, a smoothing operation using a moving average is applied 
% to the estimated signals before recalculating the cross-correlation and 
% generating plots.
%
%% Loading data

% Load the measured signals
signal_file = load("C:\Users\HeartLAB\Documents\Documents\Codes\Current codes\02 - extraction_filtering\electric\open_ephys_extraction_matlab_filtering\electric_data_E20_F01_R01_filtered.mat");

% Extract the MEA signals from the loaded data
%meas_signal_raw = signal_file.D_SYNC.EL.Data(1:80,:);
meas_signal_raw = signal_file.D_EL.Data(1:80,:);

% Load the MEA projections (electrode positions in the heart geometry)
projections = load("C:\Users\HeartLAB\Documents\Documents\Projections\projected_signals_exp20_2507\projected_signals_Rec_14_33_05_Bin=8_Cam1.mat");

% Extract the heart geometry
heart_geo = projections.geometry_HR;

% Load the estimated potentials
estimated_file_path = "C:\Users\HeartLAB\Documents\Documents\Codes\Current codes\ECGi\estimated_E20_F01_R01_20240809_145729.mat";
estimated_file = load(estimated_file_path);
estimated_signal_raw = estimated_file.estimated.Data;

% Normalization flag
norm = 'No';

clear estimated_file_path 

%% Normalizing data (optional)
% This section is useful if you want to standardize the data before performing further analyses or metrics calculations
%
% Normalizes the data along each row using the Z-score method.
% The robust option is used to make the normalization less sensitive to outliers in the data.
meas_signal_norm = normalize(meas_signal_raw, 2, 'zscore', 'robust');
estimated_signal_norm = normalize(estimated_signal_raw, 2, 'zscore', 'robust');

% Normalization flag
norm = 'Yes';

%% Choose to Use Normalized or Raw Signals
% The signals can be normalized if needed, but the raw signals are used by default.

% If normalization is selected, use the normalized signals
if strcmp(norm, 'Yes')
    meas_signal = meas_signal_norm;
    estimated_signal = estimated_signal_norm;
else
    % Otherwise, use the raw signals
    meas_signal = meas_signal_raw;
    estimated_signal = estimated_signal_raw;
end

%% Cross-Correlation Calculation for MEA and Estimated Signals
%
% This section of the script calculates the cross-correlation between the 
% MEAs and estimated potentials. The cross-correlation is performed for
% each electrode, and the maximum correlation coefficient along with the
% corresponding lag are computed.
%
% The cross-correlation results are stored 
% in the `cc` matrix, where the first column contains the correlation 
% coefficients and the second column contains the corresponding lags (in 
% milliseconds).
%

% Sampling frequency (Hz)
Fs = 4000;

% Initialize matrix for storing cross-correlation coefficients and lags
cc = zeros(80, 2); % 80 electrodes, 2 columns (correlation coefficient, lag)

% Set the resolution of the heart geometry
resolution = 'HR';

% Define the time window for the analysis (seconds)
init_time = 0.01; 
final_time = 0.05; 

% Convert to samples
init_sample = init_time * Fs + 1;
final_sample = final_time * Fs;

% Initial estimation time
init_est = estimated_file.estimated.Time(1) * Fs;

% Loop over each electrode and calculate cross-correlation
for electrode = 1:80
    
    if electrode < 33 || electrode > 64
        % Get the corresponding MEA vertex for the current electrode
        [mea_vertex] = get_electrode_position(electrode, projections, resolution);

        % Extract the measured and estimated signal data for the current electrode
        estimated = estimated_signal(mea_vertex, init_sample:final_sample); % Estimated signal
        measured = meas_signal(electrode, init_sample + init_est:final_sample + init_est); % Measured signal

        % Calculate the cross-correlation between the measured and estimated signals
        [C, lags] = xcorr(measured, estimated, 'coeff');
        
        % Find the maximum correlation coefficient and its corresponding lag
        [max_cc, max_index] = max(C);
        max_lag = lags(max_index); % Lag in samples

        % Store the maximum correlation coefficient and corresponding lag (in ms) in the matrix
        cc(electrode, 1) = max_cc;
        cc(electrode, 2) = max_lag / Fs * 1000; % Convert lag from samples to milliseconds
    end
end

clear C lags electrode init_sample init_time init_est mea_vertex final_sample final_time max_lag max_cc max_index measured estimated

%% Plotting Cross-Correlation Calculation for all MEAs and Estimated Signals After Movemean
% This section calculates and plots the cross-correlation between MEA signals 
% and the smoothed estimated potentials (using a moving mean filter). 
% It computes the mean, standard deviation, and percentage of electrodes 
% with correlation coefficients above and below a defined threshold.
%

% Setting meas electrodes
mea1 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16];
mea2 = [17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32];
mea3 = [65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80];

% Define the MEA electrodes to plot and calculate
electrodes_plot = [mea3];

% Define the window size for smoothing
window_size = 3;

% Sampling frequency (Hz)
Fs = 4000;

% Initialize matrices to store positions and signals
vertices_plot = zeros(length(electrodes_plot), 1);
measured_all = cell(length(electrodes_plot), 1);
estimated_smoothed_all = cell(length(electrodes_plot), 1);
cc = zeros(length(electrodes_plot), 1);
cc_lags = zeros(length(electrodes_plot), 1);

% Initialize matrix for storing cross-correlation coefficients and lags
cc = zeros(80, 2); % 80 electrodes, 2 columns (correlation coefficient, lag)

% Get MEA electrode positions and signals
for i = 1:length(electrodes_plot)
    electrode = electrodes_plot(i);
    [mea_vertex] = get_electrode_position(electrode, projections, resolution);
    
    % Store corresponding vertex and signals
    vertices_plot(i) = mea_vertex;
    measured_all{i} = meas_signal(electrode, :);
    estimated_all{i} = estimated_signal(vertices_plot(i), :);
    estimated_smoothed_all{i} = movmean(estimated_signal(vertices_plot(i), :), window_size); % Apply smoothing
end

% Create a figure
figure;
axis square;

% Plot the MEA and estimated signals side by side for the selected electrodes
for i = 1:length(electrodes_plot)
    electrode = electrodes_plot(i);
    
    % Define time window for plotting (in seconds)
    init_time = 0.01;
    final_time = 0.05;
    init_sample = init_time * Fs + 1;
    final_sample = final_time * Fs;
    init_est = estimated_file.estimated.Time(1) * Fs;

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
    plot(time, measured_all{i}(init_sample + init_est : final_sample + init_est), 'LineWidth', 1, 'Color', 'black');
    title(['Electrode ', num2str(electrode), ' - ', mea_region], 'FontSize', 8);
    ylabel('Amp ($\mu$V)', 'Interpreter', 'latex', 'FontSize', 8);
    xlabel('Time (s)', 'FontSize', 8);
    axis square;

    % Plot smoothed estimated signal
    subplot(8, 4, 2*i);
    plot(time, estimated_smoothed_all{i}(init_sample:final_sample), 'LineWidth', 1, 'Color', 'red');
    title(['Electrode ', num2str(electrode), ' - Vertex ', num2str(vertices_plot(i))], 'FontSize', 8);
    ylabel('Amp ($\mu$V)', 'Interpreter', 'latex', 'FontSize', 8);
    xlabel('Time (s)', 'FontSize', 8);
    axis square;

    % Calculate and plot cross-correlation for the selected time window
    [C, lags] = xcorr(measured_all{i}(init_sample + init_est : final_sample + init_est),estimated_smoothed_all{i}(init_sample:final_sample), 'coeff');
    [max_cc, max_index] = max(C);
    max_lag = lags(max_index) / Fs * 1000; % Convert lag to milliseconds

    % Store the maximum correlation coefficient and corresponding lag (in ms) in the matrix
    cc(electrode, 1) = max_cc;
    cc(electrode, 2) = max_lag / Fs * 1000; % Convert lag from samples to milliseconds

    % Set a threshold
    thresh = 0.7;

    % Calculate the percentage of electrodes with correlation coefficient
    % < threshold and >= threshold
    below_thresh = sum(cc(electrodes_plot,1) < thresh);
    above_thresh = sum(cc(electrodes_plot,1) >= thresh);
    total_electrodes = length(electrodes_plot);

    percentage_low_cc = (below_thresh / total_electrodes) * 100;
    percentage_high_cc = (above_thresh / total_electrodes) * 100;

    % Adjusted position for the correlation annotation
    subplot_pos = get(gca, 'Position');
    annotation_x = subplot_pos(1) + subplot_pos(3) * 0.65;
    annotation_y = subplot_pos(2) + subplot_pos(4) * -0.45;
    
    % Display max correlation coefficient and lag on the plot
    annotation('textbox', [annotation_x, annotation_y, 0.1, 0.1], 'String', ...
        {['Corr: ', num2str(max_cc, '%.2f')], ['Lag: ', num2str(max_lag, '%.2f'), ' ms']}, ...
        'FitBoxToText', 'on', 'FontSize', 8, 'Color','blue', 'BackgroundColor', 'none', 'EdgeColor', 'none');
end

% Calculate mean and standard deviation of the correlation coefficients
mean_cc = mean(cc(electrodes_plot,1));
std_cc = std(cc(electrodes_plot,1));

% Display the results
disp('-----------');
disp(['Mean Correlation Coefficient: ', num2str(mean_cc)]);
disp(['Standard Deviation of Correlation Coefficient: ', num2str(std_cc)]);
fprintf('Percentage of electrodes with CC < 0.7: %.2f%%\n', percentage_low_cc);
fprintf('Percentage of electrodes with CC >= 0.7: %.2f%%\n', percentage_high_cc);
disp('-----------');

clear mea1 mea2 mea3 annotation_x annotation_y inst mea_region time electrode subplot_pos C lags electrode init_sample init_time...
    init_est mea_vertex final_sample final_time max_lag max_cc max_index measured estimated electrodes_plot vertices_plot i cc_lags ...
    total_electrodes window_size

%% Plotting Cross-Correlation Calculation for a few Electrodes and Estimated Signals After Movemean
% This section calculates and plots the cross-correlation between a few MEAs electrodes 
% and the smoothed estimated potentials (using a moving mean filter). 
%

% Define the MEA electrodes
electrodes_plot = [1 2 3 4];

window_size = 3;

% Initialize matrices to store positions and signals
vertices_plot = zeros(length(electrodes_plot), 1);
measured_all = cell(length(electrodes_plot), 1);
estimated_smoothed_all = cell(length(electrodes_plot), 1); % For smoothed signals

% Get MEA electrode positions and signals
for i = 1:length(electrodes_plot)
    electrode = electrodes_plot(i);
    [mea_vertex] = get_electrode_position(electrode, projections, resolution);
    
    % Store corresponding vertex and signals
    vertices_plot(i) = mea_vertex;
    measured_all{i} = meas_signal(electrode, :);
    estimated_smoothed_all{i} = movmean(estimated_signal(vertices_plot(i), :), window_size); % Apply smoothing
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
    init_time = 0.01;
    final_time = 0.05;
    init_sample = init_time * Fs + 1;
    final_sample = final_time * Fs;
    init_est = estimated_file.estimated.Time(1) * Fs;

    time = (linspace(init_time, final_time, length(init_sample:final_sample)) - init_time) * 1000;

    % Plot measured signal in the first column
    ax1 = nexttile;
    
    plot(time, measured_all{i}(init_sample + init_est : final_sample + init_est), 'LineWidth', 1.8, 'Color', 'black');
    title(['EGM - Electrode ', num2str(electrode)], 'FontSize', 10);
%     title('EGM','FontSize', 12);
    ylabel('$\mu$V', 'Interpreter', 'latex', 'FontSize', 10);
    xlabel('Time (ms)', 'FontSize', 8);
    grid off;
    axis square;
    axis tight; % Fit axis to the data
    ylim padded;
    set(gca,'FontSize', 10, 'FontWeight', 'bold');

    xlim([0, time(end)]);
 

    % Plot smoothed estimated signal in the second column
    ax2 = nexttile;
   
    plot(time, estimated_smoothed_all{i}(init_sample:final_sample), 'LineWidth', 1.8, 'Color', 'red');
    title(['rEGM - Electrode ', num2str(electrode)], 'FontSize', 10);
%     title('rEGM', 'FontSize', 12);
    ylabel('$\mu$V', 'Interpreter', 'latex', 'FontSize', 10);
    xlabel('Time (ms)', 'FontSize', 8);
    grid off;
    axis square;
    axis tight; % Fit axis to the data
    ylim padded;
    set(gca,'FontSize', 10, 'FontWeight', 'bold');

    xlim([0, time(end)]);

    % Calculate and plot cross-correlation for the selected time window
    [C, lags] = xcorr(measured_all{i}(init_sample + init_est : final_sample + init_est), estimated_smoothed_all{i}(init_sample:final_sample), 'coeff');
    [max_cc, max_index] = max(C);
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

clear mea1 mea2 mea3 annotation_x annotation_y inst mea_region time electrode subplot_pos C lags electrode init_sample init_time...
    init_est mea_vertex final_sample final_time max_lag max_cc max_index measured estimated electrodes_plot vertices_plot i cc_lags ...
    total_electrodes window_size num_columns num_rows
