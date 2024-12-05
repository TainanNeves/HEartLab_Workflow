%% Moving average comparison
%
% Author: Angélica Quadros
% Affiliation: HeartLab, UFABC
% Year: 2024
%
% Description:
%   This script compares the original and smoothed signals (via moving average)
%   for measured signals (MEA) and estimated signals (ECGi) at selected 
%   electrodes across different regions (RA, V, LA).

% defining meas and estimated signals
meas_signal = meas_signal_raw;
estimated_signal = estimated_signal_raw;

% setting meas electrodes
mea1 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16];
mea2 = [17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32];
mea3 = [65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80];

% Define the MEA electrodes
electrodes_plot = [mea2];

% Setting the frequency sampling
fs = 4000;

% Apply the moving average filter using movmean
window_size = 3;

% Initialize matrices to store positions and signals
vertices_plot = zeros(length(electrodes_plot), 1);
pot_measured_all = cell(length(electrodes_plot), 1);
pot_measured_smoothed_all = cell(length(electrodes_plot), 1); % For smoothed signals
pot_estimated_all = cell(length(electrodes_plot), 1);
pot_estimated_smoothed_all = cell(length(electrodes_plot), 1); % For smoothed signals

% Get MEA electrode positions and signals
for i = 1:length(electrodes_plot)
    electrode = electrodes_plot(i);
    [mea, id_mea, electrode_mea] = get_electrode_position(electrode, meas_signal, file_id_meas, 20000);
    
    % Store corresponding vertex and signals
    vertices_plot(i) = id_mea(electrode_mea);
    pot_measured_all{i} = mea(electrode_mea, :);
    pot_estimated_all{i} = estimated_signal(vertices_plot(i), :);
    pot_estimated_smoothed_all{i} = movmean(estimated_signal(vertices_plot(i), :), window_size); % Apply smoothing
end

% Create a figure
figure;

% Plot the MEA and estimated signals side by side for the selected electrodes
for i = 1:length(electrodes_plot)
    electrode = electrodes_plot(i);
    
    % Define time window for plotting (in seconds)
    t_start = 1;
    t_end = 2;
    t_start_sample = t_start * fs + 1;
    t_end_sample = t_end * fs;
    t_start_sample_mea = t_start_sample + 1*fs;
    t_end_sample_mea = t_end_sample + 1*fs;
    time = linspace(t_start, t_end, length(t_start_sample:t_end_sample));

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
    
    % Plot original and smoothed MEA signal
    subplot(8, 4, 4*i-3);
    plot(time, pot_measured_all{i}(t_start_sample_mea:t_end_sample_mea), 'LineWidth', 1, 'Color', 'black');
    title(['Electrode ', num2str(electrode), ' - Original ', mea_region], 'FontSize', 8);
    ylabel('Amp ($\mu$V)', 'Interpreter', 'latex', 'FontSize', 8);
    xlabel('Time (s)', 'FontSize', 8);
    xline(inst,'blue','LineWidth', 1.5);
    
    subplot(8, 4, 4*i-2);
    plot(time, pot_measured_smoothed_all{i}(t_start_sample_mea:t_end_sample_mea), 'LineWidth', 1, 'Color', 'blue');
    title(['Electrode ', num2str(electrode), ' - Smoothed ', mea_region], 'FontSize', 8);
    ylabel('Amp ($\mu$V)', 'Interpreter', 'latex', 'FontSize', 8);
    xlabel('Time (s)', 'FontSize', 8);
    xline(inst,'blue','LineWidth', 1.5);
    
    % Plot original and smoothed estimated signal
    subplot(8, 4, 4*i-1);
    plot(time, pot_estimated_all{i}(t_start_sample:t_end_sample), 'LineWidth', 1, 'Color', 'red');
    title(['Electrode ', num2str(electrode),' - Original Vertex ', num2str(vertices_plot(i))], 'FontSize', 8);
    ylabel('Amp ($\mu$V)', 'Interpreter', 'latex', 'FontSize', 8);
    xlabel('Time (s)', 'FontSize', 8);
    xline(inst,'blue','LineWidth', 1.5);
    
    subplot(8, 4, 4*i);
    plot(time, pot_estimated_smoothed_all{i}(t_start_sample:t_end_sample), 'LineWidth', 1, 'Color', 'green');
    title(['Electrode ', num2str(electrode),' - Smoothed Vertex ', num2str(vertices_plot(i))], 'FontSize', 8);
    ylabel('Amp ($\mu$V)', 'Interpreter', 'latex', 'FontSize', 8);
    xlabel('Time (s)', 'FontSize', 8);
    xline(inst,'blue','LineWidth', 1.5);
end

%% Reconstructed signals plots original x smoothed

% setting meas electrodes
mea1 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16];
mea2 = [17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32];
mea3 = [65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80];

% Define the MEA electrodes
electrodes_plot = [1,2,3,4];

window_size = 3;

% Initialize matrices to store positions and signals
vertices_plot = zeros(length(electrodes_plot), 1);
pot_estimated_all = cell(length(electrodes_plot), 1);
pot_estimated_smoothed_all = cell(length(electrodes_plot), 1); % For smoothed signals

% Get MEA electrode positions and signals
for i = 1:length(electrodes_plot)
    electrode = electrodes_plot(i);
    [mea, id_mea, electrode_mea] = get_electrode_position(electrode, meas_signal, file_id_meas, 'HR');
    
    % Store corresponding vertex and signals
    vertices_plot(i) = id_mea(electrode_mea);
    pot_estimated_all{i} = estimated_signal(vertices_plot(i), :);
    pot_estimated_smoothed_all{i} = movmean(estimated_signal(vertices_plot(i), :), window_size); % Apply smoothing
end

% Create a tiled layout for the plots
num_plots = length(electrodes_plot);
num_columns = ceil(sqrt(num_plots)); % Number of columns
num_rows = ceil(num_plots / num_columns); % Number of rows

% Create a figure and tiled layout
figure;
tiledlayout(num_rows, num_columns, 'Padding', 'compact', 'TileSpacing', 'compact');

% Plot the original and smoothed estimated signals for the selected electrodes
for i = 1:length(electrodes_plot)
    electrode = electrodes_plot(i);

    % Define the number of samples to plot (adjust as needed)
    start_sample = 1.2*fs+1;
    end_sample = 1.25*fs;

    % Create a tile for each electrode
    ax = nexttile;
    
    % Plot original estimated signal
    h1 = plot(pot_estimated_all{i}(start_sample:end_sample), 'LineWidth', 1, 'Color', 'red');
    hold on;
    
    % Plot smoothed estimated signal
    h2 = plot(pot_estimated_smoothed_all{i}(start_sample:end_sample), 'LineWidth', 1, 'Color', 'green');
    
    % Set plot title and labels
    title(['Electrode ', num2str(electrode), ' - Vertex ', num2str(vertices_plot(i))], 'FontSize', 8);
    ylabel('Amp ($\mu$V)', 'Interpreter', 'latex', 'FontSize', 8);
    xlabel('Samples', 'FontSize', 8);
    grid on;
    
    % Make the plot square
    axis square;
end

% Add a global legend
legend([h1, h2], {'Original', 'Smoothed'}, 'Orientation', 'horizontal', 'Location', 'southoutside');



