% post processing ecgi

%% moving average

% Define the electrodes for each MEA
mea1 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]; % RA
mea2 = [17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32]; % V
mea3 = [65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80]; % LA


% Define the MEA electrodes to plot (only one MEA per time)
electrodes_plot = [mea1];

% fs = 4000; % Sampling frequency (Hz)
fs = 200;


% init_est = estimated_file.Time(1) * fs; % Starting time for the estimated signal
init_est = 2 * fs;

% Initialize matrices to store positions and signals
vertices_plot = zeros(length(electrodes_plot), 1);
pot_measured_all = cell(length(electrodes_plot), 1);
estimated_all = cell(length(electrodes_plot), 1);

resolution = 'HR';

% Define the window size for smoothing
window_size = 10;

% Loop through each selected electrode
for i = 1:length(electrodes_plot)
    electrode = electrodes_plot(i);

    % Get the MEA electrode position and the corresponding vertex in the heart geometry
    mea_vertex = get_electrode_position(electrode, projections, resolution);

    % Store the vertex index and signal data for both measured and estimated signals
    if electrodes_plot == mea3
        electrode = electrode - 32;
    end

    vertices_plot(i) = mea_vertex;
    pot_measured_all{i} = meas_signal(electrode, :);
    estimated_all{i} = estimated_signal(vertices_plot(i), :);
    estimated_smoothed_all{i} = movmean(estimated_signal(vertices_plot(i), :), window_size); % Apply smoothing

end

% Create a figure
figure;

% Plot the MEA and estimated signals side by side for the selected electrodes
for i = 1:length(electrodes_plot)
    electrode = electrodes_plot(i);
    
    % Define time window for plotting (in seconds)
    t_start = 0;
    t_end = 2;
    t_start_sample = t_start * fs + 1;
    t_end_sample = t_end * fs;
    
    time = linspace(t_start, t_end, length(t_start_sample:t_end_sample));

    inst = 0.282; % Instant of interest

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

    % Plot estimated signal
    subplot(8, 4, 2*i);
    plot(time, estimated_smoothed_all{i}(t_start_sample:t_end_sample), 'LineWidth', 1, 'Color', 'red');
    title(['Reconstructed - Electrode ', num2str(electrode)], 'FontSize', 8);
    ylabel('Amp ($\mu$V)', 'Interpreter', 'latex', 'FontSize', 8);
    xlabel('Time (s)', 'FontSize', 8);
    xline(1.2);
end


%% Median Filter

% Define the electrodes for each MEA
mea1 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]; % RA
mea2 = [17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32]; % V
mea3 = [65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80]; % LA

% Define the MEA electrodes to plot (only one MEA per time)
electrodes_plot = [mea1];

% fs = 4000; % Sampling frequency (Hz)
fs = 200;

% init_est = estimated_file.Time(1) * fs; % Starting time for the estimated signal
init_est = 2 * fs;

% Initialize matrices to store positions and signals
vertices_plot = zeros(length(electrodes_plot), 1);
pot_measured_all = cell(length(electrodes_plot), 1);
estimated_all = cell(length(electrodes_plot), 1);

resolution = 'HR';

% Define the window size for smoothing
window_size = 10;

% Loop through each selected electrode
for i = 1:length(electrodes_plot)
    electrode = electrodes_plot(i);

    % Get the MEA electrode position and the corresponding vertex in the heart geometry
    mea_vertex = get_electrode_position(electrode, projections, resolution);

    % Store the vertex index and signal data for both measured and estimated signals
    if electrodes_plot == mea3
        electrode = electrode - 32;
    end

    vertices_plot(i) = mea_vertex;
    pot_measured_all{i} = meas_signal(electrode, :);
    estimated_all{i} = estimated_signal(vertices_plot(i), :);
    estimated_filtered_all{i} = medfilt1(estimated_signal(vertices_plot(i), :), window_size); % Apply median filter
end

% Create a figure
figure;

% Plot the MEA and estimated signals side by side for the selected electrodes
for i = 1:length(electrodes_plot)
    electrode = electrodes_plot(i);
    
    % Define time window for plotting (in seconds)
    t_start = 0;
    t_end = 2;
    t_start_sample = t_start * fs + 1;
    t_end_sample = t_end * fs;
    
    time = linspace(t_start, t_end, length(t_start_sample:t_end_sample));

    inst = 0.282; % Instant of interest

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

    % Plot estimated signal
    subplot(8, 4, 2*i);
    plot(time, estimated_filtered_all{i}(t_start_sample:t_end_sample), 'LineWidth', 1, 'Color', 'red');
    title(['Reconstructed - Electrode ', num2str(electrode)], 'FontSize', 8);
    ylabel('Amp ($\mu$V)', 'Interpreter', 'latex', 'FontSize', 8);
    xlabel('Time (s)', 'FontSize', 8);
    xline(1.2);
end
