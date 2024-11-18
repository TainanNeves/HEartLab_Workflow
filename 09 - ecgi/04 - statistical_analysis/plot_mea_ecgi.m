% Define a vector of electrode numbers to plot
electrodes = [18 19 20 21]; % Example vector of electrode numbers from 1 to 10

% Define the sampling frequency (Hz)
fs = 4000;

% Define time window for plotting (in seconds)
t_start = 0;
t_end = 2;

% Convert time window to sample indices
t_start_sample = t_start * fs + 1;
t_end_sample = t_end * fs;

init_est = estimated_file.Time(1);

% Generate time array for plotting
time = linspace(t_start, t_end, length(t_start_sample:t_end_sample));

% Create a new figure for plotting
figure();

% Calculate the number of rows needed (each electrode takes 2 rows)
n_rows = ceil(length(electrodes) / 2) * 2;

% Loop through each electrode
for i = 1:length(electrodes)
    electrode = electrodes(i);
    
    % Get MEA electrode position in 3D geometry
    [mea, id_mea, electrode_mea] = get_electrode_position(electrode, meas_signal, file_id_meas, 10000);
    
    % Extract the estimated signal for the corresponding vertex
    id_est_mea = id_mea(electrode_mea); % MEA(electrode) vertex index
    pot_estimated = estimated_signal(id_est_mea, :); % Estimated signal (vertex, time)
    
    % Extract the MEA signal
    pot_measured = mea(electrode_mea, :); % Measured signal (MEA electrode, time)
    
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
    
    % Determine the column (1 or 2)
    column = mod(i-1, 2) + 1;
    
    % Determine the row (each electrode takes up 2 rows)
    row = floor((i-1) / 2) * 2 + 1;
    
    % Plot the measured signal (top subplot for this electrode)
    subplot(n_rows, 2, row*2 - (2-column));
    plot(time, pot_measured(t_start_sample+init_est:t_end_sample+init_est), 'LineWidth', 1, 'Color', 'black');
    title(['Electrode ', num2str(electrode), ' - ', mea_region], 'FontSize', 15);
    ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex');
    
    % Plot the estimated signal directly below the measured signal
    subplot(n_rows, 2, (row+1)*2 - (2-column));
    plot(time, pot_estimated(t_start_sample:t_end_sample), 'LineWidth', 1, 'Color', 'red');
    title(['Estimated Signal - Vertex ', num2str(id_est_mea)], 'FontSize', 15);
    ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex');
end

% Adjust layout for better visualization
xlabel('Time (s)', 'FontSize', 14);
