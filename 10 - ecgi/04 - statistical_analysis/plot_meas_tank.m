% Define the sampling frequency (modify as needed)
fs = 4000; % Hz 

% Define time range
start_time = 0; % Modify as needed
end_time = 8;   % Modify as needed

% Convert time to sample indices
start_sample = round(start_time * fs) + 1;
end_sample = round(end_time * fs);

% Load electric data
meas_plot = meas_signal_f;
tank_plot = tank_signal_f;

% Time vectors
To = linspace(0, length(meas_plot(1, :)) / fs, length(meas_plot(1, :)));
time = linspace(start_time, end_time, length(start_sample:end_sample));

% Define electrode mappings
el1 = 1;   % MEA1 RA
el2 = 2;   % MEA1 RA
el3 = 65;  % MEA3 LA
el4 = 66;  % MEA3 LA
el5 = 17;  % MEA2 V
el6 = 18;  % MEA2 V
el7 = 134; % Tank
el8 = 185; % Tank

% Mapping the tank electrodes numbers
el7_plot = find(el_map(:,2) == el7, 1);
el8_plot = find(el_map(:,2) == el8, 1);

% Mapping the MEA3 electrodes numbers
el3_plot = el3 - 32;
el4_plot = el4 - 32;

% Create figure
f1 = figure('color', 'white', 'Position', [40 40 600 800]);

% Font size for improved readability
fontSize = 9; 

% Subplot 1 - MEA1 RA (Electrode el1)
subplot(8, 1, 1)
plot(To, meas_plot(el1, :), 'LineWidth', 1.2, 'Color', 'black');
ylabel('$\mu$V', 'Interpreter', 'latex', 'FontSize', fontSize);
title(['MEA1 RA - Electrode ', num2str(el1)], 'FontSize', fontSize);
set(gca, 'FontSize', fontSize);
xlim([start_time end_time]);

% Subplot 2 - MEA1 RA (Electrode el2)
subplot(8, 1, 2)
plot(To, meas_plot(el2, :), 'LineWidth', 1.2, 'Color', 'black');
ylabel('$\mu$V', 'Interpreter', 'latex', 'FontSize', fontSize);
title(['MEA1 RA - Electrode ', num2str(el2)], 'FontSize', fontSize);
set(gca, 'FontSize', fontSize);
xlim([start_time end_time]);

% Subplot 3 - MEA3 LA (Electrode el3)
subplot(8, 1, 3)
plot(To, meas_plot(el3_plot, :), 'LineWidth', 1.2, 'Color', 'black');
ylabel('$\mu$V', 'Interpreter', 'latex', 'FontSize', fontSize);
title(['MEA3 LA - Electrode ', num2str(el3)], 'FontSize', fontSize);
set(gca, 'FontSize', fontSize);
xlim([start_time end_time]);

% Subplot 4 - MEA3 LA (Electrode el4)
subplot(8, 1, 4)
plot(To, meas_plot(el4_plot, :), 'LineWidth', 1.2, 'Color', 'black');
ylabel('$\mu$V', 'Interpreter', 'latex', 'FontSize', fontSize);
title(['MEA3 LA - Electrode ', num2str(el4)], 'FontSize', fontSize);
set(gca, 'FontSize', fontSize);
xlim([start_time end_time]);

% Subplot 5 - MEA2 V (Electrode el5)
subplot(8, 1, 5)
plot(To, meas_plot(el5, :), 'LineWidth', 1.2, 'Color', 'black');
ylabel('$\mu$V', 'Interpreter', 'latex', 'FontSize', fontSize);
title(['MEA2 V - Electrode ', num2str(el5)], 'FontSize', fontSize);
set(gca, 'FontSize', fontSize);
xlim([start_time end_time]);

% Subplot 6 - MEA2 V (Electrode el6)
subplot(8, 1, 6)
plot(To, meas_plot(el6, :), 'LineWidth', 1.2, 'Color', 'black');
ylabel('$\mu$V', 'Interpreter', 'latex', 'FontSize', fontSize);
title(['MEA2 V - Electrode ', num2str(el6)], 'FontSize', fontSize);
set(gca, 'FontSize', fontSize);
xlim([start_time end_time]);

% Subplot 7 - Tank (Electrode el7)
subplot(8, 1, 7)
plot(To, tank_plot(el7_plot, :), 'LineWidth', 1.2, 'Color', 'red');
ylabel('$\mu$V', 'Interpreter', 'latex', 'FontSize', fontSize);
title(['TANK - Electrode ', num2str(el7)], 'FontSize', fontSize);
set(gca, 'FontSize', fontSize);
xlim([start_time end_time]);

% Subplot 8 - Tank (Electrode el8)
subplot(8, 1, 8)
plot(To, tank_plot(el8_plot, :), 'LineWidth', 1.2, 'Color', 'red');
ylabel('$\mu$V', 'Interpreter', 'latex', 'FontSize', fontSize);
title(['TANK - Electrode ', num2str(el8)], 'FontSize', fontSize);
xlabel('Time (s)', 'FontSize', fontSize);
set(gca, 'FontSize', fontSize);
xlim([start_time end_time]);

% Link all subplots for synchronized scrolling
linkaxes(findall(gcf, 'Type', 'axes'), 'x');
