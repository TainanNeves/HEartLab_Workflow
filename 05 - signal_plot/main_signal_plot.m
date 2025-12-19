%% SIGNAL PLOTING CODE

clear; clc;

%% Loading data

% Loading variables
load("E:\Qualification\Analysis\E32F02R08\data\data_filtered_sync_E32_F02_R08.mat"); % Synchronized data
load("E:\Qualification\Analysis\E32F02R08\data\InterpolatedSignalsE32_F02_R08_filtered.mat"); % Interpolate data

%% Optical signals plot
% Define a Camera to use
Data = D_SYNC.CAM3;
Background = squeeze(Data(:,:,2000));
[x, y] = pick_up_a_trace(Background, Data,1);
% Define a pixel position
p = [x(length(x)), y(length(y))];
% Frame Sampling
Fsampling = 4000;

% Title
str_title = ['Optical Signal - Right Atria'];

% Full optical time plot
% Create a time vector
To = linspace(0, length(Data(1,1,:))/Fsampling, length(Data(1,1,:)));
% Plot the optical signal at the specified pixel position
f1 = figure('color', 'white', 'Position', [40 40 600 200]);
plot(To, squeeze(Data(p(1), p(2), :)), 'k', 'LineWidth', 2, 'DisplayName', sprintf('Pixel: (%d, %d)', p(1), p(2))); % 'k' for black
ylabel('%F');
set(gca, 'fontsize', 14);
xlim([0 8]);
legend('show', 'Location', 'northeast'); % Show the legend in the top-right
title(str_title);

% Specific optical time plot
% Define the time range for the plot
t_in = 2;
t_out = 4;
start_sample = t_in*Fsampling; % Adjust the start sample according to your data
end_sample = t_out*Fsampling;   % Adjust the end sample according to your data
% Create a time vector
To = linspace(0, length(Data(1,1,:))/Fsampling, length(Data(1,1,:)));
% Plot the optical signal within the specified time range
f2 = figure('color', 'white', 'Position', [40 40 600 200]); % Changed to f2 to open a new figure
plot(To(start_sample:end_sample), squeeze(Data(p(1), p(2), start_sample:end_sample)), 'k', 'LineWidth', 2, 'DisplayName', sprintf('Pixel: (%d, %d)', p(1), p(2))); % 'k' for black
ylabel('%F');
set(gca, 'fontsize', 14);
xlabel('Time (s)');
xlim([t_in, t_out]);
legend('show', 'Location', 'northeast'); % Show the legend in the top-right
title(str_title);


%% Electric signal plot
% Define electrode to use
el = [4 2 6 11 14 15];
Data = D_SYNC.EL(el,:); 
Fsampling = 4000;
% Title for the plots
str_title = "Electrical Signal";
To = linspace(0, length(Data)/Fsampling, size(Data, 2));

% Fultime Plot
for i = 1:length(el)
    f_full = figure('color', 'white', 'Position', [40 40 600 200]);
    plot(To, squeeze(Data(i,:)), 'k', 'LineWidth', 2, 'DisplayName', ['Electrode ' num2str(el(i))]); 
    hold on
    ylabel('Potential ($\mu$V)', 'Interpreter', 'latex');
    set(gca, 'fontsize', 14);
    xlabel('Time (s)');
    xlim([0 length(Data)/Fsampling]);
    title([str_title ' - El' num2str(el(i))]);
    legend('show', 'Location', 'northeast');
    grid on;
    hold off 
end

% Specific Time Plot
t_in = 2;
t_out = 4;
start_sample = t_in*Fsampling + 1;
end_sample = t_out*Fsampling;   
To_specific = linspace(t_in, t_out, end_sample - start_sample + 1);
for i = 1:length(el)
    f_spec = figure('color', 'white', 'Position', [40 40 600 200]);
    plot(To_specific, squeeze(Data(i, start_sample:end_sample)), 'k', 'LineWidth', 2, 'DisplayName', ['Electrode ' num2str(el(i))]); 
    hold on
    ylabel('Potential ($\mu$V)', 'Interpreter', 'latex');
    set(gca, 'fontsize', 14);
    xlabel('Time (s)');
    xlim([t_in t_out]);
    title([str_title '- El' num2str(el(i))]);
    legend('show', 'Location', 'northeast'); 
    grid on; 
    hold off
end


%% MULTIPLOT 1 (1 Cam (2 points) + 5 Electrodes)
% Define a Camera to use
Data_O = D_SYNC.CAM1;

% Optical Pixel
Background = squeeze(Data_O(:,:,2000));
[x, y] = pick_up_a_trace(Background, Data_O, 1);    
pa = [x(1) y(1)];
pv = [x(2) y(2)];

% Electric electrodes
el1 = 7; % RA
el2 = 78; % LA
el3 = 25; % V
el4 = 183; % Tank
el5 = 172; % Tank

% Frequency Sampling
Fsampling = 4000;
% Time to plot
start_time = 2;
end_time = 4;


% Loading electric data
Data_E = D_SYNC.EL;
% Ploting
To = linspace(0, length(Data_E(1, :)) / Fsampling, length(Data_E(1, :)));
f1 = figure('color', 'white', 'Position', [40 40 600 800]);

% Subplot 1 - Optic Atrium
subplot(7, 1, 1)
plot(To, squeeze(Data_O(pa(1), pa(2), :)), 'LineWidth', 1);
ylabel('%F');
title('Optical signal (Atria)');
set(gca, 'fontsize', 14);
xlim([start_time end_time])

% Subplot 2 - Optic Ventricle
subplot(7, 1, 2)
plot(To, squeeze(Data_O(pv(1), pv(2), :)), 'LineWidth', 1);
ylabel('%F');
title('Optical signal (Ventricle)');
set(gca, 'fontsize', 14);
xlim([start_time end_time])

% Subplot 3 - RA
i = el1;
subplot(7, 1, 3)
plot(To, Data_E(i, :), 'LineWidth', 1);
ylabel('$\mu$V', 'Interpreter', 'latex');
title(['el', num2str(i), ' (RA)']);
set(gca, 'fontsize', 14);
xlim([start_time end_time])

% Subplot 4 - LA
i = el2;
subplot(7, 1, 4)
plot(To, Data_E(i, :), 'LineWidth', 1);
ylabel('$\mu$V', 'Interpreter', 'latex');
title(['el', num2str(i), ' (LA)']);
set(gca, 'fontsize', 14);
xlim([start_time end_time])

% Subplot 5 - V
i = el3;
subplot(7, 1, 5)
plot(To, Data_E(i, :), 'LineWidth', 1);
ylabel('$\mu$V', 'Interpreter', 'latex');
title(['el', num2str(i), ' (V)']);
set(gca, 'fontsize', 14);
xlim([start_time end_time])

% Subplot 6 - TANK
i = el4;
subplot(7, 1, 6)
plot(To, Data_E(i, :), 'LineWidth', 1);
ylabel('$\mu$V', 'Interpreter', 'latex');
title(['TANK el', num2str(i)]);
set(gca, 'fontsize', 14);
xlim([start_time end_time])

% Subplot 7 - TANK
i = el5;
subplot(7, 1, 7)
plot(To, Data_E(i, :), 'LineWidth', 1);
ylabel('$\mu$V', 'Interpreter', 'latex');
title(['TANK el', num2str(i)]);
xlabel('Time (s)');
set(gca, 'fontsize', 14);
xlim([start_time end_time])

% Linking Axes
linkaxes([subplot(7, 1, 1), subplot(7, 1, 2), subplot(7, 1, 3), ...
    subplot(7, 1, 4), subplot(7, 1, 5), subplot(7, 1, 6), ...
    subplot(7, 1, 7)], 'x')


%% MULTIPLOT 2 (3 Cam + 5 Electrodes)
% Define Cameras to use
Data_OV = D_SYNC.CAM1;
Background = squeeze(Data_OV(:,:,2000));
[xV, yV] = pick_up_a_trace(Background, Data_OV,1);
Data_ORA = D_SYNC.CAM3;
Background = squeeze(Data_ORA(:,:,2000));
[xRA, yRA] = pick_up_a_trace(Background, Data_ORA,1);
Data_OLA = D_SYNC.CAM2;
Background = squeeze(Data_OLA(:,:,2000));
[xLA, yLA] = pick_up_a_trace(Background, Data_OLA,1);

% Optical Pixel
pV = [xV yV];
pRA = [xRA yRA];
pLA = [xLA yLA];

% Electric electrodes
el1 = 7; % RA
el2 = 22; % LA
el3 = 90; % V
el4 = 142; % Tank 1
el5 = 178; % Tank 2

% Parameters
Fsampling = 4000;
start_time = 2;
end_time = 4;

% Loading electric data
Data_E = D_SYNC.EL;
% Ploting
To = linspace(0, length(Data_E(1, :)) / Fsampling, length(Data_E(1, :)));
f1 = figure('color', 'white', 'Position', [40 40 600 800]);
% Subplot 1 - Optic Right Atrium
subplot(8, 1, 1)
plot(To, squeeze(Data_ORA(pRA(1), pRA(2), :)), 'LineWidth', 1);
ylabel('%F');
title(['Optical Signal (RA) | Points: ', num2str(pRA(1)), 'x', num2str(pRA(2))]);
set(gca, 'fontsize', 14);
xlim([start_time end_time])

% Subplot 2 - Optic Left Atrium
subplot(8, 1, 2)
plot(To, squeeze(Data_OLA(pLA(1), pLA(2), :)), 'LineWidth', 1);
ylabel('%F');
title(['Optical Signal (LA) | Points: ', num2str(pLA(1)), 'x', num2str(pLA(2))]);
set(gca, 'fontsize', 14);
xlim([start_time end_time])

% Subplot 3 - Optic Ventricle
subplot(8, 1, 3)
plot(To, squeeze(Data_OV(pV(1), pV(2), :)), 'LineWidth', 1);
ylabel('%F');
title(['Optical Signal (V) | Points: ', num2str(pV(1)), 'x', num2str(pV(2))]);
set(gca, 'fontsize', 14);
xlim([start_time end_time])

% Subplot 4 - RA
i = el1;
subplot(8, 1, 4)
plot(To, Data_E(i, :), 'LineWidth', 1);
ylabel('$\mu$V', 'Interpreter', 'latex');
title(['electrode ', num2str(i), ' (RA)']);
set(gca, 'fontsize', 14);
xlim([start_time end_time])

% Subplot 5 - LA
i = el2;
subplot(8, 1, 5)
plot(To, Data_E(i, :), 'LineWidth', 1);
ylabel('$\mu$V', 'Interpreter', 'latex');
title(['electrode ', num2str(i), ' (LA)']);
set(gca, 'fontsize', 14);
xlim([start_time end_time])

% Subplot 6 - V
i = el3;
subplot(8, 1, 6)
plot(To, Data_E(i, :), 'LineWidth', 1);
ylabel('$\mu$V', 'Interpreter', 'latex');
title(['electrode ', num2str(i), ' (V)']);
set(gca, 'fontsize', 14);
xlim([start_time end_time])

% Subplot 7 - TANK 1
i = el4;
subplot(8, 1, 7)
plot(To, Data_E(i, :), 'LineWidth', 1);
ylabel('$\mu$V', 'Interpreter', 'latex');
title(['electrode ', num2str(i), ' (TANK)']);
set(gca, 'fontsize', 14);
xlim([start_time end_time])

% Subplot 8 - TANK 2
i = el5;
subplot(8, 1, 8)
plot(To, Data_E(i, :), 'LineWidth', 1);
ylabel('$\mu$V', 'Interpreter', 'latex');
title(['electrode ', num2str(i), ' (TANK)']);
xlabel('Time (s)');
set(gca, 'fontsize', 14);
xlim([start_time end_time])

% Linking Axes (move x axes in all data together)
linkaxes([subplot(8, 1, 1), subplot(8, 1, 2), subplot(8, 1, 3), ...
    subplot(8, 1, 4), subplot(8, 1, 5), subplot(8, 1, 6), ...
    subplot(8, 1, 7), subplot(8, 1, 8)], 'x');


%% All plots (Optic + MEAs) - 3 Figures
% Configuring Electrical
case_name_e = 'MEA2';
Fsampling = 4000;
time = [2, 3]; % s
Data_E = D_SYNC.EL(:, time(1)*Fsampling:time(2)*Fsampling);
el.MEA1 = [1, 2, 3, 4, ...
            5, 6, 7, 8, ...
            9, 10, 11, 12, ...
            13, 14, 15, 16];
el.MEA2 = [17, 18, 19, 20, ...
            21, 22, 23, 24, ...
            25, 26, 27, 28, ...
            29, 30, 31, 32];
el.MEA3 = [81, 82, 83, 84, ...
            85, 86, 87, 88, ...
            89, 90, 91, 92, ...
            93, 94, 95, 96];

% Configuring Optical
case_name_o = 'CAM3';
Data_O = D_SYNC.(case_name_o)(:, :, time(1)*Fsampling:time(2)*Fsampling);
Background = squeeze(Data_O(:,:,2));
Background = single(D_SYNC.IMG.(case_name_o)); Data_O = rot90(Data_O, 1); % To use the figure with no filtering
[x, y] = pick_up_a_trace(Background, Data_O,1);

% Figure 1: Camera and electrodes identification
ni = 0; % Sum this to the number (In case you dont want to start in 1
I = single(D_SYNC.IMG.(case_name_o));
f1 = figure('color','white','Position', [40 40 500 500]);
imagesc(I); colormap('gray');
hold on;
for i = 1:length(x)
    text(y(i), x(i), num2str(ni+i), 'Color', 'red', 'FontSize', 9, 'FontWeight', 'bold'); hold on;
end
set(gca,'fontsize', 14);
title(case_name_o);
ylabel('Pixels');xlabel('Pixels');

% Figure 2: Optic plots
lo = size(Data_O,3);
To = linspace(0,lo/Fsampling,lo);
f2=figure('color','white','Position', [40 40 800 600]);
for i = 1:size(x,1)
    subplot(4,4,i); plot(To,squeeze(Data_O(x(i), y(i), :)));
    ylabel(['pixel: ', num2str(x(i)),'x', num2str(y(i))]);
    xlabel('Time [s]');
    set(gca,'fontsize', 12); xlim([0 lo/Fsampling]);
end
sgtitle(['Optic signals: ', case_name_o]);

% Figure 3: Electric Plots
ni = 0; % Sum this to the number (In case you dont want to start in 1
lo=size(Data_E,2);
To=linspace(0,lo/Fsampling,lo);
f3=figure('color','white','Position', [40 40 800 600]);
for i = 1:size(el.(case_name_e), 2)
    subplot(4,4,i); plot(To, Data_E(el.(case_name_e)(i),:));
    ylabel(['El: ', num2str(el.(case_name_e)(i))]);
    xlabel('Time [s]');
    set(gca,'fontsize', 12); xlim([0 lo/Fsampling]);
end
sgtitle(['Electric signals: ', case_name_e]);


%% All plots (Optic + MEAs) - 2 Figures
% Configuring Electrical
case_name_e = 'MEA2';
Fsampling = 4000;
time = [2, 3]; % s
Data_E = D_SYNC.EL(:, time(1)*Fsampling:time(2)*Fsampling);
el.MEA1 = [1, 2, 3, 4, ...
            5, 6, 7, 8, ...
            9, 10, 11, 12, ...
            13, 14, 15, 16];
el.MEA2 = [17, 18, 19, 20, ...
            21, 22, 23, 24, ...
            25, 26, 27, 28, ...
            29, 30, 31, 32];
el.MEA3 = [81, 82, 83, 84, ...
            85, 86, 87, 88, ...
            89, 90, 91, 92, ...
            93, 94, 95, 96];

% Configuring Optical
case_name_o = 'CAM2';
Data_O = D_SYNC.(case_name_o)(:, :, time(1)*Fsampling:time(2)*Fsampling);
Background = squeeze(Data_O(:,:,2));
Background = single(D_SYNC.IMG.(case_name_o)); Data_O = rot90(Data_O, 1); % To use the figure with no filtering
[x, y] = pick_up_a_trace(Background, Data_O,1);

% Figure 1: Camera and electrodes identification
ni = 16; % Sum this to the number in case you dont want to start in 1
I = single(D_SYNC.IMG.(case_name_o));
f1 = figure('color','white','Position', [40 40 500 500]);
imagesc(I); colormap('gray');
hold on;
for i = 1:length(x)
    text(y(i), x(i), num2str(ni+i), 'Color', 'red', 'FontSize', 9, 'FontWeight', 'bold'); hold on;
end
set(gca,'fontsize', 14);
title(case_name_o);
ylabel('Pixels');xlabel('Pixels');

% Figure 2: Optical and Electrical signals Aligned
lo = size(Data_O, 3);
To = linspace(0, lo/Fsampling, lo);
f2 = figure('color','white','Position', [100 100 1000 800]);
for i = 1:16
    subplot(4, 4, i);
    % Left Axis: Optical
    yyaxis left
    plot(To, squeeze(Data_O(x(i), y(i), :)), 'b'); % Blue for optical
    ylabel(['pixel: ', num2str(x(i)),'x', num2str(y(i))]);
    % Right Axis: Electrical
    yyaxis right
    plot(To, Data_E(el.(case_name_e)(i),:), 'r'); % Red for electrical
    ylabel(['El: ', num2str(el.(case_name_e)(i))]);
    xlabel('Time [s]');
    set(gca, 'fontsize', 10);
    xlim([0 lo/Fsampling]);
end
sgtitle(['Combined Signals: ', case_name_o, ' (Opt) & ', case_name_e, ' (Elec)']);


%% All plots (TANK)
Fsampling = 4000;
start_time = 2.25; % s
end_time = 2.85;   % s

% Ensure Data_E is loaded (assuming D_SYNC is available from previous sections)
Data_E = D_SYNC.EL;
To = linspace(0, length(Data_E(1, :)) / Fsampling, length(Data_E(1, :)));

% --- Determine Global Amplitude Limits for Tank Electrodes ---
% Collect all tank electrode numbers from both figures
all_tank_electrodes = unique([...
    145 146 155 156 165 166 129 130 139 140 181 182, ... % Figure 1, Row 1
    147 157 167 131 141 183, ... % Figure 1, Row 2
    148 149 158 159 168 169 132 133 142 143 184 185, ... % Figure 1, Row 3
    150 151 160 161 170 171 134 135 144 177 186 187, ... % Figure 2, Row 1
    152 162 172 136 178 188, ... % Figure 2, Row 2
    153 154 163 164 173 174 137 138 179 180 189 190  ... % Figure 2, Row 3
]);

% Find the time indices corresponding to start_time and end_time
time_indices = To >= start_time & To <= end_time;
% Extract data for all relevant tank electrodes within the specified time window
relevant_data = Data_E(all_tank_electrodes, time_indices);
% Calculate the global minimum and maximum amplitude
global_ymin = min(relevant_data(:));
global_ymax = max(relevant_data(:));
% Add a small margin to the limits for better visualization
amplitude_margin = 0.05 * (global_ymax - global_ymin);
global_ymin = global_ymin - amplitude_margin;
global_ymax = global_ymax + amplitude_margin;
% Define the subplot grid for consistent layout (Rows x Columns)
num_rows_grid = 3;
num_cols_grid = 12; 

% --- Figure 1: Up Tank Electrodes ---
% Create the first figure for the 'Up' set of tank electrodes
f_up = figure('color', 'white', 'Position', [100 100 1200 800], 'Name', 'Tank Electrodes (Up)');
sgtitle('Tank - Top Section', 'FontSize', 16, 'FontWeight', 'bold'); % Super title for the figure
% Define electrodes for Figure 1, grouped by the user's specified rows
tank_el_fig1_row1 = [145 146 155 156 165 166 129 130 139 140 181 182];
tank_el_fig1_row2 = [147 157 157 131 141 183];
tank_el_fig1_row3 = [148 149 158 159 168 169 132 133 142 143 184 185];
all_tank_el_fig1 = {tank_el_fig1_row1, tank_el_fig1_row2, tank_el_fig1_row3};
linked_axes_handles_fig1 = []; % Initialize an array to store axis handles for linking
current_plot_index = 0; % Keep track of the global subplot position (1 to 36)
for r = 1:num_rows_grid
    electrodes_in_current_row = all_tank_el_fig1{r};
    num_electrodes_in_this_row = length(electrodes_in_current_row);
    current_electrode_to_pick = 1; % Index to pick from electrodes_in_current_row
    for c = 1:num_cols_grid % Iterate through all 12 columns of the grid
        current_plot_index = current_plot_index + 1; % Increment global subplot index for this grid cell
        should_plot_this_cell = false;
        el_to_plot = -1; % Default value
        if r == 2 % This is the special second row where we want alternating plots
            % Plot if current column index 'c' is odd AND we still have electrodes to plot
            if mod(c, 2) == 1 && current_electrode_to_pick <= num_electrodes_in_this_row
                should_plot_this_cell = true;
                el_to_plot = electrodes_in_current_row(current_electrode_to_pick);
                current_electrode_to_pick = current_electrode_to_pick + 1; % Move to next electrode
            end
        else % For rows 1 and 3, plot from left to right as before
            if c <= num_electrodes_in_this_row
                should_plot_this_cell = true;
                el_to_plot = electrodes_in_current_row(c);
            end
        end
        if should_plot_this_cell
            % Create subplot and plot the data
            ax = subplot(num_rows_grid, num_cols_grid, current_plot_index);
            plot(To, Data_E(el_to_plot, :), 'k', 'LineWidth', 1); % Plotted in black
            % Set titles
            title(['el', num2str(el_to_plot)], 'FontSize', 8); % Smaller font for many subplots
            set(gca, 'fontsize', 7); % Adjust tick label font size
            % Apply universal y-axis limits
            ylim([global_ymin global_ymax]);
            xlim([start_time end_time]);
            % --- Label Optimization Logic ---
            % Only show Y-label on the leftmost column (c=1)
            if c == 1
                ylabel('$\mu$V', 'Interpreter', 'latex');
            else
                set(ax, 'YTickLabel', []); % Hide Y-axis tick labels for other columns
            end
            % Only show X-label on the bottom row (r=num_rows_grid)
            if r == num_rows_grid
                xlabel('Time (s)');
            else
                set(ax, 'XTickLabel', []); % Hide X-axis tick labels for non-bottom rows
            end
            % --- End Label Optimization Logic ---
            % Add the axis handle to the list for linking later
            linked_axes_handles_fig1 = [linked_axes_handles_fig1, ax];
        end
    end
end
% Link the x-axes of all subplots that were actually created in this figure
linkaxes(linked_axes_handles_fig1, 'x');

% --- Figure 2: Down Tank Electrodes ---
% Create the second figure for the 'Down' set of tank electrodes
f_down = figure('color', 'white', 'Position', [100 50 1200 800], 'Name', 'Tank Electrodes (Down)');
sgtitle('Tank - Botton Section', 'FontSize', 16, 'FontWeight', 'bold'); % Super title for the figure
% Define electrodes for Figure 2, grouped by the user's specified rows
tank_el_fig2_row1 = [150 151 160 161 170 171 134 135 144 177 186 187];
tank_el_fig2_row2 = [152 162 172 136 178 188];
tank_el_fig2_row3 = [153 154 163 164 173 174 137 138 179 180 189 190];
all_tank_el_fig2 = {tank_el_fig2_row1, tank_el_fig2_row2, tank_el_fig2_row3};
linked_axes_handles_fig2 = []; % Initialize an array for this figure's axis handles
current_plot_index = 0; % Reset global plot index for the new figure
for r = 1:num_rows_grid
    electrodes_in_current_row = all_tank_el_fig2{r};
    num_electrodes_in_this_row = length(electrodes_in_current_row);
    current_electrode_to_pick = 1; % Index to pick from electrodes_in_current_row
    for c = 1:num_cols_grid % Iterate through all 12 columns of the grid
        current_plot_index = current_plot_index + 1; % Increment global subplot index
        should_plot_this_cell = false;
        el_to_plot = -1;
        if r == 2 % This is the special second row where we want alternating plots
            % Plot if current column index 'c' is odd AND we still have electrodes to plot
            if mod(c, 2) == 1 && current_electrode_to_pick <= num_electrodes_in_this_row
                should_plot_this_cell = true;
                el_to_plot = electrodes_in_current_row(current_electrode_to_pick);
                current_electrode_to_pick = current_electrode_to_pick + 1;
            end
        else % For rows 1 and 3, plot from left to right as before
            if c <= num_electrodes_in_this_row
                should_plot_this_cell = true;
                el_to_plot = electrodes_in_current_row(c);
            end
        end
        if should_plot_this_cell
            ax = subplot(num_rows_grid, num_cols_grid, current_plot_index);
            plot(To, Data_E(el_to_plot, :), 'k', 'LineWidth', 1); % Plotted in black
            title(['el', num2str(el_to_plot)], 'FontSize', 8);
            set(gca, 'fontsize', 7);
            % Apply universal y-axis limits
            ylim([global_ymin global_ymax]);
            xlim([start_time end_time]);
            % --- Label Optimization Logic ---
            % Only show Y-label on the leftmost column (c=1)
            if c == 1
                ylabel('$\mu$V', 'Interpreter', 'latex');
            else
                set(ax, 'YTickLabel', []); % Hide Y-axis tick labels for other columns
            end
            % Only show X-label on the bottom row (r=num_rows_grid)
            if r == num_rows_grid
                xlabel('Time (s)');
            else
                set(ax, 'XTickLabel', []); % Hide X-axis tick labels for non-bottom rows
            end
            % --- End Label Optimization Logic ---
            linked_axes_handles_fig2 = [linked_axes_handles_fig2, ax];
        end
    end
end
% Link the x-axes of all subplots that were actually created in this figure
linkaxes(linked_axes_handles_fig2, 'x');





