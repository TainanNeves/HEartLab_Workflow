%% SIGNAL PLOTING CODE

clear; clc;

%% Loading data

% Loading variables
load('C:\Users\HEartLab\Downloads\Pasta de Trabalho\Subpasta 2 - Desenvolvimento\Dados Aprendizagem\1 - Analisados\data_filtered_sync_E14_F3_R4.mat'); % Synchronized data
load('C:\Users\HEartLab\Downloads\Pasta de Trabalho\Subpasta 2 - Desenvolvimento\Dados Aprendizagem\1 - Analisados\InterpolatedSignalsE18_F02_R02_filtered'); % Interpolate data

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
t_out = 6;
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
el = 183;
Data = D_SYNC.EL(el,:);
% Frame Sampling
Fsampling = 4000;
% Title for the plots
str_title = ["Electrical Signal - TANK"];

% Full electric time plot
% Create a time vector
To = linspace(0, length(Data)/Fsampling, length(Data));
% Plot the electrical signal for a specific electrode
f2 = figure('color', 'white', 'Position', [40 40 600 200]);
plot(To, squeeze(Data), 'k', 'LineWidth', 2, 'DisplayName', ['Electrode ' num2str(el)]); % Black line, size 2, legend by electrode
hold on
ylabel('Potential ($\mu$V)', 'Interpreter', 'latex');
set(gca, 'fontsize', 14);
xlim([0 8]);
title(str_title);
legend('show', 'Location', 'northeast'); % Show legend
grid on; % Adicionado para melhor visualização

% Specific electric time plot
% Define the time range for the plot
t_in = 2;
t_out = 6;
start_sample = t_in*Fsampling; % Adjust the start sample according to your data
end_sample = t_out*Fsampling;   % Adjust the end sample according to your data
% Create a time vector for the specific range
To_specific = linspace(t_in, t_out, end_sample - start_sample + 1);
% Plot the electrical signal for a specific electrode within the defined time range
f3 = figure('color', 'white', 'Position', [40 40 600 200]);
plot(To_specific, squeeze(Data(start_sample:end_sample)), 'k', 'LineWidth', 2, 'DisplayName', ['Electrode ' num2str(el)]); % Black line, size 2, legend by electrode
hold on
ylabel('Potential ($\mu$V)', 'Interpreter', 'latex');
set(gca, 'fontsize', 14);
xlabel('Time (s)');
xlim([t_in t_out]);
title(str_title);
legend('show', 'Location', 'northeast'); % Show legend
grid on; % Adicionado para melhor visualização


%% MULTIPLOT (1 Cam + 5 Electrodes)
% Define a Camera to use
Data_O = D_SYNC.CAM2;

% Optical Pixel
Background = squeeze(Data_O(:,:,2000));
[x, y] = pick_up_a_trace(Background, Data_O, 1);    
pa = [x(1) y(1)];
pv = [x(2) y(2)];

% Electric electrodes
el1 = 18; % RA
el2 = 3; % LA
el3 = 90; % V
el4 = 183; % Tank
el5 = 173; % Tank

% Frequency Sampling
Fsampling = 4000;
% Time to plot
start_time = 2;
end_time = 3;


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


%% MULTIPLOT (3 Cam + 5 Electrodes)
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
el1 = 18; % RA
el2 = 3; % LA
el3 = 90; % V
el4 = 183; % Tank 1
el5 = 172; % Tank 2

% Parameters
Fsampling = 4000;
start_time = 2;
end_time = 3;

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


%% All plots (Optic + MEAs)
%Frequency Sampling
Fsampling = 4000;
time = [2, 4]; % s
Data_E = D_SYNC.EL(:, time(1)*Fsampling:time(2)*Fsampling);


% MEA 1
Data_O = D_SYNC.CAM1(:, :, time(1)*Fsampling:time(2)*Fsampling);
Background = squeeze(Data_O(:,:,2000));
% Background = single(D_SYNC.IMG.CAM1); Data_O = rot90(Data_O, 1); % To use the figure with no filtering
[x, y] = pick_up_a_trace(Background, Data_O,1);  
p13 = [x(13), y(13)]; p14 = [x(14), y(14)]; p15 = [x(15), y(15)]; p16 = [x(16), y(16)];
p9 = [x(9), y(9)]; p10 = [x(10), y(10)]; p11 = [x(11), y(11)]; p12 = [x(12), y(12)];
p5 = [x(5), y(5)]; p6 = [x(6), y(6)]; p7 = [x(7), y(7)]; p8 = [x(8), y(8)];
p1 = [x(1), y(1)]; p2 = [x(2), y(2)]; p3 = [x(3), y(3)]; p4 = [x(4), y(4)];
plotar_pontos_1(Data_O, Data_E, Fsampling, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16);

% MEA 2
Data_O = D_SYNC.CAM2(:, :, time(1)*Fsampling:time(2)*Fsampling);
Background = squeeze(Data_O(:,:,2000));
% Background = single(D_SYNC.IMG.CAM2); Data_O = rot90(Data_O, 1); % To use the figure with no filtering
[x, y] = pick_up_a_trace(Background, Data_O,1); 
p17 = [x(1), y(1)]; p21 = [x(5), y(5)]; p25 = [x(9), y(9)]; p29 = [x(13), y(13)];
p18 = [x(2), y(2)]; p22 = [x(6), y(6)]; p26 = [x(10), y(10)]; p30 = [x(14), y(14)];
p19 = [x(3), y(3)]; p23 = [x(7), y(7)]; p27 = [x(11), y(11)]; p31 = [x(15), y(15)];
p20 = [x(4), y(4)]; p24 = [x(8), y(8)]; p28 = [x(12), y(12)]; p32 = [x(16), y(16)];
plotar_pontos_2(Data_O, Data_E, Fsampling, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32);

% MEA 3
Data_O = D_SYNC.CAM3(:, :, time(1)*Fsampling:time(2)*Fsampling);
Background = squeeze(Data_O(:,:,2000));
% Background = single(D_SYNC.IMG.CAM3); Data_O = rot90(Data_O, 1); % To use the figure with no filtering
[x, y] = pick_up_a_trace(Background, Data_O,1); 
p77 = [x(13), y(13)]; p78 = [x(14), y(14)]; p79 = [x(15), y(15)]; p80 = [x(16), y(16)];
p73 = [x(9), y(9)]; p74 = [x(10), y(10)]; p75 = [x(11), y(11)]; p76 = [x(12), y(12)];
p69 = [x(5), y(5)]; p70 = [x(6), y(6)]; p71 = [x(7), y(7)]; p72 = [x(8), y(8)];
p65 = [x(1), y(1)]; p66 = [x(2), y(2)]; p67 = [x(3), y(3)]; p68 = [x(4), y(4)];
plotar_pontos_3(Data_O, Data_E, Fsampling, p65, p66, p67, p68, p69, p70, p71, p72, p73, p74, p75, p76, p77, p78, p79, p80);


%% All plots (TANK)
Fsampling = 4000;
start_time = 2.42; % s
end_time = 2.82;   % s

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





