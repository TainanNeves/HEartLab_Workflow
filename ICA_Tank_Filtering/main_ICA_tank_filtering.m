%% Main ICA Tank Filtering
clear; clc;


%% Loading data
load(""); % Load Sync Data


%% Defining electrodes
Data = D_SYNC.EL;

% Group the indices
group_indices = {
    129:133;
    134:138;
    139:143;
    [144, 177, 178, 179, 180];
    181:185;
    186:190;
    145:149;
    150:154;
    155:159;
    160:164;
    165:169;
    170:174};

% Creating a temporary copy
Data_temp = double(Data);


%% ICA Filtering
Fsampling = 4000;
N = size(Data, 2);
To = linspace(0, N / Fsampling, N); % Time vector

numGroups = length(group_indices);

% --- START INTERACTIVE LOOP ---
for g = 1:numGroups
    % Get the indices and name for the current group
    current_indices = group_indices{g};
    current_name = ['Group ', num2str(g), ' (Electrodes: ', ...
                    strjoin(arrayfun(@num2str, current_indices, 'UniformOutput', false), ', '), ')'];
    
    disp(' '); % Spacer
    disp(['*******************************************************************']);
    disp(['*** STARTING PROCESSING FOR: ', current_name, ' ***']);
    disp(['*******************************************************************']);
    
    % 1. Extract the data subset
    Data_subset = double(Data(current_indices, :));
    
    % 2. Run the INTERACTIVE ICA analysis (The script will PAUSE here)
    % Note: The function now only takes 5 arguments (it removes componentsToUse)
    results = ica_analysis_interactive(Data_subset, current_name, current_indices, ...
                                       To, Fsampling);
    
    % 3. Store the reconstructed data
    ICA_Group_Data = results.Data_reconstructed;
    
    % 4. Consolidation
    Data_temp(current_indices, :) = ICA_Group_Data;
    
    disp(['Successfully Reconstructed with ICs: [', ...
          strjoin(arrayfun(@num2str, results.selected_components, 'UniformOutput', false), ', '), ']']);
    disp(['Processing of ', current_name, ' complete.']);
    disp(' '); % Spacer
end

disp(' ');
disp('--- All ICA groups processed and consolidated in Data_temp. ---');


%% Ploting Filtered Signals
close all; 

figure('Name', 'ICA Filtering Results: Original vs. Reconstructed', ...
       'Position', [50 50 1400 900], 'Color', 'white');

sgtitle('Comparison of Original vs. ICA-Filtered Signal (First Electrode of Each Group)', ...
        'FontSize', 16, 'FontWeight', 'bold');

numGroups = length(group_indices);

% Determine subplot layout (e.g., 4 rows by 3 columns for 12 groups)
rows = 4; 
cols = 3; 
if numGroups > 12 % Adjust for scripts with more groups
    rows = ceil(sqrt(numGroups));
    cols = ceil(numGroups / rows);
end

% Initialize an array to store the axes handles for linking
h_axes = zeros(numGroups, 1);

for g = 1:numGroups
    % Get the index of the first electrode in the current group
    original_idx = group_indices{g}(1);
    
    % Extract the signals using the variables defined in your script
    original_signal = Data(original_idx, :);       % The raw signal
    filtered_signal = Data_temp(original_idx, :);  % The ICA-filtered signal
    
    % Create subplot and store its handle
    h_axes(g) = subplot(rows, cols, g);
    
    % Plot original signal
    plot(To, original_signal, 'Color', [0.8 0.2 0.2], 'LineWidth', 1);
    hold on;
    % Plot reconstructed signal
    plot(To, filtered_signal, 'Color', [0.2 0.6 0.2], 'LineWidth', 1.5);
    hold off;
    
    % Add legend and title
    title(['Electrode ', num2str(original_idx), ' (Group ', num2str(g), ')'], ...
          'FontSize', 10, 'FontWeight', 'bold');
    xlabel('Time (s)', 'FontSize', 8);
    ylabel('Amplitude (\muV)', 'FontSize', 8);
    xlim([0 8]); 
    
    if g == 1
        legend('Original Signal', 'ICA Filtered', 'Location', 'NorthEast');
    end

    set(gca, 'FontSize', 7);
    box on;
end
linkaxes(h_axes, 'y');


%% Final Consolidation and Data Saving
D_SYNC.EL = Data_temp;

save('data_filtered_sync_EXX_FXX_RXX.mat', 'D_SYNC', '-v7.3'); 

disp(' ');
disp('*******************************************************************');
disp('*** FINAL CONSOLIDATION COMPLETE ***');
disp(['Filtered Data saved to D_SYNC.EL and written to:']);
disp('data_filtered_sync_EXX_FXX_RXX.mat');
disp('*******************************************************************');
