%% Main LAT CV Tainan
%% LAT and CV Maps - Optical & Electrical
clear; clc;

%% Loading Data
load(""); % Load Synchronized data
load(""); % Load Interpolated data (if needed for electrical)

fprintf('=== LAT and CV Analysis Pipeline ===\n');

%% Section 1: Optical LAT Analysis
fprintf('\n--- Section 1: Optical LAT Analysis ---\n');

% Parameters for optical LAT
Fs = 4000;
optical_lim1 = 9460;   % Start sample for optical
optical_lim2 = 9870;   % End sample for optical
linear_fit_length = 15;
PCL = 200;
debug_LAT = 0;

% Data preparation
Data_O = D_SYNC.CAM1;
Data_temp_O = Data_O(:,:,optical_lim1:optical_lim2);

% Preview signal
fprintf('Previewing optical signal...\n');
Background = squeeze(Data_O(:,:,2000));
[px, py] = pick_up_a_trace(Background, Data_temp_O, 1);

% Calculate LAT for each pixel
fprintf('Calculating optical LAT...\n');
LAT_O = Data_temp_O(:,:,1) * 0;

for i = 1:size(Data_temp_O,1)
    for j = 1:size(Data_temp_O,2)
        if max(squeeze(Data_temp_O(i,j,:))) ~= 0
            signal = squeeze(Data_temp_O(i,j,:));
            LAT_O(i,j) = find_LAT_linearFit_1D(signal, Fs, linear_fit_length, PCL, debug_LAT);
        end
    end
end

% Post-processing: subtract minimum
LAT_O(LAT_O < 0) = 0;
min_LAT = min(LAT_O(LAT_O ~= 0));
if ~isempty(min_LAT)
    LAT_O(LAT_O ~= 0) = LAT_O(LAT_O ~= 0) - min_LAT;
end

% Filter LAT data
fprintf('Filtering optical LAT data...\n');
LAT_O_filtered = SpatTemp_Filtering(LAT_O, 3, 0, 'GPU');

% Plot LAT map
fprintf('Plotting optical LAT map...\n');
figure('color', 'white', 'Position', [40 40 600 600]);
C = parula(256);
C(1,1:3) = [1 1 1];
C(256,1:3) = [0.5, 0.5, 0.5];
J = imrotate(LAT_O_filtered, 90);
contourf(flipud(J));
colormap(C);

% Ask user for color limits
current_limits = [min(LAT_O_filtered(LAT_O_filtered ~= 0)) max(LAT_O_filtered(LAT_O_filtered ~= 0))];
fprintf('Current LAT range: [%.2f - %.2f] ms\n', current_limits(1), current_limits(2));
user_input = input('Enter LAT color limits as [min max] (or press Enter for auto): ');
if isempty(user_input)
    caxis(current_limits);
else
    caxis(user_input);
end

hBar = colorbar('eastoutside');
ylabel(hBar, 'Local Activation Time [ms]', 'FontSize', 14);
title('Optical - Local Activation Time', 'FontSize', 16);
axis off;

% Calculate and display LAT statistics
fprintf('\n--- Optical LAT Statistics ---\n');
nonzero_LAT = LAT_O_filtered(LAT_O_filtered ~= 0);
if ~isempty(nonzero_LAT)
    fprintf('Mean LAT:      %.2f ± %.2f ms\n', mean(nonzero_LAT), std(nonzero_LAT));
    fprintf('Min LAT:       %.2f ms\n', min(nonzero_LAT));
    fprintf('Max LAT:       %.2f ms\n', max(nonzero_LAT));
    fprintf('Median LAT:    %.2f ms\n', median(nonzero_LAT));
    fprintf('Range:         %.2f ms\n', max(nonzero_LAT) - min(nonzero_LAT));
    fprintf('Number of points: %d\n', numel(nonzero_LAT));
else
    fprintf('No valid LAT data found.\n');
end

% Store results
OpticalResults.LAT = LAT_O_filtered;
OpticalResults.LAT_stats = struct(...
    'mean', mean(nonzero_LAT), ...
    'std', std(nonzero_LAT), ...
    'min', min(nonzero_LAT), ...
    'max', max(nonzero_LAT), ...
    'median', median(nonzero_LAT), ...
    'range', max(nonzero_LAT) - min(nonzero_LAT), ...
    'num_points', numel(nonzero_LAT));

fprintf('Optical LAT analysis completed.\n');


%% Section 2: Optical CV Analysis
fprintf('\n--- Section 2: Optical CV Analysis ---\n');

% Parameters for optical CV
cv_radius_optical = 5;
optical_space_scale = 0.33; % mm per pixel

% Calculate conduction velocity
fprintf('Calculating optical conduction velocity...\n');
AvgSpeed_O = zeros(size(LAT_O_filtered));
StdSpeed_O = zeros(size(LAT_O_filtered));
Angle_O = zeros(size(LAT_O_filtered));

for i = 1:size(LAT_O_filtered,1)
    for j = 1:size(LAT_O_filtered,2)
        if max(max(squeeze(LAT_O_filtered(i,j,:)))) ~= 0
            [AvgSpeed_O(i,j), StdSpeed_O(i,j), Angle_O(i,j)] = ...
                CV_CircleMethod(LAT_O_filtered, cv_radius_optical, i, j, optical_space_scale);
        end
    end
end

% Filter CV data
AvgSpeed_O_filtered = SpatTemp_Filtering(AvgSpeed_O, 5, 0, 'GPU');

% Plot CV map
fprintf('Plotting optical CV map...\n');
figure('color', 'white', 'Position', [40 40 600 600]);
C = jet(256);
C(1,1:3) = [1 1 1];
J = imrotate(AvgSpeed_O_filtered, 90);

% Ask user for color limits
current_cv_limits = [min(AvgSpeed_O_filtered(AvgSpeed_O_filtered ~= 0)) max(AvgSpeed_O_filtered(AvgSpeed_O_filtered ~= 0))];
fprintf('Current CV range: [%.2f - %.2f] cm/s\n', current_cv_limits(1), current_cv_limits(2));
user_input = input('Enter CV color limits as [min max] (or press Enter for auto): ');
if isempty(user_input)
    imagesc(J, current_cv_limits);
else
    imagesc(J, user_input);
end

colormap(C);
hBar = colorbar('eastoutside');
ylabel(hBar, 'Conduction Velocity [cm/s]', 'FontSize', 14);
title('Optical - Conduction Velocity', 'FontSize', 16);
axis off;

% Calculate and display CV statistics
fprintf('\n--- Optical CV Statistics ---\n');
nonzero_CV = AvgSpeed_O_filtered(AvgSpeed_O_filtered ~= 0);
if ~isempty(nonzero_CV)
    fprintf('Mean CV:       %.2f ± %.2f cm/s\n', mean(nonzero_CV), std(nonzero_CV));
    fprintf('Min CV:        %.2f cm/s\n', min(nonzero_CV));
    fprintf('Max CV:        %.2f cm/s\n', max(nonzero_CV));
    fprintf('Median CV:     %.2f cm/s\n', median(nonzero_CV));
    fprintf('Range:         %.2f cm/s\n', max(nonzero_CV) - min(nonzero_CV));
    fprintf('Number of points: %d\n', numel(nonzero_CV));
else
    fprintf('No valid CV data found.\n');
end

% Store results
OpticalResults.CV = AvgSpeed_O_filtered;
OpticalResults.CV_stats = struct(...
    'mean', mean(nonzero_CV), ...
    'std', std(nonzero_CV), ...
    'min', min(nonzero_CV), ...
    'max', max(nonzero_CV), ...
    'median', median(nonzero_CV), ...
    'range', max(nonzero_CV) - min(nonzero_CV), ...
    'num_points', numel(nonzero_CV));

fprintf('Optical CV analysis completed.\n');

%% Section 3: Electrical LAT Analysis
fprintf('\n--- Section 3: Electrical LAT Analysis ---\n');

% Parameters for electrical LAT
electrical_lim1 = 1;   % Start sample for electrical  
electrical_lim2 = 18000; % End sample for electrical
debug_electrical = 0;

% Data preparation
Data_E = D_SYNC.EL; % Adjust based on your data structure

% Calculate LAT for electrodes
fprintf('Calculating electrical LAT...\n');
LAT_E_raw = find_LAT_diff(Data_E, Fs, electrical_lim1, electrical_lim2, debug_electrical);

% Plot LAT maps for each configuration
cases = {'MEA1', 'MEA2', 'MEA3', 'TANK'};
ElectricalResults = struct();

for i = 1:length(cases)
    case_name = cases{i};
    fprintf('Processing %s...\n', case_name);
    
    % Plot and get LAT matrix
    LAT_matrix = plot_electric_LAT(LAT_E_raw, [0 40], 1, i);
    ElectricalResults.(case_name).LAT_matrix = LAT_matrix;
    
    % Calculate statistics
    nonzero_LAT = LAT_matrix(LAT_matrix ~= 0);
    if ~isempty(nonzero_LAT)
        ElectricalResults.(case_name).LAT_stats = struct(...
            'mean', mean(nonzero_LAT), ...
            'std', std(nonzero_LAT), ...
            'min', min(nonzero_LAT), ...
            'max', max(nonzero_LAT), ...
            'median', median(nonzero_LAT), ...
            'range', max(nonzero_LAT) - min(nonzero_LAT), ...
            'num_points', numel(nonzero_LAT));
        
        % Display statistics
        fprintf('\n--- %s LAT Statistics ---\n', case_name);
        fprintf('Mean LAT:      %.2f ± %.2f ms\n', mean(nonzero_LAT), std(nonzero_LAT));
        fprintf('Min LAT:       %.2f ms\n', min(nonzero_LAT));
        fprintf('Max LAT:       %.2f ms\n', max(nonzero_LAT));
        fprintf('Median LAT:    %.2f ms\n', median(nonzero_LAT));
        fprintf('Range:         %.2f ms\n', max(nonzero_LAT) - min(nonzero_LAT));
        fprintf('Number of points: %d\n', numel(nonzero_LAT));
    else
        fprintf('No valid LAT data found for %s.\n', case_name);
    end
end

% Find minimum LAT values for each MEA
fprintf('\n--- Minimum LAT Values ---\n');
% MEA1
mea1_id = [1:11, 14:16];
min_mea1 = min(LAT_E_raw(mea1_id));
el_min_mea1 = mea1_id(LAT_E_raw(mea1_id) == min_mea1);
fprintf('MEA1 - Min LAT: %.2f ms at electrode(s): %s\n', min_mea1, mat2str(el_min_mea1));

% MEA2
mea2_id = [17:32];
min_mea2 = min(LAT_E_raw(mea2_id));
el_min_mea2 = mea2_id(LAT_E_raw(mea2_id) == min_mea2);
fprintf('MEA2 - Min LAT: %.2f ms at electrode(s): %s\n', min_mea2, mat2str(el_min_mea2));

% MEA3
mea3_id = [65:79];
min_mea3 = min(LAT_E_raw(mea3_id));
el_min_mea3 = mea3_id(LAT_E_raw(mea3_id) == min_mea3);
fprintf('MEA3 - Min LAT: %.2f ms at electrode(s): %s\n', min_mea3, mat2str(el_min_mea3));

fprintf('Electrical LAT analysis completed.\n');

%% Section 4: Electrical CV Analysis
fprintf('\n--- Section 4: Electrical CV Analysis ---\n');

% Parameters for electrical CV
cv_radius_electrical = 3;
electrical_space_scale = 0.1;  % mm per pixel (MEA)
tank_space_scale = 0.75;       % mm per pixel (TANK)
interp_factor = 0.1;

for i = 1:length(cases)
    case_name = cases{i};
    fprintf('Calculating CV for %s...\n', case_name);
    
    LAT_matrix = ElectricalResults.(case_name).LAT_matrix;
    
    % Determine space scale based on case
    if strcmp(case_name, 'TANK')
        space_scale = tank_space_scale;
    else
        space_scale = electrical_space_scale;
    end
    
    % Interpolate for better CV calculation
    [X, Y] = meshgrid(1:size(LAT_matrix, 2), 1:size(LAT_matrix, 1));
    [Xq, Yq] = meshgrid(1:interp_factor:size(LAT_matrix, 2), 1:interp_factor:size(LAT_matrix, 1));
    LAT_interp = interp2(X, Y, LAT_matrix, Xq, Yq, 'linear');
    
    % Calculate CV
    AvgSpeed_E = zeros(size(LAT_interp));
    r = cv_radius_electrical;
    
    for ii = (r + 1):(size(LAT_interp, 1) - (r + 1))
        for jj = (r + 1):(size(LAT_interp, 2) - (r + 1))
            if max(max(squeeze(LAT_interp(ii, jj, :)))) ~= 0
                [AvgSpeed_E(ii, jj), ~, ~] = CV_CircleMethod(LAT_interp, r, ii, jj, space_scale);
            end
        end
    end
    
    % Filter CV data
    AvgSpeed_E_filtered = SpatTemp_Filtering(AvgSpeed_E, 5, 0, 'GPU');
    ElectricalResults.(case_name).CV_matrix = AvgSpeed_E_filtered;
    
    % Plot CV map
    figure('color', 'white', 'Position', [40 40 600 600]);
    C = jet(256);
    C(1,1:3) = [1 1 1];
    J = imrotate(AvgSpeed_E_filtered, 90);
    
    % Ask user for color limits
    current_cv_limits = [min(AvgSpeed_E_filtered(AvgSpeed_E_filtered ~= 0)) max(AvgSpeed_E_filtered(AvgSpeed_E_filtered ~= 0))];
    fprintf('Current CV range for %s: [%.2f - %.2f] cm/s\n', case_name, current_cv_limits(1), current_cv_limits(2));
    user_input = input(['Enter CV color limits for ' case_name ' as [min max] (or press Enter for auto): ']);
    if isempty(user_input)
        imagesc(J, current_cv_limits);
    else
        imagesc(J, user_input);
    end
    
    colormap(C);
    hBar = colorbar('eastoutside');
    ylabel(hBar, 'Conduction Velocity [cm/s]', 'FontSize', 14);
    title(['Electrical ' case_name ' - Conduction Velocity'], 'FontSize', 16);
    axis off;
    
    % Calculate and display CV statistics
    nonzero_CV = AvgSpeed_E_filtered(AvgSpeed_E_filtered ~= 0);
    if ~isempty(nonzero_CV)
        ElectricalResults.(case_name).CV_stats = struct(...
            'mean', mean(nonzero_CV), ...
            'std', std(nonzero_CV), ...
            'min', min(nonzero_CV), ...
            'max', max(nonzero_CV), ...
            'median', median(nonzero_CV), ...
            'range', max(nonzero_CV) - min(nonzero_CV), ...
            'num_points', numel(nonzero_CV));
        
        fprintf('\n--- %s CV Statistics ---\n', case_name);
        fprintf('Mean CV:       %.2f ± %.2f cm/s\n', mean(nonzero_CV), std(nonzero_CV));
        fprintf('Min CV:        %.2f cm/s\n', min(nonzero_CV));
        fprintf('Max CV:        %.2f cm/s\n', max(nonzero_CV));
        fprintf('Median CV:     %.2f cm/s\n', median(nonzero_CV));
        fprintf('Range:         %.2f cm/s\n', max(nonzero_CV) - min(nonzero_CV));
        fprintf('Number of points: %d\n', numel(nonzero_CV));
    else
        fprintf('No valid CV data found for %s.\n', case_name);
    end
end

fprintf('Electrical CV analysis completed.\n');

%% Section 5: Summary and Results Export
fprintf('\n--- Section 5: Summary and Results Export ---\n');

% Create comprehensive results table
fprintf('Creating results table...\n');

% Optical results
if exist('OpticalResults', 'var')
    OpticalTable = table(...
        {'Optical'}, ...
        OpticalResults.LAT_stats.mean, ...
        OpticalResults.LAT_stats.std, ...
        OpticalResults.LAT_stats.min, ...
        OpticalResults.LAT_stats.max, ...
        OpticalResults.LAT_stats.range, ...
        OpticalResults.CV_stats.mean, ...
        OpticalResults.CV_stats.std, ...
        OpticalResults.CV_stats.min, ...
        OpticalResults.CV_stats.max, ...
        OpticalResults.LAT_stats.num_points, ...
        'VariableNames', {'Analysis_Type', 'LAT_Mean_ms', 'LAT_Std_ms', 'LAT_Min_ms', 'LAT_Max_ms', 'LAT_Range_ms', ...
        'CV_Mean_cm_s', 'CV_Std_cm_s', 'CV_Min_cm_s', 'CV_Max_cm_s', 'Number_Points'});
end

% Electrical results
ElectricalTables = [];
for i = 1:length(cases)
    case_name = cases{i};
    if isfield(ElectricalResults, case_name) && isfield(ElectricalResults.(case_name), 'LAT_stats')
        case_table = table(...
            {['Electrical ' case_name]}, ...
            ElectricalResults.(case_name).LAT_stats.mean, ...
            ElectricalResults.(case_name).LAT_stats.std, ...
            ElectricalResults.(case_name).LAT_stats.min, ...
            ElectricalResults.(case_name).LAT_stats.max, ...
            ElectricalResults.(case_name).LAT_stats.range, ...
            ElectricalResults.(case_name).CV_stats.mean, ...
            ElectricalResults.(case_name).CV_stats.std, ...
            ElectricalResults.(case_name).CV_stats.min, ...
            ElectricalResults.(case_name).CV_stats.max, ...
            ElectricalResults.(case_name).LAT_stats.num_points, ...
            'VariableNames', {'Analysis_Type', 'LAT_Mean_ms', 'LAT_Std_ms', 'LAT_Min_ms', 'LAT_Max_ms', 'LAT_Range_ms', ...
            'CV_Mean_cm_s', 'CV_Std_cm_s', 'CV_Min_cm_s', 'CV_Max_cm_s', 'Number_Points'});
        
        if isempty(ElectricalTables)
            ElectricalTables = case_table;
        else
            ElectricalTables = [ElectricalTables; case_table];
        end
    end
end

% Combine all results
if exist('OpticalTable', 'var') && ~isempty(ElectricalTables)
    FinalResultsTable = [OpticalTable; ElectricalTables];
elseif exist('OpticalTable', 'var')
    FinalResultsTable = OpticalTable;
else
    FinalResultsTable = ElectricalTables;
end

% Display final table
fprintf('\n=== FINAL RESULTS SUMMARY ===\n');
disp(FinalResultsTable);

% Save results
save_choice = input('Save results to Excel? (y/n): ', 's');
if strcmpi(save_choice, 'y') || strcmpi(save_choice, 'yes')
    filename = 'LAT_CV_Analysis_Results.xlsx';
    writetable(FinalResultsTable, filename);
    fprintf('Results saved to: %s\n', filename);
end

% Save workspace
save_choice = input('Save workspace? (y/n): ', 's');
if strcmpi(save_choice, 'y') || strcmpi(save_choice, 'yes')
    save('LAT_CV_Analysis_Workspace.mat', 'OpticalResults', 'ElectricalResults', 'FinalResultsTable');
    fprintf('Workspace saved to: LAT_CV_Analysis_Workspace.mat\n');
end

fprintf('\n=== ANALYSIS COMPLETE ===\n');