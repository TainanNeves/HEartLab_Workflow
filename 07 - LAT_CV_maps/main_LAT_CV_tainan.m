%% LAT and CV Maps - Optical & Electrical
clear; clc;


%% Loading Data
% load(""); % Load Synchronized data
Fsampling = 4000;
% --- Optical Data & Parameters ---
Data_O = D_SYNC.CAM1;
ROI_O = D_SYNC.ROI.ROI_1;
optical_space_scale = 0.33; % mm per pixel (used for CV calculation)


%% Time Selection & Preview - Optical
optical_lim1 = 9460;
optical_lim2 = 9870; % Adjust to cover 1 or 2 activation cycles
Data_temp_O = Data_O(:,:,optical_lim1:optical_lim2);
Background = squeeze(Data_O(:,:,optical_lim1));
pick_up_a_trace(Background, Data_temp_O, 1);


%% LAT Calculation - Optical (Optical Analysis Part 1)
linear_fit_length = 15;
PCL = 200;
debug_LAT = 0;
fprintf('\nCalculating Optical LAT...\n');

% 1. Calculate LAT for each pixel
LAT_O = zeros(size(Data_temp_O, 1), size(Data_temp_O, 2));
for i = 1:size(Data_temp_O,1)
    for j = 1:size(Data_temp_O,2)
        signal = squeeze(Data_temp_O(i,j,:));
        if max(signal) ~= 0
            % Call to external function to find LAT
            LAT_O(i,j) = find_LAT_linearFit_1D(signal, Fsampling, linear_fit_length, PCL, debug_LAT);
        end
    end
end

% 2. Subtract minimum (Normalization)
LAT_O(LAT_O < 0) = 0;
min_LAT = min(LAT_O(LAT_O ~= 0));
if ~isempty(min_LAT) && ~isnan(min_LAT)
    LAT_O(LAT_O ~= 0) = LAT_O(LAT_O ~= 0) - min_LAT;
end

% 3. Filter LAT data
% Call to external function for spatial filtering
LAT_O_filtered = SpatTemp_Filtering(LAT_O, 3, 0, 'GPU');
OpticalResults.LAT = LAT_O_filtered;


%% LAT Map Plot - Optical (Optical Analysis Part 2)
levelStep = 20;
Title = 'Optical - Local Activation Time';
f1 = figure('color', 'white', 'Position', [40 40 600 600]);
C = parula(256);
C(1,1:3) = [1 1 1]; % White for background
C(256,1:3) = [0.5, 0.5, 0.5]; % Gray for the max value

J = imrotate(LAT_O_filtered, 90);
contourf(flipud(J), levelStep, 'LineStyle', 'none'); % Use contourf for smoother map
colormap(C);

% Get color limits from user
nonzero_LAT_all = LAT_O_filtered(LAT_O_filtered ~= 0);
current_limits = [min(nonzero_LAT_all) max(nonzero_LAT_all)];
fprintf('\nLAT Map Limits - Current range: [%.2f - %.2f] ms\n', current_limits(1), current_limits(2));
user_input = input('Enter LAT color limits as [min max] (or press Enter for auto): ');
if isempty(user_input)
    caxis(current_limits);
else
    caxis(user_input);
end
axis off;
hBar = colorbar('eastoutside');
ylabel(hBar, 'Local Activation Time [ms]', 'FontSize', 14);
title(Title, 'FontSize', 16);


%% In case of need to manual corrections - Optical
fprintf('\n=== Manual Corrections - Optical LAT ===\n');
% Use this section to manually find and replace erroneous LAT values.
% Corrections are applied to LAT_O_filtered, which is used for subsequent analysis.

% Example: Find zero/near-zero artifacts that shouldn't be zero
% find_value = 0.5; 
% tolerancia = 1e-10;
% indices = find(abs(LAT_O_filtered - find_value) < tolerancia);
% if ~isempty(indices)
%     fprintf('Found %d pixels with value %.2f (within tolerance). Setting them to NaN.\n', numel(indices), find_value);
%     LAT_O_filtered(indices) = NaN; % Set to NaN to exclude from statistics/CV calculation
%     OpticalResults.LAT = LAT_O_filtered; % Update the stored result matrix
% end

fprintf('Manual correction section complete. Proceeding to Statistics/CV calculation.\n');


%% LAT Statistics (ROI-Based) - Optical (Optical Analysis Part 3)
% Display the LAT map and let the user define the ROI
figure(f1); % Bring LAT Map figure to front
title('Optical - LAT: Select ROI for Statistical Analysis');
roi_lat_O = roipoly;

% Apply the ROI
LAT_O_roi = LAT_O_filtered .* roi_lat_O;
nonzero_LAT_O = LAT_O_roi(LAT_O_roi ~= 0 & ~isnan(LAT_O_roi)); % Exclude NaN values if manual correction was used
OpticalResults.roi_lat = roi_lat_O;

if ~isempty(nonzero_LAT_O)
    % Calculate and store LAT Statistics
    LAT_mean_O = mean(nonzero_LAT_O);
    LAT_std_O = std(nonzero_LAT_O);
    LAT_max_O = max(nonzero_LAT_O);
    LAT_mode_O = mode(nonzero_LAT_O); 
    LAT_var_O = var(nonzero_LAT_O); 
    LAT_num_O = numel(nonzero_LAT_O);

    OpticalResults.LAT_stats = struct(...
        'mean', LAT_mean_O, 'std', LAT_std_O, 'max', LAT_max_O, ...
        'mode', LAT_mode_O, 'var', LAT_var_O, 'num_points', LAT_num_O);

    fprintf('\n=== Optical LAT Statistics (in ROI) ===\n');
    disp(['Average LAT:      ', num2str(LAT_mean_O, '%.2f'), ' ms']);
    disp(['Mode LAT:         ', num2str(LAT_mode_O, '%.2f'), ' ms']);
    disp(['Max LAT:          ', num2str(LAT_max_O, '%.2f'), ' ms']);
    disp(['Std Deviation:    ', num2str(LAT_std_O, '%.2f'), ' ms']);
    disp(['Variance:         ', num2str(LAT_var_O, '%.2f')]);
    disp(['Number of pixels: ', num2str(LAT_num_O)]);
else
    fprintf('No valid LAT data found in the selected ROI.\n');
    LAT_mean_O = NaN; LAT_std_O = NaN; LAT_max_O = NaN; 
    LAT_mode_O = NaN; LAT_var_O = NaN; LAT_num_O = 0;
end


%% Specific ROI Plot with Min Subtraction - Optical
fprintf('\n=== Specific ROI Plot with Re-Normalization ===\n');
% This section allows plotting a specific ROI and recalculating the minimum LAT 
% inside that ROI for better visualization/comparison.

% 1. Display the filtered LAT map to select the new ROI
figure();
imagesc(imrotate(LAT_O_filtered, 90));
colormap(parula); colorbar;
title('Select SPECIFIC ROI for Subtraction/Re-Normalization Plot');
roi_specific = roipoly;

% 2. Apply the specific ROI
LAT_specific = LAT_O_filtered .* roi_specific;

% 3. Subtract the minimum LAT *within* this specific ROI
nonzero_specific = LAT_specific(LAT_specific ~= 0 & ~isnan(LAT_specific));
if ~isempty(nonzero_specific)
    min_LAT_specific = min(nonzero_specific);
    
    % Create the re-normalized map
    LAT_renormalized = LAT_specific; % Copy the specific ROI data
    LAT_renormalized(LAT_renormalized ~= 0) = LAT_renormalized(LAT_renormalized ~= 0) - min_LAT_specific;
    LAT_renormalized(LAT_renormalized < 0) = 0; % Ensure no negative values after subtraction

    % 4. Plot the re-normalized map
    figure('color', 'white', 'Position', [1250 40 600 600]);
    Title_spec = ['Optical - Specific ROI (Min Subtracted: ' num2str(min_LAT_specific, '%.2f') 'ms)'];
    C = parula(256); C(1,1:3) = [1 1 1];
    J_spec = imrotate(LAT_renormalized, 90);
    
    % Use contourf for smoothing the plot
    contourf(flipud(J_spec), levelStep, 'LineStyle', 'none'); 
    colormap(C);
    
    % Set color limits based on the new range
    new_max = max(LAT_renormalized(:));
    caxis([0 new_max]);
    
    axis off;
    hBar_spec = colorbar('eastoutside');
    ylabel(hBar_spec, 'Local Activation Time [ms]', 'FontSize', 14);
    title(Title_spec, 'FontSize', 16);
    fprintf('Re-normalized LAT map plotted. New range: [0 - %.2f] ms\n', new_max);
else
    fprintf('No valid data found in the specific ROI. Plot skipped.\n');
end

%% CV Calculation - Optical (Optical Analysis Part 4)
cv_radius_optical = 5; % in pixels
fprintf('\nCalculating Optical CV (Circle Method, radius %d)...\n', cv_radius_optical);

AvgSpeed_O = zeros(size(LAT_O_filtered));
for i = 1:size(LAT_O_filtered,1)
    for j = 1:size(LAT_O_filtered,2)
        % Only calculate if point is inside the LAT ROI and is a valid value
        if roi_lat_O(i,j) == 1 && LAT_O_filtered(i,j) ~= 0 && ~isnan(LAT_O_filtered(i,j))
            % Call to external function to calculate CV
            [AvgSpeed_O(i,j), ~, ~] = ...
                CV_CircleMethod(LAT_O_filtered, cv_radius_optical, i, j, optical_space_scale);
        end
    end
end

% Filter CV data
AvgSpeed_O_filtered = SpatTemp_Filtering(AvgSpeed_O, 5, 0, 'GPU');
OpticalResults.CV = AvgSpeed_O_filtered;

%% CV Map Plot - Optical (Optical Analysis Part 5)
Title = 'Optical - Conduction Velocity';
f2 = figure('color', 'white', 'Position', [650 40 600 600]);
C = jet(256);
C(1,1:3) = [1 1 1];

J = imrotate(AvgSpeed_O_filtered, 90);
imagesc(J);
colormap(C);

% Get color limits from user
nonzero_CV_all = AvgSpeed_O_filtered(AvgSpeed_O_filtered ~= 0 & ~isnan(AvgSpeed_O_filtered));
current_cv_limits = [min(nonzero_CV_all) max(nonzero_CV_all)];
fprintf('\nCV Map Limits - Current range: [%.2f - %.2f] cm/s\n', current_cv_limits(1), current_cv_limits(2));
user_input_cv = input('Enter CV color limits as [min max] (or press Enter for auto): ');
if isempty(user_input_cv)
    caxis(current_cv_limits);
else
    caxis(user_input_cv);
end

axis off;
hBar = colorbar('eastoutside');
ylabel(hBar, 'Conduction Velocity [cm/s]', 'FontSize', 14);
title(Title, 'FontSize', 16);

%% CV Statistics - Optical (Optical Analysis Part 6)
% Apply the same ROI used for LAT statistics
CV_O_roi = AvgSpeed_O_filtered .* roi_lat_O;
nonzero_CV_O = CV_O_roi(CV_O_roi ~= 0 & ~isnan(CV_O_roi));

if ~isempty(nonzero_CV_O)
    % Calculate and store CV Statistics
    CV_mean_O = mean(nonzero_CV_O);
    CV_std_O = std(nonzero_CV_O);
    CV_max_O = max(nonzero_CV_O);
    CV_mode_O = mode(nonzero_CV_O);
    CV_var_O = var(nonzero_CV_O); 
    CV_num_O = numel(nonzero_CV_O);

    OpticalResults.CV_stats = struct(...
        'mean', CV_mean_O, 'std', CV_std_O, 'max', CV_max_O, ...
        'mode', CV_mode_O, 'var', CV_var_O, 'num_points', CV_num_O);

    fprintf('\n=== Optical CV Statistics (in LAT ROI) ===\n');
    disp(['Average CV:      ', num2str(CV_mean_O, '%.2f'), ' cm/s']);
    disp(['Mode CV:         ', num2str(CV_mode_O, '%.2f'), ' cm/s']);
    disp(['Max CV:          ', num2str(CV_max_O, '%.2f'), ' cm/s']);
    disp(['Std Deviation:   ', num2str(CV_std_O, '%.2f'), ' cm/s']);
    disp(['Variance:        ', num2str(CV_var_O, '%.2f')]);
    disp(['Number of pixels: ', num2str(CV_num_O)]);
else
    fprintf('No valid CV data found in the selected ROI.\n');
    CV_mean_O = NaN; CV_std_O = NaN; CV_max_O = NaN; 
    CV_mode_O = NaN; CV_var_O = NaN; CV_num_O = 0;
end

%% Optical Results Table (Part 7)
% Create table with all optical statistics
Parameter_O = {'LAT (ms)'; 'CV (cm/s)'};
Max_Value_O = [LAT_max_O; CV_max_O];
Mean_Value_O = [LAT_mean_O; CV_mean_O];
Mode_Value_O = [LAT_mode_O; CV_mode_O];
Std_Deviation_O = [LAT_std_O; CV_std_O];
Variance_O = [LAT_var_O; CV_var_O];
Number_Points_O = [LAT_num_O; CV_num_O];
OpticalResults.ResultsTable = table(Parameter_O, Max_Value_O, Mean_Value_O, Mode_Value_O, Std_Deviation_O, Variance_O, Number_Points_O, ...
    'RowNames', {'LAT_O', 'CV_O'});

disp(' ');
disp('=== OPTICAL RESULTS TABLE ===');
disp(OpticalResults.ResultsTable);

%% Summary Visualizations - Optical (Part 8)
figure('Position', [100 100 1200 400]);

% Subplot 1: LAT Histogram
subplot(1,2,1);
if ~isempty(nonzero_LAT_O)
    histogram(nonzero_LAT_O, 20, 'FaceColor', 'blue', 'FaceAlpha', 0.7);
    xlabel('Local Activation Time (ms)');
    ylabel('Count');
    title('Optical LAT Distribution');
    grid on;
end

% Subplot 2: CV Histogram
subplot(1,2,2);
if ~isempty(nonzero_CV_O)
    histogram(nonzero_CV_O, 20, 'FaceColor', 'red', 'FaceAlpha', 0.7);
    xlabel('Conduction Velocity (cm/s)');
    ylabel('Count');
    title('Optical CV Distribution');
    grid on;
end
sgtitle('Optical Analysis Histograms');

%% Save workspace variables for future use - Optical
save('O_LAT_CV_CAM.mat', 'Fsampling', 'optical_lim1', 'optical_lim2', 'optical_space_scale', ...
    'LAT_O_filtered', 'AvgSpeed_O_filtered', 'OpticalResults', '-v7.3');
disp(' ');
disp('=== OPTICAL ANALYSIS COMPLETE ===');
disp('All results have been saved to: O_LAT_CV_CAM.mat');

% --- OPTICAL ANALYSIS ENDS HERE ---
% --- ELECTRICAL ANALYSIS STARTS HERE ---
% clear; % NOTE: Do NOT use clear here, as it removes the loaded data (D_SYNC) and Optical results.
% Instead, use clearvars -except if needed, but keeping the workspace is safer for final table compilation.


%% 
%% 
%%


%% Loading Data
% --- Electrical Data & Parameters ---
Data_E = D_SYNC.EL; % Assuming D_SYNC.EL holds the raw electrical signals
electrical_space_scale = 0.1;  % mm per pixel (MEA)
tank_space_scale = 0.75;       % mm per pixel (TANK)
cases = {'MEA1', 'MEA2', 'MEA3', 'TANK'};
ElectricalResults = struct();
%% Time Selection & Raw LAT Calculation - Electrical
electrical_lim1 = 1;
electrical_lim2 = 18000; % Start and end samples

% The raw LAT calculation is performed once on the entire electrical signal set.
fprintf('\nCalculating Electrical LAT from sample %d to %d...\n', electrical_lim1, electrical_lim2);
% Call to external function to calculate raw electrical LAT
LAT_E_raw = find_LAT_diff(Data_E, Fsampling, electrical_lim1, electrical_lim2, 0);
ElectricalResults.LAT_raw = LAT_E_raw;

%% LAT Map, Statistics, and CV Calculation - Electrical (Looping through cases)
AllResultsTables_E = cell(length(cases), 1);
for i = 1:length(cases)
    case_name = cases{i};
    fprintf('\n=== Analyzing Electrical Case: %s ===\n', case_name);
    
    % --- 1. LAT Map Plot ---
    % Call to external function to plot and return the LAT matrix for the case
    LAT_matrix = plot_electric_LAT(LAT_E_raw, [0 40], 1, i, case_name);
    ElectricalResults.(case_name).LAT_matrix = LAT_matrix;
    
    % --- 2. LAT Statistics (ROI-Based) ---
    figure(gcf); % Bring the current LAT map figure to front
    title([case_name ' LAT: Select ROI for Statistical Analysis']);
    roi_lat_E = roipoly;
    
    LAT_E_roi = LAT_matrix .* roi_lat_E;
    nonzero_LAT_E = LAT_E_roi(LAT_E_roi ~= 0);
    ElectricalResults.(case_name).roi_lat = roi_lat_E;
    ElectricalResults.(case_name).LAT_data = nonzero_LAT_E;

    if ~isempty(nonzero_LAT_E)
        LAT_mean_E = mean(nonzero_LAT_E);
        LAT_std_E = std(nonzero_LAT_E);
        LAT_max_E = max(nonzero_LAT_E);
        LAT_mode_E = mode(nonzero_LAT_E);
        LAT_var_E = var(nonzero_LAT_E);
        LAT_num_E = numel(nonzero_LAT_E);
        
        LAT_stats = struct('mean', LAT_mean_E, 'std', LAT_std_E, 'max', LAT_max_E, ...
                           'mode', LAT_mode_E, 'var', LAT_var_E, 'num_points', LAT_num_E);
        ElectricalResults.(case_name).LAT_stats = LAT_stats;
        
        fprintf('--- %s LAT Statistics (in ROI) ---\n', case_name);
        disp(['Average LAT:      ', num2str(LAT_mean_E, '%.2f'), ' ms']);
        disp(['Number of points: ', num2str(LAT_num_E)]);
    else
        fprintf('No valid LAT data found in ROI for %s.\n', case_name);
        LAT_stats = struct('mean', NaN, 'std', NaN, 'max', NaN, 'mode', NaN, 'var', NaN, 'num_points', 0);
        ElectricalResults.(case_name).LAT_stats = LAT_stats;
        LAT_num_E = 0; 
    end

    % --- 3. CV Calculation ---
    cv_radius_electrical = 3;
    interp_factor = 0.1;
    
    if strcmp(case_name, 'TANK')
        space_scale = tank_space_scale;
    else
        space_scale = electrical_space_scale;
    end
    
    % Interpolate LAT matrix
    [X, Y] = meshgrid(1:size(LAT_matrix, 2), 1:size(LAT_matrix, 1));
    [Xq, Yq] = meshgrid(1:interp_factor:size(LAT_matrix, 2), 1:interp_factor:size(LAT_matrix, 1));
    LAT_interp = interp2(X, Y, LAT_matrix, Xq, Yq, 'linear');
    
    % Calculate CV
    fprintf('Calculating CV for %s (interp factor 1/%g)...\n', case_name, interp_factor);
    AvgSpeed_E = zeros(size(LAT_interp));
    r = cv_radius_electrical;
    
    for ii = (r + 1):(size(LAT_interp, 1) - (r + 1))
        for jj = (r + 1):(size(LAT_interp, 2) - (r + 1))
            if LAT_interp(ii, jj) ~= 0
                [AvgSpeed_E(ii, jj), ~, ~] = CV_CircleMethod(LAT_interp, r, ii, jj, space_scale * interp_factor);
            end
        end
    end
    
    % Filter CV data
    AvgSpeed_E_filtered = SpatTemp_Filtering(AvgSpeed_E, 5, 0, 'GPU');
    ElectricalResults.(case_name).CV_matrix = AvgSpeed_E_filtered;

    % --- 4. CV Map Plot ---
    figure('color', 'white', 'Position', [40 40 600 600]);
    C = jet(256); C(1,1:3) = [1 1 1];
    if strcmp(case_name, 'TANK')
        J = AvgSpeed_E_filtered;
        imagesc(J);
        pbaspect([2 1 1]);
    else
        J = imrotate(AvgSpeed_E_filtered, 90);
        imagesc(J);
        axis equal;
    end
    colormap(C);
    axis off;
    hBar = colorbar('eastoutside');
    ylabel(hBar, 'Conduction Velocity [cm/s]', 'FontSize', 14);
    title(['Electrical ' case_name ' - Conduction Velocity'], 'FontSize', 16);
    
    % --- 5. CV Statistics ---
    % Downsample CV map to original electrode grid size for ROI application
    CV_map_original_size = imresize(AvgSpeed_E_filtered, size(LAT_matrix), 'bilinear');
    CV_E_roi = CV_map_original_size .* roi_lat_E;
    nonzero_CV_E = CV_E_roi(CV_E_roi ~= 0);
    ElectricalResults.(case_name).CV_data = nonzero_CV_E;

    if ~isempty(nonzero_CV_E)
        CV_mean_E = mean(nonzero_CV_E);
        CV_std_E = std(nonzero_CV_E);
        CV_max_E = max(nonzero_CV_E);
        CV_mode_E = mode(nonzero_CV_E);
        CV_var_E = var(nonzero_CV_E);
        CV_num_E = numel(nonzero_CV_E);
        
        CV_stats = struct('mean', CV_mean_E, 'std', CV_std_E, 'max', CV_max_E, ...
                          'mode', CV_mode_E, 'var', CV_var_E, 'num_points', CV_num_E);
        ElectricalResults.(case_name).CV_stats = CV_stats;
        
        fprintf('--- %s CV Statistics (in LAT ROI) ---\n', case_name);
        disp(['Average CV:      ', num2str(CV_mean_E, '%.2f'), ' cm/s']);
        disp(['Number of points: ', num2str(CV_num_E)]);
    else
        fprintf('No valid CV data found in ROI for %s.\n', case_name);
        CV_stats = struct('mean', NaN, 'std', NaN, 'max', NaN, 'mode', NaN, 'var', NaN, 'num_points', 0);
        ElectricalResults.(case_name).CV_stats = CV_stats;
        CV_num_E = 0;
    end
    
    % --- 6. Results Table per Case (DF-style) ---
    Parameter_E = {'LAT (ms)'; 'CV (cm/s)'};
    Max_Value_E = [LAT_stats.max; CV_stats.max];
    Mean_Value_E = [LAT_stats.mean; CV_stats.mean];
    Mode_Value_E = [LAT_stats.mode; CV_stats.mode];
    Std_Deviation_E = [LAT_stats.std; CV_stats.std];
    Variance_E = [LAT_stats.var; CV_stats.var];
    Number_Points_E = [LAT_num_E; CV_num_E];

    ResultsTable_E_Case = table(Parameter_E, Max_Value_E, Mean_Value_E, Mode_Value_E, Std_Deviation_E, Variance_E, Number_Points_E, ...
        'RowNames', {['LAT_' case_name], ['CV_' case_name]});
    
    AllResultsTables_E{i} = ResultsTable_E_Case;
    fprintf('\n=== %s - RESULTS TABLE ===\n', case_name);
    disp(ResultsTable_E_Case);
end

%% Final Results Table Compilation
% Combine Optical and Electrical Tables
FinalResultsTable = OpticalResults.ResultsTable;
for i = 1:length(cases)
    if ~isempty(AllResultsTables_E{i})
        FinalResultsTable = [FinalResultsTable; AllResultsTables_E{i}];
    end
end
ElectricalResults.AllTables = FinalResultsTable;

disp(' ');
disp('=== FINAL COMBINED RESULTS TABLE (LAT/CV) ===');
disp(FinalResultsTable);

%% Summary Visualizations - Electrical (Histograms)
figure('Position', [100 100 1200 800]);
total_cases = length(cases);
rows = 2;
cols = total_cases;

for i = 1:total_cases
    case_name = cases{i};
    if isfield(ElectricalResults, case_name) && isfield(ElectricalResults.(case_name), 'LAT_data')
        
        % Subplot: LAT Histogram
        subplot(rows, cols, i);
        if ~isempty(ElectricalResults.(case_name).LAT_data)
            histogram(ElectricalResults.(case_name).LAT_data, 10, 'FaceColor', 'blue', 'FaceAlpha', 0.7);
            title([case_name ' LAT (ms)']);
            xlabel('LAT (ms)');
            ylabel('Count');
            grid on;
        end
        
        % Subplot: CV Histogram
        subplot(rows, cols, i + cols);
        if isfield(ElectricalResults.(case_name), 'CV_data') && ~isempty(ElectricalResults.(case_name).CV_data)
            histogram(ElectricalResults.(case_name).CV_data, 10, 'FaceColor', 'red', 'FaceAlpha', 0.7);
            title([case_name ' CV (cm/s)']);
            xlabel('CV (cm/s)');
            ylabel('Count');
            grid on;
        end
    end
end
sgtitle('Electrical Analysis Histograms (LAT/CV)', 'FontSize', 16);

%% Save workspace variables for future use - Electrical
save('E_LAT_CV.mat', 'Fsampling', 'electrical_lim1', 'electrical_lim2', 'ElectricalResults', 'FinalResultsTable', '-v7.3');
disp(' ');
disp('=== ELECTRICAL ANALYSIS COMPLETE ===');
disp('All results have been saved to: E_LAT_CV.mat');