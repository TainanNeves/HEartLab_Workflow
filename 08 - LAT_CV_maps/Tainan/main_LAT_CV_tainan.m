%% LAT and CV Maps - Optical
clear; clc;


%% OPTIC ANALYSIS
%% Loading Preview Data
load("E:\Qualification\Analysis\E32F02R08\data\data_filtered_sync_E32_F02_R08.mat"); % Load Synchronized data


%% Loading Parameters Select Time Windows
Fsampling = 4000;
Data_O = D_SYNC.CAM3;
ROI_O = D_SYNC.ROI.ROI_3;
optical_space_scale = 0.33; % mm per pixel
% Preview Signal
Background = squeeze(Data_O(:,:,2000));
pick_up_a_trace(Background, Data_O,1);


%% Time Selection & Preview - Optical
optical_lim1 = 9642;
optical_lim2 = 10177;
Data_temp_O = Data_O(:,:,optical_lim1:optical_lim2);
Background = squeeze(Data_O(:,:,optical_lim1));
pick_up_a_trace(Background, Data_temp_O, 1);


%% LAT Calculation - Optical
linear_fit_length = 15;
PCL = 200;
debug_LAT = 1;

% Calculate LAT for each pixel
LAT_O = zeros(size(Data_temp_O, 1), size(Data_temp_O, 2));
for i = 1:size(Data_temp_O,1)
    for j = 1:size(Data_temp_O,2)
        signal = squeeze(Data_temp_O(i,j,:));
        if max(signal) ~= 0
            LAT_O(i,j) = find_LAT_linearFit_1D(signal, ...
                        Fsampling, linear_fit_length, ...
                        PCL, debug_LAT);
        end
    end
end

% Subtract minimum
LAT_O(LAT_O < 0) = 0;
min_LAT = min(LAT_O(LAT_O ~= 0));
if ~isempty(min_LAT) && ~isnan(min_LAT)
    LAT_O(LAT_O ~= 0) = LAT_O(LAT_O ~= 0) - min_LAT;
end

% Filter LAT data
LAT_O_filtered = SpatTemp_Filtering(LAT_O, 3, 0, 'GPU');
Results.LAT = LAT_O_filtered;


%% LAT Map Plot - Optical
levelStep = 20;
Title = 'Optical - LAT';
f1 = figure('color', 'white', 'Position', [40 40 600 600]);
C = parula(256);
% C(1,1:3) = [1 1 1]; % White for background
C(256,1:3) = [0.5, 0.5, 0.5]; % Gray for the max value
J = imrotate(LAT_O_filtered, 90);
contourf(flipud(J), levelStep);
colormap(C);

% Get color limits from non-zero LAT values
nonzero_LAT_all = LAT_O_filtered(LAT_O_filtered ~= 0);
current_limits = [min(nonzero_LAT_all) max(nonzero_LAT_all)];
% current_limits = []; % Use to specific limits
caxis(current_limits);

axis off;
hBar = colorbar('eastoutside');
ylabel(hBar, 'Local Activation Time [ms]', 'FontSize', 14);
title(Title, 'FontSize', 16);


%% In case of need to manual corrections - Optical
% Histogran to evaluate distribution
nonzero_LAT = LAT_O_filtered(LAT_O_filtered ~= 0);
figure('color', 'white');
histogram(nonzero_LAT);
title('Histogram of LAT', 'FontSize', 14);
xlabel('LAT [ms]', 'FontSize', 12);
ylabel('N of Pixels', 'FontSize', 12);
grid on;

% Configuring Parameters
find_value = 0.1; 
tolerancia = 0.1;
new_value = 0.5; % NaN to exclude

% Substitute exactly value
indices = find(abs(LAT_O_filtered - find_value) < tolerancia);
if ~isempty(indices)
    LAT_O_filtered(indices) = new_value;
    Results.LAT = LAT_O_filtered;
end

% Substitute Higher or Lower than
LAT_O_filtered(LAT_O_filtered < find_value) = new_value;

% Correct The Zero Normalization - Subtract minimum
LAT_O_filtered(LAT_O_filtered < 0) = 0;
min_LAT = min(LAT_O_filtered(LAT_O_filtered ~= 0));
if ~isempty(min_LAT) && ~isnan(min_LAT)
    LAT_O_filtered(LAT_O_filtered ~= 0) = ...
        LAT_O_filtered(LAT_O_filtered ~= 0) - min_LAT;
end


%% LAT Statistics - Optical
% Display the image and let the user define the ROI
figure();
imshow(LAT_O_filtered, C);
% imshow(Data_temp_O(:,:,1));
title('Select ROI for DF Analysis');
roi_lat_O = roipoly;

% Apply the ROI
LAT_O_roi = LAT_O_filtered .* roi_lat_O;
nonzero_LAT_O = LAT_O_roi(LAT_O_roi ~= 0 & ~isnan(LAT_O_roi));
Results.roi_lat = roi_lat_O;

% Calculate and store LAT Statistics
LAT_mean_O = mean(nonzero_LAT_O);
LAT_std_O = std(nonzero_LAT_O);
LAT_max_O = max(nonzero_LAT_O);
LAT_mode_O = mode(nonzero_LAT_O); 
LAT_var_O = var(nonzero_LAT_O); 
LAT_num_O = numel(nonzero_LAT_O);
% Displaying results
fprintf('\n=== Optical LAT Statistics (in ROI) ===\n');
disp(['Average LAT:      ', num2str(LAT_mean_O, '%.4f'), ' ms']);
disp(['Mode LAT:         ', num2str(LAT_mode_O, '%.4f'), ' ms']);
disp(['Max LAT:          ', num2str(LAT_max_O, '%.4f'), ' ms']);
disp(['Std Deviation:    ', num2str(LAT_std_O, '%.4f'), ' ms']);
disp(['Variance:         ', num2str(LAT_var_O, '%.4f')]);
disp(['Number of pixels: ', num2str(LAT_num_O)]);
% Creating Struct
Results.LAT_stats = struct(...
    'mean', LAT_mean_O, 'std', LAT_std_O, 'max', LAT_max_O, ...
    'mode', LAT_mode_O, 'var', LAT_var_O, 'num_points', LAT_num_O);


%% Specific ROI Plot with Min Subtraction - Optical
levelStep = 10;
% ROI selection
    % Same as the statistics
        roi_specific = roi_lat_O;
    % Select new ROI
        % figure();
        % imagesc(LAT_O_filtered);
        % colormap(parula); colorbar;
        % title('Select ROI');
        % roi_specific = roipoly; 

% Apply the specific ROI
LAT_specific = LAT_O_filtered .* roi_specific;
% Subtract the minimum LAT *within* this specific ROI
nonzero_specific = LAT_specific(LAT_specific ~= 0 & ~isnan(LAT_specific));
min_LAT_specific = min(nonzero_specific);
% Create the re-normalized map
LAT_renormalized = LAT_specific;
% Perform the normalization (zeroing)
LAT_renormalized(LAT_renormalized ~= 0) = LAT_renormalized(LAT_renormalized ~= 0) - min_LAT_specific;
LAT_renormalized(LAT_renormalized < 0) = 0; % Ensure no negative values after subtraction

% Plot
figure('color', 'white', 'Position', [1250 40 600 600]);
Title_spec = ['Optical - Specific ROI (Min Subtracted: ' num2str(min_LAT_specific, '%.2f') 'ms)'];
C = parula(256); C(1,1:3) = [1 1 1]; 
J_spec = imrotate(LAT_renormalized, 90);
contourf(flipud(J_spec), levelStep); 
colormap(C);
% Set color limits based on the new range
new_max = max(LAT_renormalized(:));
caxis([0 new_max]);
axis off;
hBar_spec = colorbar('eastoutside');
ylabel(hBar_spec, 'Local Activation Time [ms]', 'FontSize', 14);
title(Title_spec, 'FontSize', 16);
fprintf('Re-normalized LAT map plotted. New range: [0 - %.2f] ms\n', new_max);


%% CV Calculation - Optical
cv_radius_optical = 5; % in pixels

AvgSpeed_O = zeros(size(LAT_O_filtered));
for i = 1:size(LAT_O_filtered,1)
    for j = 1:size(LAT_O_filtered,2)
        if roi_lat_O(i,j) == 1 && LAT_O_filtered(i,j) ~= 0 && ~isnan(LAT_O_filtered(i,j))
            [AvgSpeed_O(i,j), ~, ~] = ...
                CV_CircleMethod(LAT_O_filtered, cv_radius_optical, ...
                                    i, j, optical_space_scale);
        end
    end
end

% Filter CV data
AvgSpeed_O_filtered = SpatTemp_Filtering(AvgSpeed_O, 5, 0, 'GPU');
Results.CV = AvgSpeed_O_filtered;


%% CV Map Plot - Optical
Title = 'Optical - CV';
f2 = figure('color', 'white', 'Position', [650 40 600 600]);
C = jet(256);
C(1,1:3) = [1 1 1];
J = imrotate(AvgSpeed_O_filtered, 90);
imagesc(J);
colormap(C);

% Calculate color limits from non-zero and non-NaN CV values
nonzero_CV_all = AvgSpeed_O_filtered(AvgSpeed_O_filtered ~= 0 & ~isnan(AvgSpeed_O_filtered));
current_cv_limits = [min(nonzero_CV_all) max(nonzero_CV_all)];
% current_cv_limits = []; % Specific limits
caxis(current_cv_limits);

axis off;
hBar = colorbar('eastoutside');
ylabel(hBar, 'Conduction Velocity [cm/s]', 'FontSize', 14);
title(Title, 'FontSize', 16);


%% CV Statistics - Optical
% Apply the same ROI used for LAT statistics
CV_O_roi = AvgSpeed_O_filtered .* roi_lat_O;
nonzero_CV_O = CV_O_roi(CV_O_roi ~= 0 & ~isnan(CV_O_roi) & ~isinf(CV_O_roi));

% Calculate CV Statistics
CV_mean_O = mean(nonzero_CV_O);
CV_std_O = std(nonzero_CV_O);
CV_max_O = max(nonzero_CV_O);
CV_mode_O = mode(nonzero_CV_O);
CV_var_O = var(nonzero_CV_O); 
CV_num_O = numel(nonzero_CV_O);
% Display CV Statistics
fprintf('\n=== Optical CV Statistics (in LAT ROI) ===\n');
disp(['Average CV:      ', num2str(CV_mean_O, '%.2f'), ' cm/s']);
disp(['Mode CV:         ', num2str(CV_mode_O, '%.2f'), ' cm/s']);
disp(['Max CV:          ', num2str(CV_max_O, '%.2f'), ' cm/s']);
disp(['Std Deviation:   ', num2str(CV_std_O, '%.2f'), ' cm/s']);
disp(['Variance:        ', num2str(CV_var_O, '%.2f')]);
disp(['Number of pixels: ', num2str(CV_num_O)]);
% Store CV Statistics
Results.CV_stats = struct(...
    'mean', CV_mean_O, 'std', CV_std_O, 'max', CV_max_O, ...
    'mode', CV_mode_O, 'var', CV_var_O, 'num_points', CV_num_O);


%% Optical Results Table
% Create table with all optical statistics
Parameter_O = {'LAT (ms)'; 'CV (cm/s)'};
Max_Value_O = [LAT_max_O; CV_max_O];
Mean_Value_O = [LAT_mean_O; CV_mean_O];
Mode_Value_O = [LAT_mode_O; CV_mode_O];
Std_Deviation_O = [LAT_std_O; CV_std_O];
Variance_O = [LAT_var_O; CV_var_O];
Number_Points_O = [LAT_num_O; CV_num_O];
Results.ResultsTable = table(Parameter_O, Max_Value_O, Mean_Value_O, Mode_Value_O, Std_Deviation_O, Variance_O, Number_Points_O, ...
    'RowNames', {'LAT_O', 'CV_O'});

disp(' ');
disp('=== OPTICAL RESULTS TABLE ===');
disp(Results.ResultsTable);


%% Summary Visualizations - Optical
figure('Position', [100 100 1200 400]);

% Subplot 1: LAT Histogram
subplot(1,2,1);
if ~isempty(nonzero_LAT_O)
    histogram(nonzero_LAT_O, 30, 'FaceColor', 'blue', 'FaceAlpha', 0.7);
    xlabel('Local Activation Time (ms)');
    ylabel('Count');
    title('Optical LAT Distribution');
    grid on;
end

% Subplot 2: CV Histogram
subplot(1,2,2);
if ~isempty(nonzero_CV_O)
    histogram(nonzero_CV_O, 40, 'FaceColor', 'red', 'FaceAlpha', 0.7);
    xlabel('Conduction Velocity (cm/s)');
    ylabel('Count');
    title('Optical CV Distribution');
    grid on;
end
sgtitle('Optical Analysis Histograms');


%% Save workspace variables for future use - Optical
% Joining Struct
Results.Fsampling = Fsampling;
Results.optical_lim1 = optical_lim1;
Results.optical_lim2 = optical_lim2;
Results.optical_space_scale = optical_space_scale;
Results.LAT_O_filtered = LAT_O_filtered;
Results.AvgSpeed_O_filtered = AvgSpeed_O_filtered;

% Save Struct
save('O_LAT_CV_CAM_.mat', 'Results', '-v7.3');

disp(' ');
disp('=== OPTICAL ANALYSIS COMPLETE ===');
disp('All results have been saved to "O_LAT_CV_CAM_.mat".');


%% 
%% 
%%

%% ELECTRIC ANALYSIS - INTERPOLATED SIGNALS
%% LAT and CV Maps - Electrical
clear; clc;


%% Loading Data
load("E:\Qualification\Analysis\E32F02R08\data\InterpolatedSignalsE32_F02_R08_filtered.mat"); % Load Interpolated data


%% Configuring
Data_E_Structure = InterpSignal.Sync; 
electrical_space_scale = 2;  % mm per pixel (MEA)
tank_space_scale = 3;       % mm per pixel (TANK)
cases = {'MEA1', 'MEA2', 'MEA3', 'TANK'};
ElectricalResults = struct();
Fsampling = 4000;


%% Selecting Time Window
% Preview Signal
case_name = 'MEA1';
data_temp = Data_E_Structure.(case_name);
Background = squeeze(data_temp(:,:,2000));
pick_up_a_trace(Background, data_temp,1);
clear data_temp case_name Background;


%% Defining Sample Limits
sample_limits = struct();
% MEA1
sample_limits.MEA1.lim1 = 10300; 
sample_limits.MEA1.lim2 = 10850; 
% MEA2
sample_limits.MEA2.lim1 = 9200;
sample_limits.MEA2.lim2 = 9600;
% MEA3
sample_limits.MEA3.lim1 = 9400;
sample_limits.MEA3.lim2 = 10000;
% TANK
sample_limits.TANK.lim1 = 9500;
sample_limits.TANK.lim2 = 10000;


%% Fignal Plot in the selected limits
case_name = 'MEA1';
data_temp = Data_E_Structure.(case_name)(:,:, sample_limits.(case_name).lim1:sample_limits.(case_name).lim2);
Background = squeeze(data_temp(:,:,1));
pick_up_a_trace(Background, data_temp,1);
clear data_temp case_name Background;


%% LAT Calculation - Electrical
debug_LAT = 1; % Use 0, 1 or 2
LAT_values = struct();

for i = 1:length(cases)
    case_name = cases{i};
    % --- DYNAMICALLY DETERMINE LIMITS ---
    electrical_lim1 = sample_limits.(case_name).lim1;
    electrical_lim2 = sample_limits.(case_name).lim2;
    
    % Optional: Boundary and Order Checks
    total_samples = size(Data_E_Structure.(case_name), 3);
    electrical_lim1 = max(1, min(electrical_lim1, total_samples));
    electrical_lim2 = min(total_samples, max(electrical_lim2, 1));
    % Ensure lim1 is before lim2
    if electrical_lim1 > electrical_lim2
        temp = electrical_lim1;
        electrical_lim1 = electrical_lim2;
        electrical_lim2 = temp;
    end
    
    fprintf('  -> Processing %s with window: Samples [%d, %d] | Time [%.3f, %.3f]s.\n', case_name, electrical_lim1, electrical_lim2, electrical_lim1/Fsampling, electrical_lim2/Fsampling);
    % ------------------------------------
    
    Data_current_E = Data_E_Structure.(case_name);
    Data_temp_E = Data_current_E(:,:,:);
    WinStartIdx = electrical_lim1 + 200;
    WinEndIdx = electrical_lim2 + 100;
    
    % Calculate LAT map (pixel by pixel)
    LAT_map = zeros(size(Data_temp_E, 1), size(Data_temp_E, 2));
    for x = 1:size(Data_temp_E, 1) % X-coordinate
        for y = 1:size(Data_temp_E, 2) % Y-coordinate
            signal = squeeze(Data_temp_E(x,y,:));
            if max(signal) ~= 0
                figure(88);
                plot(signal(WinStartIdx:WinEndIdx));

                signal = signal'; 
                LAT_map(x,y) = find_LAT_diff_Tainan(signal, ...
                                    Fsampling, ...
                                    WinStartIdx, ...
                                    WinEndIdx, ...
                                    case_name, x, y, ...
                                    debug_LAT);
            end
        end
    end
    
    % Normalization and Filtering
    LAT_map(LAT_map < 0) = 0;
    % Subtract minimum LAT for re-normalization
    min_LAT = min(LAT_map(LAT_map ~= 0));
    if ~isempty(min_LAT) && ~isnan(min_LAT)
        LAT_map(LAT_map ~= 0) = LAT_map(LAT_map ~= 0) - min_LAT;
    end
    % Spatial Filtering (Using S=1 | T=0 because of ...
    % the size of pixels in electrical mapping)
    LAT_map_filtered = SpatTemp_Filtering(LAT_map, 1, 0, 'GPU'); 
    
    % Store the final, filtered 2D LAT map
    LAT_values.(case_name).LAT_matrix = LAT_map_filtered;
    fprintf('  -> %s LAT matrix (Size: %dx%d) stored.\n', case_name, size(LAT_map_filtered, 1), size(LAT_map_filtered, 2));
end


%% LAT Calculation COM - Electrical
debug_COM = 0; % 1 for debugging plots 
LAT_values_COM = struct();
for i = 1:length(cases)
    case_name = cases{i};
    
    % Use the specific sample limits defined in the structure
    electrical_lim1 = sample_limits.(case_name).lim1;
    electrical_lim2 = sample_limits.(case_name).lim2;
    
    % Optional: Boundary and Order Checks
    total_samples = size(Data_E_Structure.(case_name), 3);
    electrical_lim1 = max(1, min(electrical_lim1, total_samples));
    electrical_lim2 = min(total_samples, max(electrical_lim2, 1));
    % Ensure lim1 is before lim2
    if electrical_lim1 > electrical_lim2
        temp = electrical_lim1;
        electrical_lim1 = electrical_lim2;
        electrical_lim2 = temp;
    end
    
    fprintf('  -> Processing %s (COM Method) with window: Samples [%d, %d] | Time [%.3f, %.3f]s.\n', case_name, electrical_lim1, electrical_lim2, electrical_lim1/Fsampling, electrical_lim2/Fsampling);
    % -------------------------------------------------------------------
    
    Data_current_E = Data_E_Structure.(case_name);
    Data_temp_E = Data_current_E(:,:,electrical_lim1:electrical_lim2);
    
    % Get the size of the segment
    [num_x, num_y, segment_length] = size(Data_temp_E);
    % Inicialize matrix
    LAT_map_COM = zeros(num_x, num_y);
    for x = 1:num_x % X-coordinate
        for y = 1:num_y % Y-coordinate
            signal_segment = squeeze(Data_temp_E(x,y,:))'; % Ensure it's a row vector
            if max(abs(signal_segment)) ~= 0
                % Calculate LAT map (pixel by pixel)
                LAT_map_COM(x,y) = find_LAT_com(signal_segment, ...
                                    Fsampling, x, y, ...
                                    debug_COM, case_name);
            end
        end
    end
    % Normalization and Filtering
    LAT_map_COM(LAT_map_COM < 0) = 0;
    min_LAT_COM = min(LAT_map_COM(LAT_map_COM ~= 0));
    if ~isempty(min_LAT_COM) && ~isnan(min_LAT_COM)
        LAT_map_COM(LAT_map_COM ~= 0) = LAT_map_COM(LAT_map_COM ~= 0) - min_LAT_COM;
    end
    % Spatial Filtering (Using S=1 | T=0)
    LAT_map_filtered_COM = SpatTemp_Filtering(LAT_map_COM, 1, 0, 'GPU'); 
    
    % Store
    LAT_values_COM.(case_name).LAT_matrix = LAT_map_filtered_COM;
    fprintf('  -> %s LAT matrix (COM Method) (Size: %dx%d) stored.\n', case_name, size(LAT_map_filtered_COM, 1), size(LAT_map_filtered_COM, 2));
end


%% LAT MAP PLOTTING - Electrical
C = parula(256); 
% C(1,1:3) = [1 1 1]; % White for background (0 LAT)
for i = 1:length(cases)
    case_name = cases{i};
    LAT_matrix_DIFF = LAT_values.(case_name).LAT_matrix;
    LAT_matrix_COM = LAT_values_COM.(case_name).LAT_matrix;
    % Create a figure to show both methods
    figure('color', 'white', 'Position', [40 + i*20 40 + i*20 1200 600]);
    % Determine the common color limits
    all_LATS = [LAT_matrix_DIFF(LAT_matrix_DIFF ~= 0); LAT_matrix_COM(LAT_matrix_COM ~= 0)];
    if ~isempty(all_LATS)
        max_LAT = max(all_LATS);
    else
        max_LAT = 1; % Default max
    end

    % --- Subplot 1: Derivative Method ---
    subplot(1, 2, 1);
    if strcmp(case_name, 'TANK')
        J = LAT_matrix_DIFF;
        imagesc(J);
        pbaspect([2 1 1]);
    else
        J = imrotate(LAT_matrix_DIFF, 90);
        imagesc(J);
        axis equal;
    end
    colormap(gca, C);
    caxis([0 max_LAT]);
    axis off;
    title(['Derivative Method - ' case_name], 'FontSize', 16);
    
    % --- Subplot 2: COM Method ---
    subplot(1, 2, 2);
    if strcmp(case_name, 'TANK')
        J = LAT_matrix_COM;
        imagesc(J);
        pbaspect([2 1 1]);
    else
        J = imrotate(LAT_matrix_COM, 90);
        imagesc(J);
        axis equal;
    end
    colormap(gca, C);
    caxis([0 max_LAT]);
    axis off;
    hBar = colorbar('eastoutside');
    ylabel(hBar, 'Local Activation Time [ms]', 'FontSize', 14);
    title(['COM Method - ' case_name], 'FontSize', 16);
    
    sgtitle(['Electrical ' case_name ' - LAT Map Comparison'], 'FontSize', 18, 'FontWeight', 'bold');
end


%% LAT METHOD CONSOLIDATION - Electrical
% Choose the preferred LAT detection method for each case:
% 1. 'DIFF' (Derivative Method: LAT_values)
% 2. 'COM' (Center of Mass Method: LAT_values_COM)

% --- Define Method Selection ---
MethodSelection = struct();
MethodSelection.MEA1 = 'DIFF';
MethodSelection.MEA2 = 'DIFF';
MethodSelection.MEA3 = 'DIFF';
MethodSelection.TANK = 'COM';

LAT_values_FINAL = struct();

for i = 1:length(cases)
    case_name = cases{i};
    % Get the selected method for the current case
    if isfield(MethodSelection, case_name)
        method = MethodSelection.(case_name);
    else
        method = 'DIFF'; % Default to DIFF
    end
    
    % Select the matrix based on the method
    if strcmp(method, 'COM')
        LAT_matrix_final = LAT_values_COM.(case_name).LAT_matrix;
        fprintf('  -> %s: Using COM Method (LAT_values_COM)\n', case_name);
    elseif strcmp(method, 'DIFF')
        LAT_matrix_final = LAT_values.(case_name).LAT_matrix;
        fprintf('  -> %s: Using Derivative Method (LAT_values)\n', case_name);
    else
        warning('Unknown method "%s" specified for case %s. Defaulting to DIFF.', method, case_name);
        LAT_matrix_final = LAT_values.(case_name).LAT_matrix;
    end
    
    LAT_values_FINAL.(case_name).LAT_matrix = LAT_matrix_final;
end


%% In case of need to manual corrections - Electrical
% Histogram plots
num_bins = 50;
for i = 1:length(cases)
    case_name = cases{i};
    LAT_matrix = LAT_values_FINAL.(case_name).LAT_matrix;
    valid_LATS = LAT_matrix(LAT_matrix ~= 0);
    figure();
    set(gcf, 'color', 'white', 'Position', [600 50 600 400]);
    h = histogram(valid_LATS, num_bins);
    title(['LAT Distribution: Electrical ', case_name, ' (FINAL)'], 'FontSize', 14);
    xlabel('Local Activation Time [ms]', 'FontSize', 12);
    ylabel('Count (Number of Electrodes/Pixels)', 'FontSize', 12);
    grid on;
    mean_lat = mean(valid_LATS);
    std_lat = std(valid_LATS);
    fprintf('  -> %s LAT Stats: Mean = %.2f ms, STD = %.2f ms\n', ...
            case_name, mean_lat, std_lat);
end


% Substitute specific value
case_to_correct = 'MEA1';
LAT_temp = LAT_values_FINAL.(case_to_correct).LAT_matrix;
find_value = 9; % Value to replace
tolerancia = 1;
indices = find(abs(LAT_temp - find_value) < tolerancia);
LAT_temp(indices) = 3; % Value to include
LAT_values_FINAL.(case_to_correct).LAT_matrix = LAT_temp;


% Substitute Higher or lower than
case_to_correct = 'MEA1';
find_value = 9; % Value to replace
new_value = 0;
LAT_values_FINAL.(case_to_correct).LAT_matrix( ...
            LAT_values_FINAL.(case_to_correct).LAT_matrix < ...
            find_value) = new_value;


% Correct The Zero Normalization - Subtract minimum
LAT_values_FINAL.(case_to_correct).LAT_matrix(LAT_values_FINAL.(case_to_correct).LAT_matrix < 0) = 0;
min_LAT = min(LAT_values_FINAL.(case_to_correct).LAT_matrix(LAT_values_FINAL.(case_to_correct).LAT_matrix ~= 0));
if ~isempty(min_LAT) && ~isnan(min_LAT)
    LAT_values_FINAL.(case_to_correct).LAT_matrix(LAT_values_FINAL.(case_to_correct).LAT_matrix ~= 0) = ...
        LAT_values_FINAL.(case_to_correct).LAT_matrix(LAT_values_FINAL.(case_to_correct).LAT_matrix ~= 0) - min_LAT;
end


%% FINAL LAT MAPS - Electrical
C = parula(256); 
% C(1,1:3) = [1 1 1]; % White for background

for i = 1:length(cases)
    case_name = cases{i};
    LAT_matrix_final = LAT_values_FINAL.(case_name).LAT_matrix;
    
    % Determine which method was used
    if isfield(MethodSelection, case_name)
        method = MethodSelection.(case_name);
        method_text = method;
    else
        method = 'DIFF';
        method_text = 'DIFF (default)';
    end
    
    % Create figure
    fig = figure('color', 'white', 'Position', [50 + i*30 50 + i*30 800 600]);
    fig.Name = ['LAT Map - ' case_name];
    
    % INTERPOLATION: Upsample the matrix for smoother visualization
    interpolation_factor = 5; % Increase this for smoother images
    [rows, cols] = size(LAT_matrix_final);
    % Create interpolation grid
    [X, Y] = meshgrid(1:cols, 1:rows);
    [Xq, Yq] = meshgrid(linspace(1, cols, cols*interpolation_factor), ...
                       linspace(1, rows, rows*interpolation_factor));
    % Interpolate the data
    LAT_interp = interp2(X, Y, LAT_matrix_final, Xq, Yq, 'cubic');
    
    % Handle different case geometries
    if strcmp(case_name, 'TANK')
        J = LAT_interp;
        imagesc(J);
        pbaspect([2 1 1]);
    else
        J = imrotate(LAT_interp, 90);
        imagesc(J);
        axis equal;
    end
    
    % Apply colormap
    colormap(C);
    caxis([0 max(LAT_matrix_final(:))]); % Use original max for axis limits
    axis off;
    
    % Add colorbar
    hBar = colorbar('eastoutside');
    ylabel(hBar, 'Local Activation Time [ms]', 'FontSize', 12);
    
    % Add main title
    title_str = sprintf('%s - LAT Map\n[Method: %s] (Interpolated)', case_name, method_text);
    title(title_str, 'FontSize', 16, 'FontWeight', 'bold');
    
    % Add LAT statistics as text annotation
    valid_LATS = LAT_matrix_final(LAT_matrix_final ~= 0);
    if ~isempty(valid_LATS)
        stats_str = sprintf('Min: %.1f ms\nMax: %.1f ms\nMean: %.1f ms\nStd: %.1f ms', ...
            min(valid_LATS), max(valid_LATS), mean(valid_LATS), std(valid_LATS));
        
        % Add text box with statistics
        annotation('textbox', [0.02, 0.03, 0.2, 0.15], ...
            'String', stats_str, ...
            'FontSize', 10, ...
            'BackgroundColor', 'white', ...
            'EdgeColor', 'black', ...
            'Margin', 5);
    end    
    fprintf('Created final LAT map for %s using %s method (interpolated)\n', case_name, method_text);
end


%% COMPARISON FIGURE WITH ALL CASES
if length(cases) <= 6
    fig_all = figure('color', 'white', 'Position', [100 100 1400 800]);
    fig_all.Name = 'LAT Maps - Comparison';
    
    % Determine subplot arrangement
    if length(cases) <= 4
        nrows = 2; ncols = 2;
    else
        nrows = 2; ncols = 3;
    end
    
    % Find global max for consistent color scaling
    all_LATS_global = [];
    for i = 1:length(cases)
        case_name = cases{i};
        LAT_matrix_final = LAT_values_FINAL.(case_name).LAT_matrix;
        valid_LATS = LAT_matrix_final(LAT_matrix_final ~= 0);
        all_LATS_global = [all_LATS_global; valid_LATS];
    end
    if ~isempty(all_LATS_global)
        max_LAT_global = max(all_LATS_global);
    else
        max_LAT_global = 1;
    end
    
    % Plot all cases in subplots
    for i = 1:length(cases)
        case_name = cases{i};
        LAT_matrix_final = LAT_values_FINAL.(case_name).LAT_matrix;
        
        if isfield(MethodSelection, case_name)
            method = MethodSelection.(case_name);
        else
            method = 'DIFF';
        end
        
        subplot(nrows, ncols, i);
        
        % INTERPOLATION for comparison figure too
        interpolation_factor = 3;
        [rows, cols] = size(LAT_matrix_final);
        % Create interpolation grid
        [X, Y] = meshgrid(1:cols, 1:rows);
        [Xq, Yq] = meshgrid(linspace(1, cols, cols*interpolation_factor), ...
                           linspace(1, rows, rows*interpolation_factor));
        % Interpolate the data
        LAT_interp = interp2(X, Y, LAT_matrix_final, Xq, Yq, 'cubic');
        if strcmp(case_name, 'TANK')
            J = LAT_interp;
            imagesc(J);
            pbaspect([2 1 1]);
        else
            J = imrotate(LAT_interp, 90);
            imagesc(J);
            axis equal;
        end
        
        colormap(gca, C);
        caxis([0 max_LAT_global]);
        axis off;
        
        title_str = sprintf('%s\n%s', case_name, method);
        title(title_str, 'FontSize', 12, 'FontWeight', 'bold');
    end
    
    % Add overall colorbar
    hBar = colorbar('Position', [0.93 0.15 0.02 0.7]);
    ylabel(hBar, 'Local Activation Time [ms]', 'FontSize', 12);
    
    sgtitle('LAT Maps Comparison (Interpolated)', 'FontSize', 16, 'FontWeight', 'bold');
end


%% LAT Statistics - Electrical (with ROI Selection)
C = parula(256); 
% C(1,1:3) = [1 1 1]; White Background
for i = 1:length(cases)
    case_name = cases{i};
    LAT_matrix = LAT_values_FINAL.(case_name).LAT_matrix;
    
    % ROI SELECTION
    figure('color', 'white');
    if strcmp(case_name, 'TANK')
        J = LAT_matrix;
        imagesc(J, 'Interpolation', 'bilinear');
        pbaspect([2 1 1]); % 2:1 aspect ratio
    else
        J = LAT_matrix;
        imagesc(J);
        axis equal;
    end
    colormap(C);
    axis off;
    colorbar('eastoutside');
    title(['Case: ', case_name, ' - ROI Selection (FINAL LAT)'], 'FontSize', 14);
    roi_mask = roipoly;
    close(gcf);
    
    % APPLY ROI AND ISOLATE VALID DATA
    LAT_roi = J .* roi_mask; 
    nonzero_LATS_roi = LAT_roi(LAT_roi ~= 0 & ~isnan(LAT_roi));
    % CALCULATE AND STORE STATISTICS
    LAT_mean = mean(nonzero_LATS_roi);
    LAT_std = std(nonzero_LATS_roi);
    LAT_max = max(nonzero_LATS_roi);
    LAT_mode = mode(nonzero_LATS_roi); 
    LAT_var = var(nonzero_LATS_roi); 
    LAT_num = numel(nonzero_LATS_roi); 
    % Displaying results for the current case
    fprintf('\nCase: %s (ROI Stats - FINAL LAT)\n', case_name);
    disp(['  Average LAT:      ', num2str(LAT_mean, '%.4f'), ' ms']);
    disp(['  Mode LAT:         ', num2str(LAT_mode, '%.4f'), ' ms']);
    disp(['  Max LAT:          ', num2str(LAT_max, '%.4f'), ' ms']);
    disp(['  Std Deviation:    ', num2str(LAT_std, '%.4f'), ' ms']);
    disp(['  Variance:         ', num2str(LAT_var, '%.4f')]);
    disp(['  Number of points: ', num2str(LAT_num)]);
    % Store Results
    LAT_values_FINAL.(case_name).LAT_stats = struct(...
        'mean', LAT_mean, 'std', LAT_std, 'max', LAT_max, ...
        'mode', LAT_mode, 'var', LAT_var, 'num_points', LAT_num);
    LAT_values_FINAL.(case_name).roi_mask = roi_mask;
end


%% CV Calculation - Electrical
SpaceScale.TANK = 0.75; % mm/pixel for TANK
SpaceScale.MEA = 0.1;  % mm/pixel for MEA
r = 3;                 % Radius (r) for the CV Circle Method
filter_size = 0;       % Size of the spatial filter kernel
LAT_values_CV = LAT_values_FINAL; % Start CV calculation from FINAL LAT
for i = 1:length(cases)
    case_name = cases{i};
    LAT_E_matrix = LAT_values_FINAL.(case_name).LAT_matrix;
    
    % Determine Space Scale
    if strcmp(case_name, 'TANK')
        current_space_scale = SpaceScale.TANK;
    else
        current_space_scale = SpaceScale.MEA;
    end
    % Interpolation
    [X, Y] = meshgrid(1:size(LAT_E_matrix, 2), 1:size(LAT_E_matrix, 1));
    [Xq, Yq] = meshgrid(1:0.1:size(LAT_E_matrix, 2), 1:0.1:size(LAT_E_matrix, 1));
    % Perform bilinear interpolation
    LAT_E_matrix_interp = interp2(X, Y, LAT_E_matrix, Xq, Yq, 'linear');
    
    % Initialize matrices for CV results
    [R_interp, C_interp] = size(LAT_E_matrix_interp);
    AvgSpeed_E = zeros(R_interp, C_interp);
    StdSpeed_E = zeros(R_interp, C_interp);
    Angle_E = zeros(R_interp, C_interp);
    % CV Circle Method Calculation
    for row = (r + 1):(R_interp - (r + 1))
        for col = (r + 1):(C_interp - (r + 1))
            if LAT_E_matrix_interp(row, col) ~= 0 
                [AvgSpeed_E(row, col), StdSpeed_E(row, col), Angle_E(row, col)] = ...
                    CV_CircleMethod(LAT_E_matrix_interp, r, row, col, current_space_scale);
            end
        end
    end
    % Filtering and Storage
    AvgSpeed_E_filtered = SpatTemp_Filtering(AvgSpeed_E, filter_size, 0, 'GPU');
    
    % Store the results in the FINAL structure
    LAT_values_FINAL.(case_name).AvgSpeed_E_matrix = AvgSpeed_E_filtered;
    LAT_values_FINAL.(case_name).SpaceScale_mm_px = current_space_scale;
    LAT_values_FINAL.(case_name).CV_Radius = r;
    
    fprintf('  -> %s CV map calculated and filtered (Size: %dx%d).\n', case_name, R_interp, C_interp);
end


%% CV Plot - Electric
C = jet(256);
for i = 1:length(cases)
    case_name = cases{i};
    AvgSpeed_E_matrix = LAT_values_FINAL.(case_name).AvgSpeed_E_matrix;
    J_valid = AvgSpeed_E_matrix(AvgSpeed_E_matrix ~= 0 & ~isnan(AvgSpeed_E_matrix));
    % Use percentile limits for robust visualization (e.g., 5th to 95th percentile)
    Y_limits = prctile(J_valid, [5 95], 'all');
    CV_min = Y_limits(1);
    CV_max = Y_limits(2);
    
    % --- Plotting ---
    figure('color', 'white', 'Position', [40 + i*20 40 + i*20 600 600]);
    if strcmp(case_name, 'TANK')
        % TANK case: Rectangular aspect ratio, no rotation
        J_plot = AvgSpeed_E_matrix;
        imagesc(J_plot);
        pbaspect([2 1 1]); % 2:1 aspect ratio
    else
        % MEA cases: Square aspect ratio, rotate 90 degrees
        J_plot = imrotate(AvgSpeed_E_matrix, 90);
        imagesc(J_plot);
        axis equal; % Ensures square aspect ratio
    end
    % Apply the calculated color limits
    caxis([CV_min CV_max]); 
    colormap(C);
    axis off;
    hBar = colorbar('eastoutside');
    ylabel(hBar, 'Conduction Velocity [cm/s]', 'FontSize', 14);
    title(['Electrical ' case_name ' - Conduction Velocity Map (FINAL LAT)'], 'FontSize', 16);
    fprintf('  -> %s CV Map plotted with limits [%.2f - %.2f] cm/s.\n', case_name, CV_min, CV_max);
end


%% CV Statistics - Electrical
for i = 1:length(cases)
    case_name = cases{i};
    
    % --- Retrieve the CV matrix directly from the struct ---
    if ~isfield(LAT_values_FINAL.(case_name), 'AvgSpeed_E_matrix')
        fprintf('Case %s: Warning - AvgSpeed_E_matrix field not found. Skipping CV stats.\n', case_name);
        continue; 
    end
    
    CV_matrix = LAT_values_FINAL.(case_name).AvgSpeed_E_matrix;
    
    % Check if the matrix is empty before proceeding
    if isempty(CV_matrix)
        fprintf('Case %s: Warning - CV matrix is empty. Skipping CV stats.\n', case_name);
        continue;
    end
    
    % Isolate valid CV points (non-zero and non-NaN)
    nonzero_CV = CV_matrix(CV_matrix ~= 0 & ~isnan(CV_matrix) & ~isinf(CV_matrix));
    
    if isempty(nonzero_CV)
        fprintf('Case %s: Warning - No valid non-zero CV points found. Skipping stats.\n', case_name);
        LAT_values_FINAL.(case_name).CV_stats = []; 
        continue;
    end
    
    % Calculate Statistics
    CV_mean = mean(nonzero_CV);
    CV_std = std(nonzero_CV);
    CV_max = max(nonzero_CV);
    CV_mode = mode(nonzero_CV); 
    CV_var = var(nonzero_CV); 
    CV_num = numel(nonzero_CV);
    
    % Store Results in the LAT_values_FINAL structure
    LAT_values_FINAL.(case_name).CV_stats = struct(...
        'mean', CV_mean, 'std', CV_std, 'max', CV_max, ...
        'mode', CV_mode, 'var', CV_var, 'num_points', CV_num);
    
    % Displaying results for the current case
    fprintf('\nCase: %s (FINAL LAT)\n', case_name);
    disp(['  Average CV:       ', num2str(CV_mean, '%.4f'), ' cm/s']);
    disp(['  Mode CV:          ', num2str(CV_mode, '%.4f'), ' cm/s']);
    disp(['  Max CV:           ', num2str(CV_max, '%.4f'), ' cm/s']);
    disp(['  Std Deviation:    ', num2str(CV_std, '%.4f'), ' cm/s']);
    disp(['  Variance:         ', num2str(CV_var, '%.4f')]);
    disp(['  Number of points: ', num2str(CV_num)]);
end


%% FINAL SAVING OF ALL RESULTS
% Save ALL essential variables for complete reproducibility

% 1. CORE RESULTS STRUCTURES
save_vars = struct();

% 1.1 Original Data and Configuration
save_vars.Data_E_Structure = Data_E_Structure;
save_vars.Fsampling = Fsampling;
save_vars.cases = cases;
save_vars.electrical_space_scale = electrical_space_scale;
save_vars.tank_space_scale = tank_space_scale;

% 1.2 Analysis Windows
save_vars.sample_limits = sample_limits;

% 1.3 Method Selection
save_vars.MethodSelection = MethodSelection;

% 1.4 All LAT Results (Raw, COM, Final)
save_vars.LAT_values = LAT_values;          % Derivative method results
save_vars.LAT_values_COM = LAT_values_COM;  % COM method results  
save_vars.LAT_values_FINAL = LAT_values_FINAL; % Consolidated final results

% 1.5 CV Results (already in LAT_values_FINAL, but keep separate for clarity)
if isfield(LAT_values_FINAL.(cases{1}), 'AvgSpeed_E_matrix')
    save_vars.CV_results = struct();
    for i = 1:length(cases)
        case_name = cases{i};
        save_vars.CV_results.(case_name).CV_matrix = LAT_values_FINAL.(case_name).AvgSpeed_E_matrix;
        if isfield(LAT_values_FINAL.(case_name), 'CV_stats')
            save_vars.CV_results.(case_name).CV_stats = LAT_values_FINAL.(case_name).CV_stats;
        end
        save_vars.CV_results.(case_name).SpaceScale = LAT_values_FINAL.(case_name).SpaceScale_mm_px;
        save_vars.CV_results.(case_name).CV_Radius = LAT_values_FINAL.(case_name).CV_Radius;
    end
end

% 1.6 Statistics Table - CREATE IF IT DOESN'T EXIST
if ~exist('StatsTable', 'var')
    % Create the statistics table from LAT_values_FINAL
    StatsTable = table();
    for i = 1:length(cases)
        case_name = cases{i};
        
        % Get LAT statistics
        if isfield(LAT_values_FINAL.(case_name), 'LAT_stats')
            lat_stats = LAT_values_FINAL.(case_name).LAT_stats;
            lat_mean = lat_stats.mean;
            lat_std = lat_stats.std;
            lat_max = lat_stats.max;
            lat_num = lat_stats.num_points;
        else
            % Calculate from matrix if stats don't exist
            LAT_matrix = LAT_values_FINAL.(case_name).LAT_matrix;
            valid_LATS = LAT_matrix(LAT_matrix ~= 0);
            if ~isempty(valid_LATS)
                lat_mean = mean(valid_LATS);
                lat_std = std(valid_LATS);
                lat_max = max(valid_LATS);
                lat_num = numel(valid_LATS);
            else
                lat_mean = NaN;
                lat_std = NaN;
                lat_max = NaN;
                lat_num = 0;
            end
        end
        
        % Get CV statistics
        if isfield(LAT_values_FINAL.(case_name), 'CV_stats')
            cv_stats = LAT_values_FINAL.(case_name).CV_stats;
            cv_mean = cv_stats.mean;
            cv_std = cv_stats.std;
            cv_max = cv_stats.max;
            cv_num = cv_stats.num_points;
        else
            % Set to NaN if CV stats don't exist
            cv_mean = NaN;
            cv_std = NaN;
            cv_max = NaN;
            cv_num = 0;
        end
        
        % Create table row
        NewRow = table(...
            {case_name}, ...
            lat_mean, lat_std, lat_max, lat_num, ...
            cv_mean, cv_std, cv_max, cv_num, ...
            'VariableNames', {'Case', 'LAT_Mean_ms', 'LAT_STD_ms', 'LAT_Max_ms', 'LAT_N_Pixels', ...
                             'CV_Mean_cms', 'CV_STD_cms', 'CV_Max_cms', 'CV_N_Points'});
        
        StatsTable = [StatsTable; NewRow];
    end
    fprintf('Created StatsTable from LAT_values_FINAL data.\n');
end
save_vars.StatsTable = StatsTable;

% 1.7 ROI Masks (if they exist)
save_vars.roi_masks = struct();
for i = 1:length(cases)
    case_name = cases{i};
    if isfield(LAT_values_FINAL.(case_name), 'roi_mask')
        save_vars.roi_masks.(case_name) = LAT_values_FINAL.(case_name).roi_mask;
    end
end

% 2. ADD ANALYSIS PARAMETERS (important for reproducibility)
save_vars.AnalysisParams = struct();
if exist('SpaceScale', 'var')
    save_vars.AnalysisParams.SpaceScale = SpaceScale;
end
% Add other parameters that exist in your workspace
vars_to_check = {'debug_LAT', 'debug_COM'};
for var_idx = 1:length(vars_to_check)
    var_name = vars_to_check{var_idx};
    if exist(var_name, 'var')
        save_vars.AnalysisParams.(var_name) = eval(var_name);
    end
end

% 3. ADD METADATA
save_vars.Metadata = struct();
save_vars.Metadata.analysis_date = datestr(now, 'yyyy-mm-dd HH:MM:SS');

% 4. DEFINE FILENAME
filename = sprintf('E_LAT_CV.mat');

% 5. SAVE ALL VARIABLES
save(filename, 'save_vars', '-v7.3'); % Use -v7.3 for large files
fprintf('\n✅ COMPLETE ANALYSIS SAVED TO: %s\n', filename);

% 6. EXPORT STATISTICS TO CSV
if ~isempty(StatsTable)
    csv_filename = sprintf('E_Statistics.csv');
    writetable(StatsTable, csv_filename);
    fprintf('✅ STATISTICS EXPORTED TO CSV: %s\n', csv_filename);
end










































































%% ELECTRIC ANALYSIS - REAL SIGNALS ONLY
%% LAT and CV Maps - Electrical
clear; clc;


%% Loading Data
load("E:\Qualification\Analysis\E32F02R08\data\data_filtered_sync_E32_F02_R08.mat"); % Load Synchronized data


%% Configuring
% General configurations
Data_E_Structure = D_SYNC.EL; 
electrical_space_scale = 2;  % mm per pixel (MEA)
tank_space_scale = 3;       % mm per pixel (TANK)
cases = {'MEA1', 'MEA2', 'MEA3', 'TANK'};
ElectricalResults = struct();
Fsampling = 4000;
% Electrodes distribution
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
el.TANK = [145, 146, 155, 156, 165, 166, 129, 130, 139, 140, 181, 182, ...
            147, NaN, 157, NaN, 167, NaN, 131, NaN, 141, NaN, 183, NaN, ...
            148, 149, 158, 159, 168, 169, 132, 133, 142, 143, 184, 185, ...
            150, 151, 160, 161, 170, 171, 134, 135, 144, 177, 186, 187, ...
            152, NaN, 162, NaN, 172, NaN, 136, NaN, 178, NaN, 188, NaN, ...
            153, 154, 163, 164, 173, 174, 137, 138, 179, 180, 189, 190];


%% Selecting Time Window
% Selecting electrode
el_plot = [3, 22, 90];
% Preview Signal
figure();
for i = 1:length(el_plot)
    plot(Data_E_Structure(el_plot(i), :), 'LineWidth', 1); 
    hold on;
end
hold off;
clear el_plot i;

%% Defining Sample Limits
sample_limits = struct();
% MEA1
sample_limits.MEA1.lim1 = 10300; 
sample_limits.MEA1.lim2 = 10850; 
% MEA2
sample_limits.MEA2.lim1 = 9200;
sample_limits.MEA2.lim2 = 9600;
% MEA3
sample_limits.MEA3.lim1 = 9400;
sample_limits.MEA3.lim2 = 10000;
% TANK
sample_limits.TANK.lim1 = 9500;
sample_limits.TANK.lim2 = 10000;


%% LAT Calculation - Electrical
debug_LAT = 1; % Use 0, 1 or 2
LAT_values = struct();

for i = 1:length(cases)
    case_name = cases{i};
    % --- DYNAMICALLY DETERMINE LIMITS ---
    electrical_lim1 = sample_limits.(case_name).lim1;
    electrical_lim2 = sample_limits.(case_name).lim2;
    
    % Optional: Boundary and Order Checks
    total_samples = size(Data_E_Structure.(case_name), 2);
    electrical_lim1 = max(1, min(electrical_lim1, total_samples));
    electrical_lim2 = min(total_samples, max(electrical_lim2, 1));
    % Ensure lim1 is before lim2
    if electrical_lim1 > electrical_lim2
        temp = electrical_lim1;
        electrical_lim1 = electrical_lim2;
        electrical_lim2 = temp;
    end
    
    fprintf('  -> Processing %s with window: Samples [%d, %d] | Time [%.3f, %.3f]s.\n', ...
        case_name, electrical_lim1, electrical_lim2, ...
        electrical_lim1/Fsampling, ...
        electrical_lim2/Fsampling);
    
    Data_current_E = Data_E_Structure.(case_name);
    Data_temp_E = Data_current_E(:,:,:);
    WinStartIdx = electrical_lim1 + 200;
    WinEndIdx = electrical_lim2 + 100;
    
    % Calculate LAT map (pixel by pixel)
    LAT_map = zeros(size(Data_temp_E, 1), size(Data_temp_E, 2));
    for x = 1:size(Data_temp_E, 1) % X-coordinate
        for y = 1:size(Data_temp_E, 2) % Y-coordinate
            signal = squeeze(Data_temp_E(x,y,:));
            if max(signal) ~= 0
                figure(88);
                plot(signal(WinStartIdx:WinEndIdx));

                signal = signal'; 
                LAT_map(x,y) = find_LAT_diff_Tainan(signal, ...
                                    Fsampling, ...
                                    WinStartIdx, ...
                                    WinEndIdx, ...
                                    case_name, x, y, ...
                                    debug_LAT);
            end
        end
    end
    
    % Normalization and Filtering
    LAT_map(LAT_map < 0) = 0;
    % Subtract minimum LAT for re-normalization
    min_LAT = min(LAT_map(LAT_map ~= 0));
    if ~isempty(min_LAT) && ~isnan(min_LAT)
        LAT_map(LAT_map ~= 0) = LAT_map(LAT_map ~= 0) - min_LAT;
    end
    % Spatial Filtering (Using S=1 | T=0 because of ...
    % the size of pixels in electrical mapping)
    LAT_map_filtered = SpatTemp_Filtering(LAT_map, 1, 0, 'GPU'); 
    
    % Store the final, filtered 2D LAT map
    LAT_values.(case_name).LAT_matrix = LAT_map_filtered;
    fprintf('  -> %s LAT matrix (Size: %dx%d) stored.\n', case_name, size(LAT_map_filtered, 1), size(LAT_map_filtered, 2));
end








































































%% LAT Calculation COM - Electrical
debug_COM = 0; % 1 for debugging plots 
LAT_values_COM = struct();
for i = 1:length(cases)
    case_name = cases{i};
    
    % Use the specific sample limits defined in the structure
    electrical_lim1 = sample_limits.(case_name).lim1;
    electrical_lim2 = sample_limits.(case_name).lim2;
    
    % Optional: Boundary and Order Checks
    total_samples = size(Data_E_Structure.(case_name), 3);
    electrical_lim1 = max(1, min(electrical_lim1, total_samples));
    electrical_lim2 = min(total_samples, max(electrical_lim2, 1));
    % Ensure lim1 is before lim2
    if electrical_lim1 > electrical_lim2
        temp = electrical_lim1;
        electrical_lim1 = electrical_lim2;
        electrical_lim2 = temp;
    end
    
    fprintf('  -> Processing %s (COM Method) with window: Samples [%d, %d] | Time [%.3f, %.3f]s.\n', case_name, electrical_lim1, electrical_lim2, electrical_lim1/Fsampling, electrical_lim2/Fsampling);
    % -------------------------------------------------------------------
    
    Data_current_E = Data_E_Structure.(case_name);
    Data_temp_E = Data_current_E(:,:,electrical_lim1:electrical_lim2);
    
    % Get the size of the segment
    [num_x, num_y, segment_length] = size(Data_temp_E);
    % Inicialize matrix
    LAT_map_COM = zeros(num_x, num_y);
    for x = 1:num_x % X-coordinate
        for y = 1:num_y % Y-coordinate
            signal_segment = squeeze(Data_temp_E(x,y,:))'; % Ensure it's a row vector
            if max(abs(signal_segment)) ~= 0
                % Calculate LAT map (pixel by pixel)
                LAT_map_COM(x,y) = find_LAT_com(signal_segment, ...
                                    Fsampling, x, y, ...
                                    debug_COM, case_name);
            end
        end
    end
    % Normalization and Filtering
    LAT_map_COM(LAT_map_COM < 0) = 0;
    min_LAT_COM = min(LAT_map_COM(LAT_map_COM ~= 0));
    if ~isempty(min_LAT_COM) && ~isnan(min_LAT_COM)
        LAT_map_COM(LAT_map_COM ~= 0) = LAT_map_COM(LAT_map_COM ~= 0) - min_LAT_COM;
    end
    % Spatial Filtering (Using S=1 | T=0)
    LAT_map_filtered_COM = SpatTemp_Filtering(LAT_map_COM, 1, 0, 'GPU'); 
    
    % Store
    LAT_values_COM.(case_name).LAT_matrix = LAT_map_filtered_COM;
    fprintf('  -> %s LAT matrix (COM Method) (Size: %dx%d) stored.\n', case_name, size(LAT_map_filtered_COM, 1), size(LAT_map_filtered_COM, 2));
end


%% LAT MAP PLOTTING - Electrical
C = parula(256); 
% C(1,1:3) = [1 1 1]; % White for background (0 LAT)
for i = 1:length(cases)
    case_name = cases{i};
    LAT_matrix_DIFF = LAT_values.(case_name).LAT_matrix;
    LAT_matrix_COM = LAT_values_COM.(case_name).LAT_matrix;
    % Create a figure to show both methods
    figure('color', 'white', 'Position', [40 + i*20 40 + i*20 1200 600]);
    % Determine the common color limits
    all_LATS = [LAT_matrix_DIFF(LAT_matrix_DIFF ~= 0); LAT_matrix_COM(LAT_matrix_COM ~= 0)];
    if ~isempty(all_LATS)
        max_LAT = max(all_LATS);
    else
        max_LAT = 1; % Default max
    end

    % --- Subplot 1: Derivative Method ---
    subplot(1, 2, 1);
    if strcmp(case_name, 'TANK')
        J = LAT_matrix_DIFF;
        imagesc(J);
        pbaspect([2 1 1]);
    else
        J = imrotate(LAT_matrix_DIFF, 90);
        imagesc(J);
        axis equal;
    end
    colormap(gca, C);
    caxis([0 max_LAT]);
    axis off;
    title(['Derivative Method - ' case_name], 'FontSize', 16);
    
    % --- Subplot 2: COM Method ---
    subplot(1, 2, 2);
    if strcmp(case_name, 'TANK')
        J = LAT_matrix_COM;
        imagesc(J);
        pbaspect([2 1 1]);
    else
        J = imrotate(LAT_matrix_COM, 90);
        imagesc(J);
        axis equal;
    end
    colormap(gca, C);
    caxis([0 max_LAT]);
    axis off;
    hBar = colorbar('eastoutside');
    ylabel(hBar, 'Local Activation Time [ms]', 'FontSize', 14);
    title(['COM Method - ' case_name], 'FontSize', 16);
    
    sgtitle(['Electrical ' case_name ' - LAT Map Comparison'], 'FontSize', 18, 'FontWeight', 'bold');
end


%% LAT METHOD CONSOLIDATION - Electrical
% Choose the preferred LAT detection method for each case:
% 1. 'DIFF' (Derivative Method: LAT_values)
% 2. 'COM' (Center of Mass Method: LAT_values_COM)

% --- Define Method Selection ---
MethodSelection = struct();
MethodSelection.MEA1 = 'DIFF';
MethodSelection.MEA2 = 'DIFF';
MethodSelection.MEA3 = 'DIFF';
MethodSelection.TANK = 'COM';

LAT_values_FINAL = struct();

for i = 1:length(cases)
    case_name = cases{i};
    % Get the selected method for the current case
    if isfield(MethodSelection, case_name)
        method = MethodSelection.(case_name);
    else
        method = 'DIFF'; % Default to DIFF
    end
    
    % Select the matrix based on the method
    if strcmp(method, 'COM')
        LAT_matrix_final = LAT_values_COM.(case_name).LAT_matrix;
        fprintf('  -> %s: Using COM Method (LAT_values_COM)\n', case_name);
    elseif strcmp(method, 'DIFF')
        LAT_matrix_final = LAT_values.(case_name).LAT_matrix;
        fprintf('  -> %s: Using Derivative Method (LAT_values)\n', case_name);
    else
        warning('Unknown method "%s" specified for case %s. Defaulting to DIFF.', method, case_name);
        LAT_matrix_final = LAT_values.(case_name).LAT_matrix;
    end
    
    LAT_values_FINAL.(case_name).LAT_matrix = LAT_matrix_final;
end


%% In case of need to manual corrections - Electrical
% Histogram plots
num_bins = 50;
for i = 1:length(cases)
    case_name = cases{i};
    LAT_matrix = LAT_values_FINAL.(case_name).LAT_matrix;
    valid_LATS = LAT_matrix(LAT_matrix ~= 0);
    figure();
    set(gcf, 'color', 'white', 'Position', [600 50 600 400]);
    h = histogram(valid_LATS, num_bins);
    title(['LAT Distribution: Electrical ', case_name, ' (FINAL)'], 'FontSize', 14);
    xlabel('Local Activation Time [ms]', 'FontSize', 12);
    ylabel('Count (Number of Electrodes/Pixels)', 'FontSize', 12);
    grid on;
    mean_lat = mean(valid_LATS);
    std_lat = std(valid_LATS);
    fprintf('  -> %s LAT Stats: Mean = %.2f ms, STD = %.2f ms\n', ...
            case_name, mean_lat, std_lat);
end


% Substitute specific value
case_to_correct = 'MEA1';
LAT_temp = LAT_values_FINAL.(case_to_correct).LAT_matrix;
find_value = 9; % Value to replace
tolerancia = 1;
indices = find(abs(LAT_temp - find_value) < tolerancia);
LAT_temp(indices) = 3; % Value to include
LAT_values_FINAL.(case_to_correct).LAT_matrix = LAT_temp;


% Substitute Higher or lower than
case_to_correct = 'MEA1';
find_value = 9; % Value to replace
new_value = 0;
LAT_values_FINAL.(case_to_correct).LAT_matrix( ...
            LAT_values_FINAL.(case_to_correct).LAT_matrix < ...
            find_value) = new_value;


% Correct The Zero Normalization - Subtract minimum
LAT_values_FINAL.(case_to_correct).LAT_matrix(LAT_values_FINAL.(case_to_correct).LAT_matrix < 0) = 0;
min_LAT = min(LAT_values_FINAL.(case_to_correct).LAT_matrix(LAT_values_FINAL.(case_to_correct).LAT_matrix ~= 0));
if ~isempty(min_LAT) && ~isnan(min_LAT)
    LAT_values_FINAL.(case_to_correct).LAT_matrix(LAT_values_FINAL.(case_to_correct).LAT_matrix ~= 0) = ...
        LAT_values_FINAL.(case_to_correct).LAT_matrix(LAT_values_FINAL.(case_to_correct).LAT_matrix ~= 0) - min_LAT;
end


%% FINAL LAT MAPS - Electrical
C = parula(256); 
% C(1,1:3) = [1 1 1]; % White for background

for i = 1:length(cases)
    case_name = cases{i};
    LAT_matrix_final = LAT_values_FINAL.(case_name).LAT_matrix;
    
    % Determine which method was used
    if isfield(MethodSelection, case_name)
        method = MethodSelection.(case_name);
        method_text = method;
    else
        method = 'DIFF';
        method_text = 'DIFF (default)';
    end
    
    % Create figure
    fig = figure('color', 'white', 'Position', [50 + i*30 50 + i*30 800 600]);
    fig.Name = ['LAT Map - ' case_name];
    
    % INTERPOLATION: Upsample the matrix for smoother visualization
    interpolation_factor = 5; % Increase this for smoother images
    [rows, cols] = size(LAT_matrix_final);
    % Create interpolation grid
    [X, Y] = meshgrid(1:cols, 1:rows);
    [Xq, Yq] = meshgrid(linspace(1, cols, cols*interpolation_factor), ...
                       linspace(1, rows, rows*interpolation_factor));
    % Interpolate the data
    LAT_interp = interp2(X, Y, LAT_matrix_final, Xq, Yq, 'cubic');
    
    % Handle different case geometries
    if strcmp(case_name, 'TANK')
        J = LAT_interp;
        imagesc(J);
        pbaspect([2 1 1]);
    else
        J = imrotate(LAT_interp, 90);
        imagesc(J);
        axis equal;
    end
    
    % Apply colormap
    colormap(C);
    caxis([0 max(LAT_matrix_final(:))]); % Use original max for axis limits
    axis off;
    
    % Add colorbar
    hBar = colorbar('eastoutside');
    ylabel(hBar, 'Local Activation Time [ms]', 'FontSize', 12);
    
    % Add main title
    title_str = sprintf('%s - LAT Map\n[Method: %s] (Interpolated)', case_name, method_text);
    title(title_str, 'FontSize', 16, 'FontWeight', 'bold');
    
    % Add LAT statistics as text annotation
    valid_LATS = LAT_matrix_final(LAT_matrix_final ~= 0);
    if ~isempty(valid_LATS)
        stats_str = sprintf('Min: %.1f ms\nMax: %.1f ms\nMean: %.1f ms\nStd: %.1f ms', ...
            min(valid_LATS), max(valid_LATS), mean(valid_LATS), std(valid_LATS));
        
        % Add text box with statistics
        annotation('textbox', [0.02, 0.03, 0.2, 0.15], ...
            'String', stats_str, ...
            'FontSize', 10, ...
            'BackgroundColor', 'white', ...
            'EdgeColor', 'black', ...
            'Margin', 5);
    end    
    fprintf('Created final LAT map for %s using %s method (interpolated)\n', case_name, method_text);
end


%% COMPARISON FIGURE WITH ALL CASES
if length(cases) <= 6
    fig_all = figure('color', 'white', 'Position', [100 100 1400 800]);
    fig_all.Name = 'LAT Maps - Comparison';
    
    % Determine subplot arrangement
    if length(cases) <= 4
        nrows = 2; ncols = 2;
    else
        nrows = 2; ncols = 3;
    end
    
    % Find global max for consistent color scaling
    all_LATS_global = [];
    for i = 1:length(cases)
        case_name = cases{i};
        LAT_matrix_final = LAT_values_FINAL.(case_name).LAT_matrix;
        valid_LATS = LAT_matrix_final(LAT_matrix_final ~= 0);
        all_LATS_global = [all_LATS_global; valid_LATS];
    end
    if ~isempty(all_LATS_global)
        max_LAT_global = max(all_LATS_global);
    else
        max_LAT_global = 1;
    end
    
    % Plot all cases in subplots
    for i = 1:length(cases)
        case_name = cases{i};
        LAT_matrix_final = LAT_values_FINAL.(case_name).LAT_matrix;
        
        if isfield(MethodSelection, case_name)
            method = MethodSelection.(case_name);
        else
            method = 'DIFF';
        end
        
        subplot(nrows, ncols, i);
        
        % INTERPOLATION for comparison figure too
        interpolation_factor = 3;
        [rows, cols] = size(LAT_matrix_final);
        % Create interpolation grid
        [X, Y] = meshgrid(1:cols, 1:rows);
        [Xq, Yq] = meshgrid(linspace(1, cols, cols*interpolation_factor), ...
                           linspace(1, rows, rows*interpolation_factor));
        % Interpolate the data
        LAT_interp = interp2(X, Y, LAT_matrix_final, Xq, Yq, 'cubic');
        if strcmp(case_name, 'TANK')
            J = LAT_interp;
            imagesc(J);
            pbaspect([2 1 1]);
        else
            J = imrotate(LAT_interp, 90);
            imagesc(J);
            axis equal;
        end
        
        colormap(gca, C);
        caxis([0 max_LAT_global]);
        axis off;
        
        title_str = sprintf('%s\n%s', case_name, method);
        title(title_str, 'FontSize', 12, 'FontWeight', 'bold');
    end
    
    % Add overall colorbar
    hBar = colorbar('Position', [0.93 0.15 0.02 0.7]);
    ylabel(hBar, 'Local Activation Time [ms]', 'FontSize', 12);
    
    sgtitle('LAT Maps Comparison (Interpolated)', 'FontSize', 16, 'FontWeight', 'bold');
end


%% LAT Statistics - Electrical (with ROI Selection)
C = parula(256); 
% C(1,1:3) = [1 1 1]; White Background
for i = 1:length(cases)
    case_name = cases{i};
    LAT_matrix = LAT_values_FINAL.(case_name).LAT_matrix;
    
    % ROI SELECTION
    figure('color', 'white');
    if strcmp(case_name, 'TANK')
        J = LAT_matrix;
        imagesc(J, 'Interpolation', 'bilinear');
        pbaspect([2 1 1]); % 2:1 aspect ratio
    else
        J = LAT_matrix;
        imagesc(J);
        axis equal;
    end
    colormap(C);
    axis off;
    colorbar('eastoutside');
    title(['Case: ', case_name, ' - ROI Selection (FINAL LAT)'], 'FontSize', 14);
    roi_mask = roipoly;
    close(gcf);
    
    % APPLY ROI AND ISOLATE VALID DATA
    LAT_roi = J .* roi_mask; 
    nonzero_LATS_roi = LAT_roi(LAT_roi ~= 0 & ~isnan(LAT_roi));
    % CALCULATE AND STORE STATISTICS
    LAT_mean = mean(nonzero_LATS_roi);
    LAT_std = std(nonzero_LATS_roi);
    LAT_max = max(nonzero_LATS_roi);
    LAT_mode = mode(nonzero_LATS_roi); 
    LAT_var = var(nonzero_LATS_roi); 
    LAT_num = numel(nonzero_LATS_roi); 
    % Displaying results for the current case
    fprintf('\nCase: %s (ROI Stats - FINAL LAT)\n', case_name);
    disp(['  Average LAT:      ', num2str(LAT_mean, '%.4f'), ' ms']);
    disp(['  Mode LAT:         ', num2str(LAT_mode, '%.4f'), ' ms']);
    disp(['  Max LAT:          ', num2str(LAT_max, '%.4f'), ' ms']);
    disp(['  Std Deviation:    ', num2str(LAT_std, '%.4f'), ' ms']);
    disp(['  Variance:         ', num2str(LAT_var, '%.4f')]);
    disp(['  Number of points: ', num2str(LAT_num)]);
    % Store Results
    LAT_values_FINAL.(case_name).LAT_stats = struct(...
        'mean', LAT_mean, 'std', LAT_std, 'max', LAT_max, ...
        'mode', LAT_mode, 'var', LAT_var, 'num_points', LAT_num);
    LAT_values_FINAL.(case_name).roi_mask = roi_mask;
end


%% CV Calculation - Electrical
SpaceScale.TANK = 0.75; % mm/pixel for TANK
SpaceScale.MEA = 0.1;  % mm/pixel for MEA
r = 3;                 % Radius (r) for the CV Circle Method
filter_size = 0;       % Size of the spatial filter kernel
LAT_values_CV = LAT_values_FINAL; % Start CV calculation from FINAL LAT
for i = 1:length(cases)
    case_name = cases{i};
    LAT_E_matrix = LAT_values_FINAL.(case_name).LAT_matrix;
    
    % Determine Space Scale
    if strcmp(case_name, 'TANK')
        current_space_scale = SpaceScale.TANK;
    else
        current_space_scale = SpaceScale.MEA;
    end
    % Interpolation
    [X, Y] = meshgrid(1:size(LAT_E_matrix, 2), 1:size(LAT_E_matrix, 1));
    [Xq, Yq] = meshgrid(1:0.1:size(LAT_E_matrix, 2), 1:0.1:size(LAT_E_matrix, 1));
    % Perform bilinear interpolation
    LAT_E_matrix_interp = interp2(X, Y, LAT_E_matrix, Xq, Yq, 'linear');
    
    % Initialize matrices for CV results
    [R_interp, C_interp] = size(LAT_E_matrix_interp);
    AvgSpeed_E = zeros(R_interp, C_interp);
    StdSpeed_E = zeros(R_interp, C_interp);
    Angle_E = zeros(R_interp, C_interp);
    % CV Circle Method Calculation
    for row = (r + 1):(R_interp - (r + 1))
        for col = (r + 1):(C_interp - (r + 1))
            if LAT_E_matrix_interp(row, col) ~= 0 
                [AvgSpeed_E(row, col), StdSpeed_E(row, col), Angle_E(row, col)] = ...
                    CV_CircleMethod(LAT_E_matrix_interp, r, row, col, current_space_scale);
            end
        end
    end
    % Filtering and Storage
    AvgSpeed_E_filtered = SpatTemp_Filtering(AvgSpeed_E, filter_size, 0, 'GPU');
    
    % Store the results in the FINAL structure
    LAT_values_FINAL.(case_name).AvgSpeed_E_matrix = AvgSpeed_E_filtered;
    LAT_values_FINAL.(case_name).SpaceScale_mm_px = current_space_scale;
    LAT_values_FINAL.(case_name).CV_Radius = r;
    
    fprintf('  -> %s CV map calculated and filtered (Size: %dx%d).\n', case_name, R_interp, C_interp);
end


%% CV Plot - Electric
C = jet(256);
for i = 1:length(cases)
    case_name = cases{i};
    AvgSpeed_E_matrix = LAT_values_FINAL.(case_name).AvgSpeed_E_matrix;
    J_valid = AvgSpeed_E_matrix(AvgSpeed_E_matrix ~= 0 & ~isnan(AvgSpeed_E_matrix));
    % Use percentile limits for robust visualization (e.g., 5th to 95th percentile)
    Y_limits = prctile(J_valid, [5 95], 'all');
    CV_min = Y_limits(1);
    CV_max = Y_limits(2);
    
    % --- Plotting ---
    figure('color', 'white', 'Position', [40 + i*20 40 + i*20 600 600]);
    if strcmp(case_name, 'TANK')
        % TANK case: Rectangular aspect ratio, no rotation
        J_plot = AvgSpeed_E_matrix;
        imagesc(J_plot);
        pbaspect([2 1 1]); % 2:1 aspect ratio
    else
        % MEA cases: Square aspect ratio, rotate 90 degrees
        J_plot = imrotate(AvgSpeed_E_matrix, 90);
        imagesc(J_plot);
        axis equal; % Ensures square aspect ratio
    end
    % Apply the calculated color limits
    caxis([CV_min CV_max]); 
    colormap(C);
    axis off;
    hBar = colorbar('eastoutside');
    ylabel(hBar, 'Conduction Velocity [cm/s]', 'FontSize', 14);
    title(['Electrical ' case_name ' - Conduction Velocity Map (FINAL LAT)'], 'FontSize', 16);
    fprintf('  -> %s CV Map plotted with limits [%.2f - %.2f] cm/s.\n', case_name, CV_min, CV_max);
end


%% CV Statistics - Electrical
for i = 1:length(cases)
    case_name = cases{i};
    
    % --- Retrieve the CV matrix directly from the struct ---
    if ~isfield(LAT_values_FINAL.(case_name), 'AvgSpeed_E_matrix')
        fprintf('Case %s: Warning - AvgSpeed_E_matrix field not found. Skipping CV stats.\n', case_name);
        continue; 
    end
    
    CV_matrix = LAT_values_FINAL.(case_name).AvgSpeed_E_matrix;
    
    % Check if the matrix is empty before proceeding
    if isempty(CV_matrix)
        fprintf('Case %s: Warning - CV matrix is empty. Skipping CV stats.\n', case_name);
        continue;
    end
    
    % Isolate valid CV points (non-zero and non-NaN)
    nonzero_CV = CV_matrix(CV_matrix ~= 0 & ~isnan(CV_matrix) & ~isinf(CV_matrix));
    
    if isempty(nonzero_CV)
        fprintf('Case %s: Warning - No valid non-zero CV points found. Skipping stats.\n', case_name);
        LAT_values_FINAL.(case_name).CV_stats = []; 
        continue;
    end
    
    % Calculate Statistics
    CV_mean = mean(nonzero_CV);
    CV_std = std(nonzero_CV);
    CV_max = max(nonzero_CV);
    CV_mode = mode(nonzero_CV); 
    CV_var = var(nonzero_CV); 
    CV_num = numel(nonzero_CV);
    
    % Store Results in the LAT_values_FINAL structure
    LAT_values_FINAL.(case_name).CV_stats = struct(...
        'mean', CV_mean, 'std', CV_std, 'max', CV_max, ...
        'mode', CV_mode, 'var', CV_var, 'num_points', CV_num);
    
    % Displaying results for the current case
    fprintf('\nCase: %s (FINAL LAT)\n', case_name);
    disp(['  Average CV:       ', num2str(CV_mean, '%.4f'), ' cm/s']);
    disp(['  Mode CV:          ', num2str(CV_mode, '%.4f'), ' cm/s']);
    disp(['  Max CV:           ', num2str(CV_max, '%.4f'), ' cm/s']);
    disp(['  Std Deviation:    ', num2str(CV_std, '%.4f'), ' cm/s']);
    disp(['  Variance:         ', num2str(CV_var, '%.4f')]);
    disp(['  Number of points: ', num2str(CV_num)]);
end


%% FINAL SAVING OF ALL RESULTS
% Save ALL essential variables for complete reproducibility

% 1. CORE RESULTS STRUCTURES
save_vars = struct();

% 1.1 Original Data and Configuration
save_vars.Data_E_Structure = Data_E_Structure;
save_vars.Fsampling = Fsampling;
save_vars.cases = cases;
save_vars.electrical_space_scale = electrical_space_scale;
save_vars.tank_space_scale = tank_space_scale;

% 1.2 Analysis Windows
save_vars.sample_limits = sample_limits;

% 1.3 Method Selection
save_vars.MethodSelection = MethodSelection;

% 1.4 All LAT Results (Raw, COM, Final)
save_vars.LAT_values = LAT_values;          % Derivative method results
save_vars.LAT_values_COM = LAT_values_COM;  % COM method results  
save_vars.LAT_values_FINAL = LAT_values_FINAL; % Consolidated final results

% 1.5 CV Results (already in LAT_values_FINAL, but keep separate for clarity)
if isfield(LAT_values_FINAL.(cases{1}), 'AvgSpeed_E_matrix')
    save_vars.CV_results = struct();
    for i = 1:length(cases)
        case_name = cases{i};
        save_vars.CV_results.(case_name).CV_matrix = LAT_values_FINAL.(case_name).AvgSpeed_E_matrix;
        if isfield(LAT_values_FINAL.(case_name), 'CV_stats')
            save_vars.CV_results.(case_name).CV_stats = LAT_values_FINAL.(case_name).CV_stats;
        end
        save_vars.CV_results.(case_name).SpaceScale = LAT_values_FINAL.(case_name).SpaceScale_mm_px;
        save_vars.CV_results.(case_name).CV_Radius = LAT_values_FINAL.(case_name).CV_Radius;
    end
end

% 1.6 Statistics Table - CREATE IF IT DOESN'T EXIST
if ~exist('StatsTable', 'var')
    % Create the statistics table from LAT_values_FINAL
    StatsTable = table();
    for i = 1:length(cases)
        case_name = cases{i};
        
        % Get LAT statistics
        if isfield(LAT_values_FINAL.(case_name), 'LAT_stats')
            lat_stats = LAT_values_FINAL.(case_name).LAT_stats;
            lat_mean = lat_stats.mean;
            lat_std = lat_stats.std;
            lat_max = lat_stats.max;
            lat_num = lat_stats.num_points;
        else
            % Calculate from matrix if stats don't exist
            LAT_matrix = LAT_values_FINAL.(case_name).LAT_matrix;
            valid_LATS = LAT_matrix(LAT_matrix ~= 0);
            if ~isempty(valid_LATS)
                lat_mean = mean(valid_LATS);
                lat_std = std(valid_LATS);
                lat_max = max(valid_LATS);
                lat_num = numel(valid_LATS);
            else
                lat_mean = NaN;
                lat_std = NaN;
                lat_max = NaN;
                lat_num = 0;
            end
        end
        
        % Get CV statistics
        if isfield(LAT_values_FINAL.(case_name), 'CV_stats')
            cv_stats = LAT_values_FINAL.(case_name).CV_stats;
            cv_mean = cv_stats.mean;
            cv_std = cv_stats.std;
            cv_max = cv_stats.max;
            cv_num = cv_stats.num_points;
        else
            % Set to NaN if CV stats don't exist
            cv_mean = NaN;
            cv_std = NaN;
            cv_max = NaN;
            cv_num = 0;
        end
        
        % Create table row
        NewRow = table(...
            {case_name}, ...
            lat_mean, lat_std, lat_max, lat_num, ...
            cv_mean, cv_std, cv_max, cv_num, ...
            'VariableNames', {'Case', 'LAT_Mean_ms', 'LAT_STD_ms', 'LAT_Max_ms', 'LAT_N_Pixels', ...
                             'CV_Mean_cms', 'CV_STD_cms', 'CV_Max_cms', 'CV_N_Points'});
        
        StatsTable = [StatsTable; NewRow];
    end
    fprintf('Created StatsTable from LAT_values_FINAL data.\n');
end
save_vars.StatsTable = StatsTable;

% 1.7 ROI Masks (if they exist)
save_vars.roi_masks = struct();
for i = 1:length(cases)
    case_name = cases{i};
    if isfield(LAT_values_FINAL.(case_name), 'roi_mask')
        save_vars.roi_masks.(case_name) = LAT_values_FINAL.(case_name).roi_mask;
    end
end

% 2. ADD ANALYSIS PARAMETERS (important for reproducibility)
save_vars.AnalysisParams = struct();
if exist('SpaceScale', 'var')
    save_vars.AnalysisParams.SpaceScale = SpaceScale;
end
% Add other parameters that exist in your workspace
vars_to_check = {'debug_LAT', 'debug_COM'};
for var_idx = 1:length(vars_to_check)
    var_name = vars_to_check{var_idx};
    if exist(var_name, 'var')
        save_vars.AnalysisParams.(var_name) = eval(var_name);
    end
end

% 3. ADD METADATA
save_vars.Metadata = struct();
save_vars.Metadata.analysis_date = datestr(now, 'yyyy-mm-dd HH:MM:SS');

% 4. DEFINE FILENAME
filename = sprintf('E_LAT_CV.mat');

% 5. SAVE ALL VARIABLES
save(filename, 'save_vars', '-v7.3'); % Use -v7.3 for large files
fprintf('\n✅ COMPLETE ANALYSIS SAVED TO: %s\n', filename);

% 6. EXPORT STATISTICS TO CSV
if ~isempty(StatsTable)
    csv_filename = sprintf('E_Statistics.csv');
    writetable(StatsTable, csv_filename);
    fprintf('✅ STATISTICS EXPORTED TO CSV: %s\n', csv_filename);
end