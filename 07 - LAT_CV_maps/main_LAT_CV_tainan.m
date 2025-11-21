%% LAT and CV Maps - Optical & Electrical
clear; clc;


%% Loading Preview Data
% load(""); % Load Synchronized data
Fsampling = 4000;
% --- Optical Data & Parameters ---
Data_O = D_SYNC.CAM2;
ROI_O = D_SYNC.ROI.ROI_2;
optical_space_scale = 0.33; % mm per pixel (used for CV calculation)
% Preview Signal
Background = squeeze(Data_O(:,:,2000));
pick_up_a_trace(Background, Data_O,1);


%% Time Selection & Preview - Optical
optical_lim1 = 12508;
optical_lim2 = 13340;
Data_temp_O = Data_O(:,:,optical_lim1:optical_lim2);
Background = squeeze(Data_O(:,:,optical_lim1));
pick_up_a_trace(Background, Data_temp_O, 1);


%% LAT Calculation - Optical (Optical Analysis Part 1)
linear_fit_length = 15;
PCL = 200;
debug_LAT = 0;

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
OpticalResults.LAT = LAT_O_filtered;


%% LAT Map Plot - Optical
levelStep = 20;
Title = 'Optical - LAT';
f1 = figure('color', 'white', 'Position', [40 40 600 600]);
C = parula(256);
C(1,1:3) = [1 1 1]; % White for background
C(256,1:3) = [0.5, 0.5, 0.5]; % Gray for the max value
J = imrotate(LAT_O_filtered, 90);
contourf(flipud(J), levelStep);
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
find_value = 20; 
tolerancia = 1e-10;
new_value = 20; % NaN to exclude

% Specific Value
indices = find(abs(LAT_O_filtered - find_value) < tolerancia);
if ~isempty(indices)
    LAT_O_filtered(indices) = new_value;
    OpticalResults.LAT = LAT_O_filtered;
end

% Higher or Lower than
LAT_O_filtered(LAT_O_filtered < find_value) = new_value;


%% LAT Statistics - Optical
% Display the image and let the user define the ROI
figure();
imshow(LAT_O_filtered, C);
title('Select ROI for DF Analysis');
roi_lat_O = roipoly;

% Apply the ROI
LAT_O_roi = LAT_O_filtered .* roi_lat_O;
nonzero_LAT_O = LAT_O_roi(LAT_O_roi ~= 0 & ~isnan(LAT_O_roi));
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
levelStep = 10;
figure();
imagesc(LAT_O_filtered);
colormap(parula); colorbar;
title('Select ROI');
roi_specific = roipoly; 
% Apply the specific ROI
LAT_specific = LAT_O_filtered .* roi_specific;

% Subtract the minimum LAT *within* this specific ROI
nonzero_specific = LAT_specific(LAT_specific ~= 0 & ~isnan(LAT_specific));
if ~isempty(nonzero_specific)
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
else
    fprintf('No valid data found in the specific ROI. Plot skipped.\n');
end


%% CV Calculation - Optical
cv_radius_optical = 5; % in pixels

AvgSpeed_O = zeros(size(LAT_O_filtered));
for i = 1:size(LAT_O_filtered,1)
    for j = 1:size(LAT_O_filtered,2)
        if roi_lat_O(i,j) == 1 && LAT_O_filtered(i,j) ~= 0 && ~isnan(LAT_O_filtered(i,j))
            [AvgSpeed_O(i,j), ~, ~] = ...
                CV_CircleMethod(LAT_O_filtered, cv_radius_optical, i, j, optical_space_scale);
        end
    end
end

% Filter CV data
AvgSpeed_O_filtered = SpatTemp_Filtering(AvgSpeed_O, 5, 0, 'GPU');
OpticalResults.CV = AvgSpeed_O_filtered;


%% CV Map Plot - Optical
Title = 'Optical - CV';
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


%% CV Statistics - Optical
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


%% Optical Results Table
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


%% Summary Visualizations - Optical
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


%% 
%% 
%%


%% Loading Data
load(""); % Load Interpolated data


%% Configuring
Data_E_Structure = InterpSignal.Data; 
electrical_space_scale = 2;  % mm per pixel (MEA)
tank_space_scale = 3;       % mm per pixel (TANK)
cases = {'MEA1', 'MEA2', 'MEA3', 'TANK'};
ElectricalResults = struct();
Fsampling = 4000;


%% LAT Calculation - Electrical (Applied to all cases at once, like DF)
electrical_lim1 = 12508;
electrical_lim2 = 13340;
debug_LAT = 0;
LAT_values = struct();

for i = 1:length(cases)
    case_name = cases{i};
    
    % Extract the current case's 3D data 
    Data_current_E = Data_E_Structure.(case_name);
    Data_temp_E = Data_current_E(:,:,electrical_lim1:electrical_lim2);
    relative_lim1 = 1;
    relative_lim2 = size(Data_temp_E, 3);
    
    % Calculate LAT map (pixel by pixel on the 3D data)
    LAT_map = zeros(size(Data_temp_E, 1), size(Data_temp_E, 2));
    for x = 1:size(Data_temp_E, 1)
        for y = 1:size(Data_temp_E, 2)
            signal = squeeze(Data_temp_E(x,y,:));
            if max(signal) ~= 0
                signal_2D = signal'; 
                LAT_map(x,y) = find_LAT_diff(signal_2D, ...
                            Fsampling, relative_lim1, ...
                            relative_lim2, debug_LAT);
            end
        end
    end
    
    % Post-processing (Normalization and Filtering)
    LAT_map(LAT_map < 0) = 0;
    % Subtract minimum LAT for re-normalization
    min_LAT = min(LAT_map(LAT_map ~= 0));
    if ~isempty(min_LAT) && ~isnan(min_LAT)
        LAT_map(LAT_map ~= 0) = LAT_map(LAT_map ~= 0) - min_LAT;
    end
    % Spatial Filtering
    LAT_map_filtered = SpatTemp_Filtering(LAT_map, 3, 0, 'GPU'); 
    
    % Store the final, filtered 2D LAT map
    LAT_values.(case_name).LAT_matrix = LAT_map_filtered;
    fprintf('  -> %s LAT matrix (Size: %dx%d) stored.\n', case_name, size(LAT_map_filtered, 1), size(LAT_map_filtered, 2));
end


%% LAT MAP PLOTTING - Electrical
C = parula(256); C(1,1:3) = [1 1 1]; % White for background (0 LAT)
for i = 1:length(cases)
    case_name = cases{i};
    LAT_matrix = LAT_values.(case_name).LAT_matrix;
    figure('color', 'white', 'Position', [40 40 600 600]);
    if strcmp(case_name, 'TANK')
        % TANK case: Rectangular aspect ratio
        J = LAT_matrix;
        imagesc(J);
        pbaspect([2 1 1]); % 2:1 aspect ratio
    else
        % MEA cases: Square aspect ratio (typically rotated 90 degrees)
        J = imrotate(LAT_matrix, 90);
        imagesc(J);
        axis equal;
    end
    colormap(C);
    axis off;
    hBar = colorbar('eastoutside');
    ylabel(hBar, 'Local Activation Time [ms]', 'FontSize', 14);
    title(['Electrical ' case_name ' - Local Activation Time'], 'FontSize', 16);
end


