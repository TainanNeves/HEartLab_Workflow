%% DF CL and OI Maps
clear; clc;


%% Loading Data
load(""); % Load Synchronized data


%% DF Calculation - Optical
Data = D_SYNC.CAM1;
ROI = D_SYNC.ROI.ROI_1;
Fsampling = 4000;

% Set upper and lower frequency bounds for the Dominant Frequency calculation
freq_up = 20;
freq_down = 0.5;

% Extract a subset of the optical data for analysis (adjust the sample range accordingly)
in_sample = 2*4000;
end_sample = 4*4000;
Data_temp = Data(:,:,in_sample:end_sample);

% Perform Dominant Frequency analysis on the subset
[DF_O, Sfft_O, fstep] = f_DF_optico(Data_temp, Fsampling, freq_up, freq_down);


%% Spectrum of multiple pixels - Optical
% Select points
Background = squeeze(Data(:,:,2000));
[px, py] = pick_up_a_trace(Background, Data, 1);
n_points = length(px);
% Create legend names
point_names = cell(1, n_points);
for i = 1:n_points
    point_names{i} = ['P' num2str(i)];
end

% Plot frequency spectra for all selected points
b = length(Sfft_O(1, 1, :));
figure('Position', [200 200 1000 600]);
% Define colors
colors = lines(n_points);
if n_points > 7
    colors = jet(n_points);
end
% Plot all spectra
hold on;
for i = 1:n_points
    % Check if coordinates are within valid range
    if px(i) >= 1 && px(i) <= size(Sfft_O,1) && py(i) >= 1 && py(i) <= size(Sfft_O,2)
        spectrum = squeeze(Sfft_O(px(i), py(i), :));
        plot((1:b) * fstep, spectrum, 'LineWidth', 2, 'Color', colors(i,:), 'DisplayName', point_names{i});
        % Mark dominant frequency with a star
        [max_power, max_idx] = max(spectrum);
        dominant_freq = max_idx * fstep;
        plot(dominant_freq, max_power, 'p', 'MarkerSize', 10, 'MarkerFaceColor', colors(i,:), ...
             'MarkerEdgeColor', 'k', 'LineWidth', 1, 'DisplayName', [point_names{i} ' DF=' num2str(dominant_freq, '%.2f') 'Hz']);
    else
        warning(['Point P' num2str(i) ' coordinates (' num2str(px(i)) ',' num2str(py(i)) ') are out of bounds.']);
    end
end
% Add plot formatting
legend('show', 'FontSize', 10, 'Location', 'best');
xlabel('Frequency [Hz]', 'FontSize', 12);
ylabel('Power Spectral Density [mV^2/Hz]', 'FontSize', 12);
title(['Frequency Spectrum - ' num2str(n_points) ' Selected Points'], 'FontSize', 14);
set(gca, 'fontsize', 12);
xlim([freq_down freq_up]);
grid on;
% Display point coordinates in command window
fprintf('\n=== Selected Points ===\n');
for i = 1:n_points
    if px(i) >= 1 && px(i) <= size(Sfft_O,1) && py(i) >= 1 && py(i) <= size(Sfft_O,2)
        spectrum = squeeze(Sfft_O(px(i), py(i), :));
        [max_power, max_idx] = max(spectrum);
        dominant_freq = max_idx * fstep;
        fprintf('P%d: Coordinates (%d, %d) - Dominant Frequency: %.2f Hz\n', i, px(i), py(i), dominant_freq);
    else
        fprintf('P%d: Coordinates (%d, %d) - INVALID (out of bounds)\n', i, px(i), py(i));
    end
end

% Individual Plots for better visualization
if n_points > 1
    figure('Position', [300 100 1200 800]);
    % Calculate subplot arrangement
    cols = ceil(sqrt(n_points));
    rows = ceil(n_points / cols);
    for i = 1:n_points
        if px(i) >= 1 && px(i) <= size(Sfft_O,1) && py(i) >= 1 && py(i) <= size(Sfft_O,2)
            subplot(rows, cols, i);
            spectrum = squeeze(Sfft_O(px(i), py(i), :));
            plot((1:b) * fstep, spectrum, 'LineWidth', 2, 'Color', colors(i,:));
            
            % Mark dominant frequency
            [max_power, max_idx] = max(spectrum);
            dominant_freq = max_idx * fstep;
            hold on;
            plot(dominant_freq, max_power, 'p', 'MarkerSize', 8, 'MarkerFaceColor', colors(i,:), ...
                 'MarkerEdgeColor', 'k', 'LineWidth', 1);
            
            title([point_names{i} ' (DF=' num2str(dominant_freq, '%.2f') 'Hz)'], 'FontSize', 10);
            xlabel('Frequency [Hz]');
            ylabel('Power Spectral Density [mV^2/Hz]', 'FontSize', 12);
            xlim([freq_down freq_up]);
            grid on;
        end
    end
    sgtitle('Individual Point Spectra', 'FontSize', 14);
end


%% DF Map - Optical
lim = [freq_down freq_up]; % Default [freq_down freq_up]
Title = 'Cam1 V';
C = jet(256);
C(1,1:3) = [1 1 1]; % White for background
C(256,1:3) = [0.5, 0.5, 0.5]; % Gray for the max value
d = 200;  % Define the color for the background
color = 'black';  % Specify color for text
f1 = figure('color', 'white', 'Position', [50 50 500 500]);
J = DF_O;
% Set zeros inside ROI to 0.5
% J(D_OP.ROI.ROI_2 & J == 0) = UP; %Put Max value in zeros inside the ROI
J = imrotate(J, 90);
% Plot
imagesc(J);
colormap(C);
box off;
set(gca, 'fontsize', 18);
ylabel('Pixels');
xlabel('Pixels');
axis off;
hBar1 = colorbar('eastoutside');
ylabel(hBar1, 'Dominant Frequency [Hz]', 'FontSize', 14);
caxis(lim);
title(Title);


%% In case of need to manual corrections
% Susbstitude values
find_value = 0.5; % Value to replace
tolerancia = 1e-10;
indices = find(abs(DF_O - find_value) < tolerancia);
DF_O(indices) = 0; % Value to include


%% DF CL Statistics - Optical
% Display the image and let the user define the ROI
figure();
imshow(DF_O);
title('Select ROI for DF Analysis');
roi_df = roipoly;

% Apply the ROI to your matrix (assuming DF_O is your matrix)
DF_O_roi = DF_O .* roi_df;

% Metrics within the selected ROI for non-zero values
nonzero_values_df = DF_O_roi(DF_O_roi ~= 0);

% Calculate Cycle Length (CL) in ms
CL_O_roi = (1 ./ nonzero_values_df) * 1000; % Convert to ms

% DF Statistics
HCL_roi = max(nonzero_values_df);
avg_roi = mean(nonzero_values_df);
mod_roi = mode(nonzero_values_df);
std_roi = std(nonzero_values_df);
var_roi = var(nonzero_values_df);
number_points_df = numel(nonzero_values_df);

% CL Statistics
HCL_ms = max(CL_O_roi);
avg_CL_ms = mean(CL_O_roi);
mod_CL_ms = mode(CL_O_roi);
std_CL_ms = std(CL_O_roi);
var_CL_ms = var(CL_O_roi);

disp('=== Dominant Frequency (DF) Statistics ===');
disp(['Higher DF in ROI: ', num2str(HCL_roi), ' Hz']);
disp(['Average DF in ROI: ', num2str(avg_roi), ' Hz']);
disp(['Mode DF in ROI: ', num2str(mod_roi), ' Hz']);
disp(['Standard Deviation of DF in ROI: ', num2str(std_roi), ' Hz']);
disp(['Variance of DF in ROI: ', num2str(var_roi)]);
disp(['Number of pixels in the ROI: ', num2str(number_points_df)]);

disp('=== Cycle Length (CL) Statistics ===');
disp(['Highest CL in ROI: ', num2str(HCL_ms), ' ms']);
disp(['Average CL in ROI: ', num2str(avg_CL_ms), ' ms']);
disp(['Mode CL in ROI: ', num2str(mod_CL_ms), ' ms']);
disp(['Standard Deviation of CL in ROI: ', num2str(std_CL_ms), ' ms']);


%% Optic Organization Index - OI
% Variables
MFFTi = DF_O;
Sffti = Sfft_O;
% fstep = fstep; % Already defined above
% Hzi = freq_down; % Already defined above
% Hzf = freq_up; % Already defined above
dfh_threshold_area = 0.6; % Threshold for dominant frequency harmonic area (0-1)
f_mode = 2; % Frequency mode: 1 = fundamental frequency, 2 = harmonic analysis
debug = 1; % Debug mode: 1 = show plots, 0 = no plots

% Calculating
OI = calculate_OI(MFFTi, Sffti, fstep, freq_down, freq_up, dfh_threshold_area, f_mode, debug);

% Statistics
% Select Electrodes range
figure(); imshow(OI, 'InitialMagnification', 'fit');
title('Select ROI for OI Analysis'); 
roi_oi = roipoly;

% Metrics within the selected ROI for non-zero values
nonzero_values_oi = OI(roi_oi);
nonzero_values_oi = nonzero_values_oi(nonzero_values_oi ~= 0);

% OI Statistics
HOI_roi = max(nonzero_values_oi);
avg_oi = mean(nonzero_values_oi);
mod_oi = mode(nonzero_values_oi);
std_oi = std(nonzero_values_oi);
var_oi = var(nonzero_values_oi);
number_points_oi = numel(nonzero_values_oi);

disp('=== Organization Index (OI) Statistics ===');
disp(['Higher OI in ROI: ', num2str(HOI_roi)]);
disp(['Average OI in ROI: ', num2str(avg_oi)]);
disp(['Mode OI in ROI: ', num2str(mod_oi)]);
disp(['Standard Deviation of OI in ROI: ', num2str(std_oi)]);
disp(['Variance of OI in ROI: ', num2str(var_oi)]);
disp(['Number of pixels in the ROI: ', num2str(number_points_oi)]);


%% Create Results Table
% Create table with all statistics
Parameter = {'Dominant Frequency (Hz)'; 'Cycle Length (ms)'; 'Organization Index'};
Max_Value = [HCL_roi; HCL_ms; HOI_roi];
Mean_Value = [avg_roi; avg_CL_ms; avg_oi];
Mode_Value = [mod_roi; mod_CL_ms; mod_oi];
Std_Deviation = [std_roi; std_CL_ms; std_oi];
Variance = [var_roi; var_CL_ms; var_oi];
Number_Points = [number_points_df; number_points_df; number_points_oi];

ResultsTable = table(Parameter, Max_Value, Mean_Value, Mode_Value, Std_Deviation, Variance, Number_Points);

% Display the table
disp(' ');
disp('=== RESULTS TABLE ===');
disp(ResultsTable);


%% Summary Visualizations
% Create summary figure
figure('Position', [100 100 1200 800]);

% Subplot 1: DF Map
subplot(2,3,1);
imagesc(imrotate(DF_O, 90));
colormap(jet);
colorbar;
title('Dominant Frequency (Hz)');
axis off;

% Subplot 2: CL Map
subplot(2,3,2);
CL_map = (1 ./ DF_O) * 1000;
CL_map(DF_O == 0) = 0; % Handle division by zero
imagesc(imrotate(CL_map, 90));
colormap(jet);
colorbar;
title('Cycle Length (ms)');
axis off;

% Subplot 3: OI Map
subplot(2,3,3);
imagesc(imrotate(OI, 90));
colormap(jet);
colorbar;
title('Organization Index');
axis off;

% Subplot 4: DF Histogram
subplot(2,3,4);
histogram(nonzero_values_df, 20, 'FaceColor', 'blue', 'FaceAlpha', 0.7);
xlabel('Dominant Frequency (Hz)');
ylabel('Count');
title('DF Distribution');
grid on;

% Subplot 5: CL Histogram
subplot(2,3,5);
histogram(CL_O_roi, 20, 'FaceColor', 'red', 'FaceAlpha', 0.7);
xlabel('Cycle Length (ms)');
ylabel('Count');
title('CL Distribution');
grid on;

% Subplot 6: OI Histogram
subplot(2,3,6);
histogram(nonzero_values_oi, 20, 'FaceColor', 'green', 'FaceAlpha', 0.7);
xlabel('Organization Index');
ylabel('Count');
title('OI Distribution');
grid on;

sgtitle('Optical Analysis Summary');


%% Save workspace variables for future use
save('O_DF_CL_OI_CAM.mat', 'freq_up', 'freq_down', 'in_sample', 'end_sample', ...
    'DF_O', 'Sfft_O', 'fstep', 'CL_map', 'OI', 'roi_df', 'roi_oi', 'ResultsTable');

disp(' ');
disp('=== OPTICAL ANALYSIS COMPLETE ===');
disp('All results have been saved to:');
disp('- optical_analysis_results.mat (MATLAB workspace)');


%%
%%
%%


%% DF CL and OI Maps - Electrical
clear; clc;


%% Loading Data
load(""); % Load Interpolated data


%% DF Calculation - Electrical
% Parameters
Fsampling = 4000;
freq_up = 20;
freq_down = 0.5;
in_sample = 2*4000;
end_sample = 4*4000;

% Loading Data
cases = {'MEA1', 'MEA2', 'MEA3', 'TANK'};
Data_temp = struct();
DF_values = struct();

for i = 1:length(cases)
    case_name = cases{i};
    Data = InterpSignal.Sync.(case_name);
    Data_temp.(case_name) = Data(:,:,in_sample:end_sample);
    
    % Perform Dominant Frequency analysis
    [MFFTi, Sffti, fstep] = f_DF_electric(Data_temp.(case_name), Fsampling, freq_up, freq_down);
    
    DF_values.(case_name).MFFTi = MFFTi;
    DF_values.(case_name).Sffti = Sffti;
end
DF_values.fstep = fstep;


%% Spectrum of multiple electrodes - Electrical
% Select which case to analyze
current_case = 'TANK'; % Change to 'MEA1', 'MEA2', 'MEA3', or 'TANK'
Data = InterpSignal.Sync.(current_case);
Background = squeeze(Data(:,:,2000));

% Select points using your function
[px, py] = pick_up_a_trace(Background, Data, 1);
n_points = length(px);

% Create legend names
point_names = cell(1, n_points);
for i = 1:n_points
    point_names{i} = ['P' num2str(i)];
end

% Plot frequency spectra for all selected points
source = DF_values.(current_case);
b = length(source.Sffti(1, 1, :));
figure('Position', [200 200 1000 600]);

% Define colors
colors = lines(n_points);
if n_points > 7
    colors = jet(n_points);
end

% Plot all spectra
hold on;
for i = 1:n_points
    % Check if coordinates are within valid range
    if px(i) >= 1 && px(i) <= size(source.Sffti,1) && py(i) >= 1 && py(i) <= size(source.Sffti,2)
        spectrum = squeeze(source.Sffti(px(i), py(i), :));
        plot((1:b) * fstep, spectrum, 'LineWidth', 2, 'Color', colors(i,:), 'DisplayName', point_names{i});
        
        % Mark dominant frequency with a star
        [max_power, max_idx] = max(spectrum);
        dominant_freq = max_idx * fstep;
        plot(dominant_freq, max_power, 'p', 'MarkerSize', 10, 'MarkerFaceColor', colors(i,:), ...
             'MarkerEdgeColor', 'k', 'LineWidth', 1, 'DisplayName', [point_names{i} ' DF=' num2str(dominant_freq, '%.2f') 'Hz']);
    else
        warning(['Point P' num2str(i) ' coordinates (' num2str(px(i)) ',' num2str(py(i)) ') are out of bounds.']);
    end
end

% Add plot formatting
legend('show', 'FontSize', 10, 'Location', 'best');
xlabel('Frequency [Hz]', 'FontSize', 12);
ylabel('Power Spectral Density [mV^2/Hz]', 'FontSize', 12);
title(['Frequency Spectrum - ' current_case ' - ' num2str(n_points) ' Selected Points'], 'FontSize', 14);
set(gca, 'fontsize', 12);
xlim([freq_down freq_up]);
grid on;

% Display point coordinates in command window
fprintf('\n=== Selected Points - %s ===\n', current_case);
for i = 1:n_points
    if px(i) >= 1 && px(i) <= size(source.Sffti,1) && py(i) >= 1 && py(i) <= size(source.Sffti,2)
        spectrum = squeeze(source.Sffti(px(i), py(i), :));
        [max_power, max_idx] = max(spectrum);
        dominant_freq = max_idx * fstep;
        fprintf('P%d: Coordinates (%d, %d) - Dominant Frequency: %.2f Hz\n', i, px(i), py(i), dominant_freq);
    else
        fprintf('P%d: Coordinates (%d, %d) - INVALID (out of bounds)\n', i, px(i), py(i));
    end
end

% Individual Plots for better visualization
if n_points > 1
    figure('Position', [300 100 1200 800]);
    % Calculate subplot arrangement
    cols = ceil(sqrt(n_points));
    rows = ceil(n_points / cols);
    
    for i = 1:n_points
        if px(i) >= 1 && px(i) <= size(source.Sffti,1) && py(i) >= 1 && py(i) <= size(source.Sffti,2)
            subplot(rows, cols, i);
            spectrum = squeeze(source.Sffti(px(i), py(i), :));
            plot((1:b) * fstep, spectrum, 'LineWidth', 2, 'Color', colors(i,:));
            
            % Mark dominant frequency
            [max_power, max_idx] = max(spectrum);
            dominant_freq = max_idx * fstep;
            hold on;
            plot(dominant_freq, max_power, 'p', 'MarkerSize', 8, 'MarkerFaceColor', colors(i,:), ...
                 'MarkerEdgeColor', 'k', 'LineWidth', 1);
            
            title([point_names{i} ' (DF=' num2str(dominant_freq, '%.2f') 'Hz)'], 'FontSize', 10);
            xlabel('Frequency [Hz]');
            ylabel('Power Spectral Density [mV^2/Hz]', 'FontSize', 12);
            xlim([freq_down freq_up]);
            grid on;
        end
    end
    sgtitle(['Individual Point Spectra - ' current_case], 'FontSize', 14);
end


%% In case of need to manual corrections
% Apply corrections to specific case BEFORE CL calculation
case_to_correct = 'MEA1'; % Change as needed
DF_E = DF_values.(case_to_correct).MFFTi;

% Substitute specific value in DF
find_value = 9; % Value to replace
tolerancia = 1;
indices = find(abs(DF_E - find_value) < tolerancia);
DF_E(indices) = 3; % Value to include

% Update the DF values (CL will be calculated from this corrected DF)
DF_values.(case_to_correct).MFFTi = DF_E;


%% DF Map - Electrical
% Create DF maps for all cases
lim = [freq_down freq_up];
C = jet(256);
C(1,1:3) = [1 1 1]; % White for background

for i = 1:length(cases)
    case_name = cases{i};
    DF_E = DF_values.(case_name).MFFTi;
    
    figure('color', 'white', 'Position', [50 50 500 500]);
    
    if strcmp(case_name, 'TANK')
        % Rectangular aspect ratio for TANK
        imagesc(DF_E);
        pbaspect([2 1 1]); % 2:1 aspect ratio
    else
        % Square for MEAs
        J = imrotate(DF_E, 90);
        imagesc(J);
        axis equal;
    end
    
    colormap(C);
    box off;
    set(gca, 'fontsize', 18);
    ylabel('Electrodes');
    xlabel('Electrodes');
    axis off;
    
    hBar1 = colorbar('eastoutside');
    ylabel(hBar1, 'Dominant Frequency [Hz]', 'FontSize', 14);
    caxis(lim);
    title([case_name ' - Dominant Frequency']);
end


%% CL Calculation - Electrical
% Calculate CL maps for all cases from the (corrected) DF values
CL_values = struct();

for i = 1:length(cases)
    case_name = cases{i};
    DF_E = DF_values.(case_name).MFFTi;
    
    % Calculate Cycle Length maps from DF
    CL_map = (1 ./ DF_E) * 1000; % Convert to ms
    CL_map(DF_E == 0) = 0; % Handle division by zero
    CL_values.(case_name) = CL_map; % Store CL map for each case
end


%% CL Map - Electrical
% Create CL maps for all cases using CL_values
cl_lim = [0 300]; % Typical CL range in ms
C_cl = jet(256);
C_cl(1,1:3) = [1 1 1]; % White for background

for i = 1:length(cases)
    case_name = cases{i};
    CL_E = CL_values.(case_name);
    
    figure('color', 'white', 'Position', [50 50 500 500]);
    
    if strcmp(case_name, 'TANK')
        % Rectangular aspect ratio for TANK
        imagesc(CL_E);
        pbaspect([2 1 1]); % 2:1 aspect ratio
    else
        % Square for MEAs
        J = imrotate(CL_E, 90);
        imagesc(J);
        axis equal;
    end
    
    colormap(C_cl);
    box off;
    set(gca, 'fontsize', 18);
    ylabel('Electrodes');
    xlabel('Electrodes');
    axis off;
    
    hBar1 = colorbar('eastoutside');
    ylabel(hBar1, 'Cycle Length [ms]', 'FontSize', 14);
    caxis(cl_lim);
    title([case_name ' - Cycle Length']);
end


%% DF CL Statistics - Electrical
% Analyze each case separately
Results_all = struct();

for i = 1:length(cases)
    case_name = cases{i};
    DF_E = DF_values.(case_name).MFFTi;
    CL_E = CL_values.(case_name); % Use stored CL values
    
    fprintf('\n=== Analyzing %s ===\n', case_name);
    
    % Display the image and let the user define the ROI
    figure();
    imshow(DF_E);
    title(['Select ROI for DF Analysis - ' case_name]);
    roi_df = roipoly;

    % Apply the ROI
    DF_E_roi = DF_E .* roi_df;
    CL_E_roi = CL_E .* roi_df; % Use the same ROI for CL

    % Metrics within the selected ROI for non-zero values
    nonzero_values_df = DF_E_roi(DF_E_roi ~= 0);
    nonzero_values_cl = CL_E_roi(CL_E_roi ~= 0);
    
    if isempty(nonzero_values_df)
        fprintf('No valid data in ROI for %s\n', case_name);
        continue;
    end

    % DF Statistics
    HCL_roi = max(nonzero_values_df);
    avg_roi = mean(nonzero_values_df);
    mod_roi = mode(nonzero_values_df);
    std_roi = std(nonzero_values_df);
    var_roi = var(nonzero_values_df);
    number_points_df = numel(nonzero_values_df);

    % CL Statistics
    HCL_ms = max(nonzero_values_cl);
    avg_CL_ms = mean(nonzero_values_cl);
    mod_CL_ms = mode(nonzero_values_cl);
    std_CL_ms = std(nonzero_values_cl);
    var_CL_ms = var(nonzero_values_cl);

    % Store results
    Results_all.(case_name).DF_stats = [HCL_roi, avg_roi, mod_roi, std_roi, var_roi, number_points_df];
    Results_all.(case_name).CL_stats = [HCL_ms, avg_CL_ms, mod_CL_ms, std_CL_ms, var_CL_ms];
    Results_all.(case_name).DF_data = nonzero_values_df; % Store actual data for histograms
    Results_all.(case_name).CL_data = nonzero_values_cl; % Store actual data for histograms
    
    disp(['=== ' case_name ' - Dominant Frequency (DF) Statistics ===']);
    disp(['Higher DF in ROI: ', num2str(HCL_roi), ' Hz']);
    disp(['Average DF in ROI: ', num2str(avg_roi), ' Hz']);
    disp(['Mode DF in ROI: ', num2str(mod_roi), ' Hz']);
    disp(['Standard Deviation of DF in ROI: ', num2str(std_roi), ' Hz']);
    disp(['Variance of DF in ROI: ', num2str(var_roi)]);
    disp(['Number of electrodes in the ROI: ', num2str(number_points_df)]);

    disp(['=== ' case_name ' - Cycle Length (CL) Statistics ===']);
    disp(['Highest CL in ROI: ', num2str(HCL_ms), ' ms']);
    disp(['Average CL in ROI: ', num2str(avg_CL_ms), ' ms']);
    disp(['Mode CL in ROI: ', num2str(mod_CL_ms), ' ms']);
    disp(['Standard Deviation of CL in ROI: ', num2str(std_CL_ms), ' ms']);
end

%% Electrical Organization Index - OI
% Calculate OI for each case
OI_values = struct();
dfh_threshold_area = 0.6;
f_mode = 2;
debug = 1;

for i = 1:length(cases)
    case_name = cases{i};
    MFFTi = DF_values.(case_name).MFFTi;
    Sffti = DF_values.(case_name).Sffti;
    fstep = DF_values.fstep;
    
    fprintf('Calculating OI for %s...\n', case_name);
    
    OI = calculate_OI(MFFTi, Sffti, fstep, freq_down, freq_up, dfh_threshold_area, f_mode, debug);
    OI_values.(case_name) = OI;
    
    % OI Statistics
    figure(); 
    imshow(OI, 'InitialMagnification', 'fit');
    title(['Select ROI for OI Analysis - ' case_name]); 
    roi_oi = roipoly;

    % Metrics within the selected ROI for non-zero values
    nonzero_values_oi = OI(roi_oi);
    nonzero_values_oi = nonzero_values_oi(nonzero_values_oi ~= 0);
    
    if ~isempty(nonzero_values_oi)
        HOI_roi = max(nonzero_values_oi);
        avg_oi = mean(nonzero_values_oi);
        mod_oi = mode(nonzero_values_oi);
        std_oi = std(nonzero_values_oi);
        var_oi = var(nonzero_values_oi);
        number_points_oi = numel(nonzero_values_oi);
        
        Results_all.(case_name).OI_stats = [HOI_roi, avg_oi, mod_oi, std_oi, var_oi, number_points_oi];
        Results_all.(case_name).OI_data = nonzero_values_oi; % Store actual data for histograms
        
        disp(['=== ' case_name ' - Organization Index (OI) Statistics ===']);
        disp(['Higher OI in ROI: ', num2str(HOI_roi)]);
        disp(['Average OI in ROI: ', num2str(avg_oi)]);
        disp(['Mode OI in ROI: ', num2str(mod_oi)]);
        disp(['Standard Deviation of OI in ROI: ', num2str(std_oi)]);
        disp(['Variance of OI in ROI: ', num2str(var_oi)]);
        disp(['Number of electrodes in the ROI: ', num2str(number_points_oi)]);
    end
end


%% Create Results Table
% Compile all results into tables
AllResultsTables = struct();

for i = 1:length(cases)
    case_name = cases{i};
    
    if isfield(Results_all, case_name) && isfield(Results_all.(case_name), 'DF_stats')
        Parameter = {'Dominant Frequency (Hz)'; 'Cycle Length (ms)'; 'Organization Index'};
        
        DF_stats = Results_all.(case_name).DF_stats;
        CL_stats = Results_all.(case_name).CL_stats;
        
        if isfield(Results_all.(case_name), 'OI_stats')
            OI_stats = Results_all.(case_name).OI_stats;
        else
            OI_stats = [NaN, NaN, NaN, NaN, NaN, 0];
        end
        
        Max_Value = [DF_stats(1); CL_stats(1); OI_stats(1)];
        Mean_Value = [DF_stats(2); CL_stats(2); OI_stats(2)];
        Mode_Value = [DF_stats(3); CL_stats(3); OI_stats(3)];
        Std_Deviation = [DF_stats(4); CL_stats(4); OI_stats(4)];
        Variance = [DF_stats(5); CL_stats(5); OI_stats(5)];
        Number_Points = [DF_stats(6); DF_stats(6); OI_stats(6)];
        
        ResultsTable = table(Parameter, Max_Value, Mean_Value, Mode_Value, Std_Deviation, Variance, Number_Points);
        AllResultsTables.(case_name) = ResultsTable;
        
        fprintf('\n=== %s - RESULTS TABLE ===\n', case_name);
        disp(ResultsTable);
    end
end


%% Summary Visualizations - Electrical
% Create summary figures for each case
for i = 1:length(cases)
    case_name = cases{i};
    
    if isfield(DF_values, case_name) && isfield(OI_values, case_name) && isfield(CL_values, case_name)
        DF_E = DF_values.(case_name).MFFTi;
        OI = OI_values.(case_name);
        CL_E = CL_values.(case_name); % Use stored CL values
        
        figure('Position', [100 100 1200 800]);
        
        % Subplot 1: DF Map
        subplot(2,3,1);
        if strcmp(case_name, 'TANK')
            imagesc(DF_E);
            pbaspect([2 1 1]);
        else
            imagesc(imrotate(DF_E, 90));
            axis equal;
        end
        colormap(jet);
        colorbar;
        title([case_name ' - Dominant Frequency (Hz)']);
        axis off;
        
        % Subplot 2: CL Map
        subplot(2,3,2);
        if strcmp(case_name, 'TANK')
            imagesc(CL_E);
            pbaspect([2 1 1]);
        else
            imagesc(imrotate(CL_E, 90));
            axis equal;
        end
        colormap(jet);
        colorbar;
        title([case_name ' - Cycle Length (ms)']);
        axis off;
        
        % Subplot 3: OI Map
        subplot(2,3,3);
        if strcmp(case_name, 'TANK')
            imagesc(OI);
            pbaspect([2 1 1]);
        else
            imagesc(imrotate(OI, 90));
            axis equal;
        end
        colormap(jet);
        colorbar;
        title([case_name ' - Organization Index']);
        axis off;
        
        % Add histograms if ROI data exists
        if isfield(Results_all, case_name)
            % Subplot 4: DF Histogram
            subplot(2,3,4);
            if isfield(Results_all.(case_name), 'DF_data')
                histogram(Results_all.(case_name).DF_data, 20, 'FaceColor', 'blue', 'FaceAlpha', 0.7);
                xlabel('Dominant Frequency (Hz)');
                ylabel('Count');
                title('DF Distribution');
                grid on;
            end
            
            % Subplot 5: CL Histogram
            subplot(2,3,5);
            if isfield(Results_all.(case_name), 'CL_data')
                histogram(Results_all.(case_name).CL_data, 20, 'FaceColor', 'red', 'FaceAlpha', 0.7);
                xlabel('Cycle Length (ms)');
                ylabel('Count');
                title('CL Distribution');
                grid on;
            end
            
            % Subplot 6: OI Histogram
            subplot(2,3,6);
            if isfield(Results_all.(case_name), 'OI_data')
                histogram(Results_all.(case_name).OI_data, 20, 'FaceColor', 'green', 'FaceAlpha', 0.7);
                xlabel('Organization Index');
                ylabel('Count');
                title('OI Distribution');
                grid on;
            end
        end
        
        sgtitle([case_name ' - Electrical Analysis Summary'], 'FontSize', 16);
    end
end


%% Save workspace variables for future use
save('E_DF_CL_OI.mat', 'freq_up', 'freq_down', 'in_sample', 'end_sample', ...
    'DF_values', 'CL_values', 'OI_values', 'Results_all', 'AllResultsTables', '-v7.3');

disp(' ');
disp('=== ELECTRICAL ANALYSIS COMPLETE ===');
disp('All results have been saved to:');
disp('- E_DF_CL_OI_Analysis.mat (MATLAB workspace)');
disp('- Electrical_Analysis_Results_*.xlsx (Excel tables)');