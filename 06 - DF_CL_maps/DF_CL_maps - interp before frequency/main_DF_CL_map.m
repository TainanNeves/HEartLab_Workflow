
%% DOMINANT FREQUENCY
clear; clc;


%% Loading variables

load("C:\Users\HEartLab\Documents\GitHub\HEartLab\00 - examples\data_filtered_sync_E14_F3_R4.mat"); %Filtered data

signal_file = load('C:\Users\HeartLAB\Documents\Documents\Current codes\filtered_recordings\data_filtered_sync_E14_F4_R10.mat');


%% Optic Dominant Frequency Analysis

Data = D_SYNC.CAM3;
Fsampling = 4000;

% Set upper and lower frequency bounds for the Dominant Frequency calculation
UP = 10;
DOWN = 0.5;

% Extract a subset of the optical data for analysis (adjust the sample range accordingly)
in_sample = 8366;
end_sample = 25464;
Data_temp = Data(:,:,in_sample:end_sample);

% Perform Dominant Frequency analysis on the subset
[DF_O, Sfft_O, fstep] = f_DF_optico(Data_temp, Fsampling, UP, DOWN);

% Plot the spectrum of specific pixel locations
% Select point
Background = squeeze(Data(:,:,2000));
pick_up_a_trace(Background, Data,1);    % Select a pixel in the image and shows the optical signal
                                        %Press space to stop
% Plot
pa = [51 142]; pv = [60 33];
b = length(Sfft_O(1, 1, :));
figure;
plot((1:b) * fstep, squeeze(Sfft_O(pa(1), pa(2), :)), 'LineWidth', 2);  % Atrium
hold on;
plot((1:b) * fstep, squeeze(Sfft_O(pv(1), pv(2), :)), 'LineWidth', 2);   % Ventricle
legend('Atrium', 'Ventricle', 'FontSize', 12);
xlabel('Frequency [Hz]');
title('Spectrum');
set(gca, 'fontsize', 14);
xlim([DOWN UP]);


% Susbstitude values (If needed)
find_value = 0.5; % Value to replace
tolerancia = 1e-10;
indices = find(abs(DF_O - find_value) < tolerancia);
DF_O(indices) = 0; % Value to include


% Dominant Frequency map
C = jet(256);
C(1,1:3) = [1 1 1]; % White for background
C(256,1:3) = [0.5, 0.5, 0.5]; % Gray for the max value
d = 200;  % Define the color for the background
color = 'black';  % Specify color for text
tam = 11;  % Set font size
f1 = figure('color', 'white', 'Position', [50 50 500 500]);
J = DF_O;
% Set zeros inside ROI to 0.5
J(D_OP.ROI.ROI_2 & J == 0) = UP; %Put Max value in zeros inside the ROI
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
caxis([DOWN UP]);


%% Optical Dominant Frequency Analysis - Statistics
% Display the image and let the user define the ROI
figure();
imshow(DF_O);
title('Select ROI');
roi = roipoly;

% Apply the ROI to your matrix (assuming DF_O is your matrix)
DF_O_roi = DF_O .* roi;

% Metrics within the selected ROI for non-zero values
nonzero_values = DF_O_roi(DF_O_roi ~= 0);
HDF_roi = max(nonzero_values);
disp(['Higher DF in ROI: ', num2str(HDF_roi)]);
avg_roi = mean(nonzero_values);
disp(['Average DF in ROI: ', num2str(avg_roi)]);
mod_roi = mode(nonzero_values);
disp(['Mode DF in ROI: ', num2str(mod_roi)]);
std_roi = std(nonzero_values);
disp(['Standard Deviation of DF in ROI: ', num2str(std_roi)]);
var_roi = var(nonzero_values);
disp(['Variance of DF in ROI: ', num2str(var_roi)]);


%% Electric Dominant Frequency analysis

Data = D_SYNC.EL;
Fsampling = 4000;

% Dominant Frequency calculation
freq_up = 4;
freq_down = 0.5;
in_sample = 2*4000;
end_sample = 4*4000;
% Calc
Data_temp = Data(:, in_sample:end_sample);

% Initialize the struct to store results
df_values = struct();

for i = 1:4
    V_interpolated = electric_interp(Data_temp, i);  % Interpolate data
    [MFFTi,Sffti,fstep] = f_DF_electric(V_interpolated, Fsampling, freq_up, freq_down);  % Compute frequency map
    plot_electric_DF(MFFTi, [freq_down freq_up], i);  % Plot frequency map
   
    % storing the adjusted matrices
    if i < 4
        electrodes = sprintf('MEA%d', i);
        df_values.(electrodes) = fillMatrixMEA(MFFTi);
    else
        electrodes = sprintf('TANK');
        df_values.(electrodes) = fillMatrixTANK(MFFTi);
    end
end

% Ploting frequency spectrum
el_to_plot = [149 162 164 183];
figure;
for i = 1:length(el_to_plot)
    plot(fstep:fstep:freq_up, Sffti(el_to_plot(i),:), 'DisplayName', ['Electrode ' num2str(el_to_plot(i))]);
    hold on
end
xlabel('Frequency (Hz)');
title(['Frequency spectrum']);
legend('show');


%% Electric Dominant Frequency analysis - ECGi

% Heart geometry
heart_geo_file = 'D:\3D Heart Geometry\heart_geometry_exp18.mat';
heart_data = load(heart_geo_file);
heart_geo = heart_data.(subsref(fieldnames(heart_data),substruct('{}',{1})));

% Estimated signal
estimated_signal = load('D:\Figura Frontiers\EXP18_data\estimated_signal_E18F2R2_sync_2s_4s');
estimated_signal = estimated_signal.(subsref(fieldnames(estimated_signal),substruct('{}',{1})));

Data = estimated_signal;
Fsampling = 4000;

% Parameters setting
freq_up = 4;
freq_down = 0.5;
start_time = 1;
end_time = 4;
in_sample = 1;
end_sample = 8000;

% Dominant Frequency calculation
Data_temp = Data(:, in_sample:end_sample);
[MFFTi,Sffti,fstep] = f_DF_electric(Data_temp, Fsampling, freq_up, freq_down);

% Organizing geometry
faces = heart_geo.faces;
x = heart_geo.vertices(:,1);
y = heart_geo.vertices(:,2);
z = heart_geo.vertices(:,3);

% Plotting
figure();
trisurf(faces,x,y,z,MFFTi,'facecolor','interp','LineStyle','none');
grid off; axis off;
%title(['Instante ', num2str(inst1/Fs), 's Frame', num2str(inst1)]);
colormap('jet');

% Adding colorbar with label
c = colorbar;
c.Label.String = 'Frequency (Hz)';
caxis([0 10]);


%% Electric Dominant Frequency analysis ECGi - Statistics
% Select Electrodes range
roi = [1:length(MFFTi)];

% Metrics within the selected ROI for non-zero values
nonzero_values = MFFTi(roi);
nonzero_values = nonzero_values(nonzero_values ~= 0);
HDF_roi = max(nonzero_values);
disp(['Higher DF in ROI: ', num2str(HDF_roi)]);
avg_roi = mean(nonzero_values);
disp(['Average DF in ROI: ', num2str(avg_roi)]);
mod_roi = mode(nonzero_values);
disp(['Mode DF in ROI: ', num2str(mod_roi)]);
std_roi = std(nonzero_values);
disp(['Standard Deviation of DF in ROI: ', num2str(std_roi)]);
var_roi = var(nonzero_values);
disp(['Variance of DF in ROI: ', num2str(var_roi)]);
%% Electric Dominant Frequency analysis - Statistics
% Select Electrodes range
roi = [129:174, 177:190];

% Metrics within the selected ROI for non-zero values
nonzero_values = MFFTi(roi);
nonzero_values = nonzero_values(nonzero_values ~= 0);
HDF_roi = max(nonzero_values);
disp(['Higher DF in ROI: ', num2str(HDF_roi)]);
avg_roi = mean(nonzero_values);
disp(['Average DF in ROI: ', num2str(avg_roi)]);
mod_roi = mode(nonzero_values);
disp(['Mode DF in ROI: ', num2str(mod_roi)]);
std_roi = std(nonzero_values);
disp(['Standard Deviation of DF in ROI: ', num2str(std_roi)]);
var_roi = var(nonzero_values);
disp(['Variance of DF in ROI: ', num2str(var_roi)]);


%% Electric Cicle Lenth - CL

% Loadind data
Data = D_SYNC.EL;
Fsampling = 4000;
in_sample = 9386;
end_sample = 12015;

% CL - Calculation
Data_temp = Data(:, in_sample:end_sample);
for i=1:size(Data_temp,1)
    s=Data_temp(i,:);
    [y,x]=findpeaks(s,'MinPeakHeight',max(s)*0.7,'MinPeakDistance',100); % Defaut 0.005 and 1000
        % Ploting in time with peaks (Checking if its good to select
        % automaticly = Its not)
        figure();
        plot(s);
        hold on;
        scatter(x, y, 'r','filled');
        hold off;
        title('Electrode ', num2str(i));
        xlabel('Time (samples)'); ylabel('Signal'); legend('Signal', 'Detected Peaks');
    if ~isempty(x)
        cl=x(2)-x(1);
        CL_E(i)=cl/4000*1000; % Transforming to ms
    end
end

% Ploting maps
plot_electric_CL(CL_E, [floor(min(CL_E)) ceil(max(CL_E))], 1); % MEA 1
plot_electric_CL(CL_E, [floor(min(CL_E)) ceil(max(CL_E))], 2); % MEA 2
plot_electric_CL(CL_E, [floor(min(CL_E)) ceil(max(CL_E))], 3); % MEA 3
plot_electric_CL(CL_E, [floor(min(CL_E)) ceil(max(CL_E))], 4); % TANK


%% Electric Cycle Length analysis - Statistics
% Select Electrodes range
roi = [129:174, 177:190];

% Metrics within the selected ROI for non-zero values
nonzero_values = CL_E(roi);
nonzero_values = nonzero_values(nonzero_values ~= 0);
HDF_roi = max(nonzero_values);
disp(['Higher CL in ROI: ', num2str(HDF_roi)]);
avg_roi = mean(nonzero_values);
disp(['Average CL in ROI: ', num2str(avg_roi)]);
mod_roi = mode(nonzero_values);
disp(['Mode CL in ROI: ', num2str(mod_roi)]);
std_roi = std(nonzero_values);
disp(['Standard Deviation of CL in ROI: ', num2str(std_roi)]);
var_roi = var(nonzero_values);
disp(['Variance of CL in ROI: ', num2str(var_roi)]);


%% Optic Cicle Lenth - CL

% Loadind data
Data = D_SYNC.CAM3;
Fsampling = 4000;
in_sample = 9386;
end_sample = 12015;

% CL - Calculation
Data_temp = Data(:,:,in_sample:end_sample);
[num_rows, num_cols, num_frames] = size(Data_temp);
CL_O = zeros(num_rows, num_cols); % Initialize CL matrix
for i = 1:num_rows
    for j = 1:num_cols
        s = squeeze(Data_temp(i, j, :)); % Extract time series for each pixel
        [y,x]=findpeaks(s,'MinPeakHeight',max(s)*0.5,'MinPeakDistance',100); % Defaut 0.005 and 1000
        if ~isempty(x) && length(x) >= 2 % Ensure there are at least two peaks
            cl = x(2) - x(1);
            CL_O(i, j) = cl/4000*1000; % Transforming to ms
        end
    end
end

% Susbstitude values (If needed)
find_value = 0.5; % Value to replace
tolerancia = 1e-10;
indices = find(abs(CL_O - find_value) < tolerancia);
CL_O(indices) = 0; % Value to include

% Cycle Length map
C = jet(256);
C(1,1:3) = [1 1 1]; % White for background
C(256,1:3) = [0.5, 0.5, 0.5]; % Gray for the max value
d = 200;  % Define the color for the background
color = 'black';  % Specify color for text
tam = 11;  % Set font size
f1 = figure('color', 'white', 'Position', [50 50 500 500]);
J = CL_O;
% Set zeros inside ROI to max
% J(D_OP.ROI.ROI_2 & J == 0) = ceil(max(max(CL_O)))+10; %Put Max value in zeros inside the ROI
% J = imrotate(J, 90);
% Plot
imagesc(J);
colormap(C);
box off;
set(gca, 'fontsize', 18);
ylabel('Pixels');
xlabel('Pixels');
axis off;
hBar1 = colorbar('eastoutside');
ylabel(hBar1, 'Cycle Length [ms]', 'FontSize', 14);
caxis([floor(min(CL_O(CL_O ~= 0))) ceil(max(max(CL_O)))+10]);


%% Optical Cycle Length Analysis - Statistics
% Display the image and let the user define the ROI
figure();
imshow(CL_O);
title('Select ROI');
roi = roipoly;

% Apply the ROI to your matrix (assuming DF_O is your matrix)
CL_O_roi = CL_O .* roi;

% Metrics within the selected ROI for non-zero values
nonzero_values = CL_O_roi(CL_O_roi ~= 0);
HDF_roi = max(nonzero_values);
disp(['Higher CL in ROI: ', num2str(HDF_roi)]);
avg_roi = mean(nonzero_values);
disp(['Average CL in ROI: ', num2str(avg_roi)]);
mod_roi = mode(nonzero_values);
disp(['Mode CL in ROI: ', num2str(mod_roi)]);
std_roi = std(nonzero_values);
disp(['Standard Deviation of CL in ROI: ', num2str(std_roi)]);
var_roi = var(nonzero_values);
disp(['Variance of CL in ROI: ', num2str(var_roi)]);



