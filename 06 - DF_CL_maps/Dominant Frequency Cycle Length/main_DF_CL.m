%% DOMINANT FREQUENCY AND CYCLE LENGTH

%% Reseting Setup
clear; clc;


%% Load Data
load("E:\HEartLab\experiment_analyses\exp20_analysis\03 - synchronization_optical_electric\data_filtered_sync_E20_F01_R01.mat"); % Filtered syncrhonised data
load("E:\HEartLab\experiment_analyses\exp20_analysis\04 - Interpolate signals Laplacian\InterpolatedSignalsE20_F01_R01_filtered.mat"); % Filteres interpolated data


%% Optic Dominant Frequency Analysis

Data = D_SYNC.CAM3;
ROI = D_SYNC;
Fsampling = 4000;

% Set upper and lower frequency bounds for the Dominant Frequency calculation
freq_up = 10;
freq_down = 0.5;

% Extract a subset of the optical data for analysis (adjust the sample range accordingly)
in_sample = 8366;
end_sample = 25464;
Data_temp = Data(:,:,in_sample:end_sample);

% Perform Dominant Frequency analysis on the subset
[DF_O, Sfft_O, fstep] = f_DF_optico(Data_temp, Fsampling, freq_up, freq_down);

% Plot the spectrum of specific pixel locations
% Select point
Background = squeeze(Data(:,:,2000));
[px, py] = pick_up_a_trace(Background, Data,1);    % Select a pixel in the image and shows the optical signal
                                        %Press space to stop
% Plot
pa = [px(1) py(1)]; pv = [px(2) py(2)];
b = length(Sfft_O(1, 1, :));
figure;
plot((1:b) * fstep, squeeze(Sfft_O(pa(1), pa(2), :)), 'LineWidth', 2);  % Atrium
hold on;
plot((1:b) * fstep, squeeze(Sfft_O(pv(1), pv(2), :)), 'LineWidth', 2);   % Ventricle
legend('Atrium', 'Ventricle', 'FontSize', 12);
xlabel('Frequency [Hz]');
title('Spectrum');
set(gca, 'fontsize', 14);
xlim([freq_down freq_up]);


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
caxis([freq_down freq_up]);


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
HCL_roi = max(nonzero_values);
disp(['Higher DF in ROI: ', num2str(HCL_roi)]);
avg_roi = mean(nonzero_values);
disp(['Average DF in ROI: ', num2str(avg_roi)]);
mod_roi = mode(nonzero_values);
disp(['Mode DF in ROI: ', num2str(mod_roi)]);
std_roi = std(nonzero_values);
disp(['Standard Deviation of DF in ROI: ', num2str(std_roi)]);
var_roi = var(nonzero_values);
disp(['Variance of DF in ROI: ', num2str(var_roi)]);


%% Optic Organization Index - OI

% Variables
MFFTi = DF_O;
Sffti = Sfft_O;
fstep = fstep;
Hzi = freq_down;
Hzf = freq_up;
dfh_threshold_area = 0.8;
f_mode = 1;
debug = 1;

% Calculating
OI = calculate_OI(MFFTi, Sffti, fstep, Hzi, Hzf, dfh_threshold_area, f_mode, debug);

% Statistics
% Select Electrodes range
figure(); imshow(OI, 'InitialMagnification', 'fit');
title('Select ROI'); roi = roipoly;
% Metrics within the selected ROI for non-zero values
nonzero_values = OI(roi);
nonzero_values = nonzero_values(nonzero_values ~= 0);
HOI_roi = max(nonzero_values);
disp(['Higher OI in ROI: ', num2str(HOI_roi)]);
avg_roi = mean(nonzero_values);
disp(['Average OI in ROI: ', num2str(avg_roi)]);
mod_roi = mode(nonzero_values);
disp(['Mode OI in ROI: ', num2str(mod_roi)]);
std_roi = std(nonzero_values);
disp(['Standard Deviation of OI in ROI: ', num2str(std_roi)]);
var_roi = var(nonzero_values);
disp(['Variance of OI in ROI: ', num2str(var_roi)]);


%% Electric Dominant Frequency Analysis

% Preview Signals
Data = InterpSignal.Sync.TANK;
Background = squeeze(Data(:,:,2000));
pick_up_a_trace(Background, Data,1);


% Parameters
Fsampling = 4000;
freq_up = 20;
freq_down = 0.5;
in_sample = 3*4000;
end_sample = 4*4000;


% Loading Data
Data1 = InterpSignal.Sync.MEA1;
Data2 = InterpSignal.Sync.MEA2;
Data3 = InterpSignal.Sync.MEA3;
Data4 = InterpSignal.Sync.TANK;


% Loading Temp Data into a cell array
Data_temp{1} = Data1(:,:,in_sample:end_sample);
Data_temp{2} = Data2(:,:,in_sample:end_sample);
Data_temp{3} = Data3(:,:,in_sample:end_sample);
Data_temp{4} = Data4(:,:,in_sample:end_sample);


% Calculating
DF_values = struct();
cases = {'MEA1', 'MEA2', 'MEA3', 'TANK'};
for i = 1:4
    Data = Data_temp{i};
    case_name = cases{i};

    [MFFTi, Sffti, fstep] = f_DF_electric(Data, Fsampling, freq_up, freq_down);
    
    DF_values.(case_name).MFFTi = MFFTi;
    DF_values.(case_name).Sffti = Sffti;
    DF_values.fstep = fstep;
end


% Ploting Frequency Spectrum (electrode)
el_to_plot = [10 8 74 71 23 27 131 160];
fstep = DF_values.fstep;
figure;
for i = 1:length(el_to_plot)
    [x,y,source] = getElectrodePosition(el_to_plot(i));
    case_name = cases{source};
    plot(fstep:fstep:freq_up, squeeze( DF_values.(case_name).Sffti(x,y,:) ), ...
                    'DisplayName', ['Electrode ' num2str(el_to_plot(i))]);
    hold on
end
xlabel('Frequency (Hz)');
title(['Frequency spectrum']);
legend('show');


% Ploting Frequency Spectrum (Coordinates)
source = DF_values.TANK;
points = [13 10; 13 11; 13 13; 13 14];
fstep = DF_values.fstep;
figure;
for i = 1:length(points)
    plot(fstep:fstep:freq_up, squeeze( source.Sffti(points(i,1),points(i,2),:) ), ...
                    'DisplayName', ['Point: X = ' num2str(points(i,1))...
                                                ' Y = ' num2str(points(i,2))]);
    hold on
end
xlabel('Frequency (Hz)');
title(['Frequency spectrum']);
legend('show');


% Dominant Frequencies Plot
cLimits = [0 10]; % Example colorbar limits
% Ploting
plot_MFFTi(DF_values.MEA1.MFFTi, 'Dominant Frequencies (MEA1)', cLimits, 'MEA');
plot_MFFTi(DF_values.MEA2.MFFTi, 'Dominant Frequencies (MEA2)', cLimits, 'MEA');
plot_MFFTi(DF_values.MEA3.MFFTi, 'Dominant Frequencies (MEA3)', cLimits, 'MEA');
plot_MFFTi(DF_values.TANK.MFFTi, 'Dominant Frequencies (TANK)', cLimits, 'TANK');


%% Electric Dominant Frequency analysis - Statistics
% Select MEA
DF_E = DF_values.TANK.MFFTi;

% Select Electrodes range
figure();
imshow(DF_values.MEA1.MFFTi, 'InitialMagnification', 'fit');
title('Select ROI');
roi = roipoly;

% Metrics within the selected ROI for non-zero values
nonzero_values = DF_E(roi);
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


%% Electric Organization Index - OI

% Variables
MFFTi = DF_values.MEA2.MFFTi;
Sffti = DF_values.MEA2.Sffti;
fstep = DF_values.fstep;
Hzi = freq_down;
Hzf = freq_up;
dfh_threshold_area = 0.8;
f_mode = 2;
debug = 0;

% Calculating
OI = calculate_OI(MFFTi, Sffti, fstep, Hzi, Hzf, dfh_threshold_area, f_mode, debug);

% Statistics
% Select Electrodes range
figure(); imshow(OI, 'InitialMagnification', 'fit');
title('Select ROI'); roi = roipoly;
% Metrics within the selected ROI for non-zero values
nonzero_values = OI(roi);
nonzero_values = nonzero_values(nonzero_values ~= 0);
HOI_roi = max(nonzero_values);
disp(['Higher OI in ROI: ', num2str(HOI_roi)]);
avg_roi = mean(nonzero_values);
disp(['Average OI in ROI: ', num2str(avg_roi)]);
mod_roi = mode(nonzero_values);
disp(['Mode OI in ROI: ', num2str(mod_roi)]);
std_roi = std(nonzero_values);
disp(['Standard Deviation of OI in ROI: ', num2str(std_roi)]);
var_roi = var(nonzero_values);
disp(['Variance of OI in ROI: ', num2str(var_roi)]);


%% Electric Dominant Frequency analysis - ECGi

% Heart geometry
heart_geo_file = 'C:\Users\HeartLAB\Documents\Documents\CinC 2024\ECGi\Dados\heart_geometry_exp14.mat';
heart_data = load(heart_geo_file);
heart_geo = heart_data.(subsref(fieldnames(heart_data),substruct('{}',{1})));

% Estimated signal
estimated_signal = load('C:\Users\HeartLAB\Documents\Documents\CBEB 2024\estimated_signal_E14F4R10_sync_8s.mat');
estimated_signal = estimated_signal.(subsref(fieldnames(estimated_signal),substruct('{}',{1})));

Data = estimated_signal;
Fsampling = 4000;

% Parameters setting
freq_up = 10;
freq_down = 0.5;
start_time = 1;
end_time = 2;
in_sample = start_time * Fsampling;
end_sample = end_time * Fsampling;

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


%% Optic Cicle Lenth - CL

% Preview Signals
Data = D_SYNC.CAM1;
Background = squeeze(Data(:,:,2000));
pick_up_a_trace(Background, Data,1);


% Parameters
Fsampling = 4000;
in_sample = 3*4000;
end_sample = 3.5*4000;
% Your time interval must have two peaks
Data_temp = Data(:,:,in_sample:end_sample);


% Calculating
CL_O = calculate_CL_O(Data_temp, Fsampling, 1);


% Cycle Length Plot
cLimits = [0 300]; % Example colorbar limits
% Ploting
plot_CL_O(CL_O, 'Cycle Length (CAM X)', cLimits);


% Manually correcting Wrong Values
data_modf = CL_O;
figure(); imshow(data_modf, [], 'InitialMagnification', 'fit');
title('Select ROI'); roi = roipoly;
New_Value = 300;
data_modf(roi) = New_Value;
% Preview the corrected data
plot_CL_O(data_modf, 'Cycle Length (Corrected)', cLimits);
% Update the original data with corrected values
CL_O = data_modf;


%% Optical Cycle Length Analysis - Statistics
% Display the image and let the user define the ROI
figure(); imshow(CL_O, [], 'InitialMagnification', 'fit');
title('Select ROI'); roi = roipoly;
% Apply the ROI to your matrix (assuming DF_O is your matrix)
CL_O_roi = CL_O .* roi;

% Metrics within the selected ROI for non-zero values
nonzero_values = CL_O_roi(CL_O_roi ~= 0);
HCL_roi = max(nonzero_values);
disp(['Higher CL in ROI: ', num2str(HCL_roi)]);
avg_roi = mean(nonzero_values);
disp(['Average CL in ROI: ', num2str(avg_roi)]);
mod_roi = mode(nonzero_values);
disp(['Mode CL in ROI: ', num2str(mod_roi)]);
std_roi = std(nonzero_values);
disp(['Standard Deviation of CL in ROI: ', num2str(std_roi)]);
var_roi = var(nonzero_values);
disp(['Variance of CL in ROI: ', num2str(var_roi)]);


%% Electric Cicle Length - CL

% Preview Signals
Data = InterpSignal.Sync.TANK;
Background = squeeze(Data(:,:,2000));
pick_up_a_trace(Background, Data,1);


% Parameters
Fsampling = 4000;
in_sample = 3*4000;
end_sample = 4*4000;
% Your time interval must have two peaks


% Loading Data
Data1 = InterpSignal.Sync.MEA1;
Data2 = InterpSignal.Sync.MEA2;
Data3 = InterpSignal.Sync.MEA3;
Data4 = InterpSignal.Sync.TANK;
% Loading Temp Data into a cell array
Data_temp{1} = Data1(:,:,in_sample:end_sample);
Data_temp{2} = Data2(:,:,in_sample:end_sample);
Data_temp{3} = Data3(:,:,in_sample:end_sample);
Data_temp{4} = Data4(:,:,in_sample:end_sample);


% Calculating
CL_E = calculate_CL_E(Data_temp, Fsampling, 0);


% Cycle Length Plot
cLimits = [0 300]; % Example colorbar limits
% Ploting
plot_CL_E(CL_E.MEA1, 'Cycle Length (MEA1)', cLimits, 'MEA');
plot_CL_E(CL_E.MEA2, 'Cycle Length (MEA2)', cLimits, 'MEA');
plot_CL_E(CL_E.MEA3, 'Cycle Length (MEA3)', cLimits, 'MEA');
plot_CL_E(CL_E.TANK, 'Cycle Length (TANK)', cLimits, 'TANK');


% Manually correcting Wrong Values
data_modf = CL_E.TANK;
figure(); imshow(data_modf, [],'InitialMagnification', 'fit');
title('Select ROI'); roi = roipoly;
New_Value = 300;
data_modf(roi) = New_Value;
% Preview
plot_CL_E(data_modf, 'Cycle Length', cLimits, 'TANK');
% Change in the origin
CL_E.TANK = data_modf;


%% Electric Cycle Length analysis - Statistics
% Select MEA
CL_E = CL_E.MEA1;

% Select Electrodes range
figure();
imshow(DF_values.MEA1.MFFTi, 'InitialMagnification', 'fit');
title('Select ROI');
roi = roipoly;

% Metrics within the selected ROI for non-zero values
nonzero_values = CL_E(roi);
nonzero_values = nonzero_values(nonzero_values ~= 0);
HCL_roi = max(nonzero_values);
disp(['Higher CL in ROI: ', num2str(HCL_roi)]);
avg_roi = mean(nonzero_values);
disp(['Average CL in ROI: ', num2str(avg_roi)]);
mod_roi = mode(nonzero_values);
disp(['Mode CL in ROI: ', num2str(mod_roi)]);
std_roi = std(nonzero_values);
disp(['Standard Deviation of CL in ROI: ', num2str(std_roi)]);
var_roi = var(nonzero_values);
disp(['Variance of CL in ROI: ', num2str(var_roi)]);

