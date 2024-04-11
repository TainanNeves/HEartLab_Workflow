%% DOMINANT FREQUENCY
clear; clc;


%% Loading variables

load("C:\Users\HEartLab\Documents\GitHub\HEartLab\00 - examples\data_filtered_sync_E14_F3_R4.mat"); %Filtered data


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
freq_up = 10;
freq_down = 0.5;
in_sample = 8366;
end_sample = 25464;
% Calc
Data_temp = Data(:, in_sample:end_sample);
[MFFTi,Sffti,fstep] = f_DF_electric(Data_temp, Fsampling, freq_up, freq_down);

% Ploting frequency spectrum
el_to_plot = [10 75 23 132];
figure;
for i = 1:length(el_to_plot)
    plot(fstep:fstep:freq_up, Sffti(el_to_plot(i),:), 'DisplayName', ['Electrode ' num2str(el_to_plot(i))]);
    hold on
end
xlabel('Frequency (Hz)');
title(['Frequency spectrum']);
legend('show');

% Ploting maps
plot_electric_DF(MFFTi, [freq_down freq_up], 1); % MEA 1
plot_electric_DF(MFFTi, [freq_down freq_up], 2); % MEA 2
plot_electric_DF(MFFTi, [freq_down freq_up], 3); % MEA 3
plot_electric_DF(MFFTi, [freq_down freq_up], 4); % TANK


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


%% Cicle Lenth - CL

% Loadind data
Data = D_SYNC.EL;
Fsampling = 4000;
in_sample = 8366;
end_sample = 25464;

% CL - Calculation
Data_temp = Data(:, in_sample:end_sample);
for i=1:size(Data_temp,1)
    s=Data_temp(i,:);
    [y,x]=findpeaks(s,'MinPeakHeight',0.005,'MinPeakDistance',1000);
    if ~isempty(x)
        cl=x(2)-x(1);
        CL(i)=cl/4000*1000;
    end
end

