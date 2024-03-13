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


% Display Dominant Frequency map
C = jet(256);
C(1,1:3) = [1 1 1];
d = 200;  % Define the color for the background
color = 'black';  % Specify color for text
tam = 11;  % Set font size
f1 = figure('color', 'white', 'Position', [50 50 500 500]);
J = DF_O;
J = imrotate(J, 90);
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

% Metrics
HDF = max(max(DF_O));
disp(['Higher DF: ',num2str(HDF)]);
avg = mean(mean(DF_O));
disp(['Average DF: ',num2str(avg)]);
mod = mode(mode(DF_O));
disp(['Mode DF: ',num2str(mod)]);



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

% Metrics
HDF = max(MFFTi);
disp(['Higher DF: ',num2str(HDF)]);
avg = mean(MFFTi);
disp(['Average DF: ',num2str(avg)]);
mod = mode(MFFTi);
disp(['Mode DF: ',num2str(mod)]);

% Ploting maps
plot_electric_DF(MFFTi, [freq_down freq_up], 1); % MEA 1
plot_electric_DF(MFFTi, [freq_down freq_up], 2); % MEA 2
plot_electric_DF(MFFTi, [freq_down freq_up], 3); % MEA 3
plot_electric_DF(MFFTi, [freq_down freq_up], 4); % TANK

