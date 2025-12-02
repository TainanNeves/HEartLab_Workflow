%% Extraction and filtering of electrical data

clear; clc;


%% Reads the Open Ephys files

% fulfilename is the path to the structure.oebin file contained in one of the experiments
fullfilename = "E:\experiment_data\E30\Electrical\02\noFilter\2025-07-15_12-45-13\Record Node 108\experiment1\recording23\structure.oebin"; % Put the .oebin path

% Channels to save
channels = [1:192];

DATA = load_open_ephys_binary(fullfilename, 'continuous', 1); % Load the data file
DATA.Data = DATA.Data(channels, :); % Select desired channels
EVENTS = load_open_ephys_binary(fullfilename, 'events', 1); % Load the TTL file

rectime = length(DATA.Data(1, :))/DATA.Header.sample_rate; %s
disp(['RecTime: ', num2str(rectime), ' seconds']);



%% Conversion to uV + 60dB transformation
% Convertion to uV
bitVolts = 0.1949999928474426; % Conversion factor
DATA.Data = DATA.Data*bitVolts; % Converts the data to uV
% Convert to 60dB
DATA.Data = DATA.Data / 192;
DATA.Data = DATA.Data * 1000;


%% Convert TTL to appropiate format
% Takes only every other timestamp to record the timing of each rising edge
TTL = EVENTS.Timestamps(1:2:end);

if isempty(TTL)==false
    optictime(1) = (TTL(1) - DATA.Timestamps(1)); %s
    optictime(2) = optictime(1) + 10; %s
end

% Checks if the TTL variable is empty
if isempty(TTL)
    % Prints an alert message if TTL is empty
    disp('Warning, the TTL is empty!');
else
    % Prints a confirmation message if TTL is not empty
    disp('The TTL is not empty.');
end


%% Comparison Plot
Fs = 4000;
tank_electrodes = [129:174,177:190];
% Default time
if isempty(TTL) % Without TTL
    ti = 1;
    tf = rectime;
    tline1 = ti;
    tline2 = tf;
else
    ti = optictime(1);
    tf = optictime(2);
    tline1 = (TTL(1) - DATA.Timestamps(1)); % s
    tline2 = tline1 + 10; % s
end
time = linspace(ti, tf, (tf - ti) * Fs + 1);
% Electrodes selection
el1 = 22;  % RA
el2 = 6; % LA
el3 = 87; % V
el4 = 142; % Tank
Index = {'RA', 'LA', 'V','TANK'};
% Filter range (Butterworth)
f_low = 0.5;
f_high = 200;
% Filter Configuration (Wavelet)
waveletType = {'db4'}; % Wavelet type
numLevels = 10; % Number of decomposition levels
reconstructionLevelsSets = {[7, 8, 9, 10]}; % Levels to be used on reconstruction
% Create a figure
figure;
sgtitle('EXP XX - Folder XX - REC XX');
% Plotting using filter_signal function (Butterworth)
for i = 1:4
    subplot(5, 2, 2*i - 1); % Adjusted index for Butterworth
    filteredSignal = filter_signal_2(f_low, f_high, DATA.Data(eval(['el' num2str(i)]), :), Fs);
    plot(time, filteredSignal(:, ti * Fs:tf * Fs));
    title(['Electrode: ' num2str(eval(['el' num2str(i)])) ' (' Index{i} ') - Butterworth']);
    xlabel("Time (s)");
    ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex');
    xline(tline1, 'red');
    xline(tline2, 'red');
end
% Plotting using wavelet_filter function (Wavelet)
for i = 1:4
    subplot(5, 2, 2*i); % Adjusted index for Wavelet
    filteredSignal = wavelet_filter(DATA.Data(eval(['el' num2str(i)]), :), Fs, 1, waveletType, numLevels, reconstructionLevelsSets);
    plot(time, filteredSignal(:, ti * Fs:tf * Fs));
    title(['Electrode: ' num2str(eval(['el' num2str(i)])) ' (' Index{i} ') - Wavelet']);
    xlabel("Time (s)");
    ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex');
    xline(tline1, 'red');
    xline(tline2, 'red');
end
% avg subtraction of tank
tank_signals = DATA.Data(tank_electrodes, :);
mean_tank = mean(tank_signals);
% Tank with average subtraction plot - Butterworth
subplot(5, 2, 9);
sub = DATA.Data(el4, :) - mean_tank;
filteredSignal = filter_signal_2(f_low, f_high, sub, Fs);
plot(time, filteredSignal(:, ti * Fs:tf * Fs));
title(['Electrode: ' num2str(el4) ' - ' num2str(el3) ' (avg subtracted) - Butterworth']);
xlabel("Time (s)");
ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex');
xline(tline1, 'red');
xline(tline2, 'red');
% Tank with average subtraction plot - Wavelet
subplot(5, 2, 10);
sub = DATA.Data(el4, :) - mean_tank;
filteredSignal = wavelet_filter(sub, Fs, 1, waveletType, numLevels, reconstructionLevelsSets);
plot(time, filteredSignal(:, ti * Fs:tf * Fs));
title(['Electrode: ' num2str(el4) ' - ' num2str(el3) ' (avg subtracted) - Wavelet']);
xlabel("Time (s)");
ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex');
xline(tline1, 'red');
xline(tline2, 'red');
% Link axes for all subplots
h = findobj(gcf, 'type', 'axes');
linkaxes(h,'x');


% Spectrum Frequency Plots
% Create a new figure for frequency spectrum
figure;
sgtitle('Frequency Spectrum - EXP XX - Folder XX - REC XX');
% Calculate and plot frequency spectrum for Butterworth filtered signals
for i = 1:4
    subplot(5, 2, 2*i - 1);
    filteredSignal = filter_signal_2(f_low, f_high, DATA.Data(eval(['el' num2str(i)]), :), Fs);
    Y = fft(filteredSignal);
    L = length(filteredSignal);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;
    plot(f, P1);
    title(['Electrode: ' num2str(eval(['el' num2str(i)])) ' (' Index{i} ') - Butterworth Spectrum']);
    xlabel("Frequency (Hz)");
    ylabel('Magnitude');
    xlim([0 f_high+10]);
end
% Calculate and plot frequency spectrum for Wavelet filtered signals
for i = 1:4
    subplot(5, 2, 2*i);
    filteredSignal = wavelet_filter(DATA.Data(eval(['el' num2str(i)]), :), Fs, 1, waveletType, numLevels, reconstructionLevelsSets);
    Y = fft(filteredSignal);
    L = length(filteredSignal);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;
    plot(f, P1);
    title(['Electrode: ' num2str(eval(['el' num2str(i)])) ' (' Index{i} ') - Wavelet Spectrum']);
    xlabel("Frequency (Hz)");
    ylabel('Magnitude');
    xlim([0 f_high+10]);
end
% Calculate and plot frequency spectrum for avg subtracted tank signal (Butterworth)
subplot(5, 2, 9);
sub = DATA.Data(el4, :) - mean_tank;
filteredSignal = filter_signal_2(f_low, f_high, sub, Fs);
Y = fft(filteredSignal);
L = length(filteredSignal);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
plot(f, P1);
title(['Electrode: ' num2str(el4) ' (avg subtracted) - Butterworth Spectrum']);
xlabel("Frequency (Hz)");
ylabel('Magnitude');
xlim([0 f_high+10]);
% Calculate and plot frequency spectrum for avg subtracted tank signal (Wavelet)
subplot(5, 2, 10);
sub = DATA.Data(el4, :) - mean_tank;
filteredSignal = wavelet_filter(sub, Fs, 1, waveletType, numLevels, reconstructionLevelsSets);
Y = fft(filteredSignal);
L = length(filteredSignal);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
plot(f, P1);
title(['Electrode: ' num2str(el4) ' (avg subtracted) - Wavelet Spectrum']);
xlabel("Frequency (Hz)");
ylabel('Magnitude');
xlim([0 f_high+10]);
% Link axes for all subplots in the new figure
h = findobj(gcf, 'type', 'axes');
linkaxes(h,'x');


%% Save file to .mat
% Filename to save
FileName = 'E30_F02_R23';

% Raw file export
D_EL = struct(); % Initialize structure
D_EL.Data = DATA.Data; 
D_EL.Timestamps = DATA.Timestamps;
D_EL.Header = DATA.Header;
D_EL.Header.channels = channels;
D_EL.TTL = TTL;
D_EL.opticalin = TTL(1) - DATA.Timestamps(1);
save(['electric_data_', FileName, '_raw'], 'D_EL'); 
clear D_EL;

% Defining electrode arrays
el_butter = [1:32, 81:96, 129:174, 177:190]; % for Butterworth
el_wavelet = []; % for Wavelet
el_toZero = [33:80, 97:128, 175, 176, 191, 192]; % to set to zero
tank_electrodes = [129:174,177:190]; % Tank electrodes for average subtraction
avg_subtraction_enabled = 0; % Set to 1 to enable

% Filter range (Butterworth)
f_low_butter = 0.5;
f_high_butter = 200;
% Filter Configuration (Wavelet)
waveletType = {'db4'}; % Wavelet type
numLevels = 10; % Number of decomposition levels
reconstructionLevelsSets = {[7, 8, 9, 10]}; % Levels to be used on reconstruction

% Initialize filtered data variable
filtered_data = zeros(size(DATA.Data)); % Ensure the size matches the input data
% Perform average subtraction on tank electrodes if enabled
if avg_subtraction_enabled == 1
    mean_tank = mean(DATA.Data(tank_electrodes, :), 1);
    data_subtracted = DATA.Data - mean_tank;
else
    data_subtracted = DATA.Data;
end
% Filter electrodes with Butterworth
for i = 1:length(el_butter)
    filteredSignal = filter_signal_2(f_low_butter, f_high_butter, data_subtracted(el_butter(i), :), Fs);
    filtered_data(el_butter(i), :) = filteredSignal;
end
% Filter electrodes with Wavelet
for i = 1:length(el_wavelet)
    filteredSignal = wavelet_filter(data_subtracted(el_wavelet(i), :), Fs, 1, waveletType, numLevels, reconstructionLevelsSets);
    filtered_data(el_wavelet(i), :) = filteredSignal(1:size(filteredSignal,2)-1);
end
% Set specified electrodes to zero
for i = 1:length(el_toZero)
    filtered_data(el_toZero(i), :) = zeros(1, size(filtered_data, 2));
end

% Save filtered data
D_EL = struct(); % Initialize structure
D_EL.Data = filtered_data; 
D_EL.Timestamps = DATA.Timestamps;
D_EL.Header = DATA.Header;
D_EL.Header.channels = channels;
D_EL.TTL = TTL;
D_EL.opticalin = TTL(1) - DATA.Timestamps(1);
% Adding filtering information to the Header
D_EL.Header.FilterParameters = struct();
if ~isempty(el_butter) && all(el_butter > 0) && all(el_butter < 193)
    if all(el_butter > 0) && all(el_butter < 81)
        filter_butter = ["MEA", f_low_butter, f_high_butter];
        filter_wavelet = ["Tank", string(reconstructionLevelsSets{1})];
        D_EL.Header.FilterParameters.filter_butter = filter_butter;
        D_EL.Header.FilterParameters.filter_wavelet = filter_wavelet;
    else
        if all(el_butter > 80) && all(el_butter < 193)
            filter_butter = ["Tank", f_low_butter, f_high_butter];
            filter_wavelet = ["MEA", string(reconstructionLevelsSets{1})];
            D_EL.Header.FilterParameters.filter_butter = filter_butter;
            D_EL.Header.FilterParameters.filter_wavelet = filter_wavelet;
        else
            filter_butter = ["MEA", "Tank", f_low_butter, f_high_butter];
            D_EL.Header.FilterParameters.filter_butter = filter_butter;
        end
    end
else
    filter_wavelet = ["MEA", "Tank", string(reconstructionLevelsSets{1})];
    D_EL.Header.FilterParameters.filter_wavelet = filter_wavelet;
end

% Save the filtered data
save(['electric_data_', FileName, '_filtered'], 'D_EL', '-v7.3');

% Clear unwanted variables
clear bitVolts channels DATA EVENTS D_EL FileName TTL fullfilename