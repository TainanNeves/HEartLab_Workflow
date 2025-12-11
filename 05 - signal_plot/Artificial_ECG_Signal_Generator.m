%CÓDIGO HEartLab - Artificial ECG Signal Generator - Version 3 - 19/07/2024
% This code generates artificial 12-lead ECG signals.
% The signals are generated through sinusoidal functions of different orders.
% It is possible to choose the electrical amplitudes (mV) of the P, Q, R, S, T waves of
% the 12 leads.
% It is possible to choose the heart rate in BPM and the sampling frequency.
% An Electrocardiogram Figure is plotted.
% The generated signals work with a sampling frequency of 4kHz and 8s.
% The code is simple, just run each section.
% Author - José Junior

clear; clc;

%% Values of Sine Functions

NP = 240; NQ = 9000; NR = 1800; NS = 3000; NT = 138;    NTT = 138;  % Exponents of sine functions
DP = 0.5; DQ = 0.08; DR = 0;    DS = -0.1; DT = -0.75; DTT = -0.6; % Phase shifts

FF = 3.149; % Frequency factor of the Mathematical Modeling
FA = 4000; % Sampling frequency

% Values of Amplitudes of P, Q, R, S, T, and TT waves for the 12 Leads
% ---The T wave is formed by two T and TT waves, due to its asymmetry---
Amplitudes = [
      % P        % Q       % R        % S       % T       % TT
    0.0500   -0.0500   0.7000   -0.0020   0.0010   0.0300;   % DI
    -0.080   -0.1000   -0.800   -0.1000   -0.100   -0.050;   % aVR
    0.1000   0.00000   0.1800   -0.3300   0.0500   0.0250;   % V1
    0.1000   0.00000   0.6500   -0.2700   0.1000   0.0500;   % V4
    0.1000   -0.1000   0.8000   -0.2000   0.1000   0.0500;   % DII
    0.0400   0.00000   0.2500   -0.1000   0.0500   0.0250;   % aVL
    0.1000   0.00000   0.2800   -0.5000   0.1200   0.0600;   % V2
    0.1000   -0.0500   0.8000   -0.2000   0.0500   0.0250;   % V5
    0.0700   -0.0750   0.5000   -0.1160   0.0500   0.0250;   % DIII
    0.0750   -0.0400   0.5500   -0.1160   0.1500   0.0625;   % aVF
    0.1000   0.00000   0.4500   -0.5000   0.1000   0.0500;   % V3
    0.1000   -0.0750   0.6000   -0.0020   0.0500   0.0250    % V6
];

% Time and sampling frequency
t = linspace(0,8,8*FA);

% Preallocation of array to store the ECG waves
O_ECG = cell(1, 12);

%% Determination of Heart Rate in BPM

% Choose the Heart Rate by entering a value in BPM
fc = 60;

% Assignment of Value
FF = FF * (fc/60);

%% Generation of Artificial ECG Waveforms

for i = 1:12
    % Extract the Amplitudes for the Current Lead
    AP = Amplitudes(i, 1);
    AQ = Amplitudes(i, 2);
    AR = Amplitudes(i, 3);
    AS = Amplitudes(i, 4);
    AT = Amplitudes(i, 5);
    ATT = Amplitudes(i, 6);
    
    % Calculate the Functions for the Current Lead
    FP = (AP * ((sin(t * FF + DP)).^NP));
    FQ = (AQ * ((sin(t * FF + DQ)).^NQ));
    FR = (AR * ((sin(t * FF + DR)).^NR));
    FS = (AS * ((sin(t * FF + DS)).^NS));
    FT = (AT * ((sin(t * FF + DT)).^NT)) + (ATT * ((sin(t * FF + DTT)).^NTT));
    
    % Sum the Functions to obtain the ECG Waveform
    OD = FP + FQ + FR + FS + FT;
    
    % Store the ECG Waveform in the Corresponding Cell
    O_ECG{i} = OD;
end

%% Plotting the Artificial Electrocardiogram (ECG)

% Creation of the Artificial Electrocardiogram (ECG) Figure
F = figure('color', 'white', 'Position', [40 40 1200 900]);

% List of Names of the 12 Leads
Lead = {'DI', 'aVR', 'V1', 'V4', 'DII', 'aVL', 'V2', 'V5', 'DIII', 'aVF', 'V3', 'V6'};
% List of Colors for each Lead Group
colors = {'r', 'b', 'k', 'k', 'r', 'b', 'k', 'k', 'r', 'b', 'k', 'k'};

% Define time limits
lower_time_limit = 1;  % adjust as necessary
upper_time_limit = 4;  % adjust as necessary

% Generate the Plot for each Lead
for i = 1:12
    subplot(3, 4, i);  % 3 rows and 4 columns of subplots
    plot(t, O_ECG{i}, 'LineWidth', 1.5, 'Color', colors{i});  % Plot the Waveform
    xlabel('Time (s)');  % X-axis
    ylabel('Amplitude (mV)');  % Y-axis
    title(Lead{i});  % Title of the subplot
    grid on;  % Add grid to the subplot
    ylim([-0.9, max(O_ECG{i}) + 0.5]);  % Adjust the Y-axis limits
    xlim([lower_time_limit, upper_time_limit]);  % Adjust the X-axis limits
end

% Create the Formatted Title of the ECG Figure
title = sprintf('Artificial Electrocardiogram - HEartLab - %d BPM', fc);

% Add a General Title for the Figure
sgtitle(title);  % General title of the figure

Time = t; % Time values
F_sample = FA; % Sampling Frequency
Figure_ECG = F; % Save the Figure with the Artificial Electrocardiogram (ECG) Plots
F_BPM = fc; % Heart Rate in BPM

clear AP AQ AR AS AT ATT DP DQ DR DS DT DTT NP NQ NR NS NT NTT 
clear FP FQ FR FS FT 

%% Export Files with the Generated Artificial ECG

I_R = 'Sample'; % File Name: Enter any name
FNR = ['ECG_Data_',I_R,'.mat']; % File Name

% Creating the Structure to Store the Signals and Corresponding Lead Names
ECG_Struct = struct();

% Creating the Variables
for j = 1:12
    ECG_Struct(j).LeadName = Lead{j}; % Store the lead name
    ECG_Struct(j).Signal = O_ECG{j};  % Store the corresponding signal
end

% Define the names of the columns and rows
colNames = {'P', 'Q', 'R', 'S', 'T', 'TT'};
rowNames = {'DI', 'aVR', 'V1', 'V4', 'DII', 'aVL', 'V2', 'V5', 'DIII', 'aVF', 'V3', 'V6'};

% Create the struct of the Amplitudes of the P, Q, R, S, T, and TT Waves of the 12 Leads
Amplitudes_Struct = struct();
for w = 1:length(rowNames)
    for z = 1:length(colNames)
        Amplitudes_Struct.(rowNames{w}).(colNames{z}) = Amplitudes(w, z);
    end
end

% Name of the File to Save the Figure
figureFileName = 'Artificial_Electrocardiogram.fig';

% Save the Figure
savefig(Figure_ECG, figureFileName);

% Saving the Variables
save(FNR, 'ECG_Struct', 'Amplitudes_Struct', 'Time', 'F_sample', 'F_BPM', '-v7.3');
disp(['Variables saved to ', FNR]);

clear i j w z FA FF fc t OD F
clear O_ECG Amplitudes colors Lead FNR I_R title figureFileName
clear colNames rowNames F_sample Time Figure_ECG F_BPM
clear lower_time_limit upper_time_limit

%% Close Images

close all;
