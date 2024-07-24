%CÓDIGO HEartLab - Artificial Electrogram Signal Generator - Version 3 - 20/07/2024
% This code generates artificial electrogram signals.
% The signals are generated through sinusoidal functions of different orders.
% It is possible to choose the electrical amplitudes of the Monophasic Electrogram and
% the Biphasic Electrogram.
% You can choose the heart rate in BPM and the sampling frequency.
% A figure of the Electrogram is plotted.
% The generated signals work with a sampling frequency of 4kHz and 8s.
% The code is simple, just run each section.
% Author - José Junior

% Mathematical Modeling Periods
% Monophasic Wave Period is 2.34 s -> f = 0.427 Hz -> 25.64 BPM
% Biphasic Wave Period is 2.86 s -> f = 0.349 Hz -> 20.97 BPM

clear; clc;

%% Sinusoidal Function Values

%-------------------------------------------------------------------
% Monophasic Electrogram Values
MA = 1; % Amplitude [a.u.]
MN = 50; % Order of the Sinusoidal Function - Order
Mk = 1.35; % Mathematical Modeling Frequency Factor - FMM

%-------------------------------------------------------------------
% Biphasic Electrogram Values
 %Amplitude [a.u.]    %Order           %FMM     %Phase Shift
    B1A = 0.9;      B1N = 300;      B1k = 1.1;      B1D = -3.10;
    B2A = -0.9;     B2N = 400;      B2k = 1.1;      B2D = -3.18;
    B3A = 0.2;      B3N = 256;      B3k = 1.1;      B3D = -3;
    B4A = 0.08;     B4N = 256;      B4k = 1.1;      B4D = -2.9;
    B5A = -0.2;     B5N = 256;      B5k = 1.1;      B5D = -3.2;
    B6A = -0.08;    B6N = 256;      B6k = 1.1;      B6D = -3.3;

% Sampling Time and Frequency Values
FA = 4000; % Sampling frequency
t = linspace(0,8,8*FA);

%% Determination of Heart Rate in BPM

%---Monophasic Wave------------------------------------------------------
FFM = 1; % Heart Rate - do not alter!!!

%Choose the Heart Rate by Entering a Value in BPM
fcm = 60; % Heart Rate Value in BPM

%Assign Value to Monophasic Wave
FFM = FFM * (fcm/25.641);

%---Biphasic Wave--------------------------------------------------------
FFB = 1; % Heart Rate Frequency Factor - do not alter!!!

%Choose the Heart Rate by Entering a Value in BPM
fcb = 60; % Heart Rate Value in BPM

%Assign Value to Biphasic Wave
FFB = FFB * (fcb/20.979);

%% Amplitude Modifier [a.u.]

MAOM = 1; % Constant for Calculation - Do not alter!!!
MAOB = 1; % Constant for Calculation - Do not alter!!!

%Modifier for Monophasic Wave
mao = 1; % Amplitude Value in a.u.
MAOM = MAOM * (mao/1);

%Modifier for Biphasic Wave
mab = 1; % Amplitude Value in a.u.
MAOB = MAOB * (mab/1);

%% Generation of Artificial Electrogram Waveforms

%-------------------------------------------------------------------
%Monophasic
VOM = MAOM * MA * (sin(FFM*Mk*t)).^MN; %Monophasic Wave Function

%-------------------------------------------------------------------
%Biphasic
B1 = MAOB * B1A * (sin(FFB*B1k*t+B1D).^B1N); %Function 1 of Biphasic Wave
B2 = MAOB * B2A * (sin(FFB*B2k*t+B2D).^B2N); %Function 2 of Biphasic Wave
B3 = MAOB * B3A * (sin(FFB*B3k*t+B3D).^B3N); %Function 3 of Biphasic Wave
B4 = MAOB * B4A * (sin(FFB*B4k*t+B4D).^B4N); %Function 4 of Biphasic Wave
B5 = MAOB * B5A * (sin(FFB*B5k*t+B5D).^B5N); %Function 5 of Biphasic Wave
B6 = MAOB * B6A * (sin(FFB*B6k*t+B6D).^B6N); %Function 6 of Biphasic Wave

VOB = B1 + B2 + B3 + B4 + B5 + B6; %General Function of Biphasic Wave

%-------------------------------------------------------------------
clear B1A B2A B3A B4A B5A B6A B1N B2N B3N B4N B5N B6N
clear B1D B2D B3D B4D B5D B6D B1k B2k B3k B4k B5k B6k
clear B1 B2 B3 B4 B5 B6
clear MA MN MD Mk

%% Plotting the Artificial Electrogram

% Creating Figure with both Electrograms
figure;

% Define time limits
lower_limit_time_mono = 1;  % adjust as needed for mono wave
upper_limit_time_mono = 4;  % adjust as needed for mono wave

lower_limit_time_bi = 1;  % adjust as needed for bi wave
upper_limit_time_bi = 4;  % adjust as needed for bi wave

% Monophasic Wave Subplot
subplot(2, 1, 1);
plot(t, VOM, 'r');
title(sprintf('Monophasic Morphology Signal - Heart Rate: %.2f - BPM', fcm));
xlabel('Time - [s]');
ylabel('Amplitude - [a.u.]');
grid on;
ylim([min(VOM)-0.2, max(VOM)+0.2]);
xlim([lower_limit_time_mono, upper_limit_time_mono]);  % Adjust X-axis limits

% Biphasic Wave Subplot
subplot(2, 1, 2);
plot(t, VOB, 'b');
title(sprintf('Biphasic Morphology Signal - Heart Rate: %.2f - BPM', fcb));
xlabel('Time - [s]');
ylabel('Amplitude - [a.u.]');
grid on;
ylim([min(VOB)-0.2, max(VOB)+0.2]);
xlim([lower_limit_time_bi, upper_limit_time_bi]);  % Adjust X-axis limits

% Create Formatted Title for ECG Figure
title_text = sprintf('Mono and Biphasic Signal Electrogram Artificial - HeartLab');

% Add General Title to the Figure
sgtitle(title_text);  % General Figure Title

clear lower_limit_time_bi lower_limit_time_mono
clear upper_limit_time_bi upper_limit_time_mono

%% Exporting the Signals

% Saving Generated Signals and Sampling Time and Frequency in a '.mat' File
save('Electrogram_data.mat', 'VOM', 'VOB', 't', 'FA');

clear MAOM MAOB mao mab FFM FFB fcm fcb
clear FA t title_text

%% Close Figures

close all;
