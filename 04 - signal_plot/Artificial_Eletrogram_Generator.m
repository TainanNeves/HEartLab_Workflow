%CÓDIGO HEartLab - Artificial Eletrogram Signal Generator - Versão 1 - 04/07/2024
%Author - José Junior

%-------------------------------------------------------------------
% Clear the workspace and command window
clear; clc;

%-------------------------------------------------------------------
% Define constants for signal generation
MA = 1;   % Amplitude for the Monophasic Signal
MN = 50;  % Exponent for the Monophasic Signal
Mk = 1.35; % Frequency factor for the Monophasic Signal
MD = -3.3; % Phase shift for the Monophasic Signal

% Constants for the six different biphasic components
B1A = 0.9;   % Amplitude of the 1st Biphasic Component
B1N = 300;   % Exponent of the 1st Biphasic Component
B1k = 1.1;   % Frequency factor of the 1st Biphasic Component
B1D = -3.10; % Phase shift of the 1st Biphasic Component

B2A = -0.9;  % Amplitude of the 2nd Biphasic Component
B2N = 400;   % Exponent of the 2nd Biphasic Component
B2k = 1.1;   % Frequency factor of the 2nd Biphasic Component
B2D = -3.18; % Phase shift of the 2nd Biphasic Component

B3A = 0.2;   % Amplitude of the 3rd Biphasic Component
B3N = 256;   % Exponent of the 3rd Biphasic Component
B3k = 1.1;   % Frequency factor of the 3rd Biphasic Component
B3D = -3;    % Phase shift of the 3rd Biphasic Component

B4A = 0.08;  % Amplitude of the 4th Biphasic Component
B4N = 256;   % Exponent of the 4th Biphasic Component
B4k = 1.1;   % Frequency factor of the 4th Biphasic Component
B4D = -2.9;  % Phase shift of the 4th Biphasic Component

B5A = -0.2;  % Amplitude of the 5th Biphasic Component
B5N = 256;   % Exponent of the 5th Biphasic Component
B5k = 1.1;   % Frequency factor of the 5th Biphasic Component
B5D = -3.2;  % Phase shift of the 5th Biphasic Component

B6A = -0.08; % Amplitude of the 6th Biphasic Component
B6N = 256;   % Exponent of the 6th Biphasic Component
B6k = 1.1;   % Frequency factor of the 6th Biphasic Component
B6D = -3.3;  % Phase shift of the 6th Biphasic Component

% Define the time vector for 8 seconds with a sampling rate of 4000 Hz
time = linspace(0, 8, 8*4000);

%-------------------------------------------------------------------
% Generate the Monophasic Signal using a sinusoidal function
VOM = MA * (sin(Mk*time + MD)).^MN;

% Generate each Biphasic Component
g = B1A * (sin(B1k*time+B1D).^B1N); % 1st Component
h = B2A * (sin(B2k*time+B2D).^B2N); % 2nd Component
j = B3A * (sin(B3k*time+B3D).^B3N); % 3rd Component
k = B4A * (sin(B4k*time+B4D).^B4N); % 4th Component
l = B5A * (sin(B5k*time+B5D).^B5N); % 5th Component
p = B6A * (sin(B6k*time+B6D).^B6N); % 6th Component

% Combine all Biphasic Components to generate the final Biphasic Signal
VOB = g + h + j + k + l + p;

%-------------------------------------------------------------------
% Clear variables that are no longer needed
clear B1A B2A B3A B4A B5A B6A B1N B2N B3N B4N B5N B6N
clear B1D B2D B3D B4D B5D B6D B1k B2k B3k B4k B5k B6k
clear g h j k l p
clear MA MN MD Mk

%-------------------------------------------------------------------
% Plot the generated signals

% Create a figure with two subplots
figure;

% Plot the Monophasic Signal
subplot(2, 1, 1);
plot(time, VOM); % Plot the Monophasic Signal
title('Monophasic Signal'); % Title for the Monophasic Signal Plot
xlabel('Time [s]'); % X-axis label
ylabel('Amplitude [a.u.]'); % Y-axis label
grid on; % Turn on the grid
ylim([-1, max(VOM)*2]); % Set y-axis limits

% Plot the Biphasic Signal
subplot(2, 1, 2);
plot(time, VOB); % Plot the Biphasic Signal
title('Biphasic Signal'); % Title for the Biphasic Signal Plot
xlabel('Time [s]'); % X-axis label
ylabel('Amplitude [a.u.]'); % Y-axis label
grid on; % Turn on the grid
ylim([-1, max(VOB)*1.1]); % Set y-axis limits
xlim([min(time), max(time)]); % Set x-axis limits to cover the entire time range

%-------------------------------------------------------------------
% Save the generated signals to a .mat file
save('Eletrograma_data.mat', 'VOM', 'VOB'); % Save Monophasic and Biphasic Signals