
% clear
% clc
% Define sampling frequency
Fsampling = 4000;

% ROIV=roipoly()
Data=D_SYNC.EL;

Fs = 4000;             % sampling rate in Hz (set to your value if different)
win_ms = 1;           % desired window length in ms
W = max(1, round(win_ms/1000 * Fs));
if mod(W,2)==0, W = W + 1; end
Data = movmean(Data, W, 2);

N = size(Data(:,:), 2);
% n=size(opticos_A(:,:), 2)
T=linspace(0, N / Fsampling, N);
% Time vector for plotting
To = linspace(0, N / Fsampling, N);

% Electrode indices for each array
array1_indices = 139:143;
array2_indices = 145:149;
array3_indices = 65:80;
% 
% Extract data for each array
MEA1 = double(Data(array1_indices, :));
MEA2 = double(Data(array2_indices, :));
MEA3 = double(Data(array3_indices, :));

%% 
% close all
componentsToUse = [1];
results1 = ica_analysis(MEA1, 'MEA 1', array1_indices, To, Fsampling, componentsToUse);
ICA_MEA1 = results1.Data_reconstructed;

%% results1 
% close all
componentsToUse = [139 140 141 142 143]
% Perform  ICA analysis on each array
ica_analysis(MEA2, 'MEA 2', array2_indices, To, Fsampling, componentsToUse);
results2 = ica_analysis(MEA2, 'MEA 2', array2_indices, To, Fsampling, componentsToUse);
ICA_MEA2 = results2.Data_reconstructed;
%% MEA 3
close all
componentsToUse = [1,3:13,15,16]; % Use the first three ICs;
% Perform  ICA analysis on each array
ica_analysis(MEA3, 'MEA 3', array1_indices, To, Fsampling, componentsToUse);
results3 = ica_analysis(MEA3, 'MEA 3', array1_indices, To, Fsampling, componentsToUse);
ICA_MEA3 = results3.Data_reconstructed;
%% 

% Create consolidated matrix maintaining original structure
% D_SYNC.ICA_MEAs= double(D_SYNC.EL);  % Preserve all original data
D_SYNC.ICA_MEAs(array1_indices, :) = ICA_MEA1;  % Rows 1-16
D_SYNC.ICA_MEAs(array2_indices, :) = ICA_MEA2;  % Rows 17-32
D_SYNC.ICA_MEAs(array3_indices, :) = ICA_MEA3;  % Rows 65-80
%% %% 
% % %% 
% filter_window = 40; % Window size for moving average filter (adjust based on noise level if you want)
% signals = zeros(size(D_SYNC.ICA_MEAs));
% 
% % Apply filtering to each electrode signal
% for i = 1:size(D_SYNC.ICA_MEAs, 1)
%     Data1(i,:) = movmean(D_SYNC.ICA_MEAs(i,:), filter_window);
% end
% %% % Butterworth filter parameters
% % filterOrder = 4; % 4th order filter (you can adjust this)
% % lowCutoff = 0.5; % Hz (adjust based on your signal characteristics)
% % highCutoff = 80; % Hz (adjust based on your signal characteristics)
% % Fnyquist = Fsampling/2; % Nyquist frequency
% % 
% % % Design Butterworth bandpass filter
% % [b, a] = butter(filterOrder, [lowCutoff highCutoff]/Fnyquist, 'bandpass');
% % 
% % % Apply filter to MEA1
% % MEA3_filtered = zeros(size(MEA3));
% % for channel = 1:size(MEA3, 1)
% %     MEA3_filtered(channel, :) = filtfilt(b, a, MEA1(channel, :));
% % end