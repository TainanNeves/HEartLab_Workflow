% Clear workspace and command window (optional)
% clear
% clc
% ---------------------------
% 1. Define sampling frequency and basic parameters
% ---------------------------
Fsampling = 4000; % Sampling frequency (Hz)
Data = D_SYNC.EL; % Electrical data
% ---------------------------
% 2. Electrode indices and data extraction
% ---------------------------
array1_indices = 1:16;
array2_indices = 17:32;
array3_indices = 65:80;

MEA1 = double(Data(array1_indices, :));
MEA2 = double(Data(array2_indices, :));
MEA3 = double(Data(array3_indices, :));
%. Perform ICA on MEA1, MEA2, MEA3
% Example: specify components to use for MEA1
close all
componentsToUse = [8];
results1 = ica_analysis(MEA1, 'MEA 1', array1_indices, T, Fsampling, componentsToUse);
ICA_MEA1 = results1.Data_reconstructed;

% Example: specify components for MEA2
close all
componentsToUse = [1:16]; % Example: use all components
results2 = ica_analysis(MEA2, 'MEA 2', array2_indices, T, Fsampling, componentsToUse);
ICA_MEA2 = results2.Data_reconstructed;

% Example: specify components for MEA3
close all
componentsToUse = [4];
results3 = ica_analysis(MEA3, 'MEA 3', array3_indices, T, Fsampling, componentsToUse);
ICA_MEA3 = results3.Data_reconstructed;

% ---------------------------
% 7. Consolidate ICA results back into the dataset
% ---------------------------
D_SYNC.ICA_MEAs = double(Data); % Initialize with original
D_SYNC.ICA_MEAs(array1_indices, :) = ICA_MEA1; % Insert MEA1 ICA
D_SYNC.ICA_MEAs(array2_indices, :) = ICA_MEA2; % Insert MEA2 ICA
D_SYNC.ICA_MEAs(array3_indices, :) = ICA_MEA3; % Insert MEA3 ICA