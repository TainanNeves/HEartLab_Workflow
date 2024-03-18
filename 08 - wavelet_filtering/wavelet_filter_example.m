%
%wavelet_filter function example of usage:
%
%Close all figures and clear workspace
close all;
clear all;

%Load signal
load('electric_data_E18_F2_R1_raw_allsignals.mat');

%Define parameters
electrodes = [1:16];                    % Electrodes to filter
waveletType = {'db4'};                  % Wavelet type
numLevels = 10;                         % Number of decomposition levels
reconstructionLevelsSets = {[8,9]};     % Levels to be used on reconstruction
fs  = 4000;                             % Sampling frequency

% Call the wavelet_filter function
% ensure D_EL is a struct, as exported from extraction_filtering function
wavelet_filter(D_EL, fs, electrodes, waveletType, numLevels, reconstructionLevelsSets);
