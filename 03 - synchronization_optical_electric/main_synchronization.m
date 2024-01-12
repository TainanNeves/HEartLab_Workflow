%% SYNCHRONIZATION

% Synchronize optical and electrical data
% Export the variables to be analyzed

clear all;
close all;
clc;

%% Loading variables

load('C:\Users\HEartLab\Documents\GitHub\HEartLab\00 - examples\optic_data_E14_F3_R4_filtered.mat'); % Optical data Filtered
load('C:\Users\HEartLab\Documents\GitHub\HEartLab\00 - examples\electric_data_E14_F3_R4_filtered.mat'); % Electric data filtered

%% Sample frequency from 500 to 4000 Hz

S = 8; % Multiplication factor
Fsampling = 500; % Initial Fs

% CAM1
for i = 1:size(D_OP.D_CAM1_filtered, 1)
    for j = 1:size(D_OP.D_CAM1_filtered, 2)
        optical = squeeze(D_OP.D_CAM1_filtered(i, j, (1:round(S * Fsampling))));
        xq = 0:0.00025:S;
        lo = length(optical);
        To = linspace(0, lo / Fsampling, lo);
        D_SYNC.CAM1(i, j, :) = interp1(To, optical, xq);
    end
end
% CAM2
for i = 1:size(D_OP.D_CAM2_filtered, 1)
    for j = 1:size(D_OP.D_CAM2_filtered, 2)
        optical = squeeze(D_OP.D_CAM2_filtered(i, j, (1:round(S * Fsampling))));
        xq = 0:0.00025:S;
        lo = length(optical);
        To = linspace(0, lo / Fsampling, lo);
        D_SYNC.CAM2(i, j, :) = interp1(To, optical, xq);
    end
end
% CAM3
for i = 1:size(D_OP.D_CAM3_filtered, 1)
    for j = 1:size(D_OP.D_CAM3_filtered, 2)
        optical = squeeze(D_OP.D_CAM3_filtered(i, j, (1:round(S * Fsampling))));
        xq = 0:0.00025:S;
        lo = length(optical);
        To = linspace(0, lo / Fsampling, lo);
        D_SYNC.CAM3(i, j, :) = interp1(To, optical, xq);
    end
end
% Cleaning
clear ans Fsampling i j lo optical S To xq;


%% Electric synchronization

S = 8; % Multiplication factor
Fsampling = 4000; % Electric Fs

for i = 1:size(D_EL.Data, 1)
    channel = D_EL.Data(i, :);
    [valor, ido] = min(abs(D_EL.Timestamps - D_EL.TTL(1, 1)));
    ido = round(ido + (Fsampling * 1));
    idf = round(ido + (Fsampling * S));
    x = channel(ido:idf);
    D_SYNC.EL(i, :) = x;
end
% Cleaning
clear ans Fsampling i idf ido S valor x channel;


%% Exporting work variables

FileName = 'E14_F3_R4';
save(['data_filtered_sync_', FileName, '.mat'], 'D_SYNC', '-v7.3');




