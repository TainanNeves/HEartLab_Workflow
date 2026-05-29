%% SYNCHRONIZATION

% Synchronize optical and electrical data
% Export the variables to be analyzed

clear all;
close all;
clc;

%% Loading variables

load(); % Optical data Filtered
load(); % Electric data filtered

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


%% Keep the rawcamera index
D_SYNC.IMG.CAM1 = D_OP.D_CAM1_rawimage;
D_SYNC.IMG.CAM2 = D_OP.D_CAM2_rawimage;
D_SYNC.IMG.CAM3 = D_OP.D_CAM3_rawimage;


%% Keep ROIs
D_SYNC.ROI.ROI_1 = D_OP.ROI.ROI_1;
D_SYNC.ROI.ROI_2 = D_OP.ROI.ROI_2;
D_SYNC.ROI.ROI_3 = D_OP.ROI.ROI_3;


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


%% SYNCHRONIZATION CHECK
Fsampling = 4000;
t_sync = linspace(0, size(D_SYNC.CAM1, 3)/Fsampling, size(D_SYNC.CAM1, 3));

% --- CUSTOM ELECTRODE/CAM SELECTION ---
el_idx_cam1 = 91;
el_idx_cam2 = 20;
el_idx_cam3 = 3;

configs = {D_SYNC.CAM1, el_idx_cam1, 'V', 1; ...
           D_SYNC.CAM2, el_idx_cam2, 'LA', 2; ...
           D_SYNC.CAM3, el_idx_cam3, 'RA', 3};

% Figure
figure(4); clf; 
set(gcf, 'Color', 'w', 'Units', 'normalized');

% Initialize array to store axes handles for linkaxes
h_sub = zeros(1, 3);

for i = 1:3
    Data        = configs{i, 1};
    el_idx      = configs{i, 2};
    region_name = configs{i, 3};
    cam_num     = configs{i, 4};
    
    % --- Step 1: Selection ---
    Background = squeeze(Data(:,:, 2000)); 
    [x, y] = pick_up_a_trace(Background, Data, 1);
    p = [x(end), y(end)]; 
    
    % --- Step 2: Data Extraction & Normalization ---
    opt_signal = squeeze(Data(p(1), p(2), :));
    el_signal  = D_SYNC.EL(el_idx, :);
    
    % Min-Max Normalization (avoids scale mismatch)
    opt_norm = (opt_signal - min(opt_signal)) / (max(opt_signal) - min(opt_signal));
    el_norm  = (el_signal - min(el_signal)) / (max(el_signal) - min(el_signal));
    
    % --- Step 3: Explicitly focus back on Figure 4 ---
    figure(4); 
    h_sub(i) = subplot(3, 1, i); % Store the handle
    hold on;
    
    % Plot Electrical (Blue) and Optical (Black)
    plot(t_sync, el_norm, 'b', 'LineWidth', 1.5, 'DisplayName', ['EL ' num2str(el_idx)]);
    plot(t_sync, opt_norm, 'k', 'LineWidth', 1, 'DisplayName', sprintf('Opt Pixel (%d,%d)', p(1), p(2)));
    
    % Formatting
    title(['Sync: ' region_name ' (CAM' num2str(cam_num) ' vs EL' num2str(el_idx) ')']);
    ylabel('Norm. Amp');
    if i == 3, xlabel('Time (s)'); end
    xlim([0 8]);
    grid on;
    legend('show', 'Location', 'northeast');
    set(gca, 'FontSize', 11);
    
    drawnow; 
end
linkaxes(h_sub, 'x');
sgtitle('Multi-Camera Synchronization Check', 'FontWeight', 'bold');


%% Exporting work variables
FileName = '13_32_';
save(['data_filtered_sync_', FileName, '.mat'], 'D_SYNC', '-v7.3');




