% Preprocessing Script for ECGi Data
%  
% Author: Angélica Drielly Santos de Quadros  
% Affiliation: Federal University of ABC (UFABC), Brazil  
% Year: 2025  
%  
% This script preprocesses MEAs and tank signals to enhance signal quality and prepare data for ECGi.
% The preprocessing steps include:  
%  
% 1. *Mean Subtraction*  
%    - Removes the DC component by centering each signal around zero.   
%    - Plots to compare signals before and after subtraction are generated.  
%
% 2. *Interpolation of Tank Signals*  
%    - Tank signals are interpolated to match the torso-tank geometry.   
%    - Plots to compare signals before and after interpolation are generated. 
%
% 3. *Detrending*  
%    - Baseline drift is removed using a spline-based detrending method,  
%      utilizing the function developed by Professor Oscar Barquero.  
%    - Plots to compare the signals before and after detrending are generated.  
%
% 4. *Spectral Analysis*  
%    - Frequency spectrum components of the signals are plotted.  
%    - This analysis aids in selecting the optimal cutoff frequency for the  
%      low-pass filter by inspecting the frequency spectrum.  
%
% 5. *Low-Pass Filtering*  
%    - A Butterworth filter is applied to remove high-frequency noise.
%    - The cutoff frequency is selected based on the spectrum analysis.  
%    - Plots to compare the signals before and after filtering are generated.  
%
% 6. *Downsampling*  
%    - The signals are downsampled to optimize computational efficiency.  
%    - Plots to compare the signals before and after downsampling are generated.  
%
% 7. *Upsampling*  
%    - Signals can be upsampled when alignment with the original acquisition frequency is required.  
%    - Plots to compare signals before and after upsampling are generated.  
%
% 8. *Normalization*  
%    - Signals are normalized based on the maximum value of the row.  
%    - Plots to compare the signals before and after normalization are generated.  
%
% 9. *Scaling*  
%    - Scaling is applied if adjustment of signal amplitude units is required.  
%    - Plots to compare the signals before and after scaling are generated.  
%
% 10. *Principal Component Analysis (PCA via SVD)*  
%     - PCA is applied as a project-and-reproject filter to the downsampled,
%       normalized tank signals.  
%     - SVD is used to compute principal components; explained variance is plotted.    
%     - Data is reconstructed using the selected components and compared to the  
%       original signal.  
%
% 11. *Data Storage*  
%     - The preprocessed signals are saved in .mat files for further  
%       ECGi reconstruction and analysis.  
%     - Figures generated during preprocessing can also be stored.  
%  
% The script processes both measured signals from the MEAs and tank torso separately.  
% It includes visualization steps to evaluate the impact of each preprocessing technique on the data.

%% Loading data

% MEAs and Tank raw signals
file_name =  "../../01 - data/electric_data_Exx_Fxx_Rxx_raw.mat";
signals = load(file_name);


% Extracting signals
meas_signal_raw = signals.D_EL.Data([1:32, 65:80], :);
tank_signal_raw = signals.D_EL.Data([129:174, 177:190], :);


% Tank electrodes indices
el_tank_idx = load("../../tank_el_idx.mat");
el_tank_idx = el_tank_idx.(subsref(fieldnames(el_tank_idx),substruct('{}',{1})));

% Tank geometry
tank_geo_file = ("../../LR_smoothed_tank.mat");
tank_geo = load(tank_geo_file);
tank_geo = tank_geo.(subsref(fieldnames(tank_geo),substruct('{}',{1})));


%% Subtracting the mean from the tank tank signals

mean_tank = mean(tank_signal_raw);
tank_signal_sub = tank_signal_raw - mean_tank;

%% Interpolating tank signals

y = tank_signal_sub;
idx = el_tank_idx';

[lap, edge] = mesh_laplacian(tank_geo.vertices, tank_geo.faces);
[int, keepindex, repindex] = mesh_laplacian_interp_current(lap, idx);

% Interpolate and keep measured data
tank_signal_interp = int * y;
for i = 1:60
    tank_signal_interp(idx(1, i), :) = y(i, :);
end

tank_electrodes = [129:174, 177:190];
el_map = [tank_electrodes(:) idx(:)];  
% el_map(:,1): original electrode number [129-192]
% el_map(:,2): corresponding node index in the interpolated geometry (218 nodes total;
% measured electrodes keep their original position)

%% Detrend Signals

meas = meas_signal_raw;
tank = tank_signal_interp;

fs = 4000;

% Define parameters
L_meas = size(meas,2);
t_meas = 1/fs*(0:L_meas-1);
L_tank = size(tank,2);
t_tank = 1/fs*(0:L_tank-1);

in_time = 1;
end_time = 3;
in_sample = in_time * fs + 1;
end_sample = end_time * fs;
time = linspace(in_time, end_time, length(in_sample:end_sample));

% Detrend MEA signals
meas_signal_d = zeros(size(meas));
for m = 1:size(meas,1)
    meas_signal_d(m,:) = detrendSpline(meas(m,:), t_meas, 0.2);
end

% Detrend tank signals
tank_signal_d = zeros(size(tank));
for m = 1:size(tank,1)
    tank_signal_d(m,:) = detrendSpline(tank(m,:), t_tank, 0.2);
end

%% Plot MEAs Signals before and after detrending
% Plot MEA Signals (RA, V, LA) before and after detrending

for mea_idx = 1:3
    % Set the indices for each MEA group
    if mea_idx == 1
        mea_title = 'RIGHT ATRIUM';
          = 1:16;
    elseif mea_idx == 2
        mea_title = 'VENTRICLE';
        mea_range = 17:32;
    else
        mea_title = 'LEFT ATRIUM';
        mea_range = 65:80;
        
    end
    % Loop through each set of 4 electrodes
    for i = 1:4:16
        figure(); % Create a new figure for each group of 4 electrodes
        sgtitle(['MEA ' num2str(mea_idx) ' - ' mea_title], 'FontWeight', 'bold');
        
        % Loop through the 4 electrodes in the group
        for j = 0:3
            electrode_idx = i + j;
            subplot(2, 2, j+1); % Create 2x2 subplot grid
            plot(time, meas_signal_raw(electrode_idx + (mea_idx-1)*16, in_sample:end_sample), 'color', "#00FFFF", 'LineWidth', 1.5);
            hold on;
            plot(time, meas_signal_d(electrode_idx + (mea_idx-1)*16, in_sample:end_sample), 'color', "#0000FF", 'LineWidth', 1);
            title(['Electrode ' num2str(electrode_idx)]);
            xlabel('Time (s)');
            ylabel('Amplitude (\muV)');
            legend('Original', 'Detrended');
        end
        
        % Save each figure
%         saveas(gcf, ['MEA_' num2str(mea_idx) '_Electrodes_' num2str(mea_range(i)) '-' num2str(mea_range(i+3)) '_detrend.png']);
      
    end
end

%% Plot Tank Signals before and after detrending


for i = 1:4:60
    figure(); % Create a new figure for each group of 4 tank electrodes
    sgtitle('Tank Signals', 'FontWeight', 'bold');
    
    % Loop through the 4 tank electrodes in the group
    for j = 0:3
        electrode_idx = i + j; % Tank index in the 60 rows
                
        subplot(2, 2, j+1); % Create 2x2 subplot grid
        
        % Plot Tank signals
        plot(time, tank_signal_interp(el_map(electrode_idx,1), in_sample:end_sample), 'color', "#ADD8E6");
        hold on; 
        plot(time, tank_signal_d(el_map(electrode_idx,1), in_sample:end_sample), 'color', "#00008B"); 
        title(['Electrode ' num2str(el_map(electrode_idx,2))]);
        xlabel('Time (s)');
        ylabel('Amplitude ($\mu$V)');
        legend('Original', 'Detrended');
    end
    
    % Save each figure
%     saveas(gcf, ['Tank_Electrodes_' num2str(tank_electrodes(i)) '-' num2str(tank_electrodes(i+3)) '_detrend.png']);
end

%% Low-Pass Filtering

meas = meas_signal_d;
tank = tank_signal_d;

% Define cutoff frequency
fs = 4000;
f_cut = 30 / (fs / 2); % set cutoff frequency
[b, a] = butter(4, f_cut, 'low');

% Apply butterworth filtering
meas_signal_f = filtfilt(b, a, meas')';
tank_signal_f = filtfilt(b, a, tank')';

%% Define Time and Frequency Parameters

in_time = 0;
end_time = 4;
in_sample = in_time * fs + 1;
end_sample = end_time * fs;
time = linspace(in_time, end_time, length(in_sample:end_sample));

f_up = 50; % choose frequency
f_down = 0.5; % choose frequency


%% Plot MEAs spectrum and electrical signals

for mea_idx = 1:3
    % Set the indices for each MEA group
    if mea_idx == 1
        mea_range = 1:16; 
        mea_title = 'RIGHT ATRIUM';
    elseif mea_idx == 2
        mea_range = 17:32; 
        mea_title = 'VENTRICLE';
    else
        mea_range = 33:48; 
        mea_title = 'LEFT ATRIUM';
    end

    if mea_idx == 3
        plot_range = 65:80; 
    else
        plot_range = mea_range; 
    end

    mea_temp = meas_signal_f(mea_range, in_sample:end_sample); % Extract the signals for current MEA
    [MFFTi, Sffti, fstep] = f_DF_electric(mea_temp, fs, f_up, f_down);

    % Loop through each set of 2 electrodes per plot (adjusted for indexing)
    for i = 1:2:16 % Grouping 2 electrodes per plot
        figure(); % Create a new figure for each group of 2 electrodes
        
        % Loop through the 2 electrodes in the group
        for j = 0:1
            electrode_idx = i + j;

            % Plot Frequency Spectrum for the current electrode 
            subplot(2, 2, j*2 + 1); 
            plot(fstep:fstep:f_up, Sffti(electrode_idx, :)); % Adjusted index for correct electrode
            xlabel('Frequency (Hz)');
            ylabel('Power');
            title(['Frequency Spectrum - Electrode ' num2str(plot_range(electrode_idx))],'FontSize',15);

            % Plot Electrical Signal for the current electrode 
            subplot(2, 2, j*2 + 2); 
            plot(time, meas_signal_d(mea_range(electrode_idx), in_sample:end_sample), 'r', 'LineWidth', 1);
            hold on;
            plot(time, meas_signal_f(mea_range(electrode_idx), in_sample:end_sample), 'k', 'LineWidth', 1);
            title(['Electrical Signal - Electrode ' num2str(plot_range(electrode_idx)) ' - Cutoff: ' num2str(f_cut * fs / 2) ' Hz'], 'FontSize',15);
            xlabel('Time (s)');
            ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex', 'FontSize', 12);
%           ylim([-300, 300]);
            legend('Original', 'New filtering');
        end
        
        % Add a horizontal line to separate the rows 
        annotation('line', [0.1, 0.9], [0.5, 0.5], 'LineWidth', 1, 'Color', 'k');
        
        % Save each figure
%         saveas(gcf, ['MEA_' num2str(mea_idx) '_Electrodes_' num2str(plot_range(i)) '-' num2str(plot_range(i+1)) '_spectrum.png']);
    end
end


%% Plot tank spectrum and electrical signals

tank_temp = tank_signal_f(:, in_sample:end_sample); % Extract tank signals
[MFFTi, Sffti, fstep] = f_DF_electric(tank_temp, fs, f_up, f_down);

for i = 1:2:60  % Grouping 2 electrodes per plot
    figure();
    sgtitle('Tank Signals', 'FontWeight', 'bold');

    for j = 0:1
        electrode_idx = i + j;
        
        % Spectrum 
        subplot(2, 2, 1 + j*2);
        plot(fstep:fstep:f_up, Sffti(electrode_idx, :));
        xlabel('Frequency (Hz)');
        ylabel('Power');
        title(['Frequency Spectrum - Electrode ' num2str(el_map(electrode_idx,1))]);

        % Electrical Signals 
        subplot(2, 2, 2 + j*2);
        plot(time, tank_signal_d(electrode_idx, in_sample:end_sample), 'r', 'LineWidth', 1);
        hold on;
        plot(time, tank_signal_f(electrode_idx, in_sample:end_sample), 'k', 'LineWidth', 1);
        title(['Electrical Signal - Electrode ' num2str(el_map(electrode_idx,1)) ' - Cutoff: ' num2str(f_cut * fs / 2) ' Hz']);
        xlabel('Time (s)');
        ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex', 'FontSize', 12);
        legend('Original', 'New filtering');
    end
    
    % Add a horizontal line to separate rows
    annotation('line', [0.1, 0.9], [0.5, 0.5], 'Color', 'k', 'LineWidth', 1.5);
    
    % Save each figure
%     saveas(gcf, ['TANK_Electrodes_' num2str(el_map(i,2)) '-' num2str(el_map(i+1,2)) '_spectrum.png']);

end

%% Plot MEAs and Tank

% Define the sampling frequency (modify as needed)
fs = 4000; % Hz 

% Define time range
start_time = 0; % Modify as needed
end_time = 6;   % Modify as needed

% Convert time to sample indices
start_sample = round(start_time * fs) + 1;
end_sample = round(end_time * fs);

% Load electric data
meas_plot = meas_signal_f;
tank_plot = tank_signal_f;

% Time vectors
To = linspace(0, length(meas_plot(1, :)) / fs, length(meas_plot(1, :)));
time = linspace(start_time, end_time, length(start_sample:end_sample));

% Define electrode mappings
el1 = 7;   % MEA1 RA
el2 = 10;  % MEA1 RA
el3 = 71;  % MEA3 LA
el4 = 74;  % MEA3 LA
el5 = 17;  % MEA2 V
el6 = 18;  % MEA2 V
el7 = 134; % Tank
el8 = 184; % Tank

% Mapping the tank electrodes numbers
el7_plot = find(el_map(:,2) == el7, 1);
el8_plot = find(el_map(:,2) == el8, 1);

% Mapping the MEA3 electrodes numbers (subtract offset)
el3_plot = el3 - 32;
el4_plot = el4 - 32;

% Create figure
f1 = figure('color', 'white', 'Position', [40 40 600 800]);

% Font size for improved readability
fontSize = 8;

% Subplot 1 - MEA1 RA (Electrode el1)
subplot(8, 1, 1)
plot(To, meas_plot(el1, :), 'LineWidth', 1.2, 'Color', 'black');
ylabel('$\mu$V', 'Interpreter', 'latex', 'FontSize', fontSize);
title(['MEA1 RA - Electrode ', num2str(el1)], 'FontSize', fontSize);
set(gca, 'FontSize', fontSize);
xlim([start_time end_time]);

% Subplot 2 - MEA1 RA (Electrode el2)
subplot(8, 1, 2)
plot(To, meas_plot(el2, :), 'LineWidth', 1.2, 'Color', 'black');
ylabel('$\mu$V', 'Interpreter', 'latex', 'FontSize', fontSize);
title(['MEA1 RA - Electrode ', num2str(el2)], 'FontSize', fontSize);
set(gca, 'FontSize', fontSize);
xlim([start_time end_time]);

% Subplot 3 - MEA3 LA (Electrode el3)
subplot(8, 1, 3)
plot(To, meas_plot(el3_plot, :), 'LineWidth', 1.2, 'Color', 'black');
ylabel('$\mu$V', 'Interpreter', 'latex', 'FontSize', fontSize);
title(['MEA3 LA - Electrode ', num2str(el3)], 'FontSize', fontSize);
set(gca, 'FontSize', fontSize);
xlim([start_time end_time]);

% Subplot 4 - MEA3 LA (Electrode el4)
subplot(8, 1, 4)
plot(To, meas_plot(el4_plot, :), 'LineWidth', 1.2, 'Color', 'black');
ylabel('$\mu$V', 'Interpreter', 'latex', 'FontSize', fontSize);
title(['MEA3 LA - Electrode ', num2str(el4)], 'FontSize', fontSize);
set(gca, 'FontSize', fontSize);
xlim([start_time end_time]);

% Subplot 5 - MEA2 V (Electrode el5)
subplot(8, 1, 5)
plot(To, meas_plot(el5, :), 'LineWidth', 1.2, 'Color', 'black');
ylabel('$\mu$V', 'Interpreter', 'latex', 'FontSize', fontSize);
title(['MEA2 V - Electrode ', num2str(el5)], 'FontSize', fontSize);
set(gca, 'FontSize', fontSize);
xlim([start_time end_time]);

% Subplot 6 - MEA2 V (Electrode el6)
subplot(8, 1, 6)
plot(To, meas_plot(el6, :), 'LineWidth', 1.2, 'Color', 'black');
ylabel('$\mu$V', 'Interpreter', 'latex', 'FontSize', fontSize);
title(['MEA2 V - Electrode ', num2str(el6)], 'FontSize', fontSize);
set(gca, 'FontSize', fontSize);
xlim([start_time end_time]);

% Subplot 7 - Tank (Electrode el7)
subplot(8, 1, 7)
plot(To, tank_plot(el7_plot, :), 'LineWidth', 1.2, 'Color', 'red');
ylabel('$\mu$V', 'Interpreter', 'latex', 'FontSize', fontSize);
title(['TANK - Electrode ', num2str(el7)], 'FontSize', fontSize);
set(gca, 'FontSize', fontSize);
xlim([start_time end_time]);
% 
% % Subplot 8 - Tank (Electrode el8)
subplot(8, 1, 8)
plot(To, tank_plot(el8_plot, :), 'LineWidth', 1.2, 'Color', 'red');
ylabel('$\mu$V', 'Interpreter', 'latex', 'FontSize', fontSize);
title(['TANK - Electrode ', num2str(el8)], 'FontSize', fontSize);
xlabel('Time (s)', 'FontSize', fontSize);
% set(gca, 'FontSize', fontSize_title);
xlim([start_time end_time]);

% Link all subplots for synchronized scrolling
linkaxes(findall(gcf, 'Type', 'axes'), 'x');

%% Save Filtered Signals
lp_cutoff = f_cut * fs / 2;
save('filtered_signals.mat', 'meas_signal_f', 'tank_signal_f', 'lp_cutoff');

%% Downsampling

meas = meas_signal_f;
tank = tank_signal_f;

% Downsampling factor
factor = 20;

% MEA signals downsampling
meas_signal_down = resample(meas', 1, factor)'; % Downsample MEA signals

% Tank signals downsampling
tank_signal_down = resample(tank', 1, factor)'; % Downsample Tank signals

% Update the sampling frequency
fs_down = fs / factor;

%% Upsampling

meas = meas_signal_down;
tank = tank_signal_down;

% Downsampling factor
factor = 20;

% MEA signals usampling
meas_signal_up = resample(meas', factor, 1)'; % Upsample MEA signals

% Tank signals upsampling
tank_signal_up = resample(tank', factor, 1)'; % Upsample Tank signals

% Update the sampling frequency
fs_up = fs_down * factor;

%% Plot MEAs electrical signals before and after downsampling

in_time = 2;
end_time = 6;

% Define time vectors for original and downsampled signals
t_original = (in_time*fs+1:end_time*fs) / fs; % Time vector for the original signal
t_down = (in_time*fs_down+1:end_time*fs_down) / fs_down; % Time vector for the downsampled signal

for mea_idx = 1:3
    % Set the indices for each MEA group
    if mea_idx == 1
        mea_range = 1:16;
        mea_title = 'RIGHT ATRIUM';  
    elseif mea_idx == 2
        mea_range = 17:32; 
        mea_title = 'VENTRICLE'; 
    else
        mea_range = 33:48;
        mea_title = 'LEFT ATRIUM'; 
    end

    if mea_idx == 3
        plot_range = 65:80;
    else
        plot_range = mea_range;
    end


    % Loop through each set of 2 electrodes per plot
    for i = 1:4:16 % Grouping 2 electrodes per plot
        figure();
        sgtitle(['MEA ' num2str(mea_idx) ' - ' mea_title], 'FontWeight', 'bold');


        % Loop through the 2 electrodes in the group
        for j = 0:3
            electrode_idx = i + j;
            
            % Plot Original and Downsampled Electrical Signals for the current electrode
            subplot(2, 2, j + 1);
            plot(t_original, meas_signal_f(mea_range(electrode_idx), in_time*fs+1:end_time*fs), 'r', 'LineWidth', 1);
            hold on;
            plot(t_down, meas_signal_down(mea_range(electrode_idx), in_time*fs_down+1:end_time*fs_down), 'b', 'LineWidth', 1);
            title(['Electrode ' num2str(plot_range(electrode_idx))]);
            xlabel('Time (s)');
            ylabel('Amplitude (\muV)');
            legend(['Original ' num2str(fs, '%.0f') ' Hz'], ['Downsampled ' num2str(fs_down, '%.0f') ' Hz']);

        end
         
        % Save each figure (optional)
%         saveas(gcf, ['MEA_' num2str(mea_idx) '_Electrodes_' num2str(plot_range(i)) '-' num2str(plot_range(i+1)) '_downsampling.png']);
    end
end
%% Plot Tank Signals before and after downsampling

in_time = 1;
end_time = 3;

% Define time vectors for original and downsampled signals
t_original = (in_time*fs+1:end_time*fs) / fs; % Time vector for the original signal
t_down = (in_time*fs_down+1:end_time*fs_down) / fs_down; % Time vector for the downsampled signal

for i = 1:4:60
    figure(); % Create a new figure for each group of 4 tank electrodes
    sgtitle('Tank Signals', 'FontWeight', 'bold');
    
    % Loop through the 4 tank electrodes
    for j = 0:3
        electrode_idx = i + j;
        
        % Plot Original and Downsampled Electrical Signals for the current electrode
        subplot(2, 2, j + 1);
        plot(t_original, tank_signal_f(el_map(electrode_idx,1), in_time*fs+1:end_time*fs), 'r', 'LineWidth', 1);
        hold on;
        plot(t_down, tank_signal_down(el_map(electrode_idx,1), in_time*fs_down+1:end_time*fs_down), 'b', 'LineWidth', 1);
        title(['Electrode ' num2str(el_map(electrode_idx,2))]);
        xlabel('Time (s)');
        ylabel('Amplitude (\muV)');
        legend(['Original ' num2str(fs, '%.0f') ' Hz'], ['Downsampled ' num2str(fs_down, '%.0f') ' Hz']);
        
    end
    
    % Save each figure
%   saveas(gcf, ['Tank_Electrodes_' num2str(el_map(i,2)) '-' num2str(el_map(i+3,2)) '_downsampling.png']);

end

%% Save Downsampled Signals
save('downsampled_EXX_FX_RX.mat', 'meas_signal_down', 'tank_signal_down','fs_down');

%% Normalizing signals

meas_signal_not_norm = meas_signal_down;
tank_signal_not_norm = tank_signal_down;

% Normalizing meas signal
row_max = max(abs(meas_signal_not_norm), [], 2);
meas_signal_norm = bsxfun(@rdivide, meas_signal_not_norm, row_max);

% % Normalizing tank signals
row_max = max(abs(tank_signal_not_norm), [], 2);
tank_signal_norm = bsxfun(@rdivide, tank_signal_not_norm, row_max);


%% Plotting MEAs signals before and after normalization

in_time = 1;
end_time = 3;
in_sample = in_time * fs_down + 1;
end_sample = end_time * fs_down;
time = linspace(in_time, end_time, length(in_sample:end_sample));

for mea_idx = 1:3
    % Set the indices for each MEA group
    if mea_idx == 1
        mea_range = 1:16;
        mea_title = 'RIGHT ATRIUM';  
    elseif mea_idx == 2
        mea_range = 17:32; 
        mea_title = 'VENTRICLE'; 
    else
        mea_range = 33:48;
        mea_title = 'LEFT ATRIUM'; 
    end

    if mea_idx == 3
        plot_range = 65:80;
    else
        plot_range = mea_range;
    end


    % Loop through each set of 2 electrodes per plot
    for i = 1:4:16 % Grouping 2 electrodes per plot
        figure();
        sgtitle(['MEA ' num2str(mea_idx) ' - ' mea_title], 'FontWeight', 'bold');


        % Loop through the 2 electrodes in the group
        for j = 0:3
            electrode_idx = i + j;
            
            % Plot Original and Downsampled Electrical Signals for the current electrode
            subplot(2, 2, j + 1);
            plot(time, meas_signal_down(mea_range(electrode_idx), in_sample:end_sample), 'r', 'LineWidth', 1);
            hold on;
            plot(time, meas_signal_norm(mea_range(electrode_idx), in_sample:end_sample), 'b', 'LineWidth', 1);
            title(['Electrode ' num2str(plot_range(electrode_idx))]);
            xlabel('Time (s)');
            ylabel('Amplitude (\muV)');
            legend('Original', 'Normalized');

        end
         
%         Save each figure
%         saveas(gcf, ['MEA_' num2str(mea_idx) '_Electrodes_' num2str(plot_range(i)) '-' num2str(plot_range(i+3)) '_norm.png']);

    end
end


%% Plotting Tank signals before and after normalization

in_time = 1;
end_time = 3;
in_sample = in_time * fs_down + 1;
end_sample = end_time * fs_down;
time = linspace(in_time, end_time, length(in_sample:end_sample));

tank_electrodes = [129:174, 177:190];

for i = 1:4:60
    figure(); % Create a new figure for each group of 4 tank electrodes
    sgtitle('Tank Signals - Standard normalization', 'FontWeight', 'bold');
    
    % Loop through the 4 tank electrodes
    for j = 0:3
        electrode_idx = i + j;
        electrode_idx = tank_electrodes(electrode_idx); % Map to the real electrode number
        
        % Plot Original and Downsampled Electrical Signals for the current electrode
        subplot(2, 2, j + 1);
        plot(time, tank_signal_down(electrode_idx, in_sample:end_sample), 'r', 'LineWidth', 1);
        hold on;
        plot(time, tank_signal_norm(electrode_idx, in_sample:end_sample), 'b', 'LineWidth', 1);
        title(['Electrode ' num2str(electrode_idx)]);
        xlabel('Time (s)');
        ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex');
        legend('Original', 'Normalized');

    end
%       Save each figure
%       saveas(gcf, ['Tank_Electrodes_' num2str(tank_electrodes(i)) '-' num2str(tank_electrodes(i+3)) '_norm.png']);

end

%% Scaling signals

% Scaling meas signal

meas_signal_scaled = zeros(size(meas_signal_norm)); % Preallocate
for i = 1:size(meas_signal_norm, 1)
    meas_signal_scaled(i, :) = rescale(meas_signal_norm(i, :), -1, 1);
end

% Scaling tank signals

tank_signal_scaled = zeros(size(tank_signal_norm)); % Preallocate
for i = 1:size(tank_signal_norm, 1)
    tank_signal_scaled(i, :) = rescale(tank_signal_norm(i, :), -1, 1);
end


%% Plotting MEA signals before and after normalization and scaling

in_time = 1;
end_time = 3;
in_sample = in_time * fs_down + 1;
end_sample = end_time * fs_down;
time = linspace(in_time, end_time, length(in_sample:end_sample));


for mea_idx = 1:1
    % Set the indices for each MEA group
    if mea_idx == 1
        mea_range = 1:16;
        mea_title = 'RIGHT ATRIUM';  
    elseif mea_idx == 2
        mea_range = 17:32; 
        mea_title = 'VENTRICLE'; 
    else
        mea_range = 33:48;
        mea_title = 'LEFT ATRIUM'; 
    end

    if mea_idx == 3
        plot_range = 65:80;
    else
        plot_range = mea_range;
    end

    % Loop through each set of 4 electrodes per plot
    for i = 1:4:16
        figure();
        sgtitle(['MEA ' num2str(mea_idx) ' - ' mea_title], 'FontWeight', 'bold');

        % Loop through the 4 electrodes in the group
        for j = 0:3
            electrode_idx = i + j;
            
            % Plot Original, Normalized, and Scaled Electrical Signals
            subplot(2, 2, j + 1);
            plot(time, meas_signal_norm(mea_range(electrode_idx), in_sample:end_sample), 'r', 'LineWidth', 1);
            hold on;
            plot(time, meas_signal_scaled(mea_range(electrode_idx), in_sample:end_sample), 'g', 'LineWidth', 1);
            title(['Electrode ' num2str(plot_range(electrode_idx))]);
            xlabel('Time (s)');
            ylabel('Amplitude (\muV)');
            legend('Normalized', 'Scaled');

        end

%         Save each figure
%         saveas(gcf, ['MEA_' num2str(mea_idx) '_Electrodes_' num2str(plot_range(i)) '-' num2str(plot_range(i+3)) '_norm.png']);

    end
end


%% Plotting Tank signals before and after normalization and scaling

in_time = 0;
end_time = 2;
in_sample = in_time * fs_down + 1;
end_sample = end_time * fs_down;
time = linspace(in_time, end_time, length(in_sample:end_sample));

tank_electrodes = [129:174, 177:190];

for i = 1:4:60
    figure(); % Create a new figure for each group of 4 tank electrodes
    sgtitle('Tank Signals', 'FontWeight', 'bold');
    
    % Loop through the 4 tank electrodes
    for j = 0:3
        electrode_idx = i + j;
        electrode_idx = tank_electrodes(electrode_idx); % Map to the real electrode number
        
        % Plot Original, Normalized, and Scaled Electrical Signals
        subplot(2, 2, j + 1);
        plot(time, tank_signal_norm(electrode_idx, in_sample:end_sample), 'r', 'LineWidth', 1);
        hold on;
        plot(time, tank_signal_scaled(electrode_idx, in_sample:end_sample), 'g', 'LineWidth', 1);
        title(['Electrode ' num2str(electrode_idx)]);
        xlabel('Time (s)');
        ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex');
        legend('Normalized', 'Scaled');

%       Save each figure
%       saveas(gcf, ['Tank_Electrodes_' num2str(tank_electrodes(i)) '-' num2str(tank_electrodes(i+3)) '_norm.png']);
    end
end

%% PCA
% Applies PCA using SVD to the downsampled, normalized tank signals

% Use tank downsampled signals
tank_data = tank_signal_down;

% Normalize rows
row_max = max(abs(tank_data), [], 2);
y_normalized = bsxfun(@rdivide, tank_data, row_max);

% Center the data
data_centered = y_normalized - mean(y_normalized, 2); % center each signal over time

% SVD
[U, S, V] = svd(data_centered, 'econ');

% Singular values and explained variance
singular_values = diag(S);
explained_variance = (singular_values.^2) / sum(singular_values.^2) * 100;

% Plot cumulative explained variance
figure;
plot(cumsum(explained_variance), 'o-');
title('Cumulative Explained Variance');
xlabel('Number of Components');
ylabel('Explained Variance (%)');
grid on;

% Scores (projection onto principal components)
scores = data_centered * V; % equivalent to U * S

% Visualize principal components alongside a reference MEA signal
fs_pca = fs_down; % sampling frequency of downsampled signals
tt = (0:size(y_normalized, 2)-1) / fs_pca;

% Reference signals
ref_signal_ra = meas_signal_down(7, :);  % MEA1 RA electrode 7 (adjust as needed)
ref_signal_v  = meas_signal_down(17, :); % MEA2 V  electrode 17

n_comp_plot = 10; % number of PCA components to display
n_rows = n_comp_plot + 2;

figure;
subplot(n_rows, 1, 1);
plot(tt, ref_signal_ra, 'b');
title('MEA1 RA - Reference electrode');
axis off;

subplot(n_rows, 1, 2);
plot(tt, ref_signal_v, 'b');
title('MEA2 V - Reference electrode');
axis off;

for i = 1:n_comp_plot
    subplot(n_rows, 1, i + 2);
    plot(tt, V(:, i), 'b');
    title(['Component ' num2str(i)], 'FontSize', 10);
    axis off;
end

% Select number of components to retain
k = 5; % number of components selected after inspecting the cumulative explained variance
% Reconstruct data using k components
V_k = V(:, 1:k);
scores_k = data_centered * V_k; % projection onto k components
reconstructed_data = scores_k * V_k'; % reconstruction in original space

% Verify reconstruction quality on one tank electrode
% Use el_map to find the correct row index for a given electrode number
target_electrode = 134; % adjust to electrode of interest
tank_row = find(el_map(:, 2) == target_electrode, 1);

if ~isempty(tank_row)
    figure;
    plot(y_normalized(tank_row, :) / max(abs(y_normalized(tank_row, :))));
    hold on;
    plot(reconstructed_data(tank_row, :) / max(abs(reconstructed_data(tank_row, :))));
    legend('Original (normalized)', ['Reconstructed (k=' num2str(k) ')']);
    title(['PCA Reconstruction - Tank Electrode ' num2str(target_electrode)]);
    xlabel('Sample');
    ylabel('Normalized Amplitude');
else
    warning('Electrode %d not found in el_map.', target_electrode);
end

%% Save tank signals after PCA
save('reconstructed_tank.mat', 'reconstructed_data');