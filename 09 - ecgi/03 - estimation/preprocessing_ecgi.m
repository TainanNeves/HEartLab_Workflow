% Preprocessing Script for ECGi Data  
%  
% Author: Angélica Drielly Santos de Quadros  
% Affiliation: Federal University of ABC (UFABC), Brazil  
% Year: 2025  
%  
% This script preprocesses MEAs and tank signals
% to enhance signal quality and prepare data for ECGi. The preprocessing steps include:  
%  
% 1. **Spectral Analysis**  
%    - The frequency spectrum components of the signals are plotted.  
%    - This analysis aids in selecting the optimal cutoff frequency for the low-pass filter by inspecting the frequency spectrum.  
%  
% 2. **Low-Pass Filtering**  
%    - A Butterworth filter is applied to remove high-frequency noise while preserving relevant cardiac activity.  
%    - The cutoff frequency is selected based on the spectrum analysis results to retain the essential signal components. 
%    - Plots showing the signals before and after filtering are generated.
%  
% 3. **Detrending**  
%    - Baseline drift is removed using a spline-based detrending method, utilizing the function developed by Professor Oscar Barquero.  
%    - Plots showing the signals before and after detrending are generated.
%  
% 4. **Downsampling**  
%    - The signals are downsampled to optimize computational efficiency while maintaining signal fidelity.  
%    - Plots showing the signals before and after downsampling are generated.
%  
% 5. **Normalization**  
%    - The signals are normalized, mean = 0 and std deviation = 1.  
%    - Plots showing the signals before and after normalization are generated.
%
% 6. **Data Storage**  
%    - The filtered, detrended, and downsampled signals can be stored for further processing and analysis.  
%    - Figures generated during the analysis can be saved for later use.  
%  
% The script processes both measured signals from the MEAs and tank torso separately.  
% It includes visualization steps to evaluate the impact of each preprocessing technique on the data.

%% Loading the raw signal

file_name = "C:\Users\HeartLAB\Documents\Documents\Codes\Current codes\02 - extraction_filtering\electric\open_ephys_extraction_matlab_filtering\electric_data_E20_F01_R01_raw.mat";
data = load(file_name);
meas_signal_raw = data.D_EL.Data([1:32, 65:80], :);
tank_signal_raw = data.D_EL.Data([129:174, 177:190], :);


%% Low-Pass Filtering
% Define cutoff frequency
fs = 4000;
f_cut = 50 / (fs / 2);
[b, a] = butter(6, f_cut, 'low');
% Apply zero-phase filtering
meas_signal_f = filtfilt(b, a, meas_signal_raw')';
tank_signal_f = filtfilt(b, a, tank_signal_raw')';

%% Define Time and Frequency Parameters
in_time = 1;
end_time = 3;
in_sample = in_time * fs + 1;
end_sample = end_time * fs;
time = linspace(in_time, end_time, length(in_sample:end_sample));

f_up = 60;
f_down = 0.5;

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
        sgtitle(['MEA ' num2str(mea_idx) ' - ' mea_title],'FontWeight','bold'); % Title for the whole figure

        % Loop through the 2 electrodes in the group
        for j = 0:1
            electrode_idx = i + j;

            % Plot Frequency Spectrum for the current electrode 
            subplot(2, 2, j*2 + 1); 
            plot(fstep:fstep:f_up, Sffti(electrode_idx, :)); % Adjusted index for correct electrode
            xlabel('Frequency (Hz)');
            ylabel('Power');
            title(['Frequency Spectrum - Electrode ' num2str(plot_range(electrode_idx))]);

            % Plot Electrical Signal for the current electrode 
            subplot(2, 2, j*2 + 2); 
            plot(time, meas_signal_raw(mea_range(electrode_idx), in_sample:end_sample), 'r', 'LineWidth', 1);
            hold on;
            plot(time, meas_signal_f(mea_range(electrode_idx), in_sample:end_sample), 'k', 'LineWidth', 1);
            title(['Electrical Signal - Electrode ' num2str(plot_range(electrode_idx)) ' - Cutoff: ' num2str(f_cut * fs / 2) ' Hz']);
            xlabel('Time (s)');
            ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex', 'FontSize', 12);
            legend('Original', 'New filtering');
        end
        
        % Add a horizontal line to separate the rows 
        annotation('line', [0.1, 0.9], [0.5, 0.5], 'LineWidth', 1, 'Color', 'k');
        
        % Save each figure
%         saveas(gcf, ['MEA_' num2str(mea_idx) '_Electrodes_' num2str(plot_range(i)) '-' num2str(plot_range(i+1)) '_spectrum.png']);
    end
end


%% Plot tank spectrum and electrical signals

tank_range = 1:60; 
tank_labels = [129:174, 177:190]; % Labels for plots

tank_temp = tank_signal_f(tank_range, in_sample:end_sample); % Extract tank signals
[MFFTi, Sffti, fstep] = f_DF_electric(tank_temp, fs, f_up, f_down);

for i = 1:2:60  % Grouping 2 electrodes per plot
    figure();
    sgtitle('Tank Signals', 'FontWeight', 'bold');

    for j = 0:1
        electrode_idx = i + j;
        label_idx = tank_labels(electrode_idx);

        % Spectrum 
        subplot(2, 2, 1 + j*2);
        plot(fstep:fstep:f_up, Sffti(electrode_idx, :));
        xlabel('Frequency (Hz)');
        ylabel('Power');
        title(['Frequency Spectrum - Electrode ' num2str(label_idx)]);

        % Electrical Signals 
        subplot(2, 2, 2 + j*2);
        plot(time, tank_signal(electrode_idx, in_sample:end_sample), 'r', 'LineWidth', 1);
        hold on;
        plot(time, tank_signal_f(electrode_idx, in_sample:end_sample), 'k', 'LineWidth', 1);
        title(['Electrical Signal - Electrode ' num2str(label_idx) ' - Cutoff: ' num2str(f_cut * fs / 2) ' Hz)']);
        xlabel('Time (s)');
        ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex', 'FontSize', 12);
        legend('Original', 'New filtering');
    end
    
    % Add a horizontal line to separate rows
    annotation('line', [0.1, 0.9], [0.5, 0.5], 'Color', 'k', 'LineWidth', 1.5);
    
    % Save each figure
    saveas(gcf, ['TANK_Electrodes_' num2str(tank_labels(i)) '-' num2str(tank_labels(i+1)) '_spectrum.png']);

end


%% Save Filtered Signals
lp_cutoff = f_cut * fs / 2;
save('filtered_signals.mat', 'meas_signal_f', 'tank_signal_f', 'lp_cutoff');


%% Detrend Signals

% Define parameters
L_meas = size(meas_signal_f,2);
t_meas = 1/fs*(0:L_meas-1);
L_tank = size(tank_signal_f,2);
t_tank = 1/fs*(0:L_tank-1);

in_time = 1;
end_time = 3;
in_sample = in_time * fs + 1;
end_sample = end_time * fs;
time = linspace(in_time, end_time, length(in_sample:end_sample));

% Detrend MEA signals
meas_signal_d = zeros(size(meas_signal_f));
for m = 1:size(meas_signal_f,1)
    meas_signal_d(m,:) = detrendSpline(meas_signal_f(m,:), t_meas, 0.2);
end

% Detrend tank signals
tank_signal_d = zeros(size(tank_signal_f));
for m = 1:size(tank_signal_f,1)
    tank_signal_d(m,:) = detrendSpline(tank_signal_f(m,:), t_tank, 0.2);
end

%% Plot MEAs Signals before and after detrending
% Plot MEA Signals (RA, V, LA) before and after detrending

for mea_idx = 1:3
    % Set the indices for each MEA group
    if mea_idx == 1
        mea_title = 'RIGHT ATRIUM';
        mea_range = 1:16;
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
            plot(time, meas_signal_f(electrode_idx + (mea_idx-1)*16, in_sample:end_sample), 'color', "#00FFFF", 'LineWidth', 1.5);
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

tank_electrodes = [129:174, 177:190]; % Corresponding electrode indices

for i = 1:4:60
    figure(); % Create a new figure for each group of 4 tank electrodes
    sgtitle('Tank Signals', 'FontWeight', 'bold');
    
    % Loop through the 4 tank electrodes in the group
    for j = 0:3
        tank_idx = i + j; % Tank index in the 60 rows
        electrode_idx = tank_electrodes(tank_idx); % Map to the real electrode number
        
        subplot(2, 2, j+1); % Create 2x2 subplot grid
        
        % Plot Tank signals
        plot(time, tank_signal_f(tank_idx, in_sample:end_sample), 'color', "#ADD8E6");
        hold on; 
        plot(time, tank_signal_d(tank_idx, in_sample:end_sample), 'color', "#00008B"); 
        title(['Electrode ' num2str(electrode_idx)]);
        xlabel('Time (s)');
        ylabel('Amplitude ($\mu$V)');
        legend('Original', 'Detrended');
    end
    
    % Save each figure
    saveas(gcf, ['Tank_Electrodes_' num2str(tank_electrodes(i)) '-' num2str(tank_electrodes(i+3)) '_detrend.png']);
end

%% Downsampling

% Downsampling factor
factor = 20;

% MEA signals downsampling
meas_signal_down = resample(meas_signal_d', 1, factor)'; % Downsample MEA signals
fs_meas_down = fs / factor;  % Update MEA sampling frequency

% Tank signals downsampling
tank_signal_down = resample(tank_signal_d', 1, factor)'; % Downsample Tank signals
fs_tank_down = fs / factor;  % Update Tank sampling frequency

% Update the sampling frequency
fs_down = fs / factor;

%% Plot MEAs electrical signals before and after downsampling

in_time = 1;
end_time = 3;
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
            plot(t_original, meas_signal_d(mea_range(electrode_idx), in_time*fs+1:end_time*fs), 'r', 'LineWidth', 1);
            hold on;
            plot(t_down, meas_signal_down(mea_range(electrode_idx), in_time*fs_down+1:end_time*fs_down), 'b', 'LineWidth', 1);
            title(['Electrode ' num2str(plot_range(electrode_idx))]);
            xlabel('Time (s)');
            ylabel('Amplitude (\muV)');
            legend(['Original ' num2str(fs, '%.0f') ' Hz'], ['Downsampled ' num2str(fs_down, '%.0f') ' Hz']);

        end
         
        % Save each figure (optional)
        saveas(gcf, ['MEA_' num2str(mea_idx) '_Electrodes_' num2str(plot_range(i)) '-' num2str(plot_range(i+1)) '_downsampling.png']);
    end
end
%% Plot Tank Signals before and after downsampling

tank_electrodes = [129:174, 177:190];

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
        tank_idx = i + j;
        electrode_idx = tank_electrodes(tank_idx); % Map to the real electrode number
        
        % Plot Original and Downsampled Electrical Signals for the current electrode
        subplot(2, 2, j + 1);
        plot(t_original, tank_signal_d(tank_idx, in_time*fs+1:end_time*fs), 'r', 'LineWidth', 1);
        hold on;
        plot(t_down, tank_signal_down(tank_idx, in_time*fs_down+1:end_time*fs_down), 'b', 'LineWidth', 1);
        title(['Electrode ' num2str(electrode_idx)]);
        xlabel('Time (s)');
        ylabel('Amplitude (\muV)');
        legend(['Original ' num2str(fs, '%.0f') ' Hz'], ['Downsampled ' num2str(fs_down, '%.0f') ' Hz']);
        
    end
    
    % Save each figure
  saveas(gcf, ['Tank_Electrodes_' num2str(tank_electrodes(i)) '-' num2str(tank_electrodes(i+3)) '_downsampling.png']);

end

%% Save Downsampled Signals
save('downsampled_signals.mat', 'meas_signal_down', 'tank_signal_down','fs_down');

%% Normalizing signals

% Normalizing meas signal
meas_signal_norm = normalize(meas_signal_down,2,'zscore','std');

% Normalizing tank signals
tank_signal_norm = normalize(tank_signal_down,2,'zscore','std');

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
         
        % Save each figure
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
    sgtitle('Tank Signals', 'FontWeight', 'bold');
    
    % Loop through the 4 tank electrodes
    for j = 0:3
        tank_idx = i + j;
        electrode_idx = tank_electrodes(tank_idx); % Map to the real electrode number
        
        % Plot Original and Downsampled Electrical Signals for the current electrode
        subplot(2, 2, j + 1);
        plot(time, tank_signal_down(tank_idx, in_sample:end_sample), 'r', 'LineWidth', 1);
        hold on;
        plot(time, tank_signal_norm(tank_idx, in_sample:end_sample), 'b', 'LineWidth', 1);
        title(['Electrode ' num2str(electrode_idx)]);
        xlabel('Time (s)');
        ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex');
        legend('Original', 'Normalized');
        

%       Save each figure
%       saveas(gcf, ['Tank_Electrodes_' num2str(tank_electrodes(i)) '-' num2str(tank_electrodes(i+3)) '_norm.png']);

    end
end

%% Scaling signals

% Scaling meas signal
meas_signal_scaled = rescale(meas_signal_norm, -1, 1);

% Scaling tank signals
tank_signal_scaled = rescale(tank_signal_norm, -1, 1);

%% Plotting MEA signals before and after normalization and scaling

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
            legend('Original', 'Scaled');
        end
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
        tank_idx = i + j;
        electrode_idx = tank_electrodes(tank_idx); % Map to the real electrode number
        
        % Plot Original, Normalized, and Scaled Electrical Signals
        subplot(2, 2, j + 1);
        plot(time, tank_signal_norm(tank_idx, in_sample:end_sample), 'r', 'LineWidth', 1);
        hold on;
        plot(time, tank_signal_scaled(tank_idx, in_sample:end_sample), 'g', 'LineWidth', 1);
        title(['Electrode ' num2str(electrode_idx)]);
        xlabel('Time (s)');
        ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex');
        legend('Normalized', 'Scaled', 'Location', 'best');
    end
end

