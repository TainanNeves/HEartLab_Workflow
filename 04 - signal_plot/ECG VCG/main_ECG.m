%% MAIN ECG ANALYSIS

%% Cleaning
clear; clc;


%% Loading Data
load("F:\HEartLab\experiment_analyses\exp28_analysis\data_processed\data_filtered_sync_E28_F01_R09.mat");


%% Creating leads
data = D_SYNC.EL;
Fs = 4000;
% Limb leads (Derivacoes de membros)
D1 = data(129, :) - data(145, :);
D2 = data(138, :) - data(145, :);
D3 = data(138, :) - data(129, :);
aVR = -(D1 + D2)/2;
aVL = (D1 - D2)/2;
aVF = (D2 + D3)/2;
% Precordial leads
V1 = data(158, :);
V2 = data(159, :);
V3 = data(168, :);
V4 = data(171, :);
V5 = data(134, :);
V6 = data(135, :);


%% Saving ECG leads
filename = 'E28_F01_R09';
variables = {'aVF', 'aVL', 'aVR', 'D1', 'D2', 'D3', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'Fs'};
save([filename, ' - ECG leads.mat'], variables{:});
clear filename variables;


%% Plotting clinical sheet
start_sample = 1;
time_length = 2.5;     % normally 2.5s or 5s
N = time_length * Fs;  % Number of samples to display

% Figure
Time = (0:N-1) / Fs;
figure;
tiledlayout(4, 4, 'Padding', 'compact', 'TileSpacing', 'compact');
sgtitle('ECG', 'FontWeight', 'bold', 'FontSize', 16);
% Define the lead signals
leads = {D1, aVR, V1, V4, D2, aVL, V2, V5, D3, aVF, V3, V6};
lead_names = {'I', 'aVR', 'V1', ...
            'V4', 'II', 'aVL', ...
            'V2', 'V5', 'III', ...
            'aVF', 'V3', 'V6'};
% Plot each lead
for i = 1:12
    nexttile;
    plot(Time, leads{i}(start_sample:min(end, start_sample+N-1)), 'k');
    title(lead_names{i});
    xlabel('Time (s)');
    ylabel('uV');
    % ylim([-2 2]); % Adjust limits
    grid on;
end
% Rhythm strip (Usually lead II)
nexttile([1, 4]); % Span the last row
N_long = 4 * N; % Expected length (10s)
Time_long = (0:N_long-1) / Fs; % Time axis for 10s
% Extract signal or pad with zeros if needed
if length(D2) >= start_sample + N_long - 1
    D2_long = D2(start_sample:start_sample+N_long-1);
else
    D2_long = [D2(start_sample:end), zeros(1, N_long - length(D2(start_sample:end)))];
end
plot(Time_long, D2_long, 'k');
title('Rhythm Strip (II)');
xlabel('Time (s)');
ylabel('uV');
% ylim([-2 2]); % Keep uniform scale
grid on;


%% QRS detection + Alignment
ecg_signal = D2;  % Your ECG signal lead
t = (0:length(ecg_signal)-1) / Fs;

% Plot the ECG signal for visual inspection
figure;
plot(t, ecg_signal, 'b'); hold on;
title('ECG Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Detect QRS peaks
% Ask user for threshold
threshold = input('Enter QRS detection threshold: ');
[qrs_peaks, qrs_locs] = findpeaks(ecg_signal, 'MinPeakHeight', threshold, 'MinPeakDistance', round(0.4*Fs));

% Overlay detected peaks on the ECG signal
plot(t(qrs_locs), qrs_peaks, 'ro', 'MarkerFaceColor', 'r');
legend('ECG Signal', 'Detected QRS Peaks');
hold off;

% Ask user to select two points in the plot
disp('Click twice on the ECG plot to select the QRS window (start and end).');
[x_selected, ~] = ginput(2);
sample_selected = round(x_selected * Fs);

% Determine center QRS peak in selected range
qrs_in_range = qrs_locs(qrs_locs >= sample_selected(1) & qrs_locs <= sample_selected(2));
    if isempty(qrs_in_range)
        error('No QRS peak detected in the selected range! Try again.');
    end

% Select the middle QRS peak in the range
peak_loc = qrs_in_range(round(end/2));

% Calculate window sizes before and after the peak
window_before = peak_loc - sample_selected(1);
window_after = sample_selected(2) - peak_loc;

% Calculate the expected length of each QRS segment
expected_length = window_before + window_after + 1;

% Preallocate matrix to store QRS segments
num_peaks = length(qrs_locs);
qrs_matrix = NaN(num_peaks, expected_length); % Use expected_length here

for i = 1:num_peaks
    peak = qrs_locs(i);
    start_idx = max(1, peak - window_before);
    end_idx = min(length(ecg_signal), peak + window_after);

    % Extract QRS segment
    segment_temp = ecg_signal(start_idx:end_idx); % Use a temporary variable

    % Initialize segment with NaNs to ensure correct size and padding
    segment = NaN(1, expected_length);

    % Calculate the indices to place segment_temp into segment
    % This handles cases where the segment is truncated at the start/end of the ECG signal
    if (peak - window_before) < 1
        % If window goes before start of signal, pad beginning with NaNs
        fill_start_idx = 1 - (peak - window_before);
        segment(fill_start_idx : fill_start_idx + length(segment_temp) - 1) = segment_temp;
    elseif (peak + window_after) > length(ecg_signal)
        % If window goes after end of signal, pad end with NaNs
        segment(1 : length(segment_temp)) = segment_temp;
    else
        % Normal case, segment fits perfectly
        segment = segment_temp;
    end
    
    % Ensure segment is a row vector and of the expected length
    if size(segment, 1) > 1 % If it's a column vector, transpose it
        segment = segment';
    end
    
    % If for any reason the segment is still not the expected length,
    % this is a robust way to resize and pad with NaNs.
    % This handles cases where segment_temp might be shorter than expected_length
    % due to being at the very start or end of the overall ecg_signal.
    if length(segment) < expected_length
        temp_segment = NaN(1, expected_length);
        copy_start = max(1, 1 - (peak - window_before));
        copy_end = min(expected_length, copy_start + length(segment_temp) - 1);
        temp_segment(copy_start:copy_end) = segment_temp;
        segment = temp_segment;
    end

    qrs_matrix(i, :) = segment; % Assign the properly sized segment
end

% Plot the extracted QRS complexes
figure;
x_samples = -window_before:window_after;  % Sample indices relative to peak
plot(x_samples, qrs_matrix', 'k');
title('Extracted QRS Complexes (Aligned)');
xlabel('Samples relative to QRS peak');
ylabel('Amplitude');
grid on;


%% Plot aligned beats
% Compute the peak-to-peak amplitude of each QRS complex
qrs_amplitudes = max(qrs_matrix, [], 2) - min(qrs_matrix, [], 2);
% Sort the QRS complexes in ascending order based on amplitude (lowest to highest)
[~, sort_idx] = sort(qrs_amplitudes, 'ascend');  % 'ascend' ensures it's sorted from lowest to highest
% Reorder the QRS matrix and locations
qrs_matrix_sorted = qrs_matrix(sort_idx, :);
qrs_locs_sorted = qrs_locs(sort_idx);

% Number of detected beats
num_beats = size(qrs_matrix, 1);
window_size = size(qrs_matrix, 2);

% Create figure for 3D visualization
figure('Color', 'w', 'Position', [400 400 600 450]);
for i = 1:num_beats
    x = qrs_matrix_sorted(i, :);
    % Handle missing values (NaN) in QRS segments
    x(isnan(x)) = 0;
    n = length(x);
    % 3D Plot
    plot3(i * ones(n,1), (1:n)', x', 'b', 'LineWidth', 2);
    hold on;
    fill3(i * ones(n,1), (1:n)', vertcat(-0.5, x(2:n-1)', -0.5), 'w', 'LineStyle', ':');
end

xlim([0, num_beats]);
xlabel('Beats', 'FontWeight', 'Bold');
ylabel('Time Samples', 'FontWeight', 'Bold');
zlabel('Potential', 'FontWeight', 'Bold');
set(gca, 'LineWidth', 2, 'FontWeight', 'Bold', 'CameraPosition', [-118,5942,795], 'XDir', 'reverse');
grid on;


%% Save aligned qrs
filename = 'E28_F01_R09';
variables = {'Fs', 'qrs_locs', 'qrs_matrix'};
save([filename, ' - aligned peaks.mat'], variables{:});
clear filename variables;
































