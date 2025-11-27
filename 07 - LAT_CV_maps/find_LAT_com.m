function LAT_E = find_LAT_com(Data_segment, Fs, x_coord, y_coord, debug, case_name)
    % FIND_LAT_COM Calculate Local Activation Time (LAT) for a single
    % signal segment using the Center of Mass method.
    %
    % Inputs:
    %   Data_segment: 1D vector of the electrical data segment (e.g., lim1:lim2).
    %   Fs:           Sampling frequency of the data.
    %   x_coord:      X-coordinate (row index) of the current pixel.
    %   y_coord:      Y-coordinate (column index) of the current pixel.
    %   debug:        Debug mode flag (0 or 1).
    %   case_name:    String identifying the current data case (for plotting).
    %
    % Output:
    %   LAT_E:        Local Activation Time (LAT) in milliseconds, relative to the start of the segment.
    
    cutoff_freq = 200;  % Low-pass filter cutoff frequency (Hz)
    
    % ================= FILTERING =================
    nyquist = Fs / 2;
    % Check if cutoff frequency is valid
    if cutoff_freq >= nyquist
        warning('Cutoff frequency (%d Hz) is >= Nyquist frequency (%d Hz). Skipping filter.', cutoff_freq, nyquist);
        filtered_signal = Data_segment;
    else
        [b, a] = butter(4, cutoff_freq / nyquist, 'low');
        filtered_signal = filtfilt(b, a, Data_segment);
    end
    
    % ================= LAT DETECTION (Center of Mass) =================
    % The interval is the entire segment passed to the function
    interval_signal = abs(filtered_signal);
    
    % Avoid division by zero and ensure there's a signal
    if sum(interval_signal) == 0
        LAT_E = 0;
        return;
    end
    
    cumulative_sum = cumsum(interval_signal);
    half_sum = sum(interval_signal) / 2;
    
    % find returns the index relative to the start of interval_signal
    % The '-1' is because find returns the first index where the condition is met
    com_index = find(cumulative_sum >= half_sum, 1, 'first');
    
    % Position is the index within the *segment* (1 to length(Data_segment))
    position = com_index; 
    
    % Debugging plots
    if debug == 1
        % Determine the sample indices for plotting (relative to segment start)
        sample_indices = 1:length(Data_segment);
        
        figure(101);
        clf; 
        
        subplot(2, 1, 1);
        plot(sample_indices, filtered_signal, 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1);
        hold on;
        if ~isnan(position)
            plot(position, filtered_signal(position), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
        end
        title(sprintf('Filtered Signal (Segment) - Case: %s | Pixel: (%d,%d)', case_name, x_coord, y_coord));
        xlabel('Sample Index (Relative to Segment Start)'); ylabel('Amplitude');
        xlim([1 length(Data_segment)]);
        
        subplot(2, 1, 2);
        plot(sample_indices, interval_signal, 'k', 'LineWidth', 1.2); 
        hold on;
        if ~isnan(position)
            xline(position, 'r', 'LineWidth', 2);
            % Fill the area up to the COM position
            fill([1:(position-1), position-1, 1], [interval_signal(1:com_index-1)'; 0; 0], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        end
        title('Absolute Filtered Signal (LAT Detection)');
        xlabel('Sample Index (Relative to Segment Start)'); ylabel('Abs(Amplitude)');
        xlim([1 length(Data_segment)]);
        
        sgtitle('Center of Mass LAT Detection', 'FontSize', 12, 'FontWeight', 'bold');
        drawnow;
    end
    
    % Convert the position (index in the segment) to time in ms
    if ~isnan(position)
        LAT_E = position * 1000 / Fs;
    else
        LAT_E = NaN;
    end
end