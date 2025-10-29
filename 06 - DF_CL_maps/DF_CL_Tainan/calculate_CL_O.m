function CL_O = calculate_CL_O(Data, Fsampling, debug)
    % Function to calculate cycle lengths (CL) from a given 3D array
    %
    % Inputs:
    % - Data: 3D array with dimensions [x, y, signal]
    % - Fsampling: Sampling frequency (in Hz)
    % - debug: Debug level (0 = no plots, 1 = fewer plots)
    %
    % Output:
    % - CL_O: 3D array with dimensions [x, y, 1] containing CL values
    
    % Conversion factor to milliseconds
    conversion_factor = 1000 / Fsampling;

    % Threshold and MinPeakDistance
    threshold = 0.7;
    minDist = 500;    

    % Initialize a 2D array to store CL values
    [num_x, num_y, ~] = size(Data);
    CL_O = nan(num_x, num_y); % Preallocate with NaNs
    
    for x_idx = 1:num_x
        for y_idx = 1:num_y
            s = squeeze(Data(x_idx, y_idx, :)); % Extract the signal for the given electrode
            
            % Determine if peaks are more pronounced in the positive or negative direction
            max_val = max(s);
            min_val = min(s);
            
            if abs(max_val) >= abs(min_val)
                % Detect positive peaks
                [y, x] = findpeaks(s, 'MinPeakHeight', abs(max_val) * threshold, 'MinPeakDistance', minDist);
            else
                % Detect negative peaks by inverting the signal
                [y, x] = findpeaks(-s, 'MinPeakHeight', abs(min_val) * threshold, 'MinPeakDistance', minDist);
                y = -y; % Invert the peaks back to their original negative values
            end
            
            % Plotting based on debug level
            if debug == 1 && mod(x_idx, 10) == 1 && y_idx == round(num_y / 2)
                % Plot fewer figures, one for every 10 lines
                figure();
                plot(s);
                hold on;
                scatter(x, y, 'r', 'filled');
                hold off;
                title(['Electrode (', num2str(x_idx), ', ', num2str(y_idx), ')']);
                xlabel('Time (samples)');
                ylabel('Signal');
                legend('Signal', 'Detected Peaks');
            end
            
            % Calculate CL if peaks are found
            if ~isempty(x) && length(x) > 1
                cl = x(2) - x(1);
                CL_O(x_idx, y_idx) = cl * conversion_factor; % Transforming to ms
            end
        end
    end
    
    % Reshape CL_O to include a singleton dimension for consistency
    CL_O = reshape(CL_O, [num_x, num_y, 1]);
end
