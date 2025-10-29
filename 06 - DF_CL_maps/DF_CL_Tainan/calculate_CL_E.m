function CL_E_struct = calculate_CL_E(Data_temp, Fsampling, debug)
    % Function to calculate cycle lengths (CL_E) from given data
    %
    % Inputs:
    % - Data_temp: Cell array containing data for each case
    % - Fsampling: Sampling frequency (in Hz)
    % - debug: Debug level (0 = no plots, 1 = fewer plots, 3 = all plots)
    %
    % Output:
    % - CL_E_struct: Struct containing CL_E values for each case
    
    % Initialize default cases
    cases = {'MEA1', 'MEA2', 'MEA3', 'TANK'};

    % Initialize output struct
    CL_E_struct = struct();
    
    % Conversion factor to milliseconds
    conversion_factor = 1000 / Fsampling;

    % Threshold and MinPeakDistance
    threshold = 0.7;
    minDist = 500;    
    for case_idx = 1:length(cases)
        Data = Data_temp{case_idx};
        case_name = cases{case_idx};
        
        % Initialize a 2D array to store CL_E values
        [num_x, num_y, ~] = size(Data);
        CL_E_case = nan(num_x, num_y); % Preallocate with NaNs
        
        % Figure counter for subplot
        fig_counter = 1;
        
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
                if debug == 3
                    % Plot every signal
                    figure();
                    plot(s);
                    hold on;
                    scatter(x, y, 'r', 'filled');
                    hold off;
                    title(['Electrode (', num2str(x_idx), ', ', num2str(y_idx), ') - ', case_name]);
                    xlabel('Time (samples)');
                    ylabel('Signal');
                    legend('Signal', 'Detected Peaks');
                elseif debug == 1
                    % Plot fewer figures, 2 by line
                    if mod(fig_counter, 2) == 1
                        figure();
                    end
                    subplot(1, 2, mod(fig_counter - 1, 2) + 1);
                    plot(s);
                    hold on;
                    scatter(x, y, 'r', 'filled');
                    hold off;
                    title(['Electrode (', num2str(x_idx), ', ', num2str(y_idx), ') - ', case_name]);
                    xlabel('Time (samples)');
                    ylabel('Signal');
                    legend('Signal', 'Detected Peaks');
                    fig_counter = fig_counter + 1;
                end
                
                % Calculate CL_E if peaks are found
                if ~isempty(x) && length(x) > 1
                    cl = x(2) - x(1);
                    CL_E_case(x_idx, y_idx) = cl * conversion_factor; % Transforming to ms
                end
            end
        end
        
        % Store the results in the output struct
        CL_E_struct.(case_name) = CL_E_case;
    end
end
