function LAT_E = find_LAT_diff_Tainan(Data, Fsampling, WinStartIdx, WinEndIdx, case_name, x_coord, y_coord, debug)
    % FIND_LAT_CONTEXTUAL Calculates Local Activation Time (LAT) for a small 
    % analysis window within a larger contextual data array.
    %
    % Inputs:
    %   Data:        The large, surrounding array of electrical data.
    %   Fsampling:   Sampling frequency (Hz).
    %   WinStartIdx: Start index of the Analysis Window (relative to ContextData).
    %   WinEndIdx:   End index of the Analysis Window (relative to ContextData).
    %   case_name:   String identifying the current data case.
    %   x_coord:     X-coordinate (row index) of the current pixel/location.
    %   y_coord:     Y-coordinate (column index) of the current pixel/location.
    %   debug:       Debug mode flag (0, 1 or 2).
    %
    % Output:
    %   LAT_E:       Local Activation Time (LAT) in milliseconds.

    % Ensure data is a row vector and cast to double
    if size(Data, 1) > 1 && size(Data, 2) == 1
        Data = Data.'; 
    end
    Data = double(Data);
    ContextLength = length(Data);
    
    % Basic checks
    if max(abs(Data)) == 0 || WinStartIdx < 1 || WinEndIdx > ContextLength
        LAT_E = 0; 
        return;
    end
    
    % --- Stage A: Filtering and Differentiation on ContextData (Large Window) ---
    
    % 1. Zero-Phase Butterworth Filter (Standard for EGM signals)
    f_low = 1;      % High-pass cutoff (1 Hz, removes baseline drift)
    f_high = 40;   % Low-pass cutoff (50 Hz, removes high-freq noise)
    order = 3;      
    Wn = [f_low f_high] / (Fsampling/2); 

    [b, a] = butter(order, Wn, 'bandpass');
    xs_filtered = filtfilt(b, a, Data); 
    
    % 2. Savitzky-Golay Derivative (Robust parameters for large window)
    windowLength = 15; 
    polyOrder = 3; 
    
    % FIX: Ensure windowLength is odd and less than ContextLength
    if mod(windowLength, 2) == 0
        windowLength = windowLength - 1; % Make it odd
    end
    
    if ContextLength < windowLength
         % Use simple diff as a fallback if ContextData is too short
         dvdt_context = [0, diff(xs_filtered)];
         disp(['WARNING: Context window too short (', num2str(ContextLength), '). Falling back to diff.']);
    else
         try
             % ENHANCED FIX: Use sgolay to create filter coefficients first
             % This approach gives better error handling
             if polyOrder >= windowLength
                 polyOrder = windowLength - 1;
                 disp(['WARNING: Polynomial order reduced to ', num2str(polyOrder), ...
                       ' (must be < window length)']);
             end
             
             % Create Savitzky-Golay filter coefficients
             [b, g] = sgolay(polyOrder, windowLength);
             
             % Calculate the derivative using the appropriate column
             % The (:,2) column corresponds to the first derivative
             halfWin = (windowLength - 1)/2;
             
             % Preallocate
             dvdt_context = zeros(size(xs_filtered));
             
             % Apply the filter to the interior points
             for n = (windowLength+1)/2 : ContextLength - (windowLength-1)/2
                 % Dot product of the filter coefficients with the data window
                 dvdt_context(n) = g(:,2)' * xs_filtered(n - halfWin : n + halfWin)';
             end
             
             % Handle edges with a simpler approach
             % For the first few points, use forward difference
             for n = 1:halfWin
                 dvdt_context(n) = xs_filtered(n+1) - xs_filtered(n);
             end
             
             % For the last few points, use backward difference
             for n = ContextLength-halfWin+1:ContextLength
                 dvdt_context(n) = xs_filtered(n) - xs_filtered(n-1);
             end
             
         catch ME
             % Fallback if sgolay fails
             disp(['WARNING: Savitzky-Golay filtering failed: ', ME.message]);
             disp(['Using gradient instead. ContextLength: ', num2str(ContextLength), ...
                   ', windowLength: ', num2str(windowLength)]);
             
             % Use centered difference for better results than simple diff
             dvdt_context = gradient(xs_filtered);
         end
    end
    
    % --- Stage B: Improved Peak Detection within Analysis Window ---
    
    % 3. Isolate the Derivative Trace for the Analysis Window
    dvdt_analysis = dvdt_context(WinStartIdx : WinEndIdx);
    analysisLength = length(dvdt_analysis);
    
    % 4. Find peaks using multiple strategies to avoid edge artifacts
    
    % Strategy 1: Simple global max/min (original approach)
    [dvdt_max_global, pos_max_global] = max(dvdt_analysis);
    [dvdt_min_global, pos_min_global] = min(dvdt_analysis);
    
    % Strategy 2: Find local maxima/minima (more robust)
    % Smooth the derivative slightly to reduce noise
    dvdt_smoothed = movmean(dvdt_analysis, 5);
    
    % Find local maxima (peaks)
    [pks_max, locs_max] = findpeaks(dvdt_smoothed, 'MinPeakProminence', 0.1*std(dvdt_analysis));
    
    % Find local minima (valleys)
    [pks_min, locs_min] = findpeaks(-dvdt_smoothed, 'MinPeakProminence', 0.1*std(dvdt_analysis));
    pks_min = -pks_min; % Convert back to original sign
    
    % Strategy 3: Find peaks in absolute derivative
    dvdt_abs = abs(dvdt_analysis);
    dvdt_abs_smoothed = movmean(dvdt_abs, 5);
    [pks_abs, locs_abs] = findpeaks(dvdt_abs_smoothed, 'MinPeakProminence', 0.1*std(dvdt_abs));
    
    % 5. Apply edge exclusion rules
    edgeExclusion = round(0.1 * analysisLength); % Exclude 10% from edges
    validIndices = edgeExclusion+1 : analysisLength-edgeExclusion;
    
    % Check if global max/min are in edges
    globalMaxInEdge = pos_max_global <= edgeExclusion || pos_max_global >= (analysisLength - edgeExclusion);
    globalMinInEdge = pos_min_global <= edgeExclusion || pos_min_global >= (analysisLength - edgeExclusion);
    
    % Filter peaks to exclude edge regions
    validLocsMax = locs_max(ismember(locs_max, validIndices));
    validPksMax = pks_max(ismember(locs_max, validIndices));
    
    validLocsMin = locs_min(ismember(locs_min, validIndices));
    validPksMin = pks_min(ismember(locs_min, validIndices));
    
    validLocsAbs = locs_abs(ismember(locs_abs, validIndices));
    validPksAbs = pks_abs(ismember(locs_abs, validIndices));
    
    % 6. Decision logic for selecting the best peak
    if ~isempty(validLocsAbs)
        % Prefer peaks in absolute derivative (most robust for LAT detection)
        [~, bestAbsIdx] = max(validPksAbs);
        pos_relative_to_analysis = validLocsAbs(bestAbsIdx);
    elseif ~isempty(validLocsMax) || ~isempty(validLocsMin)
        % Choose between max and min based on which has larger absolute value
        if ~isempty(validLocsMax) && ~isempty(validLocsMin)
            [~, bestMaxIdx] = max(validPksMax);
            [~, bestMinIdx] = max(abs(validPksMin));
            
            bestMaxVal = validPksMax(bestMaxIdx);
            bestMinVal = validPksMin(bestMinIdx);
            
            if abs(bestMaxVal) > abs(bestMinVal)
                pos_relative_to_analysis = validLocsMax(bestMaxIdx);
            else
                pos_relative_to_analysis = validLocsMin(bestMinIdx);
            end
        elseif ~isempty(validLocsMax)
            [~, bestIdx] = max(validPksMax);
            pos_relative_to_analysis = validLocsMax(bestIdx);
        else
            [~, bestIdx] = max(abs(validPksMin));
            pos_relative_to_analysis = validLocsMin(bestIdx);
        end
    else
        % Fallback: If no valid peaks found, use global max/min but check edges
        if globalMaxInEdge && globalMinInEdge
            % Both global extrema are in edges, this is suspicious
            % Try to find any peak in the middle 80%
            middleStart = round(0.4 * analysisLength);
            middleEnd = round(0.6 * analysisLength);
            middleIndices = middleStart:middleEnd;
            
            if ~isempty(middleIndices)
                [~, posInMiddle] = max(abs(dvdt_analysis(middleIndices)));
                pos_relative_to_analysis = middleIndices(1) + posInMiddle - 1;
            else
                % Last resort: use global max/min
                if abs(dvdt_max_global) > abs(dvdt_min_global)
                    pos_relative_to_analysis = pos_max_global;
                else
                    pos_relative_to_analysis = pos_min_global;
                end
            end
        elseif globalMaxInEdge
            % Only max is in edge, use min
            pos_relative_to_analysis = pos_min_global;
        elseif globalMinInEdge
            % Only min is in edge, use max
            pos_relative_to_analysis = pos_max_global;
        else
            % Neither is in edge, use whichever has larger absolute value
            if abs(dvdt_max_global) > abs(dvdt_min_global)
                pos_relative_to_analysis = pos_max_global;
            else
                pos_relative_to_analysis = pos_min_global;
            end
        end
    end
    
    % 7. Map the index back to the START of the Context Window (N-length index)
    % This is the final index relative to the start of the original ContextData array.
    position_final = (WinStartIdx - 1) + pos_relative_to_analysis;

    % --- Debug Modes and Plotting ---
    plotTitle = ['Case: ', case_name, ' | Pixel: (', num2str(x_coord), ',', num2str(y_coord), ')'];
    
    if debug >= 1
        figure(101);
        clf; 
        
        % Subplot 1: Full Signal Trace (Context Length)
        subplot(3,1,1);
        plot(xs_filtered);
        hold on;
        % Plot the detected point
        plot(position_final, xs_filtered(position_final), 'o', 'markersize', 10, 'color', 'red');
        % Highlight the Analysis Window region (where the search occurred)
        xline(WinStartIdx, 'k--');
        xline(WinEndIdx, 'k--');
        % Highlight edge exclusion zones
        plot([WinStartIdx+edgeExclusion, WinStartIdx+edgeExclusion], ylim, 'g--', 'LineWidth', 0.5);
        plot([WinEndIdx-edgeExclusion, WinEndIdx-edgeExclusion], ylim, 'g--', 'LineWidth', 0.5);
        % Zoom lines to show window region
        plot([WinStartIdx WinEndIdx], [xs_filtered(position_final) xs_filtered(position_final)], 'k:', 'LineWidth', 0.5);
        title('Full Signal');
        xlabel('Sample Index');
        ylabel('Amplitude');
        hold off;
        
        % Subplot 2: Signal in Analysis Window Only
        subplot(3,1,2);
        window_indices = WinStartIdx:WinEndIdx;
        window_signal = xs_filtered(window_indices);
        plot(window_indices, window_signal);
        hold on;
        % Plot detected point within window
        plot(position_final, xs_filtered(position_final), 'o', 'markersize', 10, 'color', 'red');
        % Plot edge exclusion zones
        plot([WinStartIdx+edgeExclusion, WinStartIdx+edgeExclusion], ylim, 'g--', 'LineWidth', 0.5);
        plot([WinEndIdx-edgeExclusion, WinEndIdx-edgeExclusion], ylim, 'g--', 'LineWidth', 0.5);
        title('Signal in Analysis Window');
        xlabel('Sample Index');
        ylabel('Amplitude');
        grid on;
        hold off;
        
        % Subplot 3: Derivative in Analysis Window Only
        subplot(3,1,3);
        plot(window_indices, dvdt_analysis);
        hold on;
        % Plot the derivative peak
        plot(position_final, dvdt_context(position_final), 'ro', 'markersize', 10);
        xline(position_final, 'r', 'LineWidth', 2);
        
        % Plot all detected peaks
        if ~isempty(locs_max)
            plot(window_indices(locs_max), dvdt_analysis(locs_max), 'g^', 'MarkerSize', 8);
        end
        if ~isempty(locs_min)
            plot(window_indices(locs_min), dvdt_analysis(locs_min), 'gv', 'MarkerSize', 8);
        end
        if ~isempty(locs_abs)
            plot(window_indices(locs_abs), dvdt_analysis(locs_abs), 'm*', 'MarkerSize', 8);
        end
        
        % Plot edge exclusion zones
        plot([WinStartIdx+edgeExclusion, WinStartIdx+edgeExclusion], ylim, 'g--', 'LineWidth', 0.5);
        plot([WinEndIdx-edgeExclusion, WinEndIdx-edgeExclusion], ylim, 'g--', 'LineWidth', 0.5);
        
        title('dV/dt in Analysis Window (Green: peaks, Magenta: abs peaks)');
        xlabel('Sample Index');
        ylabel('dV/dt');
        grid on;
        hold off;
        
        sgtitle(plotTitle, 'FontSize', 12, 'FontWeight', 'bold');
        drawnow; 
    end
    
    if debug == 2
        % --- Interactive Correction (Adapted for Contextual Window) ---
        figure(101);
        clf; 
        
        % Subplot 1: Full Signal Trace (Context Length)
        subplot(3,1,1);
        plot(xs_filtered);
        hold on;
        h_auto = plot(position_final, xs_filtered(position_final), 'o', 'markersize', 10, 'color', 'red');
        xline(WinStartIdx, 'k--');
        xline(WinEndIdx, 'k--');
        % Highlight edge exclusion zones
        plot([WinStartIdx+edgeExclusion, WinStartIdx+edgeExclusion], ylim, 'g--', 'LineWidth', 0.5);
        plot([WinEndIdx-edgeExclusion, WinEndIdx-edgeExclusion], ylim, 'g--', 'LineWidth', 0.5);
        title('1. Full Signal Trace (Press SPACE or Click)');
        xlabel('Sample Index');
        ylabel('Amplitude');
        
        % Subplot 2: Signal in Analysis Window Only
        subplot(3,1,2);
        window_indices = WinStartIdx:WinEndIdx;
        window_signal = xs_filtered(window_indices);
        plot(window_indices, window_signal);
        hold on;
        h_auto_win = plot(position_final, xs_filtered(position_final), 'o', 'markersize', 10, 'color', 'red');
        % Plot edge exclusion zones
        plot([WinStartIdx+edgeExclusion, WinStartIdx+edgeExclusion], ylim, 'g--', 'LineWidth', 0.5);
        plot([WinEndIdx-edgeExclusion, WinEndIdx-edgeExclusion], ylim, 'g--', 'LineWidth', 0.5);
        title('2. Signal in Analysis Window');
        xlabel('Sample Index');
        ylabel('Amplitude');
        grid on;
        
        % Subplot 3: Derivative in Analysis Window Only
        subplot(3,1,3);
        plot(window_indices, dvdt_analysis);
        hold on;
        h_auto_deriv = plot(position_final, dvdt_context(position_final), 'ro', 'markersize', 10);
        h_line = xline(position_final, 'r', 'LineWidth', 2);
        
        % Plot all detected peaks
        if ~isempty(locs_max)
            plot(window_indices(locs_max), dvdt_analysis(locs_max), 'g^', 'MarkerSize', 8);
        end
        if ~isempty(locs_min)
            plot(window_indices(locs_min), dvdt_analysis(locs_min), 'gv', 'MarkerSize', 8);
        end
        if ~isempty(locs_abs)
            plot(window_indices(locs_abs), dvdt_analysis(locs_abs), 'm*', 'MarkerSize', 8);
        end
        
        % Plot edge exclusion zones
        plot([WinStartIdx+edgeExclusion, WinStartIdx+edgeExclusion], ylim, 'g--', 'LineWidth', 0.5);
        plot([WinEndIdx-edgeExclusion, WinEndIdx-edgeExclusion], ylim, 'g--', 'LineWidth', 0.5);
        
        title('3. dV/dt in Analysis Window');
        xlabel('Sample Index');
        ylabel('dV/dt');
        grid on;
        
        sgtitle({plotTitle, '--- Press SPACE to accept, CLICK to select new point ---'}, ...
                'FontSize', 12, 'FontWeight', 'bold', 'Color', 'blue');
        
        drawnow;
        
        k = waitforbuttonpress;
        
        if k == 0 % Mouse click occurred
            % Get click from any of the subplots
            subplot(3,1,1);
            [x_click, ~] = ginput(1);
            
            % Enforce that the click must be within the Analysis Window bounds
            new_position = round(x_click);
            position_final = max(WinStartIdx, min(WinEndIdx, new_position));
            
            % Update plots
            % Subplot 1
            delete(h_auto);
            subplot(3,1,1);
            plot(position_final, xs_filtered(position_final), 's', 'markersize', 10, 'color', 'green'); 
            
            % Subplot 2
            delete(h_auto_win);
            subplot(3,1,2);
            plot(position_final, xs_filtered(position_final), 's', 'markersize', 10, 'color', 'green');
            
            % Subplot 3
            delete(h_auto_deriv);
            delete(h_line);
            subplot(3,1,3);
            plot(position_final, dvdt_context(position_final), 'gs', 'markersize', 10);
            xline(position_final, 'g', 'LineWidth', 2);
            
            disp(['Manually selected position: ', num2str(position_final)]);
        else % k == 1 (Key press, including SPACE)
            disp('Key pressed. Using automatic point.');
        end
        
        hold off; 
    end
    
    % --- Final LAT Calculation ---
    % position_final is relative to the start of ContextData.
    LAT_E = position_final * (1000 / Fsampling); 
end






% function LAT_E = find_LAT_diff_Tainan(Data, Fsampling, WinStartIdx, WinEndIdx, case_name, x_coord, y_coord, debug)
%     % FIND_LAT_CONTEXTUAL Calculates Local Activation Time (LAT) for a small 
%     % analysis window within a larger contextual data array.
%     %
%     % Inputs:
%     %   Data:        The large, surrounding array of electrical data.
%     %   Fsampling:   Sampling frequency (Hz).
%     %   WinStartIdx: Start index of the Analysis Window (relative to ContextData).
%     %   WinEndIdx:   End index of the Analysis Window (relative to ContextData).
%     %   case_name:   String identifying the current data case.
%     %   x_coord:     X-coordinate (row index) of the current pixel/location.
%     %   y_coord:     Y-coordinate (column index) of the current pixel/location.
%     %   debug:       Debug mode flag (0, 1 or 2).
%     %
%     % Output:
%     %   LAT_E:       Local Activation Time (LAT) in milliseconds.
% 
%     % Ensure data is a row vector and cast to double
%     if size(Data, 1) > 1 && size(Data, 2) == 1
%         Data = Data.'; 
%     end
%     Data = double(Data);
%     ContextLength = length(Data);
% 
%     % Basic checks
%     if max(abs(Data)) == 0 || WinStartIdx < 1 || WinEndIdx > ContextLength
%         LAT_E = 0; 
%         return;
%     end
% 
%     % --- Stage A: Filtering and Differentiation on ContextData (Large Window) ---
% 
%     % 1. Zero-Phase Butterworth Filter (Standard for EGM signals)
%     f_low = 1;      % High-pass cutoff (1 Hz, removes baseline drift)
%     f_high = 40;   % Low-pass cutoff (50 Hz, removes high-freq noise)
%     order = 3;      
%     Wn = [f_low f_high] / (Fsampling/2); 
% 
%     [b, a] = butter(order, Wn, 'bandpass');
%     xs_filtered = filtfilt(b, a, Data); 
% 
%     % 2. Savitzky-Golay Derivative (Robust parameters for large window)
%     windowLength = 15; 
%     polyOrder = 3; 
% 
%     % FIX: Ensure windowLength is odd and less than ContextLength
%     if mod(windowLength, 2) == 0
%         windowLength = windowLength - 1; % Make it odd
%     end
% 
%     if ContextLength < windowLength
%          % Use simple diff as a fallback if ContextData is too short
%          dvdt_context = [0, diff(xs_filtered)];
%          disp(['WARNING: Context window too short (', num2str(ContextLength), '). Falling back to diff.']);
%     else
%          try
%              % ENHANCED FIX: Use sgolay to create filter coefficients first
%              % This approach gives better error handling
%              if polyOrder >= windowLength
%                  polyOrder = windowLength - 1;
%                  disp(['WARNING: Polynomial order reduced to ', num2str(polyOrder), ...
%                        ' (must be < window length)']);
%              end
% 
%              % Create Savitzky-Golay filter coefficients
%              [b, g] = sgolay(polyOrder, windowLength);
% 
%              % Calculate the derivative using the appropriate column
%              % The (:,2) column corresponds to the first derivative
%              halfWin = (windowLength - 1)/2;
% 
%              % Preallocate
%              dvdt_context = zeros(size(xs_filtered));
% 
%              % Apply the filter to the interior points
%              for n = (windowLength+1)/2 : ContextLength - (windowLength-1)/2
%                  % Dot product of the filter coefficients with the data window
%                  dvdt_context(n) = g(:,2)' * xs_filtered(n - halfWin : n + halfWin)';
%              end
% 
%              % Handle edges with a simpler approach
%              % For the first few points, use forward difference
%              for n = 1:halfWin
%                  dvdt_context(n) = xs_filtered(n+1) - xs_filtered(n);
%              end
% 
%              % For the last few points, use backward difference
%              for n = ContextLength-halfWin+1:ContextLength
%                  dvdt_context(n) = xs_filtered(n) - xs_filtered(n-1);
%              end
% 
%          catch ME
%              % Fallback if sgolay fails
%              disp(['WARNING: Savitzky-Golay filtering failed: ', ME.message]);
%              disp(['Using gradient instead. ContextLength: ', num2str(ContextLength), ...
%                    ', windowLength: ', num2str(windowLength)]);
% 
%              % Use centered difference for better results than simple diff
%              dvdt_context = gradient(xs_filtered);
%          end
%     end
% 
%     % --- Stage B: Peak Detection within Analysis Window ---
% 
%     % 3. Isolate the Derivative Trace for the Analysis Window
%     dvdt_analysis = dvdt_context(WinStartIdx : WinEndIdx);
% 
%     % 4. Find the max absolute change within the isolated segment
%     [dvdt_max, pos_max_analysis] = max(dvdt_analysis);
%     [dvdt_min, pos_min_analysis] = min(dvdt_analysis);
% 
%     if abs(dvdt_max) > abs(dvdt_min)
%         pos_relative_to_analysis = pos_max_analysis;
%     else
%         pos_relative_to_analysis = pos_min_analysis;
%     end
% 
%     % 5. Map the index back to the START of the Context Window (N-length index)
%     % This is the final index relative to the start of the original ContextData array.
%     position_final = (WinStartIdx - 1) + pos_relative_to_analysis;
% 
%     % --- Debug Modes and Plotting ---
%     plotTitle = ['Case: ', case_name, ' | Pixel: (', num2str(x_coord), ',', num2str(y_coord), ')'];
% 
%     if debug >= 1
%         figure(101);
%         clf; 
% 
%         % Subplot 1: Full Signal Trace (Context Length)
%         subplot(3,1,1);
%         plot(xs_filtered);
%         hold on;
%         % Plot the detected point
%         plot(position_final, xs_filtered(position_final), 'o', 'markersize', 10, 'color', 'red');
%         % Highlight the Analysis Window region (where the search occurred)
%         xline(WinStartIdx, 'k--');
%         xline(WinEndIdx, 'k--');
%         % Zoom lines to show window region
%         plot([WinStartIdx WinEndIdx], [xs_filtered(position_final) xs_filtered(position_final)], 'k:', 'LineWidth', 0.5);
%         title('Full Signal');
%         xlabel('Sample Index');
%         ylabel('Amplitude');
%         hold off;
% 
%         % Subplot 2: Signal in Analysis Window Only
%         subplot(3,1,2);
%         window_indices = WinStartIdx:WinEndIdx;
%         window_signal = xs_filtered(window_indices);
%         plot(window_indices, window_signal);
%         hold on;
%         % Plot detected point within window
%         plot(position_final, xs_filtered(position_final), 'o', 'markersize', 10, 'color', 'red');
%         title('Signal in Analysis Window');
%         xlabel('Sample Index');
%         ylabel('Amplitude');
%         grid on;
%         hold off;
% 
%         % Subplot 3: Derivative in Analysis Window Only
%         subplot(3,1,3);
%         plot(window_indices, dvdt_analysis);
%         hold on;
%         % Plot the derivative peak
%         plot(position_final, dvdt_context(position_final), 'ro', 'markersize', 10);
%         xline(position_final, 'r', 'LineWidth', 2);
%         title('dV/dt in Analysis Window (Savitzky-Golay)');
%         xlabel('Sample Index');
%         ylabel('dV/dt');
%         grid on;
%         hold off;
% 
%         sgtitle(plotTitle, 'FontSize', 12, 'FontWeight', 'bold');
%         drawnow; 
%     end
% 
%     if debug == 2
%         % --- Interactive Correction (Adapted for Contextual Window) ---
%         figure(101);
%         clf; 
% 
%         % Subplot 1: Full Signal Trace (Context Length)
%         subplot(3,1,1);
%         plot(xs_filtered);
%         hold on;
%         h_auto = plot(position_final, xs_filtered(position_final), 'o', 'markersize', 10, 'color', 'red');
%         xline(WinStartIdx, 'k--');
%         xline(WinEndIdx, 'k--');
%         title('1. Full Signal Trace (Press SPACE or Click)');
%         xlabel('Sample Index');
%         ylabel('Amplitude');
% 
%         % Subplot 2: Signal in Analysis Window Only
%         subplot(3,1,2);
%         window_indices = WinStartIdx:WinEndIdx;
%         window_signal = xs_filtered(window_indices);
%         plot(window_indices, window_signal);
%         hold on;
%         h_auto_win = plot(position_final, xs_filtered(position_final), 'o', 'markersize', 10, 'color', 'red');
%         title('2. Signal in Analysis Window');
%         xlabel('Sample Index');
%         ylabel('Amplitude');
%         grid on;
% 
%         % Subplot 3: Derivative in Analysis Window Only
%         subplot(3,1,3);
%         plot(window_indices, dvdt_analysis);
%         hold on;
%         h_auto_deriv = plot(position_final, dvdt_context(position_final), 'ro', 'markersize', 10);
%         h_line = xline(position_final, 'r', 'LineWidth', 2);
%         title('3. dV/dt in Analysis Window');
%         xlabel('Sample Index');
%         ylabel('dV/dt');
%         grid on;
% 
%         sgtitle({plotTitle, '--- Press SPACE to accept, CLICK to select new point ---'}, ...
%                 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'blue');
% 
%         drawnow;
% 
%         k = waitforbuttonpress;
% 
%         if k == 0 % Mouse click occurred
%             % Get click from any of the subplots
%             subplot(3,1,1);
%             [x_click, ~] = ginput(1);
% 
%             % Enforce that the click must be within the Analysis Window bounds
%             new_position = round(x_click);
%             position_final = max(WinStartIdx, min(WinEndIdx, new_position));
% 
%             % Update plots
%             % Subplot 1
%             delete(h_auto);
%             subplot(3,1,1);
%             plot(position_final, xs_filtered(position_final), 's', 'markersize', 10, 'color', 'green'); 
% 
%             % Subplot 2
%             delete(h_auto_win);
%             subplot(3,1,2);
%             plot(position_final, xs_filtered(position_final), 's', 'markersize', 10, 'color', 'green');
% 
%             % Subplot 3
%             delete(h_auto_deriv);
%             delete(h_line);
%             subplot(3,1,3);
%             plot(position_final, dvdt_context(position_final), 'gs', 'markersize', 10);
%             xline(position_final, 'g', 'LineWidth', 2);
% 
%             disp(['Manually selected position: ', num2str(position_final)]);
%         else % k == 1 (Key press, including SPACE)
%             disp('Key pressed. Using automatic point.');
%         end
% 
%         hold off; 
%     end
% 
%     % --- Final LAT Calculation ---
%     % position_final is relative to the start of ContextData.
%     LAT_E = position_final * (1000 / Fsampling); 
% end
% 
% 
% 
% 
% % function LAT_E = find_LAT_diff_Tainan(Data, Fsampling, WinStartIdx, WinEndIdx, case_name, x_coord, y_coord, debug)
% %     % FIND_LAT_CONTEXTUAL Calculates Local Activation Time (LAT) for a small 
% %     % analysis window within a larger contextual data array.
% %     %
% %     % Inputs:
% %     %   Data:        The large, surrounding array of electrical data.
% %     %   Fsampling:   Sampling frequency (Hz).
% %     %   WinStartIdx: Start index of the Analysis Window (relative to ContextData).
% %     %   WinEndIdx:   End index of the Analysis Window (relative to ContextData).
% %     %   case_name:   String identifying the current data case.
% %     %   x_coord:     X-coordinate (row index) of the current pixel/location.
% %     %   y_coord:     Y-coordinate (column index) of the current pixel/location.
% %     %   debug:       Debug mode flag (0, 1 or 2).
% %     %
% %     % Output:
% %     %   LAT_E:       Local Activation Time (LAT) in milliseconds.
% % 
% %     % Ensure data is a row vector and cast to double
% %     if size(Data, 1) > 1 && size(Data, 2) == 1
% %         Data = Data.'; 
% %     end
% %     Data = double(Data);
% %     ContextLength = length(Data);
% % 
% %     % Basic checks
% %     if max(abs(Data)) == 0 || WinStartIdx < 1 || WinEndIdx > ContextLength
% %         LAT_E = 0; 
% %         return;
% %     end
% % 
% %     % --- Stage A: Filtering and Differentiation on ContextData (Large Window) ---
% % 
% %     % 1. Zero-Phase Butterworth Filter (Standard for EGM signals)
% %     f_low = 1;      % High-pass cutoff (1 Hz, removes baseline drift)
% %     f_high = 40;   % Low-pass cutoff (50 Hz, removes high-freq noise)
% %     order = 3;      
% %     Wn = [f_low f_high] / (Fsampling/2); 
% % 
% %     [b, a] = butter(order, Wn, 'bandpass');
% %     xs_filtered = filtfilt(b, a, Data); 
% % 
% %     % 2. Savitzky-Golay Derivative (Robust parameters for large window)
% %     windowLength = 15; 
% %     polyOrder = 3; 
% % 
% %     % FIX: Ensure windowLength is odd and less than ContextLength
% %     if mod(windowLength, 2) == 0
% %         windowLength = windowLength - 1; % Make it odd
% %     end
% % 
% %     if ContextLength < windowLength
% %          % Use simple diff as a fallback if ContextData is too short
% %          dvdt_context = [0, diff(xs_filtered)];
% %          disp(['WARNING: Context window too short (', num2str(ContextLength), '). Falling back to diff.']);
% %     else
% %          try
% %              % Ensure sgolayfilt parameters are valid
% %              dvdt_context = sgolayfilt(xs_filtered, polyOrder, windowLength, 1);
% %          catch ME
% %              % Fallback if sgolayfilt fails
% %              disp(['WARNING: sgolayfilt failed: ', ME.message]);
% %              disp(['Using simple diff instead. ContextLength: ', num2str(ContextLength), ...
% %                    ', windowLength: ', num2str(windowLength)]);
% %              dvdt_context = [0, diff(xs_filtered)];
% %          end
% %     end
% % 
% %     % --- Stage B: Peak Detection within Analysis Window ---
% % 
% %     % 3. Isolate the Derivative Trace for the Analysis Window
% %     dvdt_analysis = dvdt_context(WinStartIdx : WinEndIdx);
% % 
% %     % 4. Find the max absolute change within the isolated segment
% %     [dvdt_max, pos_max_analysis] = max(dvdt_analysis);
% %     [dvdt_min, pos_min_analysis] = min(dvdt_analysis);
% % 
% %     if abs(dvdt_max) > abs(dvdt_min)
% %         pos_relative_to_analysis = pos_max_analysis;
% %     else
% %         pos_relative_to_analysis = pos_min_analysis;
% %     end
% % 
% %     % 5. Map the index back to the START of the Context Window (N-length index)
% %     % This is the final index relative to the start of the original ContextData array.
% %     position_final = (WinStartIdx - 1) + pos_relative_to_analysis;
% % 
% %     % --- Debug Modes and Plotting ---
% %     plotTitle = ['Case: ', case_name, ' | Pixel: (', num2str(x_coord), ',', num2str(y_coord), ')'];
% % 
% %     if debug >= 1
% %         figure(101);
% %         clf; 
% % 
% %         % Subplot 1: Full Signal Trace (Context Length)
% %         subplot(3,1,1);
% %         plot(xs_filtered);
% %         hold on;
% %         % Plot the detected point
% %         plot(position_final, xs_filtered(position_final), 'o', 'markersize', 10, 'color', 'red');
% %         % Highlight the Analysis Window region (where the search occurred)
% %         xline(WinStartIdx, 'k--');
% %         xline(WinEndIdx, 'k--');
% %         % Zoom lines to show window region
% %         plot([WinStartIdx WinEndIdx], [xs_filtered(position_final) xs_filtered(position_final)], 'k:', 'LineWidth', 0.5);
% %         title('Full Signal');
% %         xlabel('Sample Index');
% %         ylabel('Amplitude');
% %         hold off;
% % 
% %         % Subplot 2: Signal in Analysis Window Only
% %         subplot(3,1,2);
% %         window_indices = WinStartIdx:WinEndIdx;
% %         window_signal = xs_filtered(window_indices);
% %         plot(window_indices, window_signal);
% %         hold on;
% %         % Plot detected point within window
% %         plot(position_final, xs_filtered(position_final), 'o', 'markersize', 10, 'color', 'red');
% %         title('Signal in Analysis Window');
% %         xlabel('Sample Index');
% %         ylabel('Amplitude');
% %         grid on;
% %         hold off;
% % 
% %         % Subplot 3: Derivative in Analysis Window Only
% %         subplot(3,1,3);
% %         plot(window_indices, dvdt_analysis);
% %         hold on;
% %         % Plot the derivative peak
% %         plot(position_final, dvdt_context(position_final), 'ro', 'markersize', 10);
% %         xline(position_final, 'r', 'LineWidth', 2);
% %         title('dV/dt in Analysis Window (Savitzky-Golay)');
% %         xlabel('Sample Index');
% %         ylabel('dV/dt');
% %         grid on;
% %         hold off;
% % 
% %         sgtitle(plotTitle, 'FontSize', 12, 'FontWeight', 'bold');
% %         drawnow; 
% %     end
% % 
% %     if debug == 2
% %         % --- Interactive Correction (Adapted for Contextual Window) ---
% %         figure(101);
% %         clf; 
% % 
% %         % Subplot 1: Full Signal Trace (Context Length)
% %         subplot(3,1,1);
% %         plot(xs_filtered);
% %         hold on;
% %         h_auto = plot(position_final, xs_filtered(position_final), 'o', 'markersize', 10, 'color', 'red');
% %         xline(WinStartIdx, 'k--');
% %         xline(WinEndIdx, 'k--');
% %         title('1. Full Signal Trace (Press SPACE or Click)');
% %         xlabel('Sample Index');
% %         ylabel('Amplitude');
% % 
% %         % Subplot 2: Signal in Analysis Window Only
% %         subplot(3,1,2);
% %         window_indices = WinStartIdx:WinEndIdx;
% %         window_signal = xs_filtered(window_indices);
% %         plot(window_indices, window_signal);
% %         hold on;
% %         h_auto_win = plot(position_final, xs_filtered(position_final), 'o', 'markersize', 10, 'color', 'red');
% %         title('2. Signal in Analysis Window');
% %         xlabel('Sample Index');
% %         ylabel('Amplitude');
% %         grid on;
% % 
% %         % Subplot 3: Derivative in Analysis Window Only
% %         subplot(3,1,3);
% %         plot(window_indices, dvdt_analysis);
% %         hold on;
% %         h_auto_deriv = plot(position_final, dvdt_context(position_final), 'ro', 'markersize', 10);
% %         h_line = xline(position_final, 'r', 'LineWidth', 2);
% %         title('3. dV/dt in Analysis Window');
% %         xlabel('Sample Index');
% %         ylabel('dV/dt');
% %         grid on;
% % 
% %         sgtitle({plotTitle, '--- Press SPACE to accept, CLICK to select new point ---'}, ...
% %                 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'blue');
% % 
% %         drawnow;
% % 
% %         k = waitforbuttonpress;
% % 
% %         if k == 0 % Mouse click occurred
% %             % Get click from any of the subplots
% %             subplot(3,1,1);
% %             [x_click, ~] = ginput(1);
% % 
% %             % Enforce that the click must be within the Analysis Window bounds
% %             new_position = round(x_click);
% %             position_final = max(WinStartIdx, min(WinEndIdx, new_position));
% % 
% %             % Update plots
% %             % Subplot 1
% %             delete(h_auto);
% %             subplot(3,1,1);
% %             plot(position_final, xs_filtered(position_final), 's', 'markersize', 10, 'color', 'green'); 
% % 
% %             % Subplot 2
% %             delete(h_auto_win);
% %             subplot(3,1,2);
% %             plot(position_final, xs_filtered(position_final), 's', 'markersize', 10, 'color', 'green');
% % 
% %             % Subplot 3
% %             delete(h_auto_deriv);
% %             delete(h_line);
% %             subplot(3,1,3);
% %             plot(position_final, dvdt_context(position_final), 'gs', 'markersize', 10);
% %             xline(position_final, 'g', 'LineWidth', 2);
% % 
% %             disp(['Manually selected position: ', num2str(position_final)]);
% %         else % k == 1 (Key press, including SPACE)
% %             disp('Key pressed. Using automatic point.');
% %         end
% % 
% %         hold off; 
% %     end
% % 
% %     % --- Final LAT Calculation ---
% %     % position_final is relative to the start of ContextData.
% %     LAT_E = position_final * (1000 / Fsampling); 
% % end