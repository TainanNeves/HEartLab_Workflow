function LAT_E = find_LAT_diff_Tainan(Data, Fsampling, WinStartIdx, WinEndIdx, case_name, x_coord, y_coord, debug)
    % FIND_LAT_CONTEXTUAL Calculates Local Activation Time (LAT) for a small 
    % analysis window within a larger contextual data array.
    %
    % Inputs:
    %   ContextData: The large, surrounding array of electrical data.
    %   Fs:          Sampling frequency (Hz).
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
    f_high = 50;   % Low-pass cutoff (50 Hz, removes high-freq noise)
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
             % Ensure sgolayfilt parameters are valid
             dvdt_context = sgolayfilt(xs_filtered, polyOrder, windowLength, 1);
         catch ME
             % Fallback if sgolayfilt fails
             disp(['WARNING: sgolayfilt failed: ', ME.message]);
             disp(['Using simple diff instead. ContextLength: ', num2str(ContextLength), ...
                   ', windowLength: ', num2str(windowLength)]);
             dvdt_context = [0, diff(xs_filtered)];
         end
    end
    
    % --- Stage B: Peak Detection within Analysis Window ---
    
    % 3. Isolate the Derivative Trace for the Analysis Window
    dvdt_analysis = dvdt_context(WinStartIdx : WinEndIdx);
    
    % 4. Find the max absolute change within the isolated segment
    [dvdt_max, pos_max_analysis] = max(dvdt_analysis);
    [dvdt_min, pos_min_analysis] = min(dvdt_analysis);
    
    if abs(dvdt_max) > abs(dvdt_min)
        pos_relative_to_analysis = pos_max_analysis;
    else
        pos_relative_to_analysis = pos_min_analysis;
    end
    
    % 5. Map the index back to the START of the Context Window (N-length index)
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
        % Zoom lines to show window region
        plot([WinStartIdx WinEndIdx], [xs_filtered(position_final) xs_filtered(position_final)], 'k:', 'LineWidth', 0.5);
        title('1. Full Signal Trace (Red = Detected, Dashed = Search Window)');
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
        title('2. Signal in Analysis Window Only');
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
        title('3. dV/dt in Analysis Window (Savitzky-Golay)');
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