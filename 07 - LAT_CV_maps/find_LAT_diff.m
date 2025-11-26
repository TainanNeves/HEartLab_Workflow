function LAT_E = find_LAT_diff(Data, Fs, case_name, x_coord, y_coord, debug)
    % FIND_LAT_DIFF Calculate Local Activation Time (LAT) for a single signal trace.
    %
    % Inputs:
    %   Data:      1D vector of electrical data (single trace).
    %   Fs:        Sampling frequency of the data.
    %   case_name: String identifying the current data case.
    %   x_coord:   X-coordinate (row index) of the current pixel.
    %   y_coord:   Y-coordinate (column index) of the current pixel.
    %   debug:     Debug mode flag (0, 1 or 2).
    %
    % Output:
    %   LAT_E:     Local Activation Time (LAT) in milliseconds.
    
    % Apply moving average filter
    windowSize = 100;
    b = (1/windowSize)*ones(1,windowSize);
    a = 1;
    xs = filter(b, a, Data, [], 2);
    
    % --- Core LAT Logic ---
    
    if max(xs) == 0
        LAT_E = 0;
        return;
    end
    
    % Calculate the first derivative (approximate dV/dt)
    dvdt = -diff(xs);
    [dvdt_max, position_max] = max(dvdt);
    [dvdt_min, position_min] = min(dvdt);
    
    % Assign position based on the index with highest absolute value
    if abs(dvdt_max) > abs(dvdt_min)
        position = position_max;
    else
        position = position_min;
    end
    
    % Create a clear overall title
    plotTitle = ['Case: ', case_name, ' | Pixel: (', num2str(x_coord), ',', num2str(y_coord), ')'];
    
    % --- Debug Mode 1: Plotting ---
    if debug == 1
        figure(101);
        clf; 
        
        subplot(2,1,1);
        plot(xs);
        hold on;
        plot(position, xs(position), 'o', 'markersize', 10, 'color', 'red');
        title('Signal Trace and Automatic LAT Point');
        hold off;

        subplot(2,1,2);
        plot(dvdt);
        xline(position, 'r', 'LineWidth', 2);
        title('Negative Derivative (-dV/dt)');
        hold off;
        
        sgtitle(plotTitle, 'FontSize', 12, 'FontWeight', 'bold');
        drawnow; 
    end
    
    % --- Debug Mode 2: Interactive Correction (Optimized Speed) ---
    if debug == 2
        figure(101);
        clf; % Clear the figure for the new pixel
        
        % Subplot 1: Signal Trace
        subplot(2,1,1);
        plot(xs);
        hold on;
        h_auto = plot(position, xs(position), 'o', 'markersize', 10, 'color', 'red');
        title('1. Signal Trace (Press SPACE or Click)');
        
        % Subplot 2: Derivative
        subplot(2,1,2);
        plot(dvdt);
        h_line = xline(position, 'r', 'LineWidth', 2);
        title('2. Negative Derivative (-dV/dt)');
        
        % Overall Figure Title and Instructions
        sgtitle({plotTitle, '--- Press SPACE to accept, CLICK to select new point ---'}, ...
                'FontSize', 12, 'FontWeight', 'bold', 'Color', 'blue');
        
        drawnow;
        
        % Get user input
        k = waitforbuttonpress;
        
        if k == 1 % Key press
            key = get(gcf, 'CurrentCharacter');
            if key ~= ' ' 
                % Any key other than space is pressed, do nothing (keep auto)
                disp('Key other than SPACE pressed. Defaulting to automatic point.');
            end
            % Execution proceeds immediately after key press
            
        else % Mouse click occurred
            % Prompt the user to select on the signal trace (subplot 1)
            subplot(2,1,1);
            
            % Manually select a new point
            [x_click, ~] = ginput(1);
            position = round(x_click);
            
            % Update plots
            delete(h_auto); % Remove the old signal marker
            subplot(2,1,1); % Select signal plot again for manual marker
            plot(position, xs(position), 's', 'markersize', 10, 'color', 'green'); % New green square marker
            
            delete(h_line); % Remove the old derivative line
            subplot(2,1,2); % Select derivative plot
            xline(position, 'g', 'LineWidth', 2); % New green line
            
            disp(['Manually selected position: ', num2str(position)]);
        end
        
        % The hold off and final title update will happen, and the loop moves on.
        hold off; 
    end
    
    LAT_E = position * 1000 / Fs; % Convert to time in ms
end


%% Version 2
% function LAT_E = find_LAT_diff(Data, Fs, lim1, lim2, debug)
%     % FIND_LAT_DIFF Calculate Local Activation Time (LAT) for signals from electrodes.
%     % 
%     % LAT_E = find_LAT_diff(Data, Fs, lim1, lim2, debug)
%     %
%     % Inputs:
%     %   Data:    Matrix of electrical data with electrodes in rows and time points in columns.
%     %   Fs:      Sampling frequency of the data.
%     %   lim1:    Starting index of the segment of electrical data.
%     %   lim2:    Ending index of the segment of electrical data.
%     %   debug:   Debug mode flag (0, 1 or 2). 
%     %               If 1, plots the selected point for each electrode.
%     %               If 2, plots the selected points and allow correction.
%     %
%     % Output:
%     %   LAT_E:   Local Activation Time (LAT) for each electrode.
%     %
% 
%     % Extract a segment of electrical data
%     Data_temp = Data(:, lim1:lim2);
% 
%     % Apply moving average filter along each row
%     windowSize = 100; % Choose an appropriate window size
%     b = (1/windowSize)*ones(1,windowSize);
%     a = 1;
%     Data_temp_E = filter(b, a, Data_temp, [], 2); % Apply along rows
% 
%     LAT_E = zeros(size(Data_temp_E, 1), 1); % Initialize LAT_E
%     position = zeros(size(Data_temp_E, 1), 1); % Initialize position
% 
%     % Calculate Local Activation Time (LAT) for electrical data
%     for i = 1:size(Data_temp_E, 1)
%         if max(Data_temp_E(i, :)) ~= 0
%             xs = Data_temp_E(i, :);
%             [dvdt_max, position_max] = max(-diff(xs));
%             [dvdt_min, position_min] = min(-diff(xs));
%             % Assign position based on the index with highest absolute value
%             if abs(dvdt_max) > abs(dvdt_min)
%                 position(i) = position_max;
%             else
%                 position(i) = position_min;
%             end
%             if debug == 1
%                 figure(101);
%                 title(['Electrode: ', num2str(i)]);
%                 subplot(2,1,1);
%                 plot(xs);
%                 hold on;
%                 plot(position(i), xs(position(i)), 'o', 'markersize', 10, 'color', 'red');
%                 subplot(2,1,2);
%                 plot(-diff(xs));
%                 xline(position(i), 'r', LineWidth=2);
%                 hold off;
%             end
%             if debug == 2
%                 figure(101);
%                 plot(xs);
%                 hold on;
%                 plot(position(i), xs(position(i)), 'o', 'markersize', 10, 'color', 'red');
%                 hold off;
%                 title(['Electrode: ', num2str(i)]);
% 
%                 % Wait for user input to either accept the automatically selected point or select a new one
%                 text(position(i), xs(position(i)), 'Press SPACE to accept, click to select a new point', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'blue');
%                 waitforbuttonpress;
%                 key = get(gcf, 'CurrentCharacter');
%                 if key == ' ' % If SPACE is pressed, accept the automatically selected point
%                     % Do nothing, keep the automatically selected point
%                 else
%                     % Manually select a new point
%                     [x, ~] = ginput(1);
%                     position(i) = round(x);
%                 end
%             end
%             LAT_E(i) = position(i) * 1000 / Fs; % Convert to time in ms
%         end
%     end
% 
%     % Subtract the minimum value from LAT_E
%     mini = min(LAT_E(LAT_E ~= 0));
%     LAT_E = LAT_E - mini; % Negative values will be unused electrodes or any error
% 
%     % Invert the matrix
%     LAT_E = LAT_E';
% end


%% Primeira versão Tainan Funcional
% function LAT_E = find_LAT_diff(Data, Fs, lim1, lim2, debug)
%     % FIND_LAT_DIFF Calculate Local Activation Time (LAT) for signals from electrodes.
%     % 
%     % LAT_E = find_LAT_diff(Data, Fs, lim1, lim2, debug)
%     %
%     % Inputs:
%     %   Data:    Matrix of electrical data with electrodes in rows and time points in columns.
%     %   Fs:      Sampling frequency of the data.
%     %   lim1:    Starting index of the segment of electrical data.
%     %   lim2:    Ending index of the segment of electrical data.
%     %   debug:   Debug mode flag (0 or 1). If 1, plots the selected point for each electrode.
%     %
%     % Output:
%     %   LAT_E:   Local Activation Time (LAT) for each electrode.
%     %
% 
%     % Extract a segment of electrical data
%     Data_temp = Data(:, lim1:lim2);
% 
%     % Apply moving average filter along each row
%     windowSize = 100; % Choose an appropriate window size
%     b = (1/windowSize)*ones(1,windowSize);
%     a = 1;
%     Data_temp_E = filter(b, a, Data_temp, [], 2); % Apply along rows
% 
%     LAT_E = zeros(size(Data_temp_E, 1), 1); % Initialize LAT_E
%     position = zeros(size(Data_temp_E, 1), 1); % Initialize position
%     
%     % Calculate Local Activation Time (LAT) for electrical data
%     % Plot selected point for each electrode if debug is 1
%     for i = 1:size(Data_temp_E, 1)
%         if max(Data_temp_E(i, :)) ~= 0
%             xs = Data_temp_E(i, :);
%             [dvdt_max, position_max] = max(-diff(xs));
%             [dvdt_min, position_min] = min(-diff(xs));
%             % Assign position based on the index with highest absolute value
%             if abs(dvdt_max) > abs(dvdt_min)
%                 position(i) = position_max;
%             else
%                 position(i) = position_min;
%             end
%             if debug == 1
%                 figure();
%                 title(['Electrode: ', num2str(i)]);
%                 subplot(2,1,1);
%                 plot(xs);
%                 hold on;
%                 plot(position(i), xs(position(i)), 'o', 'markersize', 10, 'color', 'red');
%                 subplot(2,1,2);
%                 plot(-diff(xs));
%                 xline(position(i), 'r', LineWidth=2);
%                 hold off;
%             end
%             LAT_E(i) = position(i) * 1000 / Fs; % Convert to time in ms
%         end
%     end
%     
%     % Subtract the minimum value from LAT_E
%     mini = min(LAT_E(LAT_E ~= 0));
%     LAT_E = LAT_E - mini; % Negative values will be unused electrodes or any error
%     
% end