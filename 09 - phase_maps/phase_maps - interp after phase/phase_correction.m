function [Pe_new] = phase_correction(A_e, Pe, F_e)
% PHASE_CORRECTION Corrects the phase of signals by interpolating between peak points.
% 
% [Pe_new] = phase_correction(A_e, Pe, F_e)
%   - A_e: Matrix of original signals (each row is a separate signal)
%   - Pe: Matrix of phase
%   - F_e: Sampling frequency
% 
% Returns:
%   - Pe_new: Matrix with corrected phase signals

% Initialize the output matrix with zeros, matching the size of Pe
Pe_new = zeros(size(Pe));

% Loop through each signal in A_e
for ele = 1:size(A_e, 1)
    % Extract the current signal and corresponding phase data
    signal = A_e(ele, :);
    phase = Pe(ele, :);

    % Find positive and negative peaks in the phase data with a minimum height of 3
    [peaks_max, peak_indices_max] = findpeaks(phase, 'MinPeakHeight', 3);
    [peaks_min, peak_indices_min] = findpeaks(-phase, 'MinPeakHeight', 3);

    % Only proceed if there are positive peaks detected
    if ~isempty(peaks_max)
        % Get the first positive peak and the last negative peak
        x1 = peak_indices_max(1);
        x2 = peak_indices_min(end);
        y1 = peaks_max(1);
        y2 = -peaks_min(end);

        % Create the x and y coordinates for linear interpolation
        x_interpol = [x1, x2];
        y_interpol = [y1, y2];

        % Generate a dense x-axis for interpolating between the peaks
        x_linea_interpol = linspace(x1, x2, x2 - x1);

        // Perform linear interpolation to get the corrected phase values
        y_linea_interpol = interp1(x_interpol, y_interpol, x_linea_interpol, 'linear');

        // Replace the phase data between x1 and x2 with the interpolated values
        corrected_phase = phase; % Make a copy of the original phase data
        corrected_phase(x1 + 1:x2) = y_linea_interpol;

        % Store the corrected phase in the output matrix
        Pe_new(ele, :) = corrected_phase;

        % Create a plot to visualize the original signal and the corrected phase
        figure('Color', 'white', 'Position', [40, 40, 300, 350]); % Set plot settings
        % Time axis based on sampling frequency
        time_axis = linspace(0, length(signal) / F_e, length(signal)); 

        % Plot the original signal
        subplot(2, 1, 1);
        plot(time_axis, signal, 'LineWidth', 1, 'Color', 'black'); 
        title('Original Signal');
        ylabel('Amplitude [mV]');
        set(gca, 'FontSize', 12); 

        % Plot the original and corrected phase
        subplot(2, 1, 2);
        plot(time_axis, phase, 'LineWidth', 1, 'Color', 'black'); % Original phase
        hold on;
        plot(time_axis, corrected_phase, 'LineWidth', 1, 'Color', 'blue'); % Corrected phase
        title('Corrected Phase');
        xlabel('Time [s]');
        ylabel('Phase [rad]');
        set(gca, 'FontSize', 12);

        % Link the x-axes of the subplots for synchronized zoom/pan
        linkaxes(findobj(gcf, 'Type', 'axes'), 'x');
    end
end
end
