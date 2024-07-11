function OI = calculate_OI(MFFTi, Sffti, fstep, Hzi, Hzf, dfh_threshold_area, f_mode, debug)
    % calculate_OI - Calculates the Organization Index (OI) for electric signal data
    %
    % Syntax: OI = calculate_OI(MFFTi, Sffti, fstep, Hzi, Hzf, dfh_threshold_area, mode, debug)
    %
    % Inputs:
    %   MFFTi - Array containing dominant frequencies (size: x x y)
    %   Sffti - Frequency spectrum (size: x x y x Hzf/fstep)
    %   fstep - Step size for the frequency spectrum
    %   Hzi - Lower limit of the frequency range of interest
    %   Hzf - Upper limit of the frequency range of interest
    %   dfh_threshold_area - Threshold around MFFTi to calculate area
    %   mode - Calculation mode (1: Only DF, 2: DF + Harmonics)
    %   debug - Debug mode flag (0: off, 1: on)
    %
    % Outputs:
    %   OI - Organization Index matrix (size: x x y)
    %
    % Description:
    %   This function calculates the Organization Index (OI) for each (x, y) coordinate
    %   in the input data. The OI is calculated as the ratio of the area around the dominant
    %   frequency (and harmonics, if mode 2 is selected) to the total area under the spectrum
    %   within a specified frequency range. The function can optionally plot the spectrum and
    %   the regions considered for OI calculation when debug mode is enabled.
    %
    % Example Usage:
    %   OI = calculate_OI(MFFTi, Sffti, fstep, 3.5, 8, 0.5, 1, 0);
    %   OI = calculate_OI(MFFTi, Sffti, fstep, 3.5, 8, 0.5, 2, 1);
    %
    % Notes:
    %   - The dominant frequency (DF) is estimated for each (x, y) coordinate.
    %   - In mode 1, OI is calculated considering only the area around the DF.
    %   - In mode 2, OI is calculated considering the area around the DF and its harmonics within the specified range.
    %   - Debug mode (when enabled) plots the spectrum and regions considered for OI calculation for visual inspection.
    %
    % Authors:
    %   Tainan Neves and Bruno 
    %   HEartLab - UFABC
    
    % Initialize OI matrix
    OI = zeros(size(MFFTi));
    
    % Loop through each electrode (x, y)
    for x = 1:size(MFFTi, 1)
        disp(['Line: ', num2str(x)]);
        for y = 1:size(MFFTi, 2)
            % Calculate total area under the spectrum in the specified range for this electrode
            Hzf_steps = size(Sffti, 3);
            isample_area = round(Hzi / fstep);
            fsample_area = round(Hzf / fstep);
            if fsample_area > Hzf_steps
                fsample_area = Hzf_steps;
            end
            area_total = sum(Sffti(x, y, isample_area:fsample_area), 'all');

            % Calculate areas around MFFTi and harmonics
            sum_area = 0;
            
            % Find indices around MFFTi within threshold
            isample = round((MFFTi(x, y) - dfh_threshold_area) / fstep);
            fsample = round((MFFTi(x, y) + dfh_threshold_area) / fstep);

            % Sum areas within the threshold
            if isample > 0 && fsample > 0 && fsample <= Hzf_steps
                area = sum(Sffti(x, y, isample:fsample), 'all');
                sum_area = sum(area, 'all');
            end

            % Mode 2: Consider harmonics in the count
            if f_mode == 2
                harmonic = MFFTi(x, y);
                while harmonic + dfh_threshold_area < Hzf
                    harmonic = harmonic + MFFTi(x, y);
                    if harmonic - dfh_threshold_area >= Hzi && harmonic + dfh_threshold_area <= Hzf
                        isample_harm = round((harmonic - dfh_threshold_area) / fstep);
                        fsample_harm = round((harmonic + dfh_threshold_area) / fstep);

                        if isample_harm > 0 && fsample_harm > 0 && fsample_harm <= Hzf_steps
                            area_harm = sum(Sffti(x, y, isample_harm:fsample_harm), 'all');
                            sum_area = sum_area + sum(area_harm, 'all');
                        end
                    end
                end
            end

            % Calculate Organization Index (OI) for this electrode
            if area_total > 0
                OI(x, y) = sum_area / area_total;
            else
                OI(x, y) = 0; % Handle division by zero case
            end
            
            % Debug mode: Plotting regions considered for OI calculation
            if debug == 1 && y == round(size(MFFTi, 2)/2)
                % Plot the spectrum and the regions considered
                figure;
                plot(fstep * (0:Hzf_steps-1), squeeze(Sffti(x, y, :)), 'b-', 'LineWidth', 1.5);
                hold on;
                plot([Hzi, Hzi], [0, max(Sffti(x, y, :))], 'r--', 'LineWidth', 1.5);
                plot([Hzf, Hzf], [0, max(Sffti(x, y, :))], 'r--', 'LineWidth', 1.5);
                plot([MFFTi(x, y) - dfh_threshold_area, MFFTi(x, y) - dfh_threshold_area], [0, max(Sffti(x, y, :))], 'g--', 'LineWidth', 1.5);
                plot([MFFTi(x, y) + dfh_threshold_area, MFFTi(x, y) + dfh_threshold_area], [0, max(Sffti(x, y, :))], 'g--', 'LineWidth', 1.5);
                
                if f_mode == 2
                    harmonic = MFFTi(x, y);
                    while harmonic + dfh_threshold_area < Hzf
                        harmonic = harmonic + MFFTi(x, y);
                        if harmonic - dfh_threshold_area >= Hzi && harmonic + dfh_threshold_area <= Hzf
                            plot([harmonic - dfh_threshold_area, harmonic - dfh_threshold_area], [0, max(Sffti(x, y, :))], 'c--', 'LineWidth', 1.5);
                            plot([harmonic + dfh_threshold_area, harmonic + dfh_threshold_area], [0, max(Sffti(x, y, :))], 'c--', 'LineWidth', 1.5);
                        end
                    end
                end
                
                title(sprintf('Electrode (%d, %d) - Spectrum and OI Calculation Regions', x, y));
                xlabel('Frequency (Hz)');
                ylabel('Amplitude');
                legend('Spectrum', 'Frequency Range [Hzi]', 'Frequency Range [Hzf]', 'DF - Threshold', 'DF + Threshold', 'Harmonics +/- threshold', 'Location', 'best');
                grid on;
                hold off;
            end
        end
    end
end
