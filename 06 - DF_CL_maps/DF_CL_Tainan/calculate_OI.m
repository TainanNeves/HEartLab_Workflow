function OI = calculate_OI(MFFTi, Sffti, fstep, Hzi, Hzf, dfh_threshold_area, f_mode, debug)
    % calculate_OI - Calculates the Organization Index (OI) for electric signal data
    %
    % Improved version with better error handling and performance
    %
    % Inputs:
    %   MFFTi - Array containing dominant frequencies (size: x x y)
    %   Sffti - Frequency spectrum (size: x x y x Hzf/fstep)
    %   fstep - Step size for the frequency spectrum
    %   Hzi - Lower limit of the frequency range of interest
    %   Hzf - Upper limit of the frequency range of interest
    %   dfh_threshold_area - Threshold around MFFTi to calculate area
    %   f_mode - Calculation mode (1: Only DF, 2: DF + Harmonics)
    %   debug - Debug mode flag (0: off, 1: on)
    %
    % Outputs:
    %   OI - Organization Index matrix (size: x x y)

    % Input validation
    if nargin < 8
        debug = 0;
    end
    if nargin < 7
        f_mode = 1;
    end
    if nargin < 6
        dfh_threshold_area = 0.5;
    end
    
    % Validate frequency range
    if Hzi >= Hzf
        error('Hzi must be less than Hzf');
    end
    if dfh_threshold_area <= 0
        error('dfh_threshold_area must be positive');
    end

    % Initialize OI matrix
    OI = zeros(size(MFFTi));
    
    % Pre-calculate frequency range indices
    Hzf_steps = size(Sffti, 3);
    isample_area = max(1, round(Hzi / fstep));
    fsample_area = min(Hzf_steps, round(Hzf / fstep));
    
    % Check if frequency range is valid
    if isample_area >= fsample_area
        warning('Invalid frequency range. Adjust Hzi and Hzf.');
        return;
    end

    % Pre-calculate frequency vector for plotting (if debug mode)
    if debug
        freq_vector = fstep * (0:Hzf_steps-1);
    end
    
    % Counter for debug plots (limit to avoid too many figures)
    debug_plot_count = 0;
    max_debug_plots = 5; % Maximum number of debug plots to show
    
    % Loop through each electrode (x, y)
    for x = 1:size(MFFTi, 1)
        fprintf('Processing line %d/%d...\n', x, size(MFFTi, 1));
        
        for y = 1:size(MFFTi, 2)
            % Skip if dominant frequency is zero or invalid
            if MFFTi(x, y) == 0 || isnan(MFFTi(x, y)) || MFFTi(x, y) < Hzi || MFFTi(x, y) > Hzf
                OI(x, y) = 0;
                continue;
            end
            
            try
                % Calculate total area under the spectrum in the specified range
                spectrum_segment = squeeze(Sffti(x, y, isample_area:fsample_area));
                area_total = sum(spectrum_segment, 'all');

                % Skip if total area is too small
                if area_total < eps
                    OI(x, y) = 0;
                    continue;
                end

                % Calculate areas around MFFTi
                sum_area = 0;
                
                % Main dominant frequency area
                isample_df = max(1, round((MFFTi(x, y) - dfh_threshold_area) / fstep));
                fsample_df = min(Hzf_steps, round((MFFTi(x, y) + dfh_threshold_area) / fstep));
                
                if isample_df <= fsample_df
                    df_area = sum(Sffti(x, y, isample_df:fsample_df), 'all');
                    sum_area = sum_area + df_area;
                end

                % Mode 2: Consider harmonics
                if f_mode == 2
                    harmonic = 2 * MFFTi(x, y); % Start with first harmonic
                    harmonic_count = 0;
                    max_harmonics = 5; % Limit number of harmonics to prevent infinite loops
                    
                    while harmonic - dfh_threshold_area <= Hzf && harmonic_count < max_harmonics
                        % Check if harmonic is within valid range
                        if harmonic + dfh_threshold_area >= Hzi
                            isample_harm = max(1, round((harmonic - dfh_threshold_area) / fstep));
                            fsample_harm = min(Hzf_steps, round((harmonic + dfh_threshold_area) / fstep));
                            
                            if isample_harm <= fsample_harm
                                harm_area = sum(Sffti(x, y, isample_harm:fsample_harm), 'all');
                                sum_area = sum_area + harm_area;
                            end
                        end
                        
                        harmonic_count = harmonic_count + 1;
                        harmonic = (harmonic_count + 2) * MFFTi(x, y); % Next harmonic
                    end
                end

                % Calculate Organization Index (OI)
                OI(x, y) = sum_area / area_total;
                
                % Ensure OI is between 0 and 1
                OI(x, y) = max(0, min(1, OI(x, y)));

                % Debug mode: Plotting (limited number of plots)
                if debug == 1 && debug_plot_count < max_debug_plots && y == round(size(MFFTi, 2)/2)
                    debug_plot_count = debug_plot_count + 1;
                    plotDebugSpectrum(x, y, Sffti, freq_vector, MFFTi, dfh_threshold_area, ...
                                     Hzi, Hzf, f_mode, area_total, sum_area, OI(x, y));
                end
                
            catch ME
                warning('Error processing electrode (%d, %d): %s', x, y, ME.message);
                OI(x, y) = 0;
                continue;
            end
        end
    end
    
    fprintf('OI calculation completed.\n');
end

function plotDebugSpectrum(x, y, Sffti, freq_vector, MFFTi, dfh_threshold_area, ...
                          Hzi, Hzf, f_mode, area_total, sum_area, oi_value)
    % Helper function for debug plotting
    
    figure;
    plot(freq_vector, squeeze(Sffti(x, y, :)), 'b-', 'LineWidth', 1.5);
    hold on;
    
    % Plot frequency range boundaries
    y_max = max(Sffti(x, y, :));
    plot([Hzi, Hzi], [0, y_max], 'r--', 'LineWidth', 1.5, 'DisplayName', 'Hzi');
    plot([Hzf, Hzf], [0, y_max], 'r--', 'LineWidth', 1.5, 'DisplayName', 'Hzf');
    
    % Plot DF region
    df = MFFTi(x, y);
    plot([df - dfh_threshold_area, df - dfh_threshold_area], [0, y_max], 'g--', 'LineWidth', 1.5, 'DisplayName', 'DF - Threshold');
    plot([df + dfh_threshold_area, df + dfh_threshold_area], [0, y_max], 'g--', 'LineWidth', 1.5, 'DisplayName', 'DF + Threshold');
    
    % Shade the DF region
    x_fill = [df - dfh_threshold_area, df + dfh_threshold_area, df + dfh_threshold_area, df - dfh_threshold_area];
    y_fill = [0, 0, y_max, y_max];
    fill(x_fill, y_fill, 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'DF Area');
    
    % Plot harmonics if mode 2
    if f_mode == 2
        harmonic = 2 * df;
        harmonic_count = 0;
        max_harmonics = 5;
        
        while harmonic - dfh_threshold_area <= Hzf && harmonic_count < max_harmonics
            if harmonic + dfh_threshold_area >= Hzi
                plot([harmonic - dfh_threshold_area, harmonic - dfh_threshold_area], [0, y_max], 'c--', 'LineWidth', 1.5);
                plot([harmonic + dfh_threshold_area, harmonic + dfh_threshold_area], [0, y_max], 'c--', 'LineWidth', 1.5);
                
                % Shade harmonic region
                x_fill_harm = [harmonic - dfh_threshold_area, harmonic + dfh_threshold_area, ...
                              harmonic + dfh_threshold_area, harmonic - dfh_threshold_area];
                fill(x_fill_harm, y_fill, 'c', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Harmonic Area');
            end
            
            harmonic_count = harmonic_count + 1;
            harmonic = (harmonic_count + 2) * df;
        end
    end
    
    title(sprintf('Electrode (%d, %d) - OI = %.3f\nTotal Area: %.3f, DF Area: %.3f', ...
          x, y, oi_value, area_total, sum_area));
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');
    legend('Location', 'best');
    grid on;
    xlim([0, Hzf * 1.1]);
    hold off;
end