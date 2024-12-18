function df_boxplot(df_values_meas_tank, df_values_ecgi, meas_signal, file_id_meas, resolution)
    % Compare DF values using boxplots by region (MEA1: RA, MEA2: V, MEA3: LA)
    
    % Initialize containers for DF values
    df_meas_mea1 = [];
    df_meas_mea2 = [];
    df_meas_mea3 = [];
    
    df_ecgi_mea1 = [];
    df_ecgi_mea2 = [];
    df_ecgi_mea3 = [];
    
    % Loop through electrodes and get DF values by region
    for electrode = 1:80
        if electrode < 33 || electrode > 64
            [mea, id_mea, electrode_mea] = get_mea_electrode(electrode, meas_signal, file_id_meas, resolution);
            
            % Append measured DF values
            if electrode < 17
                df_meas_mea1 = [df_meas_mea1; df_values_meas_tank.DF(electrode)];
            elseif electrode < 33
                df_meas_mea2 = [df_meas_mea2; df_values_meas_tank.DF(electrode)];
            else
                df_meas_mea3 = [df_meas_mea3; df_values_meas_tank.DF(electrode)];
            end
            
            % Append estimated DF values
            if electrode < 17
                df_ecgi_mea1 = [df_ecgi_mea1; df_values_ecgi.DF(electrode)];
            elseif electrode < 33
                df_ecgi_mea2 = [df_ecgi_mea2; df_values_ecgi.DF(electrode)];
            else
                df_ecgi_mea3 = [df_ecgi_mea3; df_values_ecgi.DF(electrode)];
            end
        end
    end
    
    % Prepare data for boxplot
    df_meas = [df_meas_mea1; df_meas_mea2; df_meas_mea3];
    df_ecgi = [df_ecgi_mea1; df_ecgi_mea2; df_ecgi_mea3];
    
    % Group labels for boxplot
    group_labels = [repmat({'MEA1 (RA)'}, length(df_meas_mea1), 1); ...
                    repmat({'RA - ECGi'}, length(df_ecgi_mea1), 1); ...
                    repmat({'MEA2 - V'}, length(df_meas_mea2), 1); ...
                    repmat({'V - ECGi'}, length(df_ecgi_mea2), 1); ...
                    repmat({'MEA3 - LA'}, length(df_meas_mea3), 1); ...
                    repmat({'LA - ECGi'}, length(df_ecgi_mea3), 1)];
    
    % Combined DF values for boxplot
    combined_df = [df_meas_mea1; df_ecgi_mea1; df_meas_mea2; df_ecgi_mea2; df_meas_mea3; df_ecgi_mea3];
    
    % Create the boxplot
    figure;
    boxplot(combined_df, group_labels);
    title('Dominant Frequency MEAs x ECGi');
    ylabel('Dominant Frequency (Hz)');
    xlabel('Region');
    set(gca, 'FontSize', 12);
    grid on;
    
    % Customize plot appearance
    boxplot_handle = findobj(gca, 'Tag', 'Box');
    colors = lines(2); % Define colors for measured and estimated data
    
    % Set colors for each region
    for i = 1:length(boxplot_handle)
        patch(get(boxplot_handle(i), 'XData'), get(boxplot_handle(i), 'YData'), ...
              colors(mod(i-1, 2) + 1, :), 'FaceAlpha', .5);
    end
    
    % Add legend
    legend({'ECGi', 'MEAs'}, 'Location', 'BestOutside');
end
