function [results] = ica_analysis(Data_matrix, array_name, electrode_indices, To, Fsampling, componentsToUse)
    % ica_analysis performs ICA on the given data matrix using FastICA
    %
    % Inputs:
    % - Data_matrix: [numChannels x numSamples] matrix of the data
    % - array_name: String indicating the name of the array (e.g., 'MEA 1')
    % - electrode_indices: Indices of the electrodes in the original data
    % - To: Time vector
    % - Fsampling: Sampling frequency
    % - componentsToUse: Vector of specific component indices to use (e.g., [9] or [1,2,3])
    
    % Set random seed for reproducibility
    rng(42);
    
    % Number of channels and samples
    [numChannels, numSamples] = size(Data_matrix);
    
    % Ensure data is in correct orientation (signals in rows)
    if numChannels > numSamples
        Data_matrix = Data_matrix';
        [numChannels, numSamples] = size(Data_matrix);
    end
    
    % Center the data
    Data_centered = Data_matrix - mean(Data_matrix, 2);
    
    % Determine number of components to extract from ICA
    % Always extract all possible components, then select specific ones for reconstruction
    numComponentsToExtract = numChannels;
    
    % Perform FastICA to extract ALL components
    [icasig, A, W] = fastica(Data_centered, ...
                            'numOfIC', numComponentsToExtract, ...
                            'approach', 'symm', ...
                            'g', 'tanh', ...
                            'maxNumIterations', 1000, ...
                            'epsilon', 1e-6, ...
                            'verbose', 'off');
% %     Inputs
% % X: Pre-processed (centered) data matrix, with rows as variables and columns as samples
% % Parameters:
% % 'numOfIC': Number of independent components to extract
% % 'approach': Estimation approach ('symm' for simultaneous, 'defl' for deflation)
% % 'g': Nonlinearity function ('tanh', 'pow3', 'gauss', etc.)
% % 'maxNumIterations': Maximum number of algorithm iterations
% % 'epsilon': Convergence criterion threshold
% % 'verbose': Display progress information ('on' or 'off')
% % Outputs
% % icasig: Independent components, with rows as components and columns as samples
% % A: Mixing matrix (X ≈ A × icasig)
% % W: Unmixing matrix (icasig = W × X)
    % Select only the specified components for reconstruction
    icasig_selected = icasig(componentsToUse, :);
    A_selected = A(:, componentsToUse);
    
    % Plot ALL ICA components in groups of 8 for reference
% Plot ALL ICA components in groups of 16 using a 4x4 layout
    numPerFigure = 16;
    numFigures = ceil(numComponentsToExtract / numPerFigure);
    
    for figIndex = 1:numFigures
        figure('Name', [array_name ': ALL ICA Components ' ...
                num2str((figIndex-1)*numPerFigure + 1) '-' ...
                num2str(min(figIndex*numPerFigure, numComponentsToExtract))], ...
                'Position', [100 100 1400 900]);
    
        startComp = (figIndex-1)*numPerFigure + 1;
        endComp = min(figIndex*numPerFigure, numComponentsToExtract);
    
        for i = startComp:endComp
            subplot(4, 4, i - startComp + 1);
            plot(To, icasig(i, :));
            title(['IC ' num2str(i)]);
            xlabel('Time (s)');
            ylabel('Amplitude');
            xlim([0 5]);         % Limit x-axis to 0–1 second
            grid on;
        end
    end

    
    % Highlight the selected components
    figure('Name', [array_name ': Selected ICA Components'], ...
           'Position', [100 100 1200 800]);
    
    for i = 1:length(componentsToUse)
        subplot(length(componentsToUse),1,i);
        plot(To, icasig(componentsToUse(i),:), 'LineWidth', 1.5);
        title(['Selected IC ' num2str(componentsToUse(i))]);
        xlabel('Time (s)');
        ylabel('Amplitude');
        grid on;
    end
    
    % Plot electrode contributions (mixing matrix) for selected components
%     for i = 1:length(componentsToUse)
%         figure('Name', [array_name ': Electrode Contributions to IC' num2str(componentsToUse(i))], ...
%                'Position', [100 100 800 400]);
%         bar(electrode_indices, A(:, componentsToUse(i)));
%         xlabel('Electrode Index');
%         ylabel('Weight');
%         title(['Electrode Contributions to IC' num2str(componentsToUse(i))]);
%         grid on;
%     end
    
    % Reconstruct data using only selected components
    Data_reconstructed = A_selected * icasig_selected;
    
    % Add mean back to reconstructed data
    Data_reconstructed = Data_reconstructed + mean(Data_matrix, 2);
    
    % Calculate MSE
    mse = mean((Data_matrix(:) - Data_reconstructed(:)).^2);
    disp(['Mean Squared Error of Reconstruction: ' num2str(mse)]);
    
    % Plot original vs reconstructed for first channel
    figure('Name', [array_name ': Original vs Reconstructed Signal'], ...
           'Position', [100 100 800 600]);
    
    subplot(2,1,1);
    plot(To, Data_matrix(6,:));
    title(['Original Signal - Electrode ' num2str(electrode_indices(6))]);
    xlabel('Time (s)');
    ylabel('Amplitude');
    xlim([0 1])
    grid on;
    
    subplot(2,1,2);
    plot(To, Data_reconstructed(6,:));
    title(['Reconstructed Signal using IC' strjoin(arrayfun(@num2str, componentsToUse, 'UniformOutput', false), ',')]);
    xlabel('Time (s)');
    ylabel('Amplitude');
    xlim([0 1])
    grid on;
    
%     % Power spectrum analysis of selected ICs
%     for i = 1:length(componentsToUse)
%         [pxx, f] = pwelch(icasig(componentsToUse(i),:), [], [], [], Fsampling);
%         
%         figure('Name', [array_name ': Power Spectrum of IC' num2str(componentsToUse(i))], ...
%                'Position', [100 100 800 400]);
%         plot(f, 10*log10(pxx));
%         xlabel('Frequency (Hz)');
%         ylabel('Power/Frequency (dB/Hz)');
%         title(['Power Spectrum of IC' num2str(componentsToUse(i))]);
%         grid on;
%     end
    
    % Store results in a structure
    results.icasig = icasig_selected;
    results.A = A_selected;
    results.W = W;
    results.Data_reconstructed = Data_reconstructed;
    results.mse = mse;
    results.full_icasig = icasig;  % Store full ICA results too
    results.full_A = A;
    results.selected_components = componentsToUse;
    
    % Optional: assign to base workspace if needed
    assignin('base', [strrep(array_name, ' ', '_') '_results'], results);
end