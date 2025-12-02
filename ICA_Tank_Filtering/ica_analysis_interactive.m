function [results] = ica_analysis_interactive(Data_matrix, array_name, electrode_indices, To, Fsampling)
    % ica_analysis_interactive performs ICA, prompts user for component selection
    % using a dialog box, and then reconstructs the data.
    
    % Set random seed for reproducibility
    rng(42);
    
    % --- 1. Pre-processing ---
    [numChannels, numSamples] = size(Data_matrix);
    
    % Ensure data is in correct orientation (signals in rows)
    if numChannels > numSamples
        Data_matrix = Data_matrix';
        [numChannels, numSamples] = size(Data_matrix);
    end
    
    % Center the data
    Data_centered = Data_matrix - mean(Data_matrix, 2);
    
    % --- 2. Perform FastICA (Extract ALL 5 components) ---
    numComponentsToExtract = numChannels; 
    
    [icasig, A, W] = fastica(Data_centered, ...
                            'numOfIC', numComponentsToExtract, ...
                            'approach', 'symm', ...
                            'g', 'tanh', ...
                            'maxNumIterations', 1000, ...
                            'epsilon', 1e-6, ...
                            'verbose', 'off');
    
    % --- 3. Plot ALL ICA components (for user selection) ---
    figure('Name', [array_name ': ALL ICA Components for SELECTION'], ...
            'Position', [100 100 1000 900], 'Color', 'white');
    
    for i = 1:numChannels
        subplot(numChannels, 1, i);
        plot(To, icasig(i, :), 'LineWidth', 1.2, 'Color', [0.2 0.4 0.8]);
        title(['IC ' num2str(i)], 'FontSize', 10, 'FontWeight', 'bold');
        if i == numChannels
            xlabel('Time (s)', 'FontSize', 9);
        end
        ylabel('Amplitude', 'FontSize', 9);
        xlim([0 To(end)]);
        set(gca, 'FontSize', 8, 'Color', 'white');
        box on;
    end
    sgtitle(['**' array_name '** : Select Components to KEEP (See Dialog Box)'], 'FontSize', 14, 'FontWeight', 'bold');
    
    % --- 4. Interactive Input Prompt (Using Dialog Box) ---
    
    % Display a list of all component numbers (1 to 5)
    all_components_str = strjoin(arrayfun(@num2str, 1:numChannels, 'UniformOutput', false), ', ');
    
    % Set up dialog box parameters
    prompt = {['Based on the plots, enter the IC numbers to KEEP (e.g., 1,2,5 or 1:3):']};
    title_str = ['IC Selection for ' array_name];
    dims = [1 50];
    
    % Use inputdlg to create a pop-up input box that stays on top of the figure
    user_input = inputdlg(prompt, title_str, dims, {all_components_str});
    
    if isempty(user_input)
        % User pressed Cancel
        warning('Selection cancelled. Defaulting to all components [1:5].');
        componentsToUse = 1:numChannels;
    else
        input_str = user_input{1};
        % Process the input string
        try
            componentsToUse = eval(['[' input_str ']']);
            componentsToUse = unique(componentsToUse); % Remove duplicates
            
            % Input validation (optional, but good practice)
            if any(componentsToUse < 1) || any(componentsToUse > numChannels)
                 warning('Component indices out of range [1 to %d]. Defaulting to all components [1:%d].', numChannels, numChannels);
                 componentsToUse = 1:numChannels;
            end
        catch
            warning('Invalid component selection input. Defaulting to all components [1:%d].', numChannels);
            componentsToUse = 1:numChannels;
        end
    end
    
    close(gcf); % Close the selection figure
    
    % --- 5. Reconstruction ---
    
    % Select only the specified components for reconstruction
    icasig_selected = icasig(componentsToUse, :);
    A_selected = A(:, componentsToUse);
    
    % Reconstruct data using only selected components
    Data_reconstructed = A_selected * icasig_selected;
    
    % Add mean back to reconstructed data
    Data_reconstructed = Data_reconstructed + mean(Data_matrix, 2);
    
    % --- 6. Metrics and Final Output ---
    
    % Calculate MSE
    mse = mean((Data_matrix(:) - Data_reconstructed(:)).^2);
    disp(['Mean Squared Error of Reconstruction: ' num2str(mse)]);
    
    % Store results in a structure
    results.Data_reconstructed = Data_reconstructed;
    results.selected_components = componentsToUse;
    results.mse = mse;
    
end