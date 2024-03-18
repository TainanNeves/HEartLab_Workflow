% wavelet_filter - Perform wavelet decomposition and reconstruction on signals
%
% Author: Angélica Quadros
% HeartLab
% 2024
%
% This function performs wavelet decomposition and reconstruction on signals
% from electrodes using different wavelet types and sets of reconstruction levels.
% It decomposes each signal using the specified wavelet type and number of levels,
% and then reconstructs the signal using the selected reconstruction levels.
% Optionally, it can save the MRA matrices for each electrode and configuration.
%
% Usage:
%   wavelet_filter(D_EL, fs, electrodes, waveletTypes, numLevels, reconstructionLevelsSets, saveMRA)
%
% Inputs:
%   - D_EL: Signals stored in a STRUCT, as exported from extraction_filtering function.
%   - fs: Sampling frequency of the signals.
%   - electrodes: Vector containing the electrode numbers, can be one or more.
%   - waveletTypes: Cell array of strings containing the wavelet types to be used.
%   - numLevels: Number of decomposition levels.
%   - reconstructionLevelsSets: Cell array of vectors specifying which levels will be used for reconstruction.
%   - saveMRA (optional): Boolean indicating whether to save the MRA matrices (decomposition levels variables). Default is false.
%
% Outputs:
%   - If saveMRA is true, the function creates MRA matrices and stores them as variables
%     in the workspace. The variables are named according to the electrode, wavelet type,
%     set of reconstruction levels, and the suffix "_MRA". For example, "MRA_db4_Levels_1_5_Electrode_1".
%   - Reconstructed signals named according to the wavelet type, set of reconstruction levels,
%     and the suffix "_Reconstructed". For example, "Signal_reconstructed_db4_Levels_1_5".
%
%
% This function performs wavelet decomposition and reconstruction on signals
% from electrodes using different wavelet types and sets of reconstruction levels.
% It decomposes each signal using the specified wavelet type and number of levels,
% and then reconstructs the signal using the selected reconstruction levels.
% Optionally, it can save the MRA matrices for each electrode and configuration.
%
% Example:
%   wavelet_filter(D_EL, 3600, [1:16], {'db4', 'coif1'}, 10, {[6,7], [7,8,9]});
%
% This example performs wavelet decomposition and reconstruction on signals from
% electrodes 1 to 16 using 'db4' and 'coif1' wavelet types, 10 decomposition levels,
% and using two sets of reconstruction levels: [7,8] and [7,8,9].
%
%
%
function wavelet_filter(D_EL, fs, electrodes, waveletTypes, numLevels, reconstructionLevelsSets, saveMRA)
    
    % check if saveMRA is provided, if not, set to false
    if nargin < 7
        saveMRA = false;
    end
    
    % Iterate over each set of reconstruction levels
    for i = 1:numel(reconstructionLevelsSets)
        currentLevels = reconstructionLevelsSets{i};
        
        % Get the string representation of the current reconstruction levels
        levelsString = strjoin(arrayfun(@num2str, currentLevels, 'UniformOutput', false), '_');
        
        % Iterate over each wavelet type
        for j = 1:numel(waveletTypes)
            currentWavelet = waveletTypes{j};
            
            % Initialize variables to store signals and MRA matrices for the current set and wavelet type
            signals_currentLevelsWavelet = [];
            mra_currentLevelsWavelet = [];

            % Iterate over each electrode value
            for k = 1:numel(electrodes)
                electrode = electrodes(k);
                
                % Extract the electrode number from the loop index
                electrodeName = ['Electrode_' num2str(electrode)];

                % Original signal
                Signal = D_EL.Data(electrode, :);

                % Applying conversion factor
                factor_conv = 0.1949999928474426;
                Signal = Signal * factor_conv;

                % Decomposition using modwt
                wt = modwt(Signal, currentWavelet, numLevels);

                % Construct MRA matrix using modwtmra
                mra_temp = modwtmra(wt, currentWavelet);

                % Concatenate MRA matrix with previous ones
                if isempty(mra_currentLevelsWavelet)
                    mra_currentLevelsWavelet = mra_temp;
                else
                    mra_currentLevelsWavelet = cat(1, mra_currentLevelsWavelet, mra_temp);
                end
                            
                % Assign MRA matrix to a variable named according to the electrode if saveMRA is true
                if saveMRA
                    variable_name_mra = ['MRA_' currentWavelet '_Levels_' levelsString '_Electrode_' num2str(electrode)];
                    assignin('caller', variable_name_mra, mra_temp);
                end

                % Concatenate signals from all electrodes
                if isempty(signals_currentLevelsWavelet)
                    signals_currentLevelsWavelet = Signal;
                else
                    signals_currentLevelsWavelet = cat(1, signals_currentLevelsWavelet, Signal);
                end
            end

            % Logical array for selecting reconstruction levels
            currentRange = currentLevels;
            levelForReconstruction = false(1, numLevels);
            levelForReconstruction(currentRange(1):currentRange(2)) = true;

            % Signal reconstruction by summing the selected rows
            Signal_reconstructed = sum(mra_currentLevelsWavelet(levelForReconstruction, :), 1);

            % Save the reconstructed signal for the current levels and wavelet type
            variable_name_signals = ['Signal_reconstructed_' currentWavelet '_Levels_' levelsString];
            assignin('caller', variable_name_signals, signals_currentLevelsWavelet);
        end
    end
end