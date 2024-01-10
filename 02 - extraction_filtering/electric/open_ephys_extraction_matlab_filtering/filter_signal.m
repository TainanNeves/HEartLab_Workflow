function Data_filtered = filter_signal(Lf, Hf,DATA,Fs)
    
    %% Normalization (Uncoment if you want the normalized signal)
    % Normalize 
    %Data = (DATA)/(max(DATA)-min(DATA));
    %
    
    %% PrePros
    Data = DATA;
    % Add zeros to the beggining and end to avoid distortions
    support_length = 10000;
    add = ones(1, support_length);
    data = cat(2, add*Data(1), Data, add*Data(end));
    %Timestamps = D.Timestamps;
    
    % Initial Time and Duration of the graphs
    %T0 = Timestamps(1);
    %TD = 2;
    
    %% Remove baseline wadering via butterworth HP filter
    Hw = Lf; % Cutoff frequency for the high-pass filter;
    [z,p,k] = butter(20,Hw/(Fs/2),'high'); % Design the filter
    [sos, g] = zp2sos(z,p,k); % Change the format
    
    fhdata = filtfilt(sos, g, data); % Apply the filter

    % Plot the data
    % figure;
    % plot(Timestamps, data(support_length+1: end-support_length), ...
    %     "DisplayName",'Original', LineWidth = 1);
    % hold on;
    % plot(Timestamps, fhdata(support_length+1: end-support_length), ...
    %     "DisplayName",'Filtered', LineWidth = 1);
    % title('Removing baseline wandering HP filter (Butterorth)', 'FontSize', 18);
    % set(legend, 'FontSize', 16);
    % xlim([T0 T0+TD]);
    clear z p k sos g Hw % Clear unused variables
    
    %% Remove high frequency noise via butterworth LP filter
    Lw = Hf; % Cutoff frequency for the Low-pass filter;
    [z,p,k] = butter(10,Lw/(Fs/2),'low'); % Design the filter
    [sos, g] = zp2sos(z,p,k); % Change the format
    
    fldata = filtfilt(sos, g, fhdata); % Apply the filter
    
    % % Plot the data
    % figure;
    % plot(Timestamps, fhdata(support_length+1: end-support_length), ...
    %     "DisplayName",'Original', LineWidth = 1);
    % hold on;
    % plot(Timestamps, fldata(support_length+1: end-support_length), ...
    %     "DisplayName",'Filtered', LineWidth = 1);
    % title('Removing high frequency noise with LP filter (Butterorth)', 'FontSize', 18);
    % set(legend, 'FontSize', 16);
    % xlim([T0 T0+TD]);
    clear z p k sos g Hw % Clear unused variables
    
    %% Filter harmonics via notch filter
    Q = 80; % Quality factor defined as the ratio of the center frequency to the 3 dB bandwidth.
    corr = 1.3; % Correction factor of Q; 
                % I noticed the signal got more distorted as the frequency grew
                % so corr ajusts it. HAS NO THEORETICAL BASIS.       
    frequency = 60; % Powerline frequency.
    ast = 120; % Stopband attenuation of the filter, specified as a positive scalar in dB.
    order = 6; % Order of the harmonics that will be filtered. 
               % Ex. order=3 will filter 60, 120 and 180 Hz.
    
    hdata = fldata;
    % Run the filter for every harmonic 
    for i = 1:order
        % Design the filter
        h1 = fdesign.notch('N,F0,Q,Ast', 6, frequency+frequency*(i-1), Q*i*corr, ast, Fs);
        d = design(h1);
        
        % Select desired values
        sos = d.sosMatrix;
        scale = d.ScaleValues;
        
        % Apply the filter
        hdata = filtfilt(sos, scale, hdata);
    end
    % Plot the data
    % figure;
    % plot(Timestamps, fldata(support_length+1: end-support_length), ...
    %     "DisplayName",'Original', LineWidth = 1);
    % hold on;
    % plot(Timestamps, hdata(support_length+1: end-support_length), ...
    %     "DisplayName",'Filtered', LineWidth = 1);
    % title('Removing powerline noise with notch filter', 'FontSize', 18);
    % set(legend, 'FontSize', 16);
    % xlim([T0 T0+TD]);
    clear Q corr frequency forder d order sos scale i h1 ast;
    
      %% Filter via Savitzky-Golay Filter varying frame length (This filter makes the signal smoother.)
%      % But we missed quick signals like fibrillations. 
%      % Do not use, even if you know what you are doing! 
% 
%     sgorder = 2; % Polynomial order of the filter
%     framelen = 101; % Frame length (must be odd)
%     iterations = 4; % Iterations to run the SG filter
%     if mod(framelen, 2) == 0 % Check if it's odd, add one if it isn't;
%         framelen = framelen + 1;
%     end
%     
%     % Filter the signal
%     SGdata = hdata;
%     for i = 1:iterations
%         SGdata = sgolayfilt(SGdata, sgorder, framelen);
%     end
%     
%     % Plot the data
%     % figure;
%     % plot(Timestamps, hdata(support_length+1:end-support_length), ...
%     %     'DisplayName', 'Original', 'LineWidth', 1);
%     % hold on;
%     % plot(Timestamps, SGdata(support_length+1:end-support_length), ...
%     %     'DisplayName', 'Filtered', ...
%     %     'LineWidth', 1);
%     % txt = strcat('SG-filter, frame length = ', num2str(framelen), ', order = ', num2str(sgorder), ...
%     %     ', iterations = ', iterations);
%     % title(txt, 'FontSize', 18);
%     % set(legend, 'FontSize', 14);
%     % xlim([T0 T0+TD]);
%     clear txt sgorder framelen add iterations;
    
    %% Periodogram
    % Adjusting the data to original size
    data = data(support_length+1:end-support_length);
    fldata = fldata(support_length+1:end-support_length);
    fhdata = fhdata(support_length+1:end-support_length);
    hdata = hdata(support_length+1:end-support_length);
    %SGdata = SGdata(support_length+1:end-support_length); 
    
%     % Original data
%     [pdata,fsdata]=periodogram(data,[],length(data),Fs);
%     % Filtered data
%     [sgData, fsgData]=periodogram(SGdata, [], length(SGdata), Fs);
    
    %% Plots
    
    % figure;
    % title('Periodogram');
    % hold on;
    % tiledlayout(2, 1);
    % ax1 = nexttile;
    % plot(fsdata,10*log10(abs(pdata)), 'r' );
    % title('Unfiltered data');
    % ylabel('dB')
    % 
    % ax2 = nexttile;
    % plot(fsgData,10*log10(abs(sgData)),'b');
    % title('Filtered data');
    % xlabel('Frequency (Hz)');
    % ylabel('dB');
    % 
    % linkaxes([ax1 ax2], 'xy');
    
    %%Defining Output
    Data_filtered = hdata;
end
