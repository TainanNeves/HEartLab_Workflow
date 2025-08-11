function Data_filtered = filter_signal(Lf, Hf, DATA, Fs)
    % This function filters the input signal DATA using:
    % 1. Notch filters to remove powerline noise and harmonics
    % 2. High-pass FIR filter to remove baseline wandering
    % 3. Low-pass FIR filter to remove high-frequency noise

    % Fixed parameters
    support_length = 10000; % Padding length to avoid edge effects
    powerline_freq = 60;   % Powerline frequency (Hz)
    harmonics = 12;        % Number of harmonics to filter

    % Apply notch filters
    Data = apply_notch_filters(DATA, powerline_freq, Fs, harmonics);

    % Pad signal by replicating edge samples
    Data = [repmat(Data(1), 1, support_length), Data, repmat(Data(end), 1, support_length)];

    % Apply FIR filters
    Data = apply_fir_filter(Data, Lf, Fs, 'high', 1000); % High-pass
    Data = apply_fir_filter(Data, Hf, Fs, 'low', 500);   % Low-pass

    % Remove padding and return
    Data_filtered = Data(support_length + 1 : end - support_length);
end

function filtered_data = apply_fir_filter(data, cutoff, Fs, type, order)
    % Applies FIR high-pass or low-pass filter
    % Design filter
    if strcmp(type, 'high')
        b = fir1(order, cutoff / (Fs / 2), 'high', hamming(order + 1));
    elseif strcmp(type, 'low')
        b = fir1(order, cutoff / (Fs / 2), 'low', hamming(order + 1));
    else
        error('Filter type must be "high" or "low"');
    end
    % Zero-phase filtering
    filtered_data = filtfilt(b, 1, data);
end

function filtered_data = apply_notch_filters(data, freq, Fs, order)
    % Applies IIR notch filter at fundamental frequency and harmonics
    filtered_data = data;
    for i = 1:order
        notch_freq = freq * i;
        % Calculate and normalize bandwidth
        bw = notch_freq / 35; % Fixed bandwidth ratio Q = 35
        normalized_bw = min(bw / (Fs / 2), 0.99); % Ensure < 1
        % Design and apply notch filter
        [b, a] = iirnotch(notch_freq / (Fs / 2), normalized_bw);
        filtered_data = filtfilt(b, a, filtered_data);
    end
end