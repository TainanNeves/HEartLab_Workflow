function fdata = filterdata(data, checkboxValue1_4, har_Q, har_freq, har_order,...
    har_ast, hp_freq, lp_freq, sg_order, sg_framelen, Fs)

%% Add placeholder values to the beggining and end to avoid distortions
support_length = 10000;
add = ones(1, support_length);
fdata = cat(2, add*data(1), data, add*data(end));

forder = 30;

%% Harmonic Filter
Q = har_Q; % Quality factor defined as the ratio of the center frequency to the 3 dB bandwidth.
corr = 1.3; % Correction factor of Q; 
            % I noticed the signal got more distorted as the frequency grew
            % so corr ajusts it. HAS NO THEORETICAL BASIS.
frequency = har_freq; % Base harmonic frequency
order = har_order; % Order of the harmonics that will be filtered. 
ast = har_ast; % Stopband attenuation of the filter, specified as a positive scalar in dB.

if checkboxValue1_4(1) == 1
    for i = 1:order
        % Design the filter
        h1 = fdesign.notch('N,F0,Q,Ast', forder, frequency+frequency*(i-1), Q*i*corr, ast, Fs);
        d = design(h1);

        % Select desired values
        sos = d.sosMatrix;
        scale = d.ScaleValues;
        
        % Apply the filter
        fdata = filtfilt(sos, scale, fdata);
    end
end

%% Highpass Filter
if checkboxValue1_4(3) == 1
    Hw = hp_freq; % Cutoff frequency

    % Design filter
    [z,p,k] = butter(10,Hw/(Fs/2),'high');
    [sos, g] = zp2sos(z,p,k);
    
    % Apply filter
    fdata = filtfilt(sos, g, fdata);
end

%% Lowpass Filter
if checkboxValue1_4(2) == 1
    Lw = lp_freq; % Cutoff frequency

    % Design filter 
    [b, a] = butter(10, Lw/(Fs/2));

    % Apply filter
    fdata = filtfilt(b, a, fdata);
end

% Savitzky-Golay Filter
if checkboxValue1_4(4) == 1                
    sgorder = sg_order; % Polynomial order of the filter
    framelen = sg_framelen; % Frame length (must be odd)
    if mod(framelen, 2) == 0
        framelen = framelen + 1;
    end

    % Filter the signal
    fdata = sgolayfilt(fdata, sgorder, framelen);
end

%% Adjusting the data to original size
fdata = fdata(support_length+1:end-support_length);

end