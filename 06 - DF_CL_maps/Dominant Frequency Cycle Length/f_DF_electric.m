function [MFFTi, Sffti, fstep] = f_DF_electric(DATA, fs, freq_up, freq_down)
% Function: f_DF_electric
% Authors: Angélica Quadros, Jimena Paredez, and Tainan Neves. 
% Affiliation: HEartLab - UFABC
%
% Description:
%   This function computes the dominant frequency of electric signal data
%   using a Fast Fourier Transform (FFT) approach. 
%   The dominant frequency is the frequency component with the highest 
%   amplitude in the given signal. It can be useful in various applications, 
%   such as analyzing cardiac signals or other electric signal data.
%
% Syntax:
%   [MFFTi, Sffti, fstep] = f_DF_electric(DATA, fs, freq_up, freq_down)
%
% Inputs:
%   - DATA: A 3D matrix containing the input data, where the first dimension represents the x-coordinate, 
%           the second dimension represents the y-coordinate, and the third dimension represents the time signal.
%   - fs: The sampling frequency of the data.
%   - freq_up: The upper frequency limit for dominant frequency estimation.
%   - freq_down: The lower frequency limit for dominant frequency estimation.
%
% Outputs:
%   - MFFTi: A matrix containing the estimated dominant frequencies for each (x, y) coordinate.
%   - Sffti: A 3D matrix of spectral information, where each slice corresponds to an (x, y) coordinate.
%   - fstep: The frequency step used in the FFT.
%
% Example Usage:
%   [MFFTi, Sffti, fstep] = f_DF_electric(DATA, fs, 8, 3.5);
%
% Notes:
%   - The function applies a Hamming window to the input signal and performs zero-padding to improve frequency resolution.
%   - The dominant frequency is estimated within the specified frequency range [freq_down, freq_up].
%   - The function returns a matrix of estimated dominant frequencies for each (x, y) coordinate. If no dominant frequency is found within the specified range, the value will be set to 0.
%
% Data:
%   - 24/10/2023
%   - 13/03/2024 - Filtering in 20 Hz and absolute of the peaks before fft
%   - 05/07/2024 - Code adapted to the new matrixial notation

samples = size(DATA, 3);
time = samples / fs;
factor = 5; 
sizefft = factor * fs * time; % factor * samples
fstep = fs / sizefft;
H = hamming(sizefft / factor);

[x_size, y_size, ~] = size(DATA);

% Initialize output matrices
MFFTi = zeros(x_size, y_size);
Sffti = zeros(x_size, y_size, floor(freq_up * (1 / fstep)));

% Absolute Value Calculation
DATA_2 = abs(DATA);

% Filtering
Lw = 20; % Cut off frequency
[z, p, k] = butter(5, Lw / (fs / 2), 'low'); % Design the filter
[sos, g] = zp2sos(z, p, k); % Change the format
DATA_3 = zeros(size(DATA_2));
for i = 1:x_size
    for j = 1:y_size
        DATA_3(i, j, :) = filtfilt(sos, g, squeeze(DATA_2(i, j, :))); % Apply the filter to each (x, y) coordinate
    end
end

% FFT
for i = 1:x_size
    for j = 1:y_size
        x = squeeze(DATA_3(i, j, :));
        xx = detrend(x);   
        Sfft = abs(fft(xx .* H, sizefft)); % Zero Padding
        Sfft = Sfft .* Sfft / length(Sfft);  % for power 
        Sffti(i, j, :) = Sfft(1:floor(freq_up * (1 / fstep)));     
        [~, F] = max(Sfft(floor(freq_down * (1 / fstep)):floor(freq_up * (1 / fstep))));  
        Fsam = (F + floor(freq_down * (1 / fstep)) - 1);  
        Fhz = Fsam * fstep;    
        if (Fhz >= freq_down && Fhz <= freq_up)       
            MFFTi(i, j) = Fhz;
        else
            MFFTi(i, j) = 0;
        end
    end
end
end
