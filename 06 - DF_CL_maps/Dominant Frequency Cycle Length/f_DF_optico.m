function [MFFTi, Sffti, fstep] = f_DF_optico(DATA_A, fs, UP, DOWN)
% Function: f_DF_optico
%
% Description:
%   This function computes the dominant frequency of optical signal data
%   using a Fast Fourier Transform (FFT) approach after calculating the
%   absolute value of the signal.
%
% Syntax:
%   [MFFTi, Sffti, fstep] = f_DF_optico(DATA_A, fs, UP, DOWN)
%
% Inputs:
%   - DATA_A: 3D matrix containing the input optical signal data.
%   - fs: Sampling frequency of the data.
%   - UP: Upper frequency limit for dominant frequency estimation.
%   - DOWN: Lower frequency limit for dominant frequency estimation.
%
% Outputs:
%   - MFFTi: Matrix containing the estimated dominant frequencies for each (i, j) coordinate.
%   - Sffti: 3D matrix of spectral information, where each slice corresponds to an (i, j) coordinate.
%   - fstep: Frequency step used in the FFT.
%
% Example Usage:
%   [MFFTi, Sffti, fstep] = f_DF_optico(DATA_A, 4000, 8, 3.5);
%

[a, b, ~] = size(DATA_A);

% Declare necessary values
samples = size(DATA_A, 3);
time = samples / fs;
factor = 5; 
sizefft = factor * fs * time;
fstep = fs / sizefft;
H = hamming(sizefft / factor);
freq_up_AF = UP;
freq_down_AF = DOWN;

% Initialize output matrices
MFFTi = zeros(a, b);
Sffti = zeros(a, b, floor(freq_up_AF * (1 / fstep)));

% Calculate dominant frequency
for i = 1:a
    for j = 1:b
        x = squeeze(DATA_A(i, j, :));
        
        % Calculate absolute value of the signal
        x = abs(x);
        
        % FFT
        Sfft = abs(fft(x .* H, sizefft)); % Zero Padding
        Sfft = Sfft .* Sfft / length(Sfft);  % for power 
        Sffti(i, j, :) = Sfft(1:floor(freq_up_AF * (1 / fstep)));
        
        % Find dominant frequency
        [~, F] = max(Sfft(floor(freq_down_AF * (1 / fstep)):floor(freq_up_AF * (1 / fstep))));  
        Fsam = (F + floor(freq_down_AF * (1 / fstep)) - 1);  
        Fhz = Fsam * fstep;
        
        if (Fhz >= freq_down_AF && Fhz <= freq_up_AF)       
            MFFTi(i, j) = Fhz;
        else
            MFFTi(i, j) = 0;
        end   
    end 
end
end
