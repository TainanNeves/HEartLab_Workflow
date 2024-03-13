function [MFFTi,Sffti,fstep]=f_DF_electric(DATA,fs, freq_up, freq_down)
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
%   - DATA: A matrix containing the input data, where each row represents a different signal.
%   - fs: The sampling frequency of the data.
%   - freq_up: The upper frequency limit for dominant frequency estimation.
%   - freq_down: The lower frequency limit for dominant frequency estimation.
%
% Outputs:
%   - MFFTi: An array containing the estimated dominant frequencies for each input signal.
%   - Sffti: A matrix of spectral information, where each row corresponds to an input signal.
%   - fstep: The frequency step used in the FFT.
%
% Example Usage:
%   [MFFTi, Sffti, fstep] = f_DF_electric(DATA, fs, 8, 3.5);
%
% Notes:
%   - The function applies a Hamming window to the input signal and performs zero-padding to improve frequency resolution.
%   - The dominant frequency is estimated within the specified frequency range [freq_down, freq_up].
%   - The function returns an array of estimated dominant frequencies for each input signal. If no dominant frequency is found within the specified range, the value will be set to 0.
%
% Data:
%   - 24/10/2023
%   - 13/03/2024 - Filtering in 20 Hz and absolute of the peaks before fft

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%To AF
% freq_up=100;
% freq_down=1;
%To Sinus
% freq_up = 8;
% freq_down = 2.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%To plot the frequency map from an electrode use this script:
% electrode = ;
% plot(fstep:fstep:freq_up,Sffti(electrode,:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

samples = length(DATA(1, :));
time = samples/fs;
factor = 5; 
sizefft = factor * fs * time; % factor * samples
fstep = fs / sizefft;
H = hamming(sizefft/factor);
[a, ~] = size(DATA);

% Analyse over absolute values
DATA_2=abs(DATA);
% Filtering
Lw=20; %Cut off frequency
[z,p,k] = butter(5,Lw/(fs/2),'low'); % Design the filter
[sos, g] = zp2sos(z,p,k); % Change the format
DATA_3 = zeros(size(DATA_2));
for i = 1:size(DATA_2, 1)
    DATA_3(i, :) = filtfilt(sos, g, DATA_2(i, :)); % Apply the filter to each row
end

% FFT
for i = 1:a
    x = DATA_3(i,:)';
    xx = detrend(x);   
    Sfft = abs(fft(xx .* H, sizefft)); % Zero Padding
    Sfft = Sfft .* Sfft / length(Sfft);  % for power 
    Sffti(i, :) = Sfft(1:floor(freq_up*(1/fstep)));     
    [~, F] = max(Sfft(floor(freq_down*(1/fstep)):floor(freq_up*(1/fstep))));  
    Fsam = (F + floor(freq_down*(1/fstep)) - 1);  
    Fhz = Fsam * fstep;    
    if (Fhz >= freq_down && Fhz <= freq_up)       
        MFFTi(i) = Fhz;
    else
        MFFTi(i) = 0;
    end
end
