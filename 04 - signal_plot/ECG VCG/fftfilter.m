function filtersignal = fftfilter(ori_signal,lowband,highband)
% This is a program for bandpass filtering using FFT 
% Author: Hui Yang
% Affiliation: 
       %The Pennsylvania State University
       %310 Leohard Building, University Park, PA
       %Email: yanghui@gmail.com
% ori_signal: input time series
% lowband: low band 
% highband: high band
% Note: It will be better for the length of input time series to be 2 power
% This demo shows an animated plot of vectorcardiogram. The color values are associated
% with the moving speed of VCG vectors. 
% If you find this demo useful, please cite the following paper:
% [1] H. Yang, S. T. S. Bukkapatnam, and R. Komanduri, �Spatiotemporal Representation of 
% Cardiac Vectorcardiogram (VCG) Signal,� Biomedical Engineering Online, Vol.11, No. 16,
% 2012, DOI: 10.1186/1475-925X-11-16
% [2] H. Yang*, S. T. S. Bukkapatnam, L. Trung and R. Komanduri, �Identification of 
% myocardial infarction (MI) using spatio-temporal heart dynamics,� Medical Engineering and Physics, 
% Vol. 34, No. 4, p485-497, 2011,  DOI: 10.1016/j.medengphy.2011.08.009

fs=4000;
passband(1) = lowband;
passband(2) = highband;

N = length(ori_signal);
y = fft(ori_signal);

lowicut = round(passband(1)*N/fs);
lowmirror = N-lowicut+2;
highicut =  round(passband(2)*N/fs);
highmirror = N-highicut+2;

y([1:(lowicut-1) (lowmirror+1):end])=0;
y((highicut+1):(highmirror-1))=0;

filtersignal = ifft(y);
