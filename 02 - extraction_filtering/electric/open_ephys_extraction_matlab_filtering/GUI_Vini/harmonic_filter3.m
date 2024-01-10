function fData = harmonic_filter3(data, frequency, order, Fs)

forder = 20;
Q = 25;
corr = 1.3;
fData = data;

for i = 1:order
    h1 = fdesign.notch(forder, frequency+frequency*(i-1), Q*i*corr, Fs);
    d = design(h1);
    
    sos = d.sosMatrix;
    scale = d.ScaleValues;

    fData = filtfilt(sos, scale, fData);
end

