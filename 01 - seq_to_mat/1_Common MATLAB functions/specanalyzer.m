function [f,Y,M]=specanalyzer(X,Fs,Min,Max,type,Display)

% Version 2.0
% specanalyzer(X,Fs,Min, Max,type, Display)
% This function plots the frequency spectrum of a signal X without 0Hz value, that is a DC. 
% Fs is the sampling frequency 
% Min and Max specify the frequency range that will be plotted
% type is either 'normal' or 'log' and defined the scale of the y axis

X = double(X);

if nargin<2
    display('Not enough input arguments')
    return
end
if nargin>6
    display('Too many input arguments')
    return
end

Y = abs(fft(X));
Y=Y(1:round(length(Y)/2));

f = linspace(0,Fs/2,length(Y));  %increment is the smallest detectable frequency

f(1)=[];Y(1)=[]; % % removing DC
Y=Y(find(f>=Min));
f=f(find(f>=Min));
Y=Y(find(f<=Max));
f=f(find(f<=Max));

if Display == 1
    if strcmp(type,'normal')
        plot( f, Y);title('Single-Sided Amplitude Spectrum of y(t)')
        xlabel('Frequency (Hz)');ylabel('|Y(f)|')
    elseif strcmp(type,'log')
        plot( f,20*log10( Y) )
        title('Single-Sided Amplitude Spectrum of y(t)')
        xlabel('Frequency (Hz)')
        ylabel('|Y(f)| (dB)')
    else
    display('wrong type was specified')
    end
end
M=f(find( Y == max(Y) ));
display(['peak freq = ', num2str( M )] )