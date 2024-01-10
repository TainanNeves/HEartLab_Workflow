%% Baseline filter
%This filter works with a high pass filter 
%DATA= data in an array of x*y*z (example: 125*125*5000)
%fc1= cut frecuency high pass filter 
%fs= sample frequency 
%n= filter's order
%less= number of samples to be eliminated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [R]=f_fil(DATA,fc1,fs,n,less)
[a b c] = size(DATA);[z,p] = butter(n,(fc1/(fs/2)),'high');
ii=0;jj=0;
for ii=1:a
        for jj=1:b
           signal=squeeze(DATA(ii,jj,:));
           R(ii,jj,:)=filtfilt(z,p,double(signal));
        end
end
[a b c]=size(R);%1 segundo es 500 samples 
R=R(:,:,less:(c-less));%(c-less)
