function [MFFTi,Sffti,fstep]=f_DF_optico(DATA_A,fs,UP,DOWN)
% DATA_A=O_SIN_LA;
% Fs=4000;
% i=62;
% j=133;

[a b c]=size(DATA_A);

%declaro os valores que eu precisso
samples = length(DATA_A(1,1, :));
time = samples/fs;
factor=5; 
sizefft=factor*fs*time;
fstep=fs/sizefft;
H=hamming (sizefft/factor);
freq_up_AF=UP;
freq_down_AF=DOWN;


%calculou a frequencua
for i=1:a
    for j=1:b
        x=squeeze(DATA_A(i,j,:));
        xx=detrend(x);
        avg_value = mean (xx);
        if  (avg_value==0)
            MFFTi(i,j)=0;
           %Sffti (i,j,:)= 0;
           %Sfftf (i,j,:)= 0;
            Sfft=abs(fft(xx.*H, sizefft)); % Zero Padding
            Sfft=Sfft.*Sfft/length(Sfft);  % for power 
            Sffti (i,j,:)=Sfft(1:floor(freq_up_AF*(1/fstep))); 
        else
%             keyboard
            Sfft=abs(fft(xx.*H, sizefft)); % Zero Padding
            Sfft=Sfft.*Sfft/length(Sfft);  % for power 
            Sffti (i,j,:)=Sfft(1:floor(freq_up_AF*(1/fstep))); 
%           b = length(Sffti (i,j,:));
            %figure;plot((1:b)*fstep, squeeze(Sfft_O_8S(62,133,:)));hold on;
            [P,F] = max(Sfft(floor(freq_down_AF*(1/fstep)):floor(freq_up_AF*(1/fstep))));  
            Fsam=(F+floor(freq_down_AF*(1/fstep))-1);  
             Fhz=Fsam*fstep;
                if  (Fhz>=freq_down_AF&&Fhz<=freq_up_AF)       
                    MFFTi(i,j)=Fhz;
                else
                    MFFTi(i,j)=0;
                end   
        end 
    end
end

%