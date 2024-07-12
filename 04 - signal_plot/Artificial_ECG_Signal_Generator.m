%CÓDIGO HEartLab - Artificial ECG Signal Generator - Versão 1 - 04/07/2024
%Author - José Junior

%Clear the Command Window and Workspace
clc; clear;

%Constants for Offset and Exponent Values
NP = 240;  NQ = 9000;  NR = 1800;  NS = 3000;  NT = 138;    NTT = 138;
DP = 0.5;  DQ = 0.08;  DR = 0;     DS = -0.1;  DT = -0.75;  DTT = -0.6;

%Constants for Frequency Values and Time Conversion Factor
k = 1.25;  CD = 0.4;

%Test Vector to Traverse 12 ECG Leads
VTT = 1:12;

%Vectors to Store Function Point Values
VZAs = cell(1, 12); %Vector to Store Amplitude
VZTs = cell(1, 12); %Vector to Store Time of Amplitude Recording

%Main FOR Loop - Traverses the 12 ECG Leads
for i = 1:12
    
    vt = VTT(i); %Variable that Automatically Selects One of the 12 ECG Leads
    
    switch vt %Test the Variable Value and Assign Amplitude Values (A)
        case 1 %ECG Lead DI
           AP = 0.05;  AQ = -0.05;   AR = 0.7;   AS = -0.002; AT = 0.001; ATT = 0.03;
        case 2 %ECG Lead DII
           AP = 0.1;   AQ = -0.1;    AR = 0.8;   AS = -0.2;   AT = 0.1;   ATT = 0.05;
        case 3 %ECG Lead DIII
           AP = 0.07;  AQ = -0.075;  AR = 0.5;   AS = -0.116; AT = 0.05;  ATT = 0.025;
        case 4 %ECG Lead aVR 
           AP = -0.08; AQ = 0.1;     AR = -0.8;  AS = -0.1;   AT = -0.1;  ATT = -0.05;
        case 5 %ECG Lead aVL 
           AP = 0.04;  AQ = 0;       AR = 0.25;  AS = -0.1;   AT = 0.05;  ATT = 0.025;
        case 6 %ECG Lead aVF
           AP = 0.075; AQ = -0.0416; AR = 0.55;  AS = -0.116; AT = 0.15;  ATT = 0.0625;
        case 7 %ECG Lead V1 
           AP = 0.1;   AQ = 0;       AR = 0.18;  AS = -0.33;  AT = 0.05;  ATT = 0.025;
        case 8 %ECG Lead V2
           AP = 0.1;   AQ = 0;       AR = 0.28;  AS = -0.5;   AT = 0.12;  ATT = 0.06;
        case 9 %ECG Lead V3
           AP = 0.1;   AQ = 0;       AR = 0.453; AS = -0.5;   AT = 0.1;   ATT = 0.05;
        case 10 %ECG Lead V4
           AP = 0.1;   AQ = 0;       AR = 0.653; AS = -0.27;  AT = 0.1;   ATT = 0.05;
        case 11 %ECG Lead V5
           AP = 0.1;   AQ = -0.05;   AR = 0.8;   AS = -0.2;   AT = 0.05;  ATT = 0.025;
        case 12 %ECG Lead V6
           AP = 0.1;   AQ = -0.075;  AR = 0.6;   AS = -0.002; AT = 0.05;  ATT = 0.025;
    end %End of Switch

%Define the Time Variable and Command Auxiliary Variable
    t = 0;  vc = 0;

%Define the Vectors to Store Amplitude and Time
    VZA = zeros(1,10000);  VZT = zeros(1,10000);

%FOR Loop for Calculating Function Values

    for j = 0:0.01:30 %Increment of 0.01 corresponds to the Sampling Time Value in seconds.
        
        t = t + j; %Assign Sampling Value
        t = t * CD; %Convert Time Value
        vc = vc + 1; %Command Auxiliary Variable

        FP = sin((k * j) + DP);   FP = (FP)^(NP);    FP = FP * AP; %Calculation of P Wave
        FQ = sin((k * j) + DQ);   FQ = (FQ)^(NQ);    FQ = FQ * AQ; %Calculation of Q Wave
        FR = sin((k * j) + DR);   FR = (FR)^(NR);    FR = FR * AR; %Calculation of R Wave
        FS = sin((k * j) + DS);   FS = (FS)^(NS);    FS = FS * AS; %Calculation of S Wave
        FT = sin((k * j) + DT);   FT = (FT)^(NT);    FT = FT * AT; %Calculation of T Wave
        FTT = sin((k * j) + DTT); FTT = (FTT)^(NTT); FTT = FTT * ATT; %Calculation of TT Wave

        FF = FP + FQ + FR + FS + FT + FTT; %Sum of Amplitudes P, Q, R, S, T, and TT

%Store the Calculated Amplitude and Time Values

        VZA(vc) = FF; %Vector for Amplitude    
        VZT(vc) = t;  %Vector for Time

    end %End of FOR Loop for Calculation

%Store Amplitude and Time Values for the 12 ECG Leads

    VZAs{i} = VZA; %Store Amplitude Values in Cell i
    VZTs{i} = VZT; %Store Time Values in Cell i

end %End of Main FOR Loop

%Clear Variables that Will No Longer Be Used - Part 1
clear AP AQ AR AS AT ATT DP DQ DR DS DT DTT NP NQ NR NS NT NTT
clear FP FQ FR FS FT FTT FF CD vc vt VTT  i  j  t  k VZA VZT

%Plot the Electrocardiogram
figure;

%Vector with Lead Names
name_lead = {'DI','aVR','V1','V4','DII','aVL','V2','V5','DIII','aVF','V3','V6'};

%FOR Loop to Traverse All Cells to Plot the 12 ECGs
for i = 1:12
    subplot(3,4,i); %Arrangement of 12 Figures - 3 rows and 4 columns

    %Choose the Time Range for Plotting
    ti = 5;   %Initial Time
    tf = 8.5;  %Final Time
    
    %Determine the Index for the Specified Time Interval
    idx = find(VZTs{i} >= ti & VZTs{i} <= tf);
    
    %Plot Amplitude and Time Values
    plot(VZTs{i}(idx), VZAs{i}(idx));
    
    %Graph Parameters
    title(sprintf('ECG Lead %s', name_lead{i}));
    xlabel('Time [s]');
    ylabel('Amplitude [mV]');
    grid on;
end %End of FOR Loop for Plotting

%Clear Variables that Will No Longer Be Used - Part 2
clear i ti tf idx

%Export Data
save('ECG_data.mat', 'VZTs', 'VZAs', 'name_lead'); %Export Data in .mat Format