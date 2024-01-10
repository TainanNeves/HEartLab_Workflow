%% Extraction and filtering of electrical data

clear; clc;


%% Reads the Open Ephys files

% fulfilename is the path to the structure.oebin file contained in one of
% the experiments
fullfilename = 'E:\HEartLab\TAINAN WORKFLOW\00 - examples\electric_E14_F3_recording4\structure.oebin'; % Put the .oebin path

% Channels to save
channels = [1:192];

DATA = load_open_ephys_binary(fullfilename, 'continuous', 1); % Load the data file
DATA.Data = DATA.Data(channels, :); % Select desired channels
EVENTS = load_open_ephys_binary(fullfilename, 'events', 1); % Load the TTL file

rectime = length(DATA.Data(1, :))/DATA.Header.sample_rate; %s



%% Comment this section if conversion to uV is not wanted
bitVolts = 0.1949999928474426; % Conversion factor
DATA.Data = DATA.Data*bitVolts; % Converts the data to uV


%% Convert TTL to appropiate format
% Takes only every other timestamp to record the timing of each rising edge
TTL = EVENTS.Timestamps(1:2:end);

if isempty(TTL)==false
    optictime(1) = (TTL(1) - DATA.Timestamps(1)); %s
    optictime(2) = optictime(1) + 10; %s
end

%% Plot novo
Fs = 4000;
%Default time
if isempty(TTL)         %Sem TTL
    ti = 1;
    tf = rectime;
    tline1 = ti;
    tline2 = tf;
else
    ti = optictime(1);
    tf = optictime(2);
    tline1 = (TTL(1) - DATA.Timestamps(1)); %s
    tline2 = tline1 + 10; %s
end

%Specific Time
% ti = 10;
% tf = 15;
% tline1 = ti;
% tline2 = tf;

time = linspace(ti, tf, (tf - ti) * Fs + 1);


%Electrodes selection
el1 = 6;       %RA
el2 = 70;      %LA
el3 = 26;      %V
el4 = 132;      %Tank



% Create a figure
figure;
sgtitle('EXP XX - Folder XX - REC XX');

subplot(5, 1, 1);
filteredSignal1 = filter_signal(0.5, 250, DATA.Data(el1, :), Fs);
plot(time, filteredSignal1(:, ti * Fs:tf * Fs));
title(['Eletrodo: ' num2str(el1) ' (RA)']);
xlabel("Time (s)");
ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex');
xline(tline1, 'red');
xline(tline2, 'red');

subplot(5, 1, 2);
filteredSignal1 = filter_signal(0.5, 250, DATA.Data(el2, :), Fs);
plot(time, filteredSignal1(:, ti * Fs:tf * Fs));
title(['Eletrodo: ' num2str(el2) ' (LA)']);
xlabel("Time (s)");
ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex');
xline(tline1, 'red');
xline(tline2, 'red');

subplot(5, 1, 3);
filteredSignal2 = filter_signal(0.5, 250, DATA.Data(el3, :), Fs);
plot(time, filteredSignal2(:, ti * Fs:tf * Fs));
title(['Eletrodo: ' num2str(el3) ' (V)']);
xlabel("Time (s)");
ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex');
xline(tline1, 'red');
xline(tline2, 'red');

subplot(5, 1, 4);
filteredSignal3 = filter_signal(0.5, 250, DATA.Data(el4, :), Fs);
plot(time, filteredSignal3(:, ti * Fs:tf * Fs));
title(['Eletrodo: ' num2str(el4) ' (Tanque)']);
xlabel("Time (s)");
ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex');
xline(tline1, 'red');
xline(tline2, 'red');




% subplot(5, 1, 5);
% sub = DATA.Data(el1,:) - DATA.Data(el3, :);
% filteredSignal4 = filter_signal(0.5, 250, sub, Fs);
% plot(time, filteredSignal4(:, ti * Fs:tf * Fs));
% title(['Subtração: ' num2str(el1) ' - ' num2str(el3) ' (RA - V)']);
% xlabel("Time (s)");
% ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex');
% xline(tline1, 'red');
% xline(tline2, 'red');

subplot(5, 1, 5);
sub = DATA.Data(el2,:) - DATA.Data(el3, :);
filteredSignal4 = filter_signal(0.5, 250, sub, Fs);
plot(time, filteredSignal4(:, ti * Fs:tf * Fs));
title(['Subtração: ' num2str(el2) ' - ' num2str(el3) ' (LA - V)']);
xlabel("Time (s)");
ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex');
xline(tline1, 'red');
xline(tline2, 'red');



%% Save file to .mat

% Filename to save
FileName = 'EXX_FXX_RXX';

% Raw file export
D_EL.Data = DATA.Data; 
D_EL.Timestamps = DATA.Timestamps;
D_EL.Header = DATA.Header;
D_EL.Header.channels = channels;
D_EL.TTL = TTL;
D_EL.opticalin = TTL(1) - DATA.Timestamps(1);
save(['electric_data_', FileName, '_raw'], 'D_EL'); 

clear D_EL;

% Filtered file export
dataRaw = DATA.Data;
dataFiltered = zeros(size(dataRaw));
for i = 1:size(dataRaw, 1)
    Data = dataRaw(i, :);
    dataFiltered(i, :) = filter_signal(0.5, 250, Data, DATA.Header.sample_rate);
end

% Subtraindo o sinal - média (apenas no tanque)
tank_signals = dataFiltered([129:174,177:190], :);
mean_tank = mean(tank_signals);
dataFiltered([129:174,177:190], :) = tank_signals - mean_tank;

D_EL.Data = dataFiltered; 
D_EL.Timestamps = DATA.Timestamps;
D_EL.Header = DATA.Header;
D_EL.Header.channels = channels;
D_EL.TTL = TTL;
D_EL.opticalin = TTL(1) - DATA.Timestamps(1);

save(['electric_data_', FileName, '_filtered'], 'D_EL', '-v7.3');

% Limpar variáveis indesejadas
clear bitVolts channels DATA EVENTS D_EL FileName TTL fullfilename

