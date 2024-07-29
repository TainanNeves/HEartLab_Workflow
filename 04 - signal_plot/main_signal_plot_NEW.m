%% SIGNAL PLOTING CODE

clear; clc;

%% Loading data

% Loading variables
load('C:\Users\HEartLab\Downloads\Pasta de Trabalho\Subpasta 2 - Desenvolvimento\Dados Aprendizagem\1 - Analisados\data_filtered_sync_E14_F3_R4.mat'); % Synchronized data
load('C:\Users\HEartLab\Downloads\InterpolatedSignalsE18_F02_R02_filtered'); % Interpolate data


%% Optical signals plot

% Define a Camera to use
Data = D_SYNC.CAM1;
Background = squeeze(Data(:,:,2000)); % Select a pixel in the image and shows the optical signal
[x, y] = pick_up_a_trace(Background, Data,1);    %Press space to stop
                                        
% Define a pixel position
p = [x(length(x)), y(length(y))];
%Frame Sampling
Fsampling = 4000;

% Full optic time plot
% Create a time vector
To = linspace(0, length(Data(1,1,:))/Fsampling, length(Data(1,1,:)));
% Plot the optical signal at the specified pixel position
f1 = figure('color', 'white', 'Position', [40 40 600 200]);
plot(To, squeeze(Data(p(1), p(2), :)), 'LineWidth', 1);
ylabel('%F');
set(gca, 'fontsize', 14);
xlim([0 8]);
title('Optical Signal Time Plot');

% Specific optic time plot
% Define the time range for the plot
t_in = 4;
t_out = 6;
start_sample = t_in*Fsampling; % Adjust the start sample according to your data
end_sample = t_out*Fsampling;   % Adjust the end sample according to your data
% Create a time vector
To = linspace(0, length(Data(1,1,:))/Fsampling, length(Data(1,1,:)));
% Plot the optical signal within the specified time range
f1 = figure('color', 'white', 'Position', [40 40 600 200]);
plot(To(start_sample:end_sample), squeeze(Data(p(1), p(2), start_sample:end_sample)), 'LineWidth', 1);
ylabel('%F');
set(gca, 'fontsize', 14);
xlabel('Time (s)');
xlim([t_in, t_out]);
title('Optical Signal Time Plot');


%% Electric signal plot
%Run the code section by section using F9

Data_ = InterpSignal.Data.MEA1;
Background = squeeze(Data(:,:,2000)); % Select a pixel in the image and shows the optical signal
[x, y] = pick_up_a_trace(Background, Data,1); 

% Define electrode to use
el = 129;
Data = D_SYNC.EL(el,:);

%Frame Sampling
Fsampling = 4000;

cases = {'MEA1', 'MEA2', 'MEA3', 'TANK'};
[x, y, source] = getElectrodePosition(el);
case_name = cases{source};

% Full electric time plot
% Create a time vector
To = linspace(0, length(Data)/Fsampling, length(Data));
% Plot the oelectrical signal for an specific electrode
f1 = figure('color', 'white', 'Position', [40 40 600 200]);
plot(To, squeeze(Data),'DisplayName', ['Electrode ' num2str(el)]);
hold on
ylabel('Potential ($\mu$V)', 'Interpreter', 'latex');
set(gca, 'fontsize', 14);
xlim([0 8]);
title('Electric Signal Time Plot');
legend('show');

% Specific electric time plot
% Define the time range for the plot
t_in = 4;
t_out = 6;
start_sample = t_in*Fsampling; % Adjust the start sample according to your data
end_sample = t_out*Fsampling;   % Adjust the end sample according to your data
% Create a time vector
To = linspace(0, length(Data)/Fsampling, length(Data));
% Plot the electricl signal for an specific electrode
f1 = figure('color', 'white', 'Position', [40 40 600 200]);
plot(To(start_sample:end_sample), Data(start_sample:end_sample), 'DisplayName', ['Electrode ' num2str(el)]);
hold on
ylabel('Potential ($\mu$V)', 'Interpreter', 'latex');
set(gca, 'fontsize', 14);
xlabel('Time (s)');
xlim([t_in t_out]);
title('Electric Signal Time Plot');
legend('show');

%% Mixed Plot

% Define a Camera to use
Data_O = D_SYNC.CAM1;
Background = squeeze(Data_O(:,:,2000)); % Select a pixel in the image and shows the optical signal
[x, y] = pick_up_a_trace(Background, Data_O, 1);     %Press space to stop
                                                                             
% Optical Pixel
pa = [x(1) y(1)];
pv = [x(2) y(2)];

% Electric electrodes
el1 = 14; % RA MEA1
el2 = 74; % LA MEA3
el3 = 27; % V MEA2
el4 = 134; % Tank

% Frequency Sampling
Fsampling = 4000;

% Time to plot
start_time = 4;
end_time = 6;

% Loading electric data
Data_E = D_SYNC.EL;

% Ploting
To = linspace(0, length(Data_E(1, :)) / Fsampling, length(Data_E(1, :)));
f1 = figure('color', 'white', 'Position', [40 40 600 700]);
% Subplot 1 - Optic Atrium
subplot(6, 1, 1)
plot(To, squeeze(Data_O(pa(1), pa(2), :)), 'LineWidth', 1);
ylabel('%F');
title('Optical signal (Atria)');
set(gca, 'fontsize', 14);
xlim([start_time end_time])
% Subplot 2 - Optic Ventricle
subplot(6, 1, 2)
plot(To, squeeze(Data_O(pv(1), pv(2), :)), 'LineWidth', 1);
ylabel('%F');
title('Optical signal (Ventricle)');
set(gca, 'fontsize', 14);
xlim([start_time end_time])
% Subplot 3 - MEA 1 (RA)
i = el1;
subplot(6, 1, 3)
plot(To, Data_E(i, :), 'LineWidth', 1);
ylabel('$\mu$V', 'Interpreter', 'latex');
title(['MEA1 el', num2str(i), ' (RA)']);
set(gca, 'fontsize', 14);
xlim([start_time end_time])
% Subplot 4 - MEA 3 (LA)
i = el2;
subplot(6, 1, 4)
plot(To, Data_E(i, :), 'LineWidth', 1);
ylabel('$\mu$V', 'Interpreter', 'latex');
title(['MEA3 el', num2str(i), ' (LA)']);
set(gca, 'fontsize', 14);
xlim([start_time end_time])
% Subplot 5 - MEA 2 (V)
i = el3;
subplot(6, 1, 5)
plot(To, Data_E(i, :), 'LineWidth', 1);
ylabel('$\mu$V', 'Interpreter', 'latex');
title(['MEA2 el', num2str(i), ' (V)']);
set(gca, 'fontsize', 14);
xlim([start_time end_time])
% Subplot 6 - TANK
i = el4;
subplot(6, 1, 6)
plot(To, Data_E(i, :), 'LineWidth', 1);
ylabel('$\mu$V', 'Interpreter', 'latex');
title(['TANK el', num2str(i)]);
xlabel('Time (s)');
set(gca, 'fontsize', 14);
xlim([start_time end_time])
% Linking Axes (move x axes in all data together)
linkaxes([subplot(6, 1, 1), subplot(6, 1, 2), subplot(6, 1, 3), subplot(6, 1, 4), subplot(6, 1, 5), subplot(6, 1, 6)], 'x')


%% MULTIMIXED PLOT

% Selecting eletrodes and pixels
% Define Cameras to use
Data_OV = D_SYNC.CAM1;
Background = squeeze(Data_OV(:,:,2000));
[xV, yV] = pick_up_a_trace(Background, Data_OV,1);

Data_ORA = D_SYNC.CAM2;
Background = squeeze(Data_ORA(:,:,2000));
[xRA, yRA] = pick_up_a_trace(Background, Data_ORA,1);

Data_OLA = D_SYNC.CAM3;
Background = squeeze(Data_OLA(:,:,2000));
[xLA, yLA] = pick_up_a_trace(Background, Data_OLA,1);

% Optical Pixel
pV = [xV yV];
pRA = [xRA yRA];
pLA = [xLA yLA];

% Electric electrodes
el1 = 6; % RA MEA1
el2 = 70; % LA MEA3
el3 = 26; % V MEA2
el4 = 132; % Tank

%Ploting
% Frequency Sampling
Fsampling = 4000;

% Time to plot
start_time = 2.4;
end_time = 6.4;

% Loading electric data
Data_E = D_SYNC.EL;

% Ploting
To = linspace(0, length(Data_E(1, :)) / Fsampling, length(Data_E(1, :)));
f1 = figure('color', 'white', 'Position', [40 40 600 700]);
% Subplot 1 - Optic Right Atrium
subplot(7, 1, 1)
plot(To, squeeze(Data_ORA(pRA(1), pRA(2), :)), 'LineWidth', 1);
ylabel('%F');
title(['Optical Signal (RA) | Points: ', num2str(pRA(1)), 'x', num2str(pRA(2))]);
set(gca, 'fontsize', 14);
xlim([start_time end_time])
% Subplot 2 - Optic Left Atrium
subplot(7, 1, 2)
plot(To, squeeze(Data_OLA(pLA(1), pLA(2), :)), 'LineWidth', 1);
ylabel('%F');
title(['Optical Signal (LA) | Points: ', num2str(pLA(1)), 'x', num2str(pLA(2))]);
set(gca, 'fontsize', 14);
xlim([start_time end_time])
% Subplot 3 - Optic Ventricle
subplot(7, 1, 3)
plot(To, squeeze(Data_OV(pV(1), pV(2), :)), 'LineWidth', 1);
ylabel('%F');
title(['Optical Signal (V) | Points: ', num2str(pV(1)), 'x', num2str(pV(2))]);
set(gca, 'fontsize', 14);
xlim([start_time end_time])
% Subplot 4 - MEA 1 (RA)
i = el1;
subplot(7, 1, 4)
plot(To, Data_E(i, :), 'LineWidth', 1);
ylabel('$\mu$V', 'Interpreter', 'latex');
title(['MEA1 el', num2str(i), ' (RA)']);
set(gca, 'fontsize', 14);
xlim([start_time end_time])
% Subplot 5 - MEA 3 (LA)
i = el2;
subplot(7, 1, 5)
plot(To, Data_E(i, :), 'LineWidth', 1);
ylabel('$\mu$V', 'Interpreter', 'latex');
title(['MEA3 el', num2str(i), ' (LA)']);
set(gca, 'fontsize', 14);
xlim([start_time end_time])
% Subplot 6 - MEA 2 (V)
i = el3;
subplot(7, 1, 6)
plot(To, Data_E(i, :), 'LineWidth', 1);
ylabel('$\mu$V', 'Interpreter', 'latex');
title(['MEA2 el', num2str(i), ' (V)']);
set(gca, 'fontsize', 14);
xlim([start_time end_time])
% Subplot 7 - TANK
i = el4;
subplot(7, 1, 7)
plot(To, Data_E(i, :), 'LineWidth', 1);
ylabel('$\mu$V', 'Interpreter', 'latex');
title(['TANK el', num2str(i)]);
xlabel('Time (s)');
set(gca, 'fontsize', 14);
xlim([start_time end_time])
% Linking Axes (move x axes in all data together)
linkaxes([subplot(7, 1, 1), subplot(7, 1, 2), subplot(7, 1, 3), subplot(7, 1, 4), subplot(7, 1, 5), subplot(7, 1, 6), subplot(7, 1, 7)], 'x');


%% Multielectrode optical signal plot

%Frequency Sampling
Fsampling = 4000;
time = [4, 5]; % s
Data_E = D_SYNC.EL(:, time(1)*Fsampling:time(2)*Fsampling);


% MEA 1
Data_O = D_SYNC.CAM1(:, :, time(1)*Fsampling:time(2)*Fsampling);
Background = squeeze(Data_O(:,:,2000));
% Background = single(D_SYNC.IMG.CAM1); Data_O = rot90(Data_O, 1); % To use the figure with no filtering
[x, y] = pick_up_a_trace(Background, Data_O,1);  
p13 = [x(13), y(13)]; p14 = [x(14), y(14)]; p15 = [x(15), y(15)]; p16 = [x(16), y(16)];
p9 = [x(9), y(9)]; p10 = [x(10), y(10)]; p11 = [x(11), y(11)]; p12 = [x(12), y(12)];
p5 = [x(5), y(5)]; p6 = [x(6), y(6)]; p7 = [x(7), y(7)]; p8 = [x(8), y(8)];
p1 = [x(1), y(1)]; p2 = [x(2), y(2)]; p3 = [x(3), y(3)]; p4 = [x(4), y(4)];
plotar_pontos_1(Data_O, Data_E, Fsampling, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16);

% MEA 2
Data_O = D_SYNC.CAM2(:, :, time(1)*Fsampling:time(2)*Fsampling);
Background = squeeze(Data_O(:,:,2000));
% Background = single(D_SYNC.IMG.CAM2); Data_O = rot90(Data_O, 1); % To use the figure with no filtering
[x, y] = pick_up_a_trace(Background, Data_O,1); 
p17 = [x(1), y(1)]; p21 = [x(5), y(5)]; p25 = [x(9), y(9)]; p29 = [x(13), y(13)];
p18 = [x(2), y(2)]; p22 = [x(6), y(6)]; p26 = [x(10), y(10)]; p30 = [x(14), y(14)];
p19 = [x(3), y(3)]; p23 = [x(7), y(7)]; p27 = [x(11), y(11)]; p31 = [x(15), y(15)];
p20 = [x(4), y(4)]; p24 = [x(8), y(8)]; p28 = [x(12), y(12)]; p32 = [x(16), y(16)];
plotar_pontos_2(Data_O, Data_E, Fsampling, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32);

% MEA 3
Data_O = D_SYNC.CAM3(:, :, time(1)*Fsampling:time(2)*Fsampling);
Background = squeeze(Data_O(:,:,2000));
% Background = single(D_SYNC.IMG.CAM3); Data_O = rot90(Data_O, 1); % To use the figure with no filtering
[x, y] = pick_up_a_trace(Background, Data_O,1); 
p77 = [x(13), y(13)]; p78 = [x(14), y(14)]; p79 = [x(15), y(15)]; p80 = [x(16), y(16)];
p73 = [x(9), y(9)]; p74 = [x(10), y(10)]; p75 = [x(11), y(11)]; p76 = [x(12), y(12)];
p69 = [x(5), y(5)]; p70 = [x(6), y(6)]; p71 = [x(7), y(7)]; p72 = [x(8), y(8)];
p65 = [x(1), y(1)]; p66 = [x(2), y(2)]; p67 = [x(3), y(3)]; p68 = [x(4), y(4)];
plotar_pontos_3(Data_O, Data_E, Fsampling, p65, p66, p67, p68, p69, p70, p71, p72, p73, p74, p75, p76, p77, p78, p79, p80);































