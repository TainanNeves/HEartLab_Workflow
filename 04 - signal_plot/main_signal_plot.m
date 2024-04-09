%% SIGNAL PLOTING CODE

clear; clc;


%% Loading data

% Loading variables
load('C:\Users\HEartLab\Documents\GitHub\HEartLab\00 - examples\data_filtered_sync_E14_F3_R4.mat'); % Synchronized data


%% Optical signals plot

% Define a Camera to use
Data = D_SYNC.CAM1;
Background = squeeze(Data(:,:,2000));
pick_up_a_trace(Background, Data,1);    % Select a pixel in the image and shows the optical signal
                                        %Press space to stop

% Define a pixel position
p = [68, 58];
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
start_sample = 4*Fsampling; % Adjust the start sample according to your data
end_sample = 6*Fsampling;   % Adjust the end sample according to your data
% Create a time vector
To = linspace(0, length(Data(1,1,:))/Fsampling, length(Data(1,1,:)));
% Plot the optical signal within the specified time range
f1 = figure('color', 'white', 'Position', [40 40 600 200]);
plot(To(start_sample:end_sample), squeeze(Data(p(1), p(2), start_sample:end_sample)), 'LineWidth', 1);
ylabel('%F');
set(gca, 'fontsize', 14);
xlabel('Time (s)');
title('Optical Signal Time Plot');


%% Electric signal plot

% Define electrode to use
el = 75;
Data = D_SYNC.EL(el,:);

%Frame Sampling
Fsampling = 4000;

% Full electric time plot
% Create a time vector
To = linspace(0, length(Data)/Fsampling, length(Data));
% Plot the oelectrical signal for an specific electrode
f1 = figure('color', 'white', 'Position', [40 40 600 200]);
plot(To, Data, 'LineWidth', 1);
ylabel('Potential ($\mu$V)', 'Interpreter', 'latex');
set(gca, 'fontsize', 14);
xlim([0 8]);
title('Electric Signal Time Plot');

% Specific electric time plot
% Define the time range for the plot
start_sample = 4*Fsampling; % Adjust the start sample according to your data
end_sample = 6*Fsampling;   % Adjust the end sample according to your data
% Create a time vector
To = linspace(0, length(Data)/Fsampling, length(Data));
% Plot the electricl signal for an specific electrode
f1 = figure('color', 'white', 'Position', [40 40 600 200]);
plot(To(start_sample:end_sample), Data(start_sample:end_sample), 'LineWidth', 1);
ylabel('Potential ($\mu$V)', 'Interpreter', 'latex');
set(gca, 'fontsize', 14);
xlabel('Time (s)');
title('Electric Signal Time Plot');


%% Mixed Plot

% Define a Camera to use
Data_O = D_SYNC.CAM1;
Background = squeeze(Data_O(:,:,2000));
pick_up_a_trace(Background, Data_O,1);    % Select a pixel in the image and shows the optical signal
                                        %Press space to stop                                      
% Optical Pixel
pa = [95 130];
pv = [75 49];

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
pick_up_a_trace(Background, Data_OV,1);

Data_ORA = D_SYNC.CAM2;
Background = squeeze(Data_ORA(:,:,2000));
pick_up_a_trace(Background, Data_ORA,1);

Data_OLA = D_SYNC.CAM3;
Background = squeeze(Data_OLA(:,:,2000));
pick_up_a_trace(Background, Data_OLA,1);

% Optical Pixel
pV = [37 51];
pRA = [28 141];
pLA = [43 162];

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

% Define a Camera to use
Data_O = D_SYNC.CAM1;
Background = squeeze(Data_O(:,:,2000));
pick_up_a_trace(Background, Data_O,1);    % Select a pixel in the image and shows the optical signal
                                        %Press space to stop

%Frequency Sampling
Fsampling = 4000;

% Positioning end ploting electrodes (In optical signal)
% MEA 1
p13 = [68,142]; p14 = [78,144]; p15 = [87,146]; p16 = [97,149];
p9 = [72,131]; p10 = [81,134]; p11 = [92,138]; p12 = [101,138];
p5 = [75,120]; p6 = [86,125]; p7 = [96,126]; p8 = [104,128];
p1 = [78,112]; p2 = [88,115]; p3 = [97,118]; p4 = [106,118];
plotar_pontos_1(Data_O, Data_O, Fsampling, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16);
% MEA 2
p17 = [68,142]; p21 = [78,144]; p25 = [87,146]; p29 = [97,149];
p18 = [72,131]; p22 = [81,134]; p26 = [92,138]; p30 = [101,138];
p19 = [75,120]; p23 = [86,125]; p27 = [96,126]; p31 = [104,128];
p20 = [78,112]; p24 = [88,115]; p28 = [97,118]; p32 = [106,118];
plotar_pontos_2(Data_O, Data_O, Fsampling, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32);
% MEA 3
p77 = [68,142]; p78 = [78,144]; p79 = [87,146]; p80 = [97,149];
p73 = [72,131]; p74 = [81,134]; p75 = [92,138]; p76 = [101,138];
p69 = [75,120]; p70 = [86,125]; p71 = [96,126]; p72 = [104,128];
p65 = [78,112]; p66 = [88,115]; p67 = [97,118]; p68 = [106,118];
plotar_pontos_3(Data_O, Data_O, Fsampling, p65, p66, p67, p68, p69, p70, p71, p72, p73, p74, p75, p76, p77, p78, p79, p80);






























