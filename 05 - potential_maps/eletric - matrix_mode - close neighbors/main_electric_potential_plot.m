%% MAIN ELECTRIC POTENTIAL PLOTS
% Use the electric data filtered
clear; 
clc;

%% Loading Data If using matlab extraction
load("D:\Qualification\Analysis\E14F03R04\data\electric_data_E14F03R04_filtered.mat");
data = D_EL.Data;

%%Transformar para mV (Codigo de export tras em uV)
data = data/1000;


%% 
% Plot 3
el1 = 28; % Ventricle MEA
el2 = 75; % Atria MEA
el3 = 142; % Tank electrode
tin = 2; % Initial time (s)
tfin = 4; % Final time (s)
time = linspace(tin, tfin, length(tin * D_EL.Header.sample_rate : tfin * D_EL.Header.sample_rate));
tPlot = 3.1; % Time at which the red line will appear

figure();
subplot(3, 1, 1); 
plot(time, data(el2, tin * D_EL.Header.sample_rate : tfin * D_EL.Header.sample_rate));
title(['Electrode ', num2str(el2), ' (Left Atria MEA)']);
ylabel('mV');
xline(tPlot, 'r');
text(tPlot, max(data(el2, tin * D_EL.Header.sample_rate : tfin * D_EL.Header.sample_rate)), [num2str(tPlot), ' s'], 'Color', 'r');

subplot(3, 1, 2); 
plot(time, data(el1, tin * D_EL.Header.sample_rate : tfin * D_EL.Header.sample_rate));
title(['Electrode ', num2str(el1), ' (Ventricle MEA)']);
ylabel('mV');
xline(tPlot, 'r');

subplot(3, 1, 3); 
plot(time, data(el3, tin * D_EL.Header.sample_rate : tfin * D_EL.Header.sample_rate));
title(['Electrode ', num2str(el3), ' (Tank Electrode)']);
ylabel('mV');
xlabel('Time (s)');
xline(tPlot, 'r');


%% MEAs Isopotential Map

%Loading electrodes
MEA1 = zeros(4, 4);
MEA1 = [0, data(14,fix(tPlot*4000)), data(15,fix(tPlot*4000)), data(16,fix(tPlot*4000));     %Bad electrode
    data(9,fix(tPlot*4000)), data(10,fix(tPlot*4000)), data(11,fix(tPlot*4000)), 0;          %Bad electrode
    data(5,fix(tPlot*4000)), data(6,fix(tPlot*4000)), data(7,fix(tPlot*4000)), data(8,fix(tPlot*4000));
    data(1,fix(tPlot*4000)), data(2,fix(tPlot*4000)), data(3,fix(tPlot*4000)),data(4,fix(tPlot*4000))];

MEA2 = zeros(4, 4);
MEA2 = [data(17,fix(tPlot*4000)), data(21,fix(tPlot*4000)), data(25,fix(tPlot*4000)), data(29,fix(tPlot*4000));
    data(18,fix(tPlot*4000)), data(22,fix(tPlot*4000)), data(26,fix(tPlot*4000)), data(30,fix(tPlot*4000));
    data(19,fix(tPlot*4000)), data(23,fix(tPlot*4000)), data(27,fix(tPlot*4000)), data(31,fix(tPlot*4000));
    data(20,fix(tPlot*4000)), data(24,fix(tPlot*4000)),data(28,fix(tPlot*4000)), data(32,fix(tPlot*4000))];

MEA3 = zeros(4, 4);
MEA3 = [data(77,fix(tPlot*4000)), data(78,fix(tPlot*4000)), data(79,fix(tPlot*4000)), 0;     %Bad electrode
    data(73,fix(tPlot*4000)), data(74,fix(tPlot*4000)), data(75,fix(tPlot*4000)), data(76,fix(tPlot*4000));
    data(69,fix(tPlot*4000)), data(70,fix(tPlot*4000)), data(71,fix(tPlot*4000)), data(72,fix(tPlot*4000));
    data(65,fix(tPlot*4000)), data(66,fix(tPlot*4000)), data(67,fix(tPlot*4000)), data(68,fix(tPlot*4000))];

%Interpolatin zeros in the matrix
interpolatedMatrix1 = interpolateZeros(MEA1);
disp(interpolatedMatrix1);

interpolatedMatrix2 = interpolateZeros(MEA2);
disp(interpolatedMatrix2);

interpolatedMatrix3 = interpolateZeros(MEA3);
disp(interpolatedMatrix3);

%Creating figure
figure;

% First Subplot
subplot(1, 3, 1);
pcolor(interpolatedMatrix1);     % Plot the pcolor plot
shading interp;                 % Smooth shading between cells
colorbar;                       % Display color bar
c = colorbar;
c.Label.String = 'Potential (mV)'; % Set the colorbar label
axis tight;                     % Set the axis limits to fit the data
title(['MEA 1 (RA)']);

% Second Subplot
subplot(1, 3, 2);
pcolor(interpolatedMatrix2);
shading interp;
colorbar;
c = colorbar;
c.Label.String = 'Potential (mV)';
axis tight;
title(['MEA 2 (V)']);

% Third Subplot
subplot(1, 3, 3);
pcolor(interpolatedMatrix3);
shading interp;
colorbar;
c = colorbar;
c.Label.String = 'Potential (mV)';
axis tight;
title(['MEA 3 (LA)']);

% % Define the common colormap for all subplots
% c = jet(256);
% colormap(c);

% Set the color axis limits (if needed)
% colorbarLimits = [-0.1, 0.1]; % Replace with your desired range
% caxis(colorbarLimits);

% Remove the axes box for all subplots
for i = 1:3
    subplot(1, 3, i);
    box off;
end


%% Tank Isopotential Map

% carregando eletrodos do tanque em momento tPlot desejado:
%Preenchendo sinais dos eletrodos
%Linha 1
el145 = data(145,fix(tPlot*4000));
el146 = data(146,fix(tPlot*4000));
el155 = data(155,fix(tPlot*4000));
el156 = data(156,fix(tPlot*4000));
el165 = data(165,fix(tPlot*4000));
el166 = data(166,fix(tPlot*4000));
el129 = data(129,fix(tPlot*4000));
el130 = data(130,fix(tPlot*4000));
el139 = data(139,fix(tPlot*4000));
el140 = data(140,fix(tPlot*4000));
el181 = data(181,fix(tPlot*4000));
el182 = data(182,fix(tPlot*4000));
%Linha 2
el147 = data(147,fix(tPlot*4000));
el157 = data(157,fix(tPlot*4000));
el167 = data(167,fix(tPlot*4000));
el131 = data(131,fix(tPlot*4000));
el141 = data(141,fix(tPlot*4000));
el183 = data(183,fix(tPlot*4000));
%Linha 3
el148 = data(148,fix(tPlot*4000));
el149 = data(149,fix(tPlot*4000));
el158 = data(158,fix(tPlot*4000));
el159 = data(159,fix(tPlot*4000));
el168 = data(168,fix(tPlot*4000));
el169 = data(169,fix(tPlot*4000));
el132 = data(132,fix(tPlot*4000));
el133 = data(133,fix(tPlot*4000));
el142 = data(142,fix(tPlot*4000));
el143 = data(143,fix(tPlot*4000));
el184 = data(184,fix(tPlot*4000));
el185 = data(185,fix(tPlot*4000));
%Linha 4
el150 = data(150,fix(tPlot*4000));
el151 = data(151,fix(tPlot*4000));
el160 = data(160,fix(tPlot*4000));
el161 = data(161,fix(tPlot*4000));
el170 = data(170,fix(tPlot*4000));
el171 = data(171,fix(tPlot*4000));
el134 = data(134,fix(tPlot*4000));
el135 = data(135,fix(tPlot*4000));
el144 = data(144,fix(tPlot*4000));
el177 = data(177,fix(tPlot*4000));
el186 = data(186,fix(tPlot*4000));
el187 = data(187,fix(tPlot*4000));
%Linha 5
el152 = data(152,fix(tPlot*4000));
el162 = data(162,fix(tPlot*4000));
el172 = data(172,fix(tPlot*4000));
el136 = data(136,fix(tPlot*4000));
el178 = data(178,fix(tPlot*4000));
el188 = data(188,fix(tPlot*4000));
%Linha 6
el153 = data(153,fix(tPlot*4000));
el154 = data(154,fix(tPlot*4000));
el163 = data(163,fix(tPlot*4000));
el164 = data(164,fix(tPlot*4000));
el173 = data(173,fix(tPlot*4000));
el174 = data(174,fix(tPlot*4000));
el137 = data(137,fix(tPlot*4000));
el138 = data(138,fix(tPlot*4000));
el179 = data(179,fix(tPlot*4000));
el180 = data(180,fix(tPlot*4000));
el189 = data(189,fix(tPlot*4000));
el190 = data(190,fix(tPlot*4000));

%Criando Matrix do tanque com zeros
M = zeros(6, 18);

% Define the input values for each line of tank
line1 = [el145, el146, el155, el156, el165, el166, el129, el130, el139, el140, el181, el182];
line2 = [el147, el157, el167, el131, el141, el183];
line3 = [el148, el149, el158, el159, el168, el169, el132, el133, el142, el143, el184, el185];
line4 = [el150, el151, el160, el161, el170, el171, el134, el135, el144, el177, el186, el187];
line5 = [el152, el162, el172, el136, el178, el188];
line6 = [el153, el154, el163, el164, el173, el174, el137, el138, el179, el180, el189, el190];

% Call the completeMatrix function with the values for each line
completedMatrix = completeMatrix(line1, line2, line3, line4, line5, line6);
disp(completedMatrix);

%Interpolatin zeros in the matrix
interpolatedMatrix = interpolateZeros(completedMatrix);
disp(interpolatedMatrix);


%Gerando figura do mapa isopotencial
figure;
pcolor(interpolatedMatrix);     % Plot the pcolor plot
shading interp;                 % Smooth shading between cells
colorbar;                       % Display color bar
c = colorbar;
c.Label.String = 'mV';          % Set the colorbar label to 'mV'
axis tight;                     % Set the axis limits to fit the data
title(['Tank Isopotential Map (Moment: ', num2str(tPlot),' s)']);

% Set the color axis limits
% colorbarLimits = [-0.1, 0.1];   % Replace with your desired range
%                                 %Max = max(interpolatedMatrix(:))
%                                 %min = min(interpolatedMatrix(:))
% caxis(colorbarLimits);

% Remove the axes box
box off;

% Define the x-values where you want to add vertical lines
x_values = [3.5, 6.5, 9.5, 12.5, 15.5];

% Add vertical black lines
for x = x_values
    line([x, x], ylim, 'Color', 'k', 'LineWidth', 1);
end


%% Plot submetido para o SIIM
el1 = 22; % Ventricle MEA
el2 = 75; % Atria MEA
el3 = 142; % Tank electrode
tin = 27; % Initial time (s)
tfin = 29; % Final time (s)
time = linspace(tin, tfin, length(tin * D_EL.Header.sample_rate : tfin * D_EL.Header.sample_rate));
tPlot = 7:0.25:8; % Time at which the red line will appear

figure();

% Subplot 1
subplot(4, 1, 1); 
plot(time, data(el2, tin * D_EL.Header.sample_rate : tfin * D_EL.Header.sample_rate));
title(['Electrode ', num2str(el2), ' (Left Atria MEA)']);
ylabel('mV');
xline(tPlot, 'r');
text(tPlot, max(data(el2, tin * D_EL.Header.sample_rate : tfin * D_EL.Header.sample_rate)), [num2str(tPlot), ' s'], 'Color', 'r');

% Subplot 2
subplot(4, 1, 2); 
plot(time, data(el3, tin * D_EL.Header.sample_rate : tfin * D_EL.Header.sample_rate));
title(['Electrode ', num2str(el3), ' (Tank Electrode)']);
ylabel('mV');
xlabel('Time (s)');
xline(tPlot, 'r');

% Subplot 3 (occupies two rows)
subplot(4, 1, [3, 4]); 
pcolor(interpolatedMatrix);     % Plot the pcolor plot
shading interp;                 % Smooth shading between cells
colorbar;                       % Display color bar
c = colorbar;
c.Label.String = 'mV';          % Set the colorbar label to 'mV'
axis tight;                     % Set the axis limits to fit the data
title(['Tank Isopotential Map (Moment: ', num2str(tPlot),' s)']);

% Set the color axis limits
% colorbarLimits = [-0.3, 0.05];   % Replace with your desired range
%                                 %Max = max(interpolatedMatrix(:))
%                                 %min = min(interpolatedMatrix(:))
% caxis(colorbarLimits);

% Remove the axes box
box off;

% Define the x-values where you want to add vertical lines
x_values = [3.5, 6.5, 9.5, 12.5, 15.5];

% Add vertical black lines
for x = x_values
    line([x, x], ylim, 'Color', 'k', 'LineWidth', 1);
end

% Plot a green circle at (3, 13)
hold on; % This allows you to overlay additional elements on the plot
scatter(13, 4, 100, 'magenta', 'filled'); % Adjust the size (100) and color ('g' for green) as needed
hold off; % Release the hold on the plot


%% Generating video
% Define the time range
tplotin = 2;
tplotfin = 2.3;
tin = 8.2;               % Starting time
tfin = 10.2;             % Ending time
frameRate = 60;          % Frame rate for the video
durationMultiplier = 10; % Multiply the duration
numFrames = frameRate * (tfin - tin) * durationMultiplier; % Number of frames

% Create a VideoWriter object
videoFile = 'outputVideo.mp4';  % Specify the filename and format (e.g., '.mp4')
videoObj = VideoWriter(videoFile, 'MPEG-4');
videoObj.Quality = 100;  % Set video quality (1-100)
videoObj.VideoCompressionMethod = 'H.264';
videoObj.FrameRate = frameRate;
open(videoObj);

% Set the figure size to full HD (1920x1080)
figure('Position', [0, 0, 1920, 1080]);

for t = linspace(tin, tfin, numFrames)
    %Preenchendo sinais dos eletrodos
    %Linha 1
    el145 = data(145, fix(t * 4000));
    el146 = data(146, fix(t * 4000));
    el155 = data(155, fix(t * 4000));
    el156 = data(156, fix(t * 4000));
    el165 = data(165, fix(t * 4000));
    el166 = data(166, fix(t * 4000));
    el129 = data(129, fix(t * 4000));
    el130 = data(130, fix(t * 4000));
    el139 = data(139, fix(t * 4000));
    el140 = data(140, fix(t * 4000));
    el181 = data(181, fix(t * 4000));
    el182 = data(182, fix(t * 4000));
    %Linha 2
    el147 = data(147, fix(t * 4000));
    el157 = data(157, fix(t * 4000));
    el167 = data(167, fix(t * 4000));
    el131 = data(131, fix(t * 4000));
    el141 = data(141, fix(t * 4000));
    el183 = data(183, fix(t * 4000));
    %Linha 3
    el148 = data(148, fix(t * 4000));
    el149 = data(149, fix(t * 4000));
    el158 = data(158, fix(t * 4000));
    el159 = data(159, fix(t * 4000));
    el168 = data(168, fix(t * 4000));
    el169 = data(169, fix(t * 4000));
    el132 = data(132, fix(t * 4000));
    el133 = data(133, fix(t * 4000));
    el142 = data(142, fix(t * 4000));
    el143 = data(143, fix(t * 4000));
    el184 = data(184, fix(t * 4000));
    el185 = data(185, fix(t * 4000));
    %Linha 4
    el150 = data(150, fix(t * 4000));
    el151 = data(151, fix(t * 4000));
    el160 = data(160, fix(t * 4000));
    el161 = data(161, fix(t * 4000));
    el170 = data(170, fix(t * 4000));
    el171 = data(171, fix(t * 4000));
    el134 = data(134, fix(t * 4000));
    el135 = data(135, fix(t * 4000));
    el144 = data(144, fix(t * 4000));
    el177 = data(177, fix(t * 4000));
    el186 = data(186, fix(t * 4000));
    el187 = data(187, fix(t * 4000));
    %Linha 5
    el152 = data(152, fix(t * 4000));
    el162 = data(162, fix(t * 4000));
    el172 = data(172, fix(t * 4000));
    el136 = data(136, fix(t * 4000));
    el178 = data(178, fix(t * 4000));
    el188 = data(188, fix(t * 4000));
    %Linha 6
    el153 = data(153, fix(t * 4000));
    el154 = data(154, fix(t * 4000));
    el163 = data(163, fix(t * 4000));
    el164 = data(164, fix(t * 4000));
    el173 = data(173, fix(t * 4000));
    el174 = data(174, fix(t * 4000));
    el137 = data(137, fix(t * 4000));
    el138 = data(138, fix(t * 4000));
    el179 = data(179, fix(t * 4000));
    el180 = data(180, fix(t * 4000));
    el189 = data(189, fix(t * 4000));
    el190 = data(190, fix(t * 4000));

    %Criando Matrix do tanque com zeros
    M = zeros(6, 18);

    % Define the input values for each line of the tank
    line1 = [el145, el146, el155, el156, el165, el166, el129, el130, el139, el140, el181, el182];
    line2 = [el147, el157, el167, el131, el141, el183];
    line3 = [el148, el149, el158, el159, el168, el169, el132, el133, el142, el143, el184, el185];
    line4 = [el150, el151, el160, el161, el170, el171, el134, el135, el144, el177, el186, el187];
    line5 = [el152, el162, el172, el136, el178, el188];
    line6 = [el153, el154, el163, el164, el173, el174, el137, el138, el179, el180, el189, el190];

    % Call the completeMatrix function with the values for each line
    completedMatrix = completeMatrix(line1, line2, line3, line4, line5, line6);
    
    %Interpolating zeros in the matrix
    interpolatedMatrix = interpolateZeros(completedMatrix);
    

    time = linspace(tplotin, tplotfin, length(tplotin * D_EL.Header.sample_rate : tplotfin * D_EL.Header.sample_rate));
    % Subplot 1
    subplot(4, 1, 1); 
    plot(time, data(el2, tplotin * D_EL.Header.sample_rate : tplotfin * D_EL.Header.sample_rate));
    title(['Electrode ', num2str(el2), ' (Left Atria MEA)']);
    ylabel('mV');
    xline(t, 'r');
    
    % Subplot 2
    subplot(4, 1, 2); 
    plot(time, data(el3, tplotin * D_EL.Header.sample_rate : tplotfin * D_EL.Header.sample_rate));
    title(['Electrode ', num2str(el3), ' (Tank Electrode)']);
    ylabel('mV');
    xlabel('Time (s)');
    xline(t, 'r');
    
    %Gerando figura do mapa isopotencial
    % Subplot 3 (occupies two rows)
    subplot(4, 1, [3, 4]); 
    pcolor(interpolatedMatrix);     % Plot the pcolor plot
    shading interp;                 % Smooth shading between cells
    colorbar;                       % Display color bar
    c = colorbar;
    c.Label.String = 'mV';          % Set the colorbar label to 'mV'
    axis tight;                     % Set the axis limits to fit the data
    %title(['Tank Isopotential Map (Moment: ', num2str(t),' s)']);

    % Set the color axis limits
    colorbarLimits = [min(min(data(129:190, round(tin*D_EL.Header.sample_rate):round(tfin*D_EL.Header.sample_rate)))), max(max(data(129:190, (tin*D_EL.Header.sample_rate):(tfin*D_EL.Header.sample_rate))))];  % Replace with your desired range
    caxis(colorbarLimits);

    % Remove the axes box
    box off;

    % Define the x-values where you want to add vertical lines
    x_values = [3.5, 6.5, 9.5, 12.5, 15.5];

    % Add vertical black lines
    for x = x_values
        line([x, x], ylim, 'Color', 'k', 'LineWidth', 1);
    end

    % Capture the current figure as a frame
    frame = getframe(gcf);  % 'gcf' gets the current figure

    % Write the frame to the video file
    writeVideo(videoObj, frame);

    % Close the current figure to avoid accumulating plots
    close(gcf);
end

% Close the video file
close(videoObj);


