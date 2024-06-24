%% Main - Potential Map Plots
%this code allows you to plot a 2d isopotential maps from TANK and MEAs
%using laplacian interpolator
%Authors: Angélica Quadros, Ismael Romero, Tainan Neves.

clear; clc;


%% Load data

%Load the filtered data
load('E:\HEartLab\TAINAN WORKFLOW\00 - examples\data_filtered_sync_E14_F3_R4.mat'); % Synchronized data
Data = D_SYNC.EL;

% Transform to mV (Exported default data are in uV)
% data = data/1000;


%% Previsualizing data from electric time view

% Defining parameters
el1 = 10; % Right Atria MEA (add this line)
el2 = 74; % Left Atria MEA
el3 = 23; % Ventricle MEA
el4 = 142; % Tank electrode
tin = 4; % Initial time (s)
tfin = 7; % Final time (s)
Fsampling = 4000;
time = linspace(tin, tfin, length(tin*Fsampling:tfin*Fsampling));

figure();
subplot(4, 1, 1); % Change the number of subplots to 5
plot(time, Data(el1, tin*Fsampling:tfin*Fsampling));
title(['Electrode ', num2str(el1), ' (Right Atria MEA)']);
ylabel('$\mu$V', 'Interpreter', 'latex');
xlabel('Time (s)');


subplot(4, 1, 2); % Adjust subplot index
plot(time, Data(el2, tin*Fsampling:tfin*Fsampling));
title(['Electrode ', num2str(el2), ' (Left Atria MEA)']);
ylabel('$\mu$V', 'Interpreter', 'latex');

subplot(4, 1, 3); % Adjust subplot index
plot(time, Data(el3, tin*Fsampling:tfin*Fsampling));
title(['Electrode ', num2str(el3), ' (Ventricle MEA)']);
ylabel('$\mu$V', 'Interpreter', 'latex');

subplot(4, 1, 4); % Adjust subplot index
plot(time, Data(el4, tin*Fsampling:tfin*Fsampling));
title(['Electrode ', num2str(el4), ' (Tank Electrode)']);
ylabel('$\mu$V', 'Interpreter', 'latex');

clear el1 el2 el3 el4; 


%% POTENTIAL MAPS
lim = [-50 50];
sample = 4.5*4000;
plot_electric_pot(Data, lim, sample, 1); % MEA 1
plot_electric_pot(Data, lim, sample, 2); % MEA 2
plot_electric_pot(Data, lim, sample, 3); % MEA 3
plot_electric_pot(Data, lim, sample, 4); % TANK

%% SAVING STRUCT WITH INTERPOLATED DATA

% Initialize the struct to store results
pot_values = struct();

for i = 1:4
    V_interpolated = electric_interp(Data, i);  % Interpolate data
   
    % storing the adjusted matrices
    if i < 4
        electrodes = sprintf('MEA%d', i);
        pot_values.(electrodes) = fillMatrixMEA(phase_el);
    else
        electrodes = sprintf('TANK');
        pot_values.(electrodes) = fillMatrixTANK(phase_el);
    end
end

%% VIDEO - calculating time of video
isample = 4.0*4000;
fsample = 4.5*4000;
step = 4; % sample step between plots
videoFrameRate = 30; % fps

Tvideo = ((fsample - isample)/step)/videoFrameRate;

fprintf('Video final com %.2f s de duração, levando em conta um step de %d e um export de video em %d fps.\n', Tvideo, step, videoFrameRate);

clear   Tvideo;

%% VIDEO - TANK CREATING VIDEO
% Create the VideoWriter object
writerObj = VideoWriter(videoFileName, 'MPEG-4');
writerObj.FrameRate = 30;

% Open the video writer
open(writerObj);
% Write frames to the video
for i = 1:length(F)
    % Extract the current frame
    frame = F(i);
    % Write the current frame to the video
    writeVideo(writerObj, frame);
end
% Close the video writer
close(writerObj);

% Display a message indicating successful video creation
disp(['Video successfully saved as: ', videoFileName]);
%choosing the electrodes to plot
el1 = 10; %(RA)
el2 = 75; %(LA)
el3 = 28; %(V)
el4 = 129; %(Tank - High)
el5 = 160; %(Tank - Mid)
el6 = 190; %(Tank - Low)

% General Information
videoFileName = 'teste';
dataFiltered = Data;
lim = [-50 50];
[V_TANK, Plane, tank_plane_indx] = potentialMatrix(dataFiltered, 4);


%Creating Video
figure;
set(gcf, 'color', 'white');
colormap(jet)
fr=1;
for ii=isample:step:fsample

    % Subplot 1
    subplot(5, 2, 1); 
    plot(dataFiltered(el1, isample:fsample));
    title(['Electrode ', num2str(el1), ' (Right Atrium)']);
    ylabel('uV');
    xlim([0 fsample-isample]);
    xline((ii - isample), 'r');
    

    % Subplot 2
    subplot(5, 2, 3); 
    plot(dataFiltered(el2,isample:fsample));
    title(['Electrode ', num2str(el2), ' (Left Atrium)']);
    ylabel('uV');
    xlim([0 fsample-isample]);
    xline((ii - isample), 'r');

    % Subplot 3
    subplot(5, 2, 5); 
    plot(dataFiltered(el3,isample:fsample));
    title(['Electrode ', num2str(el3), ' (Ventricle)']);
    ylabel('uV');
    xlim([0 fsample-isample]);
    xlabel('Time (s)');
    xline((ii - isample), 'r');

    % Subplot 4
    subplot(5, 2, 2); 
    plot(dataFiltered(el4,isample:fsample)); %time,dataFiltered(el4,tplotin*Fs:tplotfin*Fs));
    title(['Electrode ', num2str(el4), ' (Tank - High)']);
    ylabel('uV');
    xlim([0 fsample-isample]);
    xlabel('Time (s)');
    xline((ii - isample), 'r');

    % Subplot 5
    subplot(5, 2, 4); 
    plot(dataFiltered(el5,isample:fsample)); %time,dataFiltered(el4,tplotin*Fs:tplotfin*Fs));
    title(['Electrode ', num2str(el5), ' (Tank - Mid)']);
    ylabel('uV');
    xlim([0 fsample-isample]);
    xlabel('Time (s)');
    xline((ii - isample), 'r');

    % Subplot 6
    subplot(5, 2, 6); 
    plot(dataFiltered(el6,isample:fsample)); %time,dataFiltered(el4,tplotin*Fs:tplotfin*Fs));
    title(['Electrode ', num2str(el6), ' (Tank - Low)']);
    ylabel('uV');
    xlim([0 fsample-isample]);
    xlabel('Time (s)');
    xline((ii - isample), 'r');

    % subplot 7
    subplot(5, 2, [7:10]);
    plotTank_video(V_TANK, ii, Plane, tank_plane_indx, lim);
    

    F(fr) = getframe(gcf) ;
    pause(0.01);
    fr=fr+1;%
end

% Create the VideoWriter object
writerObj = VideoWriter(videoFileName, 'MPEG-4');
writerObj.FrameRate = videoFrameRate;

% Open the video writer
open(writerObj);
% Write frames to the video
for i = 1:length(F)
    % Extract the current frame
    frame = F(i);
    % Write the current frame to the video
    writeVideo(writerObj, frame);
end
% Close the video writer
close(writerObj);

% Display a message indicating successful video creation
disp(['Video successfully saved as: ', videoFileName]);


%% VIDEO - MEAS CREATING VIDEO

%electrodes Selection
el1 = 10; %(RA)
el2 = 75; %(LA)
el3 = 28; %(V)
el4 = 129; %(Tank)

% general info
videoFileName = 'teste';
dataFiltered = Data;
lim = [-800 800];
[V_MEA1, MEA, MEA1_plane_indx] = potentialMatrix(dataFiltered, 1);
[V_MEA2, ~, MEA2_plane_indx] = potentialMatrix(dataFiltered, 2);
[V_MEA3, ~, MEA3_plane_indx] = potentialMatrix(dataFiltered, 3);

%plot
figure();
set(gcf, 'color', 'white');
colormap(jet)
fr=1;
for ii=isample:step:fsample
    
    % Subplot 1
    subplot(4, 6, [1,2,3]); 
    plot(dataFiltered(el1, isample:fsample));
    title(['Electrode ', num2str(el1), ' (Right Atrium)']);
    ylabel('uV');
    xlim([0 fsample-isample]);
    xline(ii-isample, 'r');
    
    % Subplot 2
    subplot(4, 6, [4,5,6]); 
    plot(dataFiltered(el2,isample:fsample));
    title(['Electrode ', num2str(el2), ' (Left Atrium)']);
    ylabel('uV');
    xlim([0 fsample-isample]);
    xline(ii-isample, 'r');

    % Subplot 3
    subplot(4, 6, [7,8,9]); 
    plot(dataFiltered(el3,isample:fsample));
    title(['Electrode ', num2str(el3), ' (Ventricle)']);
    ylabel('uV');
    xlim([0 fsample-isample]);
    xlabel('Time (s)');
    xline(ii-isample, 'r');

    % Subplot 4
    subplot(4, 6, [10,11,12]); 
    plot(dataFiltered(el4,isample:fsample));
    title(['Electrode ', num2str(el4), ' (Tank)']);
    ylabel('uV');
    xlim([0 fsample-isample]);
    xlabel('Time (s)');
    xline(ii-isample, 'r');

    %subplot 5
    subplot(4, 6, [13,14,19,20]);
    plotMEA_video(V_MEA1, ii, MEA, MEA1_plane_indx, lim);

    %subplot 6
    subplot(4, 6, [15,16,21,22]);
    plotMEA_video(V_MEA2, ii, MEA, MEA2_plane_indx, lim);

    %subplot 7
    subplot(4, 6, [17,18,23,24]);
    plotMEA_video(V_MEA3, ii, MEA, MEA3_plane_indx, lim);

    %config
    F(fr) = getframe(gcf) ;
    pause(0.01)
    fr=fr+1;
end

% Create the VideoWriter object
writerObj = VideoWriter(videoFileName, 'MPEG-4');
writerObj.FrameRate = videoFrameRate;

% Open the video writer
open(writerObj);
% Write frames to the video
for i = 1:length(F)
    % Extract the current frame
    frame = F(i);
    % Write the current frame to the video
    writeVideo(writerObj, frame);
end
% Close the video writer
close(writerObj);

% Display a message indicating successful video creation
disp(['Video successfully saved as: ', videoFileName]);

