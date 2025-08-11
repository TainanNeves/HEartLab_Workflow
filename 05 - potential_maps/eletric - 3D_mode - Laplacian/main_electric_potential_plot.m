%% Main - Potential Map Plots
%this code allows you to plot a 2d isopotential maps from TANK and MEAs
%using laplacian interpolator
%Authors: Angélica Quadros, Ismael Romero, Tainan Neves.

clear; clc;


%% Load data

%Load the filtered data
load('E:\HEartLab\TAINAN WORKFLOW\00 - examples\data_filtered_sync_E14_F3_R4.mat'); % Synchronized data


% Transform to mV (Exported default data are in uV)
% data = data/1000;


%% Previsualizing data from electric time view
% Previsualizing data from electric time view
Data = D_SYNC.EL;

% Defining parameters
el1 = 10; pos1 = 'RA';
el2 = 74; pos2 = 'LA';
el3 = 23; pos3 = 'V';
el4 = 142; pos4 = 'TANK';
tin = 3; % Initial time (s)
tfin = 4; % Final time (s)
Fsampling = 4000;
time = linspace(tin, tfin, length(tin*Fsampling:tfin*Fsampling));
redX = [4000];

figure();

% Plotting each subplot
for i = 1:4
    subplot(4, 1, i);
    plot(time, Data(eval(['el', num2str(i)]), tin*Fsampling:tfin*Fsampling));
    title(['Electrode ', num2str(eval(['el', num2str(i)])), ' (', eval(['pos', num2str(i)]), ')']);
    ylabel('$\mu$V', 'Interpreter', 'latex');
    xlabel('Time (s)');
    hold on;
    
    % Adding red line at redX points
    for x = redX
        xline(x,'--r',{[num2str(x)]});
    end
    hold off;
end

% Adjusting overall plot
sgtitle('Electrode Data');

clear el1 el2 el3 el4 pos1 pos2 pos3 pos4; 


%% POTENTIAL MAPS
lim1 = [-100 65];
lim2 = [-50 50];
lim3 = [-30 20];
lim4 = [-10 10];
for sample = [4000:10:4150]
    plot_electric_pot(Data, lim1, sample, 1); % MEA 1
    plot_electric_pot(Data, lim2, sample, 2); % MEA 2
    plot_electric_pot(Data, lim3, sample, 3); % MEA 3
    plot_electric_pot(Data, lim4, sample, 4); % TANK
end


%% VIDEO - calculating time of video
isample = 3.22*4000;
fsample = 3.35*4000;
step = 1; % sample step between plots
videoFrameRate = 30; % fps

Tvideo = ((fsample - isample)/step)/videoFrameRate;

fprintf('Video final com %.2f s de duração, levando em conta um step de %d e um export de video em %d fps.\n', Tvideo, step, videoFrameRate);

clear   Tvideo;

%% VIDEO - TANK CREATING VIDEO
%Electrodes Selection
el1 = 10; %(RA)
el2 = 74; %(LA)
el3 = 23; %(V)
el4 = 131; %(Tank - High)
el5 = 136; %(Tank - Mid)
el6 = 172; %(Tank - Low)

% General Information
videoFileName = 'EXXFXXRXX - tank_video';
dataFiltered = Data;
lim = [-10 10];
[V_TANK, Plane, tank_plane_indx] = potentialMatrix(dataFiltered, 4);

% Plot
figure('color','white','Position', [40 40 1000 600]);
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
    xlabel(['Frame: ' num2str(ii)]);
    xline((ii - isample), 'r');

    % Subplot 4
    subplot(5, 2, 2); 
    plot(dataFiltered(el4,isample:fsample)); %time,dataFiltered(el4,tplotin*Fs:tplotfin*Fs));
    title(['Electrode ', num2str(el4), ' (Tank - High)']);
    ylabel('uV');
    xlim([0 fsample-isample]);
    xline((ii - isample), 'r');

    % Subplot 5
    subplot(5, 2, 4); 
    plot(dataFiltered(el5,isample:fsample)); %time,dataFiltered(el4,tplotin*Fs:tplotfin*Fs));
    title(['Electrode ', num2str(el5), ' (Tank - Mid)']);
    ylabel('uV');
    xlim([0 fsample-isample]);
    xline((ii - isample), 'r');

    % Subplot 6
    subplot(5, 2, 6); 
    plot(dataFiltered(el6,isample:fsample)); %time,dataFiltered(el4,tplotin*Fs:tplotfin*Fs));
    title(['Electrode ', num2str(el6), ' (Tank - Low)']);
    ylabel('uV');
    xlim([0 fsample-isample]);
    xlabel(['Frame: ' num2str(ii)]);
    xline((ii - isample), 'r');

    % subplot 7
    subplot(5, 2, [7:10]);
    plotTank_video(V_TANK, ii, Plane, tank_plane_indx, lim);
    
    %config
    F(fr) = getframe(gcf) ;
    pause(0.01);
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


%% VIDEO - MEAS CREATING VIDEO

%electrodes Selection
el1 = 10; %(RA)
el2 = 74; %(LA)
el3 = 23; %(V)
el4 = 136; %(Tank)

% general info
videoFileName = 'EXXFXXRXX - MEAs_video';
dataFiltered = Data;
lim1 = [-100 60];
lim2 = [-40 40];
lim3 = [-30 20];
[V_MEA1, MEA, MEA1_plane_indx] = potentialMatrix(dataFiltered, 1);
[V_MEA2, ~, MEA2_plane_indx] = potentialMatrix(dataFiltered, 2);
[V_MEA3, ~, MEA3_plane_indx] = potentialMatrix(dataFiltered, 3);

%plot
figure('color','white','Position', [40 40 1000 600]);
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
    xlabel(['Frame: ' num2str(ii)]);
    xline(ii-isample, 'r');

    % Subplot 4
    subplot(4, 6, [10,11,12]); 
    plot(dataFiltered(el4,isample:fsample));
    title(['Electrode ', num2str(el4), ' (Tank)']);
    ylabel('uV');
    xlim([0 fsample-isample]);
    xlabel(['Frame: ' num2str(ii)]);
    xline(ii-isample, 'r');

    % Subplot 5
    subplot(4, 6, [13,14,19,20]);
    plotMEA_video(V_MEA1, ii, MEA, MEA1_plane_indx, lim1);
    text(0.5, -0.1, 'MEA 1 - RA', 'Units', 'normalized', 'FontSize', 10, 'HorizontalAlignment', 'center');
    
    % Subplot 6
    subplot(4, 6, [15,16,21,22]);
    plotMEA_video(V_MEA2, ii, MEA, MEA2_plane_indx, lim2);
    text(0.5, -0.1, 'MEA 2 - V', 'Units', 'normalized', 'FontSize', 10, 'HorizontalAlignment', 'center');
    
    % Subplot 7
    subplot(4, 6, [17,18,23,24]);
    plotMEA_video(V_MEA3, ii, MEA, MEA3_plane_indx, lim3);
    text(0.5, -0.1, 'MEA 3 - LA', 'Units', 'normalized', 'FontSize', 10, 'HorizontalAlignment', 'center');



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
