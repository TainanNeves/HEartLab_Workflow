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


%% TANK


%% TANK - Loading Data

%loading the 3d tank geometry
load('Tank_geometry.mat');

%choosing the number of the tank electrodes
tank_list_indx =[145 146 155 156 165 166 129 130 139 140 181 182, ...
            147 157 167 131 141 183, ...
            148 149 158 159 168 169 132 133 142 143 184 185, ...
            150 151 160 161 170 171 134 135 144 177 186 187, ...
            152 162 172 136 178 188, ...
            153 154 163 164 173 174 137 138 179 180 189 190];
%choosing the positions that represents each electrode
tank_plane_indx =[77:2:99, ...
            153:4:173, ...
            227:2:249, ...
            402:2:424, ...
            478:4:498, ...
            552:2:574];


%% TANK - Plane Patch
%Tank
figure;
set(gcf, 'color', 'white');
patch('faces', Plane.faces, 'vertices', Plane.vertices,'FaceVertexCData', ones(1201,1), 'FaceColor', 'interp', ...
        'FaceAlpha', 1, 'EdgeAlpha', 1, ...
        'FaceLighting','gouraud');
axis equal
axis off
hold on
title('Tank electrodes plane patch');
plot3(Plane.vertices(tank_plane_indx,1),Plane.vertices(tank_plane_indx,2),Plane.vertices(tank_plane_indx,3),'or','LineWidth',2)


%% TANK - Laplacian interpolator
[lap,edge] = mesh_laplacian(Plane.vertices,Plane.faces);
[int, keepindex, repindex] = mesh_laplacian_interp(lap,tank_plane_indx);


%% TANK - Interpolation
dataFiltered = Data;
V_tank = int * sgolayfilt(dataFiltered(tank_list_indx,:),2,201,[],2);
V_tank = int * sgolayfilt(dataFiltered(tank_list_indx,:),2,131,[],2);


%% TANK - Individual Plot
%Time selection
tPlot = 4.5;%s
Fsampling = 4000;%Hz

%Plot
figure;
set(gcf, 'color', 'white');
colormap(jet);
patch('faces', Plane.faces, 'vertices', Plane.vertices,'FaceVertexCData', V_tank(:,tPlot*Fsampling), 'FaceColor', 'interp', ...
            'FaceAlpha', 1, 'EdgeAlpha', 0, ...
            'FaceLighting','gouraud');
title(['TANK', ' - ',num2str(tPlot), ' s']);
axis equal
axis off
clim([-50 50]);
colorbar;


%% TANK - Calculating time of video
tin = 4; % s
tfin = 4.5; % s
Fsampling = 4000; % Hz
step_plot = 4; % sample step between plots
videoFrameRate = 30; % fps

Tvideo = (((tfin - tin) * Fsampling) / step_plot) / videoFrameRate;

mensagem = sprintf('Video final com %.2f s de duração, levando em conta um step de %d e um export de video em %d fps.', Tvideo, step_plot, videoFrameRate);
disp(mensagem);
clear   Fsampling step_plot videoFrameRate Tvideo mensagem;


%% TANK - Generating Tank Video

% Specify time window
tplotin = tin;
tplotfin = tfin;

%choosing the electrodes to plot
el1 = 10; %(RA)
el2 = 75; %(LA)
el3 = 28; %(V)
el4 = 129; %(Tank - High)
el5 = 160; %(Tank - Mid)
el6 = 190; %(Tank - Low)
Fsampling = 4000;

%Creating Video
figure;
set(gcf, 'color', 'white');
colormap(jet)
fr=1;
for ii=tplotin*Fsampling:2:tplotfin*Fsampling
    time = linspace(tplotin, tplotfin, length(tplotin*Fsampling : tplotfin*Fsampling));

    % Subplot 1
    subplot(5, 2, 1); 
    plot(time, dataFiltered(el1, tplotin*Fsampling:tplotfin*Fsampling));
    title(['Electrode ', num2str(el1), ' (Right Atrium)']);
    ylabel('uV');
    xline(ii/Fsampling, 'r');

    % Subplot 2
    subplot(5, 2, 3); 
    plot(time, dataFiltered(el2,tplotin*Fsampling:tplotfin*Fsampling));
    title(['Electrode ', num2str(el2), ' (Left Atrium)']);
    ylabel('uV');
    xline(ii/Fsampling, 'r');

    % Subplot 3
    subplot(5, 2, 5); 
    plot(time, dataFiltered(el3,tplotin*Fsampling:tplotfin*Fsampling));
    title(['Electrode ', num2str(el3), ' (Ventricle)']);
    ylabel('uV');
    xlabel('Time (s)');
    xline(ii/Fsampling, 'r');

    % Subplot 4
    subplot(5, 2, 2); 
    plot( time,V_tank(7,tplotin*Fsampling:tplotfin*Fsampling)); %time,dataFiltered(el4,tplotin*Fs:tplotfin*Fs));
    title(['Electrode ', num2str(el4), ' (Tank - High)']);
    ylabel('uV');
    xlabel('Time (s)');
    xline(ii/Fsampling, 'r');

    % Subplot 5
    subplot(5, 2, 4); 
    plot( time,V_tank(33,tplotin*Fsampling:tplotfin*Fsampling)); %time,dataFiltered(el4,tplotin*Fs:tplotfin*Fs));
    title(['Electrode ', num2str(el5), ' (Tank - Mid)']);
    ylabel('uV');
    xlabel('Time (s)');
    xline(ii/Fsampling, 'r');

    % Subplot 6
    subplot(5, 2, 6); 
    plot( time,V_tank(60,tplotin*Fsampling:tplotfin*Fsampling)); %time,dataFiltered(el4,tplotin*Fs:tplotfin*Fs));
    title(['Electrode ', num2str(el6), ' (Tank - Low)']);
    ylabel('uV');
    xlabel('Time (s)');
    xline(ii/Fsampling, 'r');

    % subplot 7
    subplot(5, 2, [7:10]);
    patch('faces', Plane.faces, 'vertices', Plane.vertices,'FaceVertexCData', V_tank(:,ii), 'FaceColor', 'interp', ...
            'FaceAlpha', 1, 'EdgeAlpha', 0, ...
            'FaceLighting','gouraud');
    axis equal
    axis off
    clim([-50 100]);
    colorbar;

    F(fr) = getframe(gcf) ;
    pause(0.01);
    fr=fr+1;
end


%% TANK - Video Export

% Specify the video file name
videoFileName = 'teste';
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


%% MEAs


%% Loading data

%loading the 3d MEAs geometry
load('MEA.mat');

%Choosing the number of the MEAs electrodes
%MEA1
MEA1_list_indx = [13 14 15 16, ...
                9 10 11 12, ...
                5 6 7 8, ...
                1 2 3 4];
%MEA2
MEA2_list_indx = [17 21 25 29, ...
                18 22 26 30, ...
                19 23 27 31, ...
                20 24 28 32];
%MEA3
MEA3_list_indx = [77 78 79 80, ...
                73 74 75 76, ...
                69 70 71 72, ...
                65 66 67 68];
%Choosing th    e position that represents each electrode
MEA_plane_indx = [25:2:31, ...
                47:2:53, ...
                69:2:75, ...
                91:2:97];


%% MEAs - Plane Patch

figure;
set(gcf, 'color', 'white');
patch('faces', MEA.faces, 'vertices', MEA.vertices,'FaceVertexCData', ones(221,1), 'FaceColor', 'interp', ...
        'FaceAlpha', 1, 'EdgeAlpha', 1, ...
        'FaceLighting','gouraud');
axis equal
axis off
hold on
title('MEAs electrodes plane patch');
plot3(MEA.vertices(MEA_plane_indx,1),MEA.vertices(MEA_plane_indx,2),MEA.vertices(MEA_plane_indx,3),'or','LineWidth',2)


%% MEAs - Laplacian interpolator
[lap,edge] = mesh_laplacian(MEA.vertices,MEA.faces);

%MEA1
MEA1_plane_indx = MEA_plane_indx;
MEA1_plane_indx([1 8]) = [];%bad electrodes
[int1, keepindex1, repindex1] = mesh_laplacian_interp(lap,MEA1_plane_indx);

%MEA2
MEA2_plane_indx = MEA_plane_indx;
[int2, keepindex2, repindex2] = mesh_laplacian_interp(lap,MEA2_plane_indx);

%MEA3
MEA3_plane_indx = MEA_plane_indx;
MEA3_plane_indx([4]) = [];%bad electrodes
[int3, keepindex3, repindex3] = mesh_laplacian_interp(lap,MEA3_plane_indx);


%% MEAs - Interpolation
dataFiltered = Data;

%MEA1
MEA1_list_indx([1 8]) = [];%bad electrodes
V_MEA1 = int1 * sgolayfilt(dataFiltered(MEA1_list_indx,:),2,201,[],2);

%MEA2
V_MEA2 = int2 * sgolayfilt(dataFiltered(MEA2_list_indx,:),2,201,[],2);

%MEA3
MEA3_list_indx([4]) = [];%bad electrodes
V_MEA3 = int3 * sgolayfilt(dataFiltered(MEA3_list_indx,:),2,201,[],2);


%% MEAs - Individual Plot
%MEA selection
% V_MEA = V_MEA1; t = 'MEA 1 (Right Atrium)';
% V_MEA = V_MEA2; t = 'MEA 2 (Ventricle)';
V_MEA = V_MEA3; t = 'MEA 3 (Left Atrium)';

%Time selection
tPlot = 4.5;%s
Fsampling = 4000;%Hz

%plot
figure;
set(gcf, 'color', 'white');
colormap(jet);
patch('faces', MEA.faces, 'vertices', MEA.vertices,'FaceVertexCData', V_MEA(:,tPlot*Fsampling), 'FaceColor', 'interp', ...
        'FaceAlpha', 1, 'EdgeAlpha', 0, ...
        'FaceLighting','gouraud');
title([t, ' - ',num2str(tPlot), ' s']);
axis equal
axis off
caxis([-800 800])
colorbar


%% MEAs - Calculating time of video

tin = 4; % s
tfin = 4.5; % s
Fsampling = 4000; % Hz
step_plot = 5; % sample step between plots
videoFrameRate = 30; % fps

Tvideo = (((tfin - tin) * Fsampling) / step_plot) / videoFrameRate;

mensagem = sprintf('Video final com %.2f s de duração, levando em conta um step de %d e um export de video em %d fps.', Tvideo, step_plot, videoFrameRate);
disp(mensagem);
clear Fsampling step_plot videoFrameRate Tvideo mensagem;


%% MEAs - Generation Video

%specific time
tplotin = tin;
tplotfin = tfin;

%electrodes
el1 = 10; %(RA)
el2 = 75; %(LA)
el3 = 28; %(V)
el4 = 129; %(Tank)
Fsampling = 4000;

%Limits
iLin = -800;
fLin = 800;

%plot
figure(5);
set(gcf, 'color', 'white');
colormap(jet)
fr=1;
for ii=tplotin*Fsampling:2:tplotfin*Fsampling
    time = linspace(tplotin, tplotfin, length(tplotin*Fsampling : tplotfin*Fsampling));
    % Subplot 1
    subplot(4, 6, [1,2,3]); 
    plot(time, dataFiltered(el1, tplotin*Fsampling:tplotfin*Fsampling));
    title(['Electrode ', num2str(el1), ' (Right Atrium)']);
    ylabel('uV');
    
    xline(ii/Fsampling, 'r');
    
    % Subplot 2
    subplot(4, 6, [4,5,6]); 
    plot(time, dataFiltered(el2,tplotin*Fsampling:tplotfin*Fsampling));
    title(['Electrode ', num2str(el2), ' (Left Atrium)']);
    ylabel('uV');
    xline(ii/Fsampling, 'r');

    % Subplot 3
    subplot(4, 6, [7,8,9]); 
    plot(time, dataFiltered(el3,tplotin*Fsampling:tplotfin*Fsampling));
    title(['Electrode ', num2str(el3), ' (Ventricle)']);
    ylabel('uV');
    xlabel('Time (s)');
    xline(ii/Fsampling, 'r');

    % Subplot 4
    subplot(4, 6, [10,11,12]); 
    plot(time,dataFiltered(el4,tplotin*Fsampling:tplotfin*Fsampling));
    title(['Electrode ', num2str(el4), ' (Tank)']);
    ylabel('uV');
    xlabel('Time (s)');
    xline(ii/Fsampling, 'r');

    %subplot 5
    subplot(4, 6, [13,14,19,20]);
    %clf
    patch('faces', MEA.faces, 'vertices', MEA.vertices,'FaceVertexCData', V_MEA1(:,ii), 'FaceColor', 'interp', ...
            'FaceAlpha', 1, 'EdgeAlpha', 0, ...
            'FaceLighting','gouraud');
    title('MEA1 (Right Atrium) ');
    axis equal
    axis off
    caxis([iLin fLin])
    colorbar

    %subplot 6
    subplot(4, 6, [15,16,21,22]);
    %clf
    patch('faces', MEA.faces, 'vertices', MEA.vertices,'FaceVertexCData', V_MEA2(:,ii), 'FaceColor', 'interp', ...
            'FaceAlpha', 1, 'EdgeAlpha', 0, ...
            'FaceLighting','gouraud');
    title('MEA2 (Ventricle) ');
    axis equal
    axis off
    caxis([iLin fLin])
    colorbar

    %subplot 7
    subplot(4, 6, [17,18,23,24]);
    %clf
    patch('faces', MEA.faces, 'vertices', MEA.vertices,'FaceVertexCData', V_MEA3(:,ii), 'FaceColor', 'interp', ...
            'FaceAlpha', 1, 'EdgeAlpha', 0, ...
            'FaceLighting','gouraud');
    title('MEA3 (Left Atrium) ');
    axis equal
    axis off
    caxis([iLin fLin])
    colorbar

    %config
    F(fr) = getframe(gcf) ;
    pause(0.01)
    fr=fr+1;
end


%% MEAs - Video Export
% Specify the video file name
videoFileName = 'teste';
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




















