%% Optic Potential Map

clear; clc;


%% Loading Variables

load('E:\HEartLab\TAINAN WORKFLOW\00 - examples\data_filtered_sync_E14_F3_R4.mat'); % Synchronized data
% Use the IMPLAY to analyse frames and pic a sample or sample interval to
% use un your plots.


%% Potential map
% Select sample and camera to plot
sample = 5400;
Data_O = D_SYNC.CAM1;

% Initialize color map and adjust first color to white
C = jet(256);
C(1,1:3) = [1 1 1];
% Extract a single frame from the processed data and rotate it 90 degrees
I = squeeze(Data_O(:,:,sample));
J = imrotate(I, 90);
% Create a figure and display the rotated image with a color map
f1 = figure('color', 'white', 'Position', [50 50 500 500]);
imagesc(J, [0 10]);
colormap(C);
% Add a colorbar and adjust labels
hBar1 = colorbar('eastoutside');
ylabel(hBar1, 'Fluorescence [%]', 'FontSize', 18);
set(gca, 'fontsize', 18);
title('Cam ');
ylabel('Pixels');
xlabel('Pixels');
axis off


%% Potential video map

% Define the start and end samples for the video
start_sample = 1;  % Adjust the start sample according to your data
end_sample = 4000;  % Adjust the end sample according to your data
Fsampling = 4000;

% Define the frames per second (fps) for the video
fps = 30;  % 30 fps is a standard value to video codec
% Define the step for plot
step = 2;

% Define the title of the video
video_title = 'optic_potential_map_video';

% Calculate the time duration of the video
duration = length(start_sample:step:end_sample) / fps;
disp(['Predicted time: ', num2str(duration), ' s']);
disp(['Predicted time: ', num2str(round(duration/60, 1)), ' min']);

% Create a VideoWriter object
video_file = VideoWriter([video_title, '.mp4'], 'MPEG-4');
video_file.FrameRate = fps;
% Open the VideoWriter object
open(video_file);
% Loop through the selected samples and create frames for the video
for sample = start_sample:step:end_sample
    % Extract a single frame from the processed data and rotate it 90 degrees
    I = squeeze(Data_O(:,:,sample));
    J = imrotate(I, 90);
    % Create a figure and display the rotated image with a color map
    f1 = figure('color', 'white', 'Position', [50 50 500 500]);
    imagesc(J, [0 10]);
    colormap(C);
    % Add a colorbar and adjust labels
    hBar1 = colorbar('eastoutside');
    ylabel(hBar1, 'Fluorescence [%]', 'FontSize', 18);
    set(gca, 'fontsize', 18);
    title(['Cam 3 - Frame: ', num2str(sample)]);
    ylabel('Pixels');
    xlabel('Pixels');
    axis off
    % Convert the figure to a frame
    frame = getframe(f1); 
    % Write the frame to the video file
    writeVideo(video_file, frame);
    % Close the figure to prevent unnecessary figures from being created
    close(f1);
end
% Close the VideoWriter object
close(video_file);

disp(['Video exported: ', video_title, '.mp4']);
disp(['Video duration: ', num2str(duration), ' seconds']);