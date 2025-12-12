%% Optic Potential Map
clear; clc;


%% Loading Variables
load("E:\Qualification\Analysis\E32F02R01\data\data_filtered_sync_E32_F02_R01.mat"); % Synchronized data


%% Overview Optical plot
% Loading variables
data_o1 = D_SYNC.CAM1;
data_o2 = D_SYNC.CAM2;
data_o3 = D_SYNC.CAM3;
data_e = D_SYNC.EL;

% Plot time
timepoint_seconds = 2.46; 
% Plot window
plot_window_seconds = 1; 

% Select 6 electrodes
electrodes = [7 22 90 142 178 172]; 

% Selecting Colorbar limits
vmin_cam1 = 0; 
vmax_cam1 = 10; 
vmin_cam2 = 0; 
vmax_cam2 = 10;
vmin_cam3 = 0;
vmax_cam3 = 10;

% Defining Fixed parameters
Fs = 4000;
sample = round(timepoint_seconds * Fs);
total_time = size(data_o1, 3) / Fs;

% --- CALCULATE PLOT WINDOW LIMITS ---
% Calculate the desired start and end times
t_start_desired = timepoint_seconds - (plot_window_seconds / 2);
t_end_desired = timepoint_seconds + (plot_window_seconds / 2);
% Apply boundary checks
t_start_plot = max(0, t_start_desired);
t_end_plot = min(total_time, t_end_desired);

% Selecting Optical points
Background_CAM1 = squeeze(imrotate(D_SYNC.IMG.CAM1(:,:,1), -90));
Background_CAM2 = squeeze(imrotate(D_SYNC.IMG.CAM2(:,:,1), -90));
Background_CAM3 = squeeze(imrotate(D_SYNC.IMG.CAM3(:,:,1), -90));
disp('Select 2 points for CAM1');
[x1, y1] = pick_up_a_trace(Background_CAM1, data_o1, 1);
disp('Select 2 points for CAM2');
[x2, y2] = pick_up_a_trace(Background_CAM2, data_o2, 1);
disp('Select 2 points for CAM3');
[x3, y3] = pick_up_a_trace(Background_CAM3, data_o3, 1);

% Define the subplot grid
rows = 6;
cols = 6;
% Create a figure
figure('color', 'white');
time_axis = (0:size(data_o1, 3)-1) / Fs;

% --- TIME PLOTS ---

% CAM1 Point 1
subplot(rows, cols, 1);
plot(time_axis, squeeze(data_o1(x1(1), y1(1), :))); 
hold on; 
xline(timepoint_seconds, 'r', 'LineWidth', 2);
xlim([t_start_plot t_end_plot]); % <<< APPLIED WINDOW
title(sprintf('CAM1 Pt 1\n(x=%d, y=%d)', x1(1), y1(1))); 
xlabel('Time (s)'); ylabel('%F'); 

% CAM1 Point 2
subplot(rows, cols, 2);
plot(time_axis, squeeze(data_o1(x1(2), y1(2), :))); 
hold on; 
xline(timepoint_seconds, 'r', 'LineWidth', 2);
xlim([t_start_plot t_end_plot]); % <<< APPLIED WINDOW
title(sprintf('CAM1 Pt 2\n(x=%d, y=%d)', x1(2), y1(2))); 
xlabel('Time (s)'); ylabel('%F'); 

% CAM2 Point 1
subplot(rows, cols, 3);
plot(time_axis, squeeze(data_o2(x2(1), y2(1), :))); 
hold on; 
xline(timepoint_seconds, 'r', 'LineWidth', 2);
xlim([t_start_plot t_end_plot]); % <<< APPLIED WINDOW
title(sprintf('CAM2 Pt 1\n(x=%d, y=%d)', x2(1), y2(1))); 
xlabel('Time (s)'); ylabel('%F'); 

% CAM2 Point 2
subplot(rows, cols, 4);
plot(time_axis, squeeze(data_o2(x2(2), y2(2), :))); 
hold on; 
xline(timepoint_seconds, 'r', 'LineWidth', 2);
xlim([t_start_plot t_end_plot]); % <<< APPLIED WINDOW
title(sprintf('CAM2 Pt 2\n(x=%d, y=%d)', x2(2), y2(2))); 
xlabel('Time (s)'); ylabel('%F'); 

% CAM3 Point 1
subplot(rows, cols, 5); 
plot(time_axis, squeeze(data_o3(x3(1), y3(1), :))); 
hold on; 
xline(timepoint_seconds, 'r', 'LineWidth', 2);
xlim([t_start_plot t_end_plot]); % <<< APPLIED WINDOW
title(sprintf('CAM3 Pt 1\n(x=%d, y=%d)', x3(1), y3(1))); 
xlabel('Time (s)'); ylabel('%F'); 

% CAM3 Point 2
subplot(rows, cols, 6); 
plot(time_axis, squeeze(data_o3(x3(2), y3(2), :))); 
hold on; 
xline(timepoint_seconds, 'r', 'LineWidth', 2);
xlim([t_start_plot t_end_plot]); % <<< APPLIED WINDOW
title(sprintf('CAM3 Pt 2\n(x=%d, y=%d)', x3(2), y3(2))); 
xlabel('Time (s)'); ylabel('%F'); 

% Electrical Plots
for k = 1:6 % Loop through all 6 electrodes
    subplot_index = k + 6; 
    subplot(rows, cols, subplot_index);
    plot(time_axis, data_e(electrodes(k), :), 'k');
    hold on; 
    xline(timepoint_seconds, 'r', 'LineWidth', 2);
    xlim([t_start_plot t_end_plot]); % <<< APPLIED WINDOW
    title(sprintf('Electrode %d', electrodes(k))); 
    xlabel('Time (s)'); ylabel('\muV'); 
end

% --- OPTICAL FRAME PLOTS ---

% CAM1 Frame
subplot(rows, cols, [13 14 19 20 25 26 31 32]);
data_o1 = imrotate(data_o1, 90);
imagesc(data_o1(:, :, sample), [vmin_cam1 vmax_cam1]);
colormap(jet(256));
hold on;
% Rotatinf the markers position
x1_temp = x1;
y1_temp = y1;
y1 = x1_temp;
x1 = size(data_o1, 1) - y1_temp + 1;
clear x1_temp y1_temp;
% Ploting markers
plot(y1, x1, 'ro', 'MarkerSize', 10, 'LineWidth', 2); 
text(y1(1)+1, x1(1)-1, ' Point 1', 'Color', 'r', 'FontSize', 10, 'VerticalAlignment', 'bottom'); 
text(y1(2)+1, x1(2)-1, ' Point 2', 'Color', 'r', 'FontSize', 10, 'VerticalAlignment', 'bottom'); 
title('CAM1');
colorbar;

% CAM2 Frame
subplot(rows, cols, [15 16 21 22 27 28 33 34]);
data_o2 = imrotate(data_o2, 90);
imagesc(data_o2(:, :, sample), [vmin_cam2 vmax_cam2]);
colormap(jet(256));
hold on;
% Rotatinf the markers position
x2_temp = x2;
y2_temp = y2;
y2 = x2_temp;
x2 = size(data_o2, 1) - y2_temp + 1;
clear x2_temp y2_temp;
% Ploting markers
plot(y2, x2, 'ro', 'MarkerSize', 10, 'LineWidth', 2); 
text(y2(1)+1, x2(1)-1, ' Point 1', 'Color', 'r', 'FontSize', 10, 'VerticalAlignment', 'bottom'); 
text(y2(2)+1, x2(2)-1, ' Point 2', 'Color', 'r', 'FontSize', 10, 'VerticalAlignment', 'bottom'); 
title('CAM2');
colorbar;

% CAM3 Frame
subplot(rows, cols, [17 18 23 24 29 30 35 36]);
data_o3 = imrotate(data_o3, 90);
imagesc(data_o3(:, :, sample), [vmin_cam3 vmax_cam3]);
colormap(jet(256));
hold on;
% Rotatinf the markers position
x3_temp = x3;
y3_temp = y3;
y3 = x3_temp;
x3 = size(data_o3, 1) - y3_temp + 1;
clear x3_temp y3_temp;
% Ploting markers
plot(y3, x3, 'ro', 'MarkerSize', 10, 'LineWidth', 2); 
text(y3(1)+1, x3(1)-1, ' Point 1', 'Color', 'r', 'FontSize', 10, 'VerticalAlignment', 'bottom'); 
text(y3(2)+1, x3(2)-1, ' Point 2', 'Color', 'r', 'FontSize', 10, 'VerticalAlignment', 'bottom'); 
title('CAM3');
colorbar;

% Adjust overall figure appearance
sgtitle(['Synchronized Signals at t = ' num2str(timepoint_seconds) ' s '], ...
            'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'latex');


%% Overview Optical video
data_o1 = D_SYNC.CAM1;
data_o2 = D_SYNC.CAM2;
data_o3 = D_SYNC.CAM3;
data_e = D_SYNC.EL;

% Select 6 electrodes
electrodes = [7 22 90 142 178 172]; 

% Time parameters
ti = 2.274; % Start time [s]
to = 2.635; % End time [s]
step = 0.001; % Step size [s]
frame_rate = 30; % Video frame rate (FPS)

% --- WINDOW PARAMETER ---
plot_window_seconds = 1; % The total duration shown in the time plots

% Video Output filename
videoFile = 'O_overview_video.mp4';

% Colorbar limits for image frames
clim_cam1 = [0 10]; 
clim_cam2 = [0 10];
clim_cam3 = [0 10]; 

% Fixed parameters
Fs = 4000;
time_axis = (0:size(data_o1, 3)-1) / Fs;
total_time = size(data_o1, 3) / Fs;
rows = 6;
cols = 6;

% --- Calculate and Confirm Duration ---
total_frames = length(ti:step:to);
predicted_video_duration = total_frames / frame_rate;
disp(['Total Frames to process: ', num2str(total_frames)]);
disp(['Predicted Video Duration at ' num2str(frame_rate) ' FPS: ', num2str(predicted_video_duration), ' s']);
% Prompt user to confirm
user_input = input('Do you want to proceed with these values? (y/n): ', 's');
if lower(user_input) ~= 'y'
    fprintf('Please modify the values for ti, to, step, and frame_rate in the code and rerun the section.\n');
    return;
end

% --- Select Optical Points and Initialize Figure ---
Background_CAM1 = squeeze(imrotate(D_SYNC.IMG.CAM1(:,:,1), -90));
Background_CAM2 = squeeze(imrotate(D_SYNC.IMG.CAM2(:,:,1), -90));
Background_CAM3 = squeeze(imrotate(D_SYNC.IMG.CAM3(:,:,1), -90));
disp('Select 2 points for CAM1');
[x1, y1] = pick_up_a_trace(Background_CAM1, data_o1, 1);
disp('Select 2 points for CAM2');
[x2, y2] = pick_up_a_trace(Background_CAM2, data_o2, 1);
disp('Select 2 points for CAM3');
[x3, y3] = pick_up_a_trace(Background_CAM3, data_o3, 1);

figure('color', 'white', 'Name', 'Optical/Electric Overview Video');
set(gcf, 'Position', get(0, 'ScreenSize'));

% --- Initialize Video Writer (MP4 Format) ---
v = VideoWriter(videoFile, 'MPEG-4'); 
v.FrameRate = frame_rate; 
open(v);

% --- Main Video Generation Loop ---
for timepoint_seconds = ti:step:to
    sample = round(timepoint_seconds * Fs);
    
    % --- DYNAMICALLY CALCULATE PLOT WINDOW LIMITS ---
    t_start_desired = timepoint_seconds - (plot_window_seconds / 2);
    t_end_desired = timepoint_seconds + (plot_window_seconds / 2);
    
    % Apply boundary checks
    t_start_plot = max(0, t_start_desired);
    t_end_plot = min(total_time, t_end_desired);
    
    % Rotate the current image frame for plotting (0 for no rotation)
    I1 = imrotate(data_o1(:, :, sample), 0); 
    I2 = imrotate(data_o2(:, :, sample), 0);
    I3 = imrotate(data_o3(:, :, sample), 0);
    
    % --- ROW 1: Optical Traces (APPLYING XLIM) ---
    % CAM1 Point 1
    subplot(rows, cols, 1);
    plot(time_axis, squeeze(data_o1(x1(1), y1(1), :)));
    hold on; xline(timepoint_seconds, 'r', 'LineWidth', 2); hold off;
    xlim([t_start_plot t_end_plot]);
    title(sprintf('CAM1 Pt 1\n(x=%d, y=%d)', x1(1), y1(1))); 
    xlabel('Time (s)'); ylabel('%F');
    
    % CAM1 Point 2
    subplot(rows, cols, 2);
    plot(time_axis, squeeze(data_o1(x1(2), y1(2), :)));
    hold on; xline(timepoint_seconds, 'r', 'LineWidth', 2); hold off;
    xlim([t_start_plot t_end_plot]);
    title(sprintf('CAM1 Pt 2\n(x=%d, y=%d)', x1(2), y1(2))); 
    xlabel('Time (s)'); ylabel('%F');
    
    % CAM2 Point 1
    subplot(rows, cols, 3);
    plot(time_axis, squeeze(data_o2(x2(1), y2(1), :)));
    hold on; xline(timepoint_seconds, 'r', 'LineWidth', 2); hold off;
    xlim([t_start_plot t_end_plot]); 
    title(sprintf('CAM2 Pt 1\n(x=%d, y=%d)', x2(1), y2(1))); 
    xlabel('Time (s)'); ylabel('%F');
    
    % CAM2 Point 2
    subplot(rows, cols, 4);
    plot(time_axis, squeeze(data_o2(x2(2), y2(2), :)));
    hold on; xline(timepoint_seconds, 'r', 'LineWidth', 2); hold off;
    xlim([t_start_plot t_end_plot]);
    title(sprintf('CAM2 Pt 2\n(x=%d, y=%d)', x2(2), y2(2))); 
    xlabel('Time (s)'); ylabel('%F');
    
    % CAM3 Point 1
    subplot(rows, cols, 5);
    plot(time_axis, squeeze(data_o3(x3(1), y3(1), :)));
    hold on; xline(timepoint_seconds, 'r', 'LineWidth', 2); hold off;
    xlim([t_start_plot t_end_plot]);
    title(sprintf('CAM3 Pt 1\n(x=%d, y=%d)', x3(1), y3(1))); 
    xlabel('Time (s)'); ylabel('%F');
    
    % CAM3 Point 2
    subplot(rows, cols, 6);
    plot(time_axis, squeeze(data_o3(x3(2), y3(2), :)));
    hold on; xline(timepoint_seconds, 'r', 'LineWidth', 2); hold off;
    xlim([t_start_plot t_end_plot]);
    title(sprintf('CAM3 Pt 2\n(x=%d, y=%d)', x3(2), y3(2))); 
    xlabel('Time (s)'); ylabel('%F');
    
    % --- ROW 2: Electrical Traces (APPLYING XLIM) ---
    for k = 1:6 
        subplot_index = k + 6; 
        subplot(rows, cols, subplot_index);
        plot(time_axis, data_e(electrodes(k), :), 'k');
        hold on; xline(timepoint_seconds, 'r', 'LineWidth', 2); hold off;
        xlim([t_start_plot t_end_plot]);
        title(sprintf('Electrode %d', electrodes(k))); 
        xlabel('Time (s)'); ylabel('\muV'); 
    end
    
    % --- Optical Frames ---
    % CAM1 Frame
    subplot(rows, cols, [13 14 19 20 25 26 31 32]);
    data_o1 = imrotate(data_o1, 90);
    imagesc(data_o1(:, :, sample), [vmin_cam1 vmax_cam1]);
    colormap(jet(256));
    hold on;
    % Rotatinf the markers position
    x1_temp = x1;
    y1_temp = y1;
    y1 = x1_temp;
    x1 = size(data_o1, 1) - y1_temp + 1;
    clear x1_temp y1_temp;
    % Ploting markers
    plot(y1, x1, 'ro', 'MarkerSize', 10, 'LineWidth', 2); 
    text(y1(1)+1, x1(1)-1, ' Point 1', 'Color', 'r', 'FontSize', 10, 'VerticalAlignment', 'bottom'); 
    text(y1(2)+1, x1(2)-1, ' Point 2', 'Color', 'r', 'FontSize', 10, 'VerticalAlignment', 'bottom'); 
    title('CAM1');
    colorbar;
    hold off;
    
    % CAM2 Frame
    subplot(rows, cols, [15 16 21 22 27 28 33 34]);
    data_o2 = imrotate(data_o2, 90);
    imagesc(data_o2(:, :, sample), [vmin_cam2 vmax_cam2]);
    colormap(jet(256));
    hold on;
    % Rotatinf the markers position
    x2_temp = x2;
    y2_temp = y2;
    y2 = x2_temp;
    x2 = size(data_o2, 1) - y2_temp + 1;
    clear x2_temp y2_temp;
    % Ploting markers
    plot(y2, x2, 'ro', 'MarkerSize', 10, 'LineWidth', 2); 
    text(y2(1)+1, x2(1)-1, ' Point 1', 'Color', 'r', 'FontSize', 10, 'VerticalAlignment', 'bottom'); 
    text(y2(2)+1, x2(2)-1, ' Point 2', 'Color', 'r', 'FontSize', 10, 'VerticalAlignment', 'bottom'); 
    title('CAM2');
    colorbar;
    hold off;
    
    % CAM3 Frame
    subplot(rows, cols, [17 18 23 24 29 30 35 36]);
    data_o3 = imrotate(data_o3, 90);
    imagesc(data_o3(:, :, sample), [vmin_cam3 vmax_cam3]);
    colormap(jet(256));
    hold on;
    % Rotatinf the markers position
    x3_temp = x3;
    y3_temp = y3;
    y3 = x3_temp;
    x3 = size(data_o3, 1) - y3_temp + 1;
    clear x3_temp y3_temp;
    % Ploting markers
    plot(y3, x3, 'ro', 'MarkerSize', 10, 'LineWidth', 2); 
    text(y3(1)+1, x3(1)-1, ' Point 1', 'Color', 'r', 'FontSize', 10, 'VerticalAlignment', 'bottom'); 
    text(y3(2)+1, x3(2)-1, ' Point 2', 'Color', 'r', 'FontSize', 10, 'VerticalAlignment', 'bottom'); 
    title('CAM3');
    colorbar;
    hold off;
    
    sgtitle(['Overview Optical Video | t=' num2str(timepoint_seconds, '%.3f') ' s'], ...
                    'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'latex');
    
    % Capture the frame
    frame = getframe(gcf);
    writeVideo(v, frame);
end

% Close the video writer
close(v);
close(gcf); % Close the figure after video generation
fprintf('MP4 video saved as: %s\n', videoFile);


%% Potential Map - Individual Figures
% Select sample and camera to plot
sample = round(2.3*4000):100:round(2.6*4000);
Data_O = D_SYNC.CAM3;

% Initialize color map and adjust first color to white
C = jet(256);
C(1,1:3) = [1 1 1];

% Loop through each sample point
for i = 1:length(sample)
    % Extract current sample index
    current_sample = sample(i);
    
    % Extract a single frame from the processed data and rotate it 90 degrees
    I = squeeze(Data_O(:,:,current_sample));
    J = imrotate(I, 90);
    
    % Create a figure for each sample
    f1 = figure('color', 'white', 'Position', [50 50 500 500]);
    imagesc(J, [0 15]);
    colormap(C);
    
    % Add a colorbar and adjust labels
    hBar1 = colorbar('eastoutside');
    ylabel(hBar1, 'Fluorescence [%]', 'FontSize', 18);
    set(gca, 'fontsize', 18);
    title(['Cam - Sample ' num2str(current_sample)]);
    ylabel('Pixels');
    xlabel('Pixels');
    axis off
    
    % Optional: Add a pause to see each figure if plotting many
    % pause(0.5);
end


%% Potential Map - Multiplot
% Potential map - Subplot version
sample = round(2.3*4000):80:round(2.6*4000);
Data_O = D_SYNC.CAM3;

% Initialize color map
C = jet(256);
C(1,1:3) = [1 1 1];

% Calculate subplot dimensions - always use 8 columns
num_samples = length(sample);
cols = 8; % Num of coll
rows = ceil(num_samples / cols); % Calculate rows based on 8 columns

% Create single figure with subplots
f1 = figure('color', 'white');

for i = 1:num_samples
    current_sample = sample(i);
    
    I = squeeze(Data_O(:,:,current_sample));
    J = imrotate(I, 90);
    
    subplot(rows, cols, i);
    imagesc(J, [0 10]);
    colormap(C);
    title(['Sample ' num2str(current_sample)], 'FontSize', 8);
    axis off
end

% Add one colorbar for the entire figure
hBar1 = colorbar('eastoutside');
ylabel(hBar1, 'Fluorescence [%]', 'FontSize', 18);


%% Potential Map - Video
Data_O = D_SYNC.CAM1;
% Define the start and end samples for the video
start_sample = round(2.3*4000);  % Adjust the start sample according to your data
end_sample = round(2.6*4000);  % Adjust the end sample according to your data
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


%%
%%
%%
%%
%% Electric Potential Plot
clear; clc;


%% Loading Variables
load("E:\Qualification\Analysis\E32F02R01\data\data_filtered_sync_E32_F02_R01.mat"); % Load Sync data
load("E:\Qualification\Analysis\E32F02R01\data\InterpolatedSignalsE32_F02_R01_filtered.mat"); % Load the interpolated data


%% Overview electrical plot
% Loading variables
data_o1 = D_SYNC.CAM1;
data_o2 = D_SYNC.CAM2;
data_o3 = D_SYNC.CAM3;
data_e = D_SYNC.EL;
data_e1 = InterpSignal.Sync.MEA1;
data_e2 = InterpSignal.Sync.MEA2;
data_e3 = InterpSignal.Sync.MEA3;
data_e4 = InterpSignal.Sync.TANK;

% Plot time
timepoint_seconds = 2.46; 
% Plot window
plot_window_seconds = 1; 

% Select 6 electrodes
electrodes = [7 22 90 142 178 172]; 

% Selecting Colorbar limits
vmin_MEA1 = 0;
vmax_MEA1 = 5;
vmin_MEA2 = 0;
vmax_MEA2 = 5;
vmin_MEA3 = 0;
vmax_MEA3 = 5;
vmin_TANK = 0;
vmax_TANK = 5;

% Defining Fixed parameters
Fs = 4000;
sample = round(timepoint_seconds * Fs);
total_time = size(data_o1, 3) / Fs;

% --- CALCULATE PLOT WINDOW LIMITS ---
% Calculate the desired start and end times
t_start_desired = timepoint_seconds - (plot_window_seconds / 2);
t_end_desired = timepoint_seconds + (plot_window_seconds / 2);
% Apply boundary checks
t_start_plot = max(0, t_start_desired);
t_end_plot = min(total_time, t_end_desired);

% Selecting Optical points
Background_CAM1 = squeeze(imrotate(D_SYNC.IMG.CAM1(:,:,1), -90));
Background_CAM2 = squeeze(imrotate(D_SYNC.IMG.CAM2(:,:,1), -90));
Background_CAM3 = squeeze(imrotate(D_SYNC.IMG.CAM3(:,:,1), -90));
disp('Select 2 points for CAM1');
[x1, y1] = pick_up_a_trace(Background_CAM1, data_o1, 1);
disp('Select 2 points for CAM2');
[x2, y2] = pick_up_a_trace(Background_CAM2, data_o2, 1);
disp('Select 2 points for CAM3');
[x3, y3] = pick_up_a_trace(Background_CAM3, data_o3, 1);

% Define the subplot grid
rows = 6;
cols = 6;
% Create a figure
figure('color', 'white');
time_axis = (0:size(data_o1, 3)-1) / Fs;

% --- TIME PLOTS ---

% CAM1 Point 1
subplot(rows, cols, 1);
plot(time_axis, squeeze(data_o1(x1(1), y1(1), :))); 
hold on; 
xline(timepoint_seconds, 'r', 'LineWidth', 2);
xlim([t_start_plot t_end_plot]);
title(sprintf('CAM1 Pt 1\n(x=%d, y=%d)', x1(1), y1(1))); 
xlabel('Time (s)'); ylabel('%F'); 

% CAM1 Point 2
subplot(rows, cols, 2);
plot(time_axis, squeeze(data_o1(x1(2), y1(2), :))); 
hold on; 
xline(timepoint_seconds, 'r', 'LineWidth', 2);
xlim([t_start_plot t_end_plot]);
title(sprintf('CAM1 Pt 2\n(x=%d, y=%d)', x1(2), y1(2))); 
xlabel('Time (s)'); ylabel('%F'); 

% CAM2 Point 1
subplot(rows, cols, 3);
plot(time_axis, squeeze(data_o2(x2(1), y2(1), :))); 
hold on; 
xline(timepoint_seconds, 'r', 'LineWidth', 2);
xlim([t_start_plot t_end_plot]);
title(sprintf('CAM2 Pt 1\n(x=%d, y=%d)', x2(1), y2(1))); 
xlabel('Time (s)'); ylabel('%F'); 

% CAM2 Point 2
subplot(rows, cols, 4);
plot(time_axis, squeeze(data_o2(x2(2), y2(2), :))); 
hold on; 
xline(timepoint_seconds, 'r', 'LineWidth', 2);
xlim([t_start_plot t_end_plot]);
title(sprintf('CAM2 Pt 2\n(x=%d, y=%d)', x2(2), y2(2))); 
xlabel('Time (s)'); ylabel('%F'); 

% CAM3 Point 1
subplot(rows, cols, 5); 
plot(time_axis, squeeze(data_o3(x3(1), y3(1), :))); 
hold on; 
xline(timepoint_seconds, 'r', 'LineWidth', 2);
xlim([t_start_plot t_end_plot]);
title(sprintf('CAM3 Pt 1\n(x=%d, y=%d)', x3(1), y3(1))); 
xlabel('Time (s)'); ylabel('%F'); 

% CAM3 Point 2
subplot(rows, cols, 6); 
plot(time_axis, squeeze(data_o3(x3(2), y3(2), :))); 
hold on; 
xline(timepoint_seconds, 'r', 'LineWidth', 2);
xlim([t_start_plot t_end_plot]);
title(sprintf('CAM3 Pt 2\n(x=%d, y=%d)', x3(2), y3(2))); 
xlabel('Time (s)'); ylabel('%F'); 

% Electrical Plots
for k = 1:6 % Loop through all 6 electrodes
    subplot_index = k + 6; 
    subplot(rows, cols, subplot_index);
    plot(time_axis, data_e(electrodes(k), :), 'k');
    hold on; 
    xline(timepoint_seconds, 'r', 'LineWidth', 2);
    xlim([t_start_plot t_end_plot]);
    title(sprintf('Electrode %d', electrodes(k))); 
    xlabel('Time (s)'); ylabel('\muV'); 
end

% --- ELECTRICAL FRAME PLOTS ---

% MEA1 Frame
subplot(rows, cols, [13 14 19 20]); 
imagesc(data_e1(:, :, sample), [vmin_MEA1 vmax_MEA1]);
colormap(jet(256));
hold on;
title('MEA1');
colorbar;

% MEA2 Frame
subplot(rows, cols, [15 16 21 22]); 
imagesc(data_e2(:, :, sample), [vmin_MEA2 vmax_MEA2]);
colormap(jet(256));
hold on;
title('MEA2');
colorbar;

% MEA3 Frame
subplot(rows, cols, [17 18 23 24]);
imagesc(data_e3(:, :, sample), [vmin_MEA3 vmax_MEA3]);
colormap(jet(256));
hold on;
title('MEA3');
colorbar;

% TANK Frame
subplot(rows, cols, [26 27 28 29 32 33 34 35]);
imagesc(data_e4(:, :, sample), [vmin_TANK vmax_TANK]);
colormap(jet(256));
hold on;
title('TANK');
colorbar;

% Adjust overall figure appearance
sgtitle(['Synchronized Signals at t = ' num2str(timepoint_seconds) ' s '], ...
            'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'latex');


%% Overview Electrical Video
% Loading variables
data_o1 = D_SYNC.CAM1;
data_o2 = D_SYNC.CAM2;
data_o3 = D_SYNC.CAM3;
data_e = D_SYNC.EL;
data_e1 = InterpSignal.Sync.MEA1;
data_e2 = InterpSignal.Sync.MEA2;
data_e3 = InterpSignal.Sync.MEA3;
data_e4 = InterpSignal.Sync.TANK;
% Parâmetros Elétricos/Ópticos
electrodes = [7 22 90 142 178 172]; % 6 Eletrodos

% Parâmetros de Tempo do Vídeo
ti = 2.28; % Tempo inicial [s]
to = 2.62; % Tempo final [s]
step = 0.002; % [s]
frame_rate = 30; % Taxa de quadros (FPS)
% Parâmetro de Janela
plot_window_seconds = 0.5; 

% Limites da Colorbar para os Mapas Elétricos (do seu código)
vmin_MEA1 = 0; 
vmax_MEA1 = 5; 
vmin_MEA2 = 0; 
vmax_MEA2 = 5;
vmin_MEA3 = 0; 
vmax_MEA3 = 5;
vmin_TANK = 0; 
vmax_TANK = 5;

% Saída de Vídeo
videoFile = 'E_overview_video.mp4';

% Parâmetros Fixos
Fs = 4000;
time_axis = (0:size(data_o1, 3)-1) / Fs;
total_time = size(data_o1, 3) / Fs;
rows = 6;
cols = 6;

% --- Calcular e Confirmar Duração ---
total_frames = length(ti:step:to);
predicted_video_duration = total_frames / frame_rate;
disp(['Total de Quadros a processar: ', num2str(total_frames)]);
disp(['Duração Prevista do Vídeo a ' num2str(frame_rate) ' FPS: ', num2str(predicted_video_duration), ' s']);
% Prompt user to confirm
user_input = input('Deseja prosseguir com estes valores? (y/n): ', 's');
if lower(user_input) ~= 's' && lower(user_input) ~= 'y' 
    fprintf('Por favor, modifique os valores de ti, to, step e frame_rate no código e execute a seção novamente.\n');
    return;
end

% --- Selecionar Pontos Ópticos e Inicializar Figura ---
Background_CAM1 = squeeze(imrotate(D_SYNC.IMG.CAM1(:,:,1), -90));
Background_CAM2 = squeeze(imrotate(D_SYNC.IMG.CAM2(:,:,1), -90));
Background_CAM3 = squeeze(imrotate(D_SYNC.IMG.CAM3(:,:,1), -90));
disp('Selecione 2 pontos para CAM1');
[x1, y1] = pick_up_a_trace(Background_CAM1, data_o1, 1);
disp('Selecione 2 pontos para CAM2');
[x2, y2] = pick_up_a_trace(Background_CAM2, data_o2, 1);
disp('Selecione 2 pontos para CAM3');
[x3, y3] = pick_up_a_trace(Background_CAM3, data_o3, 1);

figure('color', 'white', 'Name', 'Optical/Electric Maps Video');
set(gcf, 'Position', get(0, 'Screensize')); % Maximizar figura

% --- Inicializar Gravador de Vídeo (MP4) ---
v = VideoWriter(videoFile, 'MPEG-4'); 
v.FrameRate = frame_rate; 
open(v);

% --- Loop Principal de Geração de Vídeo ---
for timepoint_seconds = ti:step:to
    sample = round(timepoint_seconds * Fs);
    
    % --- CALCULAR LIMITES DA JANELA DE PLOTAGEM ---
    t_start_desired = timepoint_seconds - (plot_window_seconds / 2);
    t_end_desired = timepoint_seconds + (plot_window_seconds / 2);
    
    % Aplicar verificação de limites
    t_start_plot = max(0, t_start_desired);
    t_end_plot = min(total_time, t_end_desired);
    

    % CAM1 Point 1
    subplot(rows, cols, 1);
    plot(time_axis, squeeze(data_o1(x1(1), y1(1), :)));
    hold on; xline(timepoint_seconds, 'r', 'LineWidth', 2); hold off;
    xlim([t_start_plot t_end_plot]); 
    title(sprintf('CAM1 Pt 1\n(x=%d, y=%d)', x1(1), y1(1))); 
    xlabel('Time (s)'); ylabel('%F');
    
    % CAM1 Point 2
    subplot(rows, cols, 2);
    plot(time_axis, squeeze(data_o1(x1(2), y1(2), :)));
    hold on; xline(timepoint_seconds, 'r', 'LineWidth', 2); hold off;
    xlim([t_start_plot t_end_plot]); 
    title(sprintf('CAM1 Pt 2\n(x=%d, y=%d)', x1(2), y1(2))); 
    xlabel('Time (s)'); ylabel('%F');
    
    % CAM2 Point 1
    subplot(rows, cols, 3);
    plot(time_axis, squeeze(data_o2(x2(1), y2(1), :)));
    hold on; xline(timepoint_seconds, 'r', 'LineWidth', 2); hold off;
    xlim([t_start_plot t_end_plot]); 
    title(sprintf('CAM2 Pt 1\n(x=%d, y=%d)', x2(1), y2(1))); 
    xlabel('Time (s)'); ylabel('%F');
    
    % CAM2 Point 2
    subplot(rows, cols, 4);
    plot(time_axis, squeeze(data_o2(x2(2), y2(2), :)));
    hold on; xline(timepoint_seconds, 'r', 'LineWidth', 2); hold off;
    xlim([t_start_plot t_end_plot]); 
    title(sprintf('CAM2 Pt 2\n(x=%d, y=%d)', x2(2), y2(2))); 
    xlabel('Time (s)'); ylabel('%F');
    
    % CAM3 Point 1
    subplot(rows, cols, 5);
    plot(time_axis, squeeze(data_o3(x3(1), y3(1), :)));
    hold on; xline(timepoint_seconds, 'r', 'LineWidth', 2); hold off;
    xlim([t_start_plot t_end_plot]); 
    title(sprintf('CAM3 Pt 1\n(x=%d, y=%d)', x3(1), y3(1))); 
    xlabel('Time (s)'); ylabel('%F');
    
    % CAM3 Point 2
    subplot(rows, cols, 6);
    plot(time_axis, squeeze(data_o3(x3(2), y3(2), :)));
    hold on; xline(timepoint_seconds, 'r', 'LineWidth', 2); hold off;
    xlim([t_start_plot t_end_plot]); 
    title(sprintf('CAM3 Pt 2\n(x=%d, y=%d)', x3(2), y3(2))); 
    xlabel('Time (s)'); ylabel('%F');
    
    % Electric Electrodes
    for k = 1:6 
        subplot_index = k + 6; 
        subplot(rows, cols, subplot_index);
        plot(time_axis, data_e(electrodes(k), :), 'k');
        hold on; xline(timepoint_seconds, 'r', 'LineWidth', 2); hold off;
        xlim([t_start_plot t_end_plot]); 
        title(sprintf('Electrode %d', electrodes(k))); 
        xlabel('Time (s)'); ylabel('\muV'); 
    end
    
    
    % MEA1 Frame
    subplot(rows, cols, [13 14 19 20]); 
    I_e1 = squeeze(data_e1(:, :, sample)); 
    imagesc(I_e1, [vmin_MEA1 vmax_MEA1]);
    colormap(jet(256));
    title(sprintf('MEA1'));
    colorbar;
    axis off; 
    
    % MEA2 Frame
    subplot(rows, cols, [15 16 21 22]); 
    I_e2 = squeeze(data_e2(:, :, sample)); 
    imagesc(I_e2, [vmin_MEA2 vmax_MEA2]);
    colormap(jet(256));
    title(sprintf('MEA2'));
    colorbar;
    axis off; 
    
    % MEA3 Frame
    subplot(rows, cols, [17 18 23 24]);
    I_e3 = squeeze(data_e3(:, :, sample)); 
    imagesc(I_e3, [vmin_MEA3 vmax_MEA3]);
    colormap(jet(256));
    title(sprintf('MEA3'));
    colorbar;
    axis off;
    
    % TANK Frame
    subplot(rows, cols, [26 27 28 29 32 33 34 35]);
    I_e4 = squeeze(data_e4(:, :, sample)); 
    imagesc(I_e4, [vmin_TANK vmax_TANK]);
    colormap(jet(256));
    title(sprintf('TANK'));
    colorbar;
    axis off;
    
    sgtitle(['Overview Electric Video | t=' num2str(timepoint_seconds, '%.3f') ' s'], ...
                'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'latex');
    
    % Capturar o frame
    frame = getframe(gcf);
    writeVideo(v, frame);
end

% Fechar o gravador de vídeo
close(v);
close(gcf); 
fprintf('Vídeo MP4 salvo como: %s\n', videoFile);


%% Potential Map - Indivudual Figures
% Select parameters to plot
sample = round(2.3*4000):40:round(2.6*4000);
Data_E = InterpSignal.Sync.TANK;
lim = [0 5];
Title = 'Potential Plot';

% Loop through each sample point
for i = 1:length(sample)
    % Extract current sample index
    current_sample = sample(i);
    
    % Plot
    plot_potential_map(Data_E, current_sample, ...
        'lim', lim, 'title', Title);
end


%% Potential Map - Multiplot
% Select parameters to plot
sample = round(2.3*4000):26:round(2.4*4000);
Data_E = InterpSignal.Sync.TANK;
lim = [0 5];
Title = 'Potential Plot - TANK';

% Plot
plot_potential_map(Data_E, sample, ...
    'lim', lim, 'title', Title);


%% Potential Map - Video
Data_E = InterpSignal.Data.TANK;
% Define plot parameters
lim = [0 5]; % Set to [] for auto-scaling or specific limits like [-100 100]
Title = 'Potential Map'; % Custom title for the video frames
video_title = 'E_potential_map_video_';
% Define the start and end samples for the video
start_sample = round(2.3*4000);
end_sample = round(2.4*4000);
Fsampling = 4000;

% Define the frames per second (fps) for the video
fps = 30;  % 30 fps is a standard value for video codec
% Define the step for plot
step = 1; 

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
for current_sample = start_sample:step:end_sample
    % Create figure using our plot_potential_map function
    f1 = figure('color', 'white', 'Position', [50 50 500 400], 'Visible', 'off');
    
    % Extract single frame
    I = squeeze(Data_E(:,:,current_sample));
    
    % Create 2D surface plot
    I = flipud(I); % Vertical Flip because surf() invert the y-axis
    surf(I, 'EdgeColor', 'none', 'FaceColor', 'interp');
    view(2); % 2D view from above
    axis equal;
    axis tight;
    
    % Set color limits if specified
    if ~isempty(lim)
        caxis(lim);
    end
    
    % Set title with sample number
    title([Title, ' | Sample: ', num2str(current_sample)], 'FontSize', 14, 'Interpreter', 'none');
    
    % Add colorbar
    hBar1 = colorbar('Location', 'eastoutside', 'Position', [0.92 0.1 0.02 0.8]);
    colormap(jet(256));
    ylabel(hBar1, 'Potential [$\mu$V]', 'FontSize', 12, 'Interpreter', 'latex');
    
    % Remove axis ticks
    set(gca, 'XTick', [], 'YTick', []);
    set(gca,'fontsize', 12);
    
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



