%% Main - Potential Plot
% uses the interpolated signal to plot
clear; clc;


%% Loading Variables
load(""); % Load the interpolated data


%% Potential Map - Indivudual Figures
% Select parameters to plot
sample = round(2.23*4000):20:round(2.62*4000);
Data_E = InterpSignal.Sync.TANK;
lim = [-400 400];
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
sample = round(2.23*4000):20:round(2.62*4000);
Data_E = InterpSignal.Sync.MEA3;
lim = [-3000 3000];
Title = 'Potential Plot - MEA3';

% Plot
plot_potential_map(Data_E, sample, ...
    'lim', lim, 'title', Title);


%% Potential Map - Video
Data_E = InterpSignal.Data.MEA1;
% Define plot parameters
lim = [-3000 3000]; % Set to [] for auto-scaling
Title = 'Potential Map'; % Custom title for the video frames
video_title = 'E_potential_map_video_';
% Define the start and end samples for the video
start_sample = round(2.23*4000);
end_sample = round(2.62*4000);
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
    
    % Extract single frame and transpose for correct orientation
    I = squeeze(Data_E(:,:,current_sample))';
    
    % Create 2D surface plot
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