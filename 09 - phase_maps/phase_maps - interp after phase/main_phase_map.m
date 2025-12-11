%% PHASE MAPS

clear; clc;


%% Loading variables

load("C:\Users\HEartLab\Documents\GitHub\HEartLab\00 - examples\data_filtered_sync_E14_F3_R4.mat"); %Filtered data


%% Optic Phase Map

% Loading variables
Data = D_SYNC.CAM3;
Fsampling = 4000;
% Extract a subset of the optical data for analysis (adjust the sample range accordingly)
start_sample = 5092;
end_sample = 6092;
temp_o = Data(:,:,start_sample:end_sample);

% Calculate phase map for selected window
[S, phase_o] = phasecalc_opticos(temp_o);

% Selecting Point
Background = squeeze(temp_o(:,:,500));
pick_up_a_trace(Background, temp_o,1);
point = [41 121];

% Signal Plot
frame = 500;
figure;
sgtitle("CAM 03: V | Ventricular Signal");
subplot(3,1,1);
plot(squeeze(temp_o(point(1),point(2),:)));
title("Original Signal");
subplot(3,1,2);
plot(squeeze(phase_o(point(1),point(2),:)));
title("Signal Phase");
subplot(3,1,3);
plot(squeeze(S(point(1),point(2),:)));
title("Filtered Signal (Internally in the function)");
for i = [1:3]
    subplot(3,1,i);
    xline(frame, '-r', {num2str(frame)});
end


% Optical Map
frame = 500;
ROI = D_SYNC.ROI.ROI_3;
% Load colormap
load("PS_colormaps.mat");
% Create a figure for the phase map
f1 = figure('Color', 'white', 'Position', [500 500 400 500]);
I = squeeze(phase_o(:, :, frame));
J = imrotate(I, 90);
ROI = imrotate(ROI, 90);
% Color correction
J2 = J;
J2(1,1) = +pi; 
J2(1,2) = -pi;
J2 = mat2gray(J2); 
J2 = gray2ind(J2);
J2 = ind2rgb(J2,colormap((mycmap))); % flipud
% Make background white
for i=1:size(J2,1)
     for j=1:size(J2,2)
         if ROI(i,j) == 0
            J2(i,j,1:3) = [1,1,1];
         end
     end
 end
imagesc(J2,[-pi pi]); % Display image with a specified data range
h1 = colorbar('southoutside');
title(['Phase Map | Frame: ', num2str(frame), ' | ']);
axis equal;
axis off;
ylabel(h1, ['Radians [' char(960) ']'],'FontSize',14);
set(gca,'fontsize', 14);


%% Electric Phase Map

% Loading variables
Data = D_SYNC.EL;
Fsampling = 4000;

% Selecting sample interval
start_sample = 5092;
end_sample = 6092;
% Loading temporary variables
temp_el = Data(:, start_sample:end_sample);


% Calculating Phase
[S, phase_el] = phasecalc_electrico(temp_el);

% Signal Plot
% Selecting Electrode
el = 188;
frame = 500;
% Plot
figure;
sgtitle(['electrode: ', num2str(el), ' | TANK']);
subplot(3,1,1);
plot(temp_el(el,:));
title("Original Signal");
subplot(3,1,2);
plot(phase_el(el,:));
title("Signal Phase");
subplot(3,1,3);
plot(S(el,:));
title("Filtered Signal (Internally in the function)");
for i = [1:3]
    subplot(3,1,i);
    xline(frame, '-r', {'*'});
end


% Ploting Maps
frame = 500;
phase_e_matrix1 = plot_electric_phasemap(phase_el, [-3 3], frame, 1); % MEA 1
phase_e_matrix2 = plot_electric_phasemap(phase_el, [-3 3], frame, 2); % MEA 2
phase_e_matrix3 = plot_electric_phasemap(phase_el, [-3 3], frame, 3); % MEA 3
phase_e_matrix4 = plot_electric_phasemap(phase_el, [-3 3], frame, 4); % TANK


%% Colorbar

% HORIZONTAL
% Create a new figure for the colorbar
f1 = figure('color', 'white', 'Position', [800, 150, 800, 300]);
% Set the title for the plot
title('Phase Colorbar', 'FontSize', 16);
% Adjust axis properties
axis equal; % Ensures a uniform scaling
axis off;   % Hides the axes
caxis([-3, 3]); % Set the color scale range
% Load the colormap and add the colorbar
load("PS_colormaps.mat"); % Ensure your custom colormap is loaded
colormap(mycmap);         % Apply the colormap
% Add the colorbar and set its position
hBar1 = colorbar('southoutside'); % Start with default position
ylabel(hBar1, ['Radians [' char(960) ']'], 'FontSize', 14); % Add a label
% Numbers
hBar1.Ticks = [-2, 0, 2]; % Define the tick marks
hBar1.TickLabels = {'-2', '0', '2'}; % Set the tick labels
% States
hBar1.Ticks = [-2, 0, 2, 3]; % Define the tick marks
hBar1.TickLabels = {'Repolarization', 'Plateau', 'Depolarization','Resting'}; % Set the tick labels
% Adjust the colorbar's position within the figure
set(hBar1, 'Position', [0.05, 0.4, 0.9, 0.2]); % Adjust the position: [left, bottom, width, height]
% Set font size for the axis and colorbar
set(gca, 'FontSize', 14); % Font size for the plot


% VERTICAL
% Create a new figure for the colorbar
f1 = figure('color', 'white', 'Position', [800, 150, 400, 800]); % Taller figure for vertical colorbar
% Set the title for the plot
title('Phase Colorbar', 'FontSize', 16);
% Adjust axis properties
axis equal; % Ensures uniform scaling
axis off;   % Hides the axes
caxis([-3, 3]); % Set the color scale range
% Load the colormap and add the colorbar
load("PS_colormaps.mat"); % Ensure your custom colormap is loaded
colormap(mycmap);         % Apply the colormap
% Add the colorbar with a vertical orientation ('eastoutside' for vertical on the right)
hBar1 = colorbar('eastoutside'); % Place colorbar on the right side of the figure
% Set the label for the colorbar
ylabel(hBar1, ['Radians [' char(960) ']'], 'FontSize', 14, 'Rotation', 270); % Rotate label to align with vertical orientation
% Numbers
hBar1.Ticks = [-2, 0, 2]; % Define the tick marks
hBar1.TickLabels = {'-2', '0', '2'}; % Set the tick labels
% States
hBar1.Ticks = [-2, 0, 2, 3]; % Define the tick marks
hBar1.TickLabels = {'Repolarization', 'Plateau', 'Depolarization', 'Resting'}; % Custom tick labels
% Adjust the position and size of the vertical colorbar
set(hBar1, 'Position', [0.5, 0.1, 0.15, 0.8]); % Adjust the position: [left, bottom, width, height]
% Set font size for the axis and colorbar
set(gca, 'FontSize', 14); % Font size for the plot





