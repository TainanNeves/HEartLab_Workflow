%% PHASE MAPS

clear; clc;


%% Loading variables

load("E:\HEartLab\TAINAN WORKFLOW\00 - examples\data_filtered_sync_E14_F3_R4.mat"); %Filtered data


%% Optic Phase Map

% Loading variables
Data = D_SYNC.CAM1;
Fsampling = 4000;

% Loading colormap
load("PS_colormaps.mat");
mycmap_2 = mycmap;
mycmap_2(1, 1:3) = [1 1 1];

% Extract a subset of the optical data for analysis (adjust the sample range accordingly)
temp_o = Data(:,:,1:1000); % Usually XX samples?????????????

% Calculate phase map for selected window
[phase_op] = phasemap_opticos(temp_o,1);

% Ploting
% Click space to see next figure. See all and select one.
for i = 1:1000
    f1 = figure('color', 'white', 'Position', [40 40 450 350]);
    I = squeeze(phase_op(:, :, i));
    J = imrotate(I, 90);
    imagesc(J, [-3 3]);
    colormap(mycmap);
    colorbar;
    box off;
    set(gca, 'fontsize', 18);
    ylabel('Pixels');
    xlabel('Pixels');
    axis off;
    pause;
end
% close all (To chlose figures)

% Saving Variable
CAM_1_PHASE = J;

% Visualizing unic frame
i = 30; % Put yor selected frame (do not use the first 300 and last 300)
f1 = figure('color', 'white', 'Position', [40 40 450 350]);
I = squeeze(phase_op(:, :, i));
J = imrotate(I, 90);
imagesc(J, [-3 3]);
colormap(mycmap);
colorbar;
box off;
set(gca, 'fontsize', 18);
ylabel('Pixels');
xlabel('Pixels');
titulo = ['Frame: ', num2str(i)];
title(titulo);
axis off;


%% Electric Phase Map

% Loading variables
Data = D_SYNC.EL;
Fsampling = 4000;

% Selecting sample interval
start_sample = 1;
end_sample = 1000;
% Loading temporary variables
temp_e = Data(:, start_sample:end_sample);

% Loading colormap
load("PS_colormaps.mat");
mycmap_2 = mycmap;
mycmap_2(1, 1:3) = [1 1 1];

% Calculating Phase
[phase_el] = phasemap_electrico(temp_e);

% imput variables
frame=30;

% Ploting Maps
plot_electric_phasemap(phase_el, [-3 3], frame, 1); % MEA 1
plot_electric_phasemap(phase_el, [-3 3], frame, 2); % MEA 2
plot_electric_phasemap(phase_el, [-3 3], frame, 3); % MEA 3
plot_electric_phasemap(phase_el, [-3 3], frame, 4); % TANK

