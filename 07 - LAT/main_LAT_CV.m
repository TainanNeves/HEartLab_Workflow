%% LOCAL ACTIVATION TIME and CONDUCTION VELOCITY
clear; clc;


%% Loading Variables

load("C:\Users\HEartLab\Documents\GitHub\HEartLab\00 - examples\data_filtered_sync_E14_F3_R4.mat"); %Filtered data


%% LAT Optic Analysis

% Parameters
Data = D_SYNC.CAM2;
Fs = 4000;

% Preview Signal
Background = squeeze(Data(:,:,2000));
pick_up_a_trace(Background, Data,1);
% Selecting interval
lim1 = 7710; % Start sample index
lim2 = 9220; % End sample index

% Extract a segment of optical data
Data_temp_O = Data(:,:,lim1:lim2);
% Initialize a matrix for Local Activation Times (LAT) for each pixel
LAT_O = Data_temp_O(:,:,1) * 0;
% Calculate LAT for each pixel using find_LAT_linearFit_1D function
for i = 1:size(Data_temp_O,1)
    for j = 1:size(Data_temp_O,2)
        if max(max(squeeze(Data_temp_O(i,j,:)))) ~= 0
            % (y = 1D array, fr = Frame Rate, L = length of the linear fit line, PCL = 200 (not Used), debug = 1 or 0 (Plot or not point and trace))
            [LAT_O(i,j)] = find_LAT_linearFit_1D(squeeze(Data_temp_O(i,j,:)), Fs, 15, 200, 0); 
        end
    end
end

% Plot the LAT data on a colormap
f1 = figure('color', 'white', 'Position', [40 40 600 600]);
C = parula(256); %hsv é um colormap de arcoiris, mas não termina em purple
C(1,1:3) = [1 1 1]; % Make 0 = white
C(256,1:3) = [0.5, 0.5, 0.5];
J = LAT_O(:);
% Set zeros inside ROI
% J(D_OP.ROI.ROI_2 & J == 0) = max(max(LAT_O))+10;
J = imrotate(LAT_O,90);
Y = prctile(nonzeros(J),[1 95],'all');
imagesc(J-Y(1), [-1, 100]);
colormap(C);
hBar1 = colorbar('eastoutside');
ylabel(hBar1, 'Local Activation Time [ms]', 'FontSize', 14);
box off
set(gca, 'fontsize', 18);
ylabel('Pixels');
xlabel('Pixels');
axis off;
title('Local Activation Time | Cam ');


%% CV Optic Analysis

% Apply spatial and temporal filtering to LAT data
LAT_O2 = SpatTemp_Filtering(LAT_O, 3, 0, 'GPU');
% Initialize matrices for Average Speed, Standard Deviation, and Angle
AvgSpeed_O = zeros(size(Data,1), size(Data,2));
StdSpeed_O = zeros(size(Data,1), size(Data,2));
Angle_O = zeros(size(Data,1), size(Data,2));
% Calculate conduction velocity using CV_CircleMethod
for i = 1:size(Data,1)
    for j = 1:size(Data,2)
        if max(max(squeeze(Data(i,j,:)))) ~= 0
            % (LAT matrix, radious to calc, row, col, SapaceScale in millimeters per pixel)
            [AvgSpeed_O(i,j),  StdSpeed_O(i,j), Angle_O(i,j)] = CV_CircleMethod(LAT_O2, 5, i, j, 0.33); 
        end
    end
end

% Plot the conduction velocity data on a colormap
f1 = figure('color', 'white', 'Position', [40 40 600 600]);
C = jet(256);
C(1,1:3) = [1 1 1];
J = AvgSpeed_O(:);
Y = prctile(nonzeros(J),[5 95],'all');
J = imrotate(AvgSpeed_O,90);
teto = Y(2) - Y(1);
imagesc(J - Y(1), [0, Y(2) - Y(1)]);
colormap(C);
hBar1 = colorbar('eastoutside');
ylabel(hBar1, 'Conduction velocity [cm/s]', 'FontSize', 18);
hold on;
d = 200; 
color = 'black'; 
tam = 14;
set(gca, 'fontsize', 18);
ylabel('Pixels');
xlabel('Pixels');
axis off
title('Conduction Velocity | Cam ');


%% LAT Electric Analysis

% Parameters
Data = D_SYNC.EL;
Fs = 4000;
lim1 = 7710; % Start sample index
lim2 = 9220; % End sample index

% Calculate LAT for electrodes using find_LAT_diff_1D function
LAT_E = find_LAT_diff(Data, Fs, lim1, lim2, 0);

% Ploting maps
LAT_E_matrix1 = plot_electric_LAT(LAT_E, [0 30], 1, 1); % MEA 1
LAT_E_matrix2 = plot_electric_LAT(LAT_E, [0 30], 1, 2); % MEA 2
LAT_E_matrix3 = plot_electric_LAT(LAT_E, [0 30], 1, 3); % MEA 3
LAT_E_matrix4 = plot_electric_LAT(LAT_E, [0 30], 1, 4); % TANK


%% CV Electric Analysis

% Select Matrix
LAT_E_matrix = LAT_E_matrix4;

% Create a meshgrid of coordinates for activation times
[X, Y] = meshgrid(1:size(LAT_E_matrix, 2), 1:size(LAT_E_matrix, 1));
% Create a meshgrid of coordinates for interpolation
[Xq, Yq] = meshgrid(1:0.1:size(LAT_E_matrix, 2), 1:0.1:size(LAT_E_matrix, 1));
% Perform bilinear interpolation
LAT_E_matrix_interp = interp2(X, Y, LAT_E_matrix, Xq, Yq, 'linear');
% Perform CV cicle method for CV calculation
r = 3;
for i = (r + 1):(size(LAT_E_matrix_interp, 1) - (r + 1))
    for j = (r + 1):(size(LAT_E_matrix_interp, 2) - (r + 1))
        if max(max(squeeze(LAT_E_matrix_interp(i, j, :)))) ~= 0
            % (LAT matrix, radious to calc, row, col, SapaceScale in millimeters per pixel (TANK = 0.75 | MEA = 0.1)
            [AvgSpeed_E(i, j), StdSpeed_E(i, j), Angle_E(i, j)] = CV_CircleMethod(LAT_E_matrix_interp, r, i, j, 0.1);
        end
    end
end

% Spatially filter the average speed
AvgSpeed_E2 = SpatTemp_Filtering(AvgSpeed_E, 5, 0, 'GPU');

% Plot the results using jet colormap
f1 = figure('color', 'white', 'Position', [40 40 230 210]);
C = jet(256);
J = AvgSpeed_E2(:);
Y = prctile(nonzeros(J), [5 95], 'all');
J = imrotate(AvgSpeed_E2, 90);
imagesc(J - Y(1), [0, Y(2) - Y(1)]); colormap(C);
hBar1 = colorbar('eastoutside'); ylabel(hBar1, 'Conduction velocity [cm/s]', 'FontSize', 18);
set(gca, 'fontsize', 18);
title('Cam 2');
ylabel('Pixels'); xlabel('Pixels');
axis off




















