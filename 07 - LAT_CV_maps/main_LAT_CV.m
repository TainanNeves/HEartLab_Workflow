%% LOCAL ACTIVATION TIME and CONDUCTION VELOCITY
clear; clc;


%% Loading Variables

load("F:\HEartLab\Activities\CCC06 - CBEB\analises 02\signals\data_filtered_sync_E14_F03_R09.mat"); %Filtered data


%% LAT Optic Analysis

% Parameters
Data = D_SYNC.CAM1;
Fs = 4000;

% Preview Signal
Background = squeeze(Data(:,:,2000));
pick_up_a_trace(Background, Data,1);
% Selecting interval
lim1 = 8863; % Start sample index
lim2 = 9508; % End sample index

% Visualize in the interval
Background = squeeze(Data(:,:,2000));
pick_up_a_trace(Background, Data(:,:,lim1:lim2),1);

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
contourf(flipud(J));
% imagesc(J-Y(1), [-1 20]);
colormap(C);
hBar1 = colorbar('eastoutside');
ylabel(hBar1, 'Local Activation Time [ms]', 'FontSize', 14);
box off
set(gca, 'fontsize', 18);
ylabel('Pixels');
xlabel('Pixels');
axis off;
title('Local Activation Time | Cam 1');


%% subtractin Minimum
LAT_O(LAT_O < 0) = 0;
min_LAT = min(LAT_O(LAT_O ~= 0));



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
% Data = signal_file.D_SYNC.EL;
Fs = 4000;
lim1 = start_index; % Start sample index
lim2 = 18000; % End sample index

% Calculate LAT for electrodes using find_LAT_diff_1D function
LAT_E = find_LAT_diff(Data, Fs, lim1, lim2, 0);

% % Ploting maps
% LAT_E_matrix1 = plot_electric_LAT(LAT_E, [0 30], 1, 1); % MEA 1
% LAT_E_matrix2 = plot_electric_LAT(LAT_E, [0 30], 1, 2); % MEA 2
% LAT_E_matrix3 = plot_electric_LAT(LAT_E, [0 30], 1, 3); % MEA 3
% LAT_E_matrix4 = plot_electric_LAT(LAT_E, [0 40], 1, 4); % TANK

%% Minimum LAT values for MEAs

% MEA1
mea1_id = [1:11, 14:16]; % Define indices for MEA1 electrodes
min_mea1 = min(LAT_E(mea1_id)); % Find the minimum LAT value for MEA1
el_min_mea1 = mea1_id(LAT_E(mea1_id) == min_mea1); % Find the electrode with the minimum LAT value for MEA1

% MEA2
mea2_id = [17:32]; % Define indices for MEA2 electrodes
min_mea2 = min(LAT_E(mea2_id)); % Find the minimum LAT value for MEA2
el_min_mea2 = mea2_id(LAT_E(mea2_id) == min_mea2); % Find the electrode with the minimum LAT value for MEA2

% MEA3
mea3_id = [65:79]; % Define indices for MEA3 electrodes
min_mea3 = min(LAT_E(mea3_id)); % Find the minimum LAT value for MEA3
el_min_mea3 = mea3_id(LAT_E(mea3_id) == min_mea3); % Find the electrode with the minimum LAT value for MEA3

% Correcting indices for MEA2 and MEA3
el_min_mea2 = el_min_mea2; % correcting the electrode number
el_min_mea3 = el_min_mea3 + 64; % correcting the electrode number


%% LAT from ECGi

% Load the estimated potentials
% estimated_signal = load('estimated_signal_E18F2R2_sync.mat');
% estimated_signal = estimated_signal.(subsref(fieldnames(estimated_signal),substruct('{}',{1})));
% 
% % Load the positions of the MEAs
% files_id_meas = load('C:\Users\HeartLAB\Documents\Documents\CinC 2024\ECGi\Dados\projected_signals_exp14.mat');

% Combine the indices of MEA1, MEA2, and MEA3 to obtain all vertices
% vertices = cat(2, [file_id_meas.MEAS_HR_IDX.MEA1], [file_id_meas.MEAS_HR_IDX.MEA2], [file_id_meas.MEAS_HR_IDX.MEA3]);

% Extract the estimated potentials corresponding to the vertices
% estimated_meas = estimated_signal(vertices,:);

Fs = 4000; % Sampling frequency (Hz)
start_time = 0.270; % Start time in seconds
end_time = 0.289;
% Calculate the start and end sample indices
start_index = start_time * Fs + 1; % Start sample index
end_index = end_time * Fs; % End sample index

% Calculate LAT for electrodes using the find_LAT_diff function
LAT_ECGI = find_LAT_diff(estimated_signal, Fs, start_index, end_index, 0);

%% Plot LAT ECGi map

% heart_geo = estimated_file.estimated.Heart_geometry;
LAT_ECGI = LAT_ECGI';
% Sampling frequency and time details (adjust based on LAT if necessary)
n_frames = size(LAT_ECGI, 2);  % Adjust as necessary based on LAT dimensions

% Organizing geometry
faces = heart_geo.faces;
x = heart_geo.vertices(:, 1);
y = heart_geo.vertices(:, 2);
z = heart_geo.vertices(:, 3);

% Create a figure for the LAT map
figure();
set(gcf, 'Color', 'w');  % Set background color to white

% Prepare the trisurf plot for LAT data
trisurf_plot = trisurf(faces, x, y, z, LAT_ECGI(:, 1), 'FaceColor', 'interp', 'EdgeColor', 'none');
title('Local Activation Time Map');
c = colorbar('Location', 'eastoutside');  % Position the colorbar on the right
c.Label.String = 'LAT (ms)';
c.Label.FontWeight = 'bold';
c.Label.FontSize = 14;
c.TickLabelInterpreter = 'tex';
set(c, 'FontWeight', 'bold', 'FontSize', 14);

colormap('parula'); 

% Adjust color axis settings based on LAT data range
caxis([0, max(LAT_ECGI(:))]);
grid off; axis off;

% Optional: Adjust view angle, lighting, and material properties for better visualization
% view(3); % Set view for a 3D perspective
% lighting gouraud;  % Smooth lighting
% camlight('headlight');  % Light from the camera direction
% material dull;  % Set surface material properties for LAT visualization



%% Mininum ECGi LAT values in MEAs positions

% Correspondence for MEA1
min_mea1_ecgi = min(LAT_ECGI([1:11,14:16])); % Find the minimum LAT value in MEA1
el_min_mea1_ecgi = find(LAT_ECGI([1:11,14:16]) == min_mea1_ecgi); % Find the electrode index with the minimum value for MEA1

% Correspondence for MEA2
min_mea2_ecgi = min(LAT_ECGI([17:32])); % Find the minimum LAT value in MEA2
el_min_mea2_ecgi = find(LAT_ECGI([17:32]) == min_mea2_ecgi); % Find the electrode index with the minimum value for MEA2

% Correspondence for MEA3
min_mea3_ecgi = min(LAT_ECGI([33:47])); % Find the minimum LAT value in MEA3
el_min_mea3_ecgi = find(LAT_ECGI([33:47]) == min_mea3_ecgi); % Find the electrode index with the minimum value for MEA3

% Correcting indices for MEA2 and MEA3
el_min_mea2_ecgi = el_min_mea2_ecgi + 16; % correcting the electrode number
el_min_mea3_ecgi = el_min_mea3_ecgi + 64; % correcting the electrode number

%% Distance between electrodes
% calculate the geodesic distance between electrodes in the 3D heart
% geometry
%inputs: electrodes, projections_file (with all geometries), resolution

el1 = 1;
el2 = 2;
geodesicDistance(el1,el2,file_id_meas,'HR');


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




















