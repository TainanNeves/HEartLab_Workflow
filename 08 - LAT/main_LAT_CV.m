%% LOCAL ACTIVATION TIME and CONDUCTION VELOCITY
clear; clc;


%% Loading Variables

load("C:\Users\HEartLab\Documents\GitHub\HEartLab\00 - examples\data_filtered_sync_E14_F3_R4.mat"); %Filtered data


%% LAT Optic Analysis

% Parameters
Data = D_SYNC.CAM2;
Fs = 4000;
lim1 = 240; % Start sample index
lim2 = 400; % End sample index

% Extract a segment of optical data
Data_temp_O = Data(:,:,lim1:lim2);
% Initialize a matrix for Local Activation Times (LAT) for each pixel
T_LAT_ms = Data_temp_O(:,:,1) * 0;
% Calculate LAT for each pixel using find_LAT_linearFit_1D function
for i = 1:size(Data_temp_O,1)
    for j = 1:size(Data_temp_O,2)
        if max(max(squeeze(Data_temp_O(i,j,:)))) ~= 0
            % (y = 1D array, fr = Frame Rate, L = length of the linear fit line, PCL = 200 (not Used), debug = 1 or 0 (Plot or not point and trace))
            [T_LAT_ms(i,j)] = find_LAT_linearFit_1D(squeeze(Data_temp_O(i,j,:)), Fs, 15, 200, 0); 
        end
    end
end

% Plot the LAT data on a colormap
f1 = figure('color', 'white', 'Position', [40 40 400 300]);
C = parula(256); %hsv é um colormap de arcoiris, mas não termina em purple
C(1,1:3) = [1 1 1]; % Make 0 = white
J = T_LAT_ms(:);
J = imrotate(T_LAT_ms,90);
Y = prctile(nonzeros(J),[1 95],'all');
imagesc(J-Y(1), [-1, 100]);colormap(C);
hBar1 = colorbar('eastoutside');
ylabel(hBar1, 'Local Activation Time [ms]', 'FontSize', 14);
box off
set(gca, 'fontsize', 18);
ylabel('Pixels');
xlabel('Pixels');
axis off;
% title('Title');
% xlim([30, 115]);
% ylim([30, 120]);


%% CV Optic Analysis

% Apply spatial and temporal filtering to LAT data
T_LAT_ms2 = SpatTemp_Filtering(T_LAT_ms, 3, 0, 'GPU');
% Initialize matrices for Average Speed, Standard Deviation, and Angle
AvgSpeed_O = zeros(size(Data,1), size(Data,2));
StdSpeed_O = zeros(size(Data,1), size(Data,2));
Angle_O = zeros(size(Data,1), size(Data,2));
% Calculate conduction velocity using CV_CircleMethod
for i = 1:size(Data,1)
    for j = 1:size(Data,2)
        if max(max(squeeze(Data(i,j,:)))) ~= 0
            % (LAT matrix, radious to calc, row, col, SapaceScale in millimeters per pixel)
            [AvgSpeed_O(i,j),  StdSpeed_O(i,j), Angle_O(i,j)] = CV_CircleMethod(T_LAT_ms2, 5, i, j, 0.33); 
        end
    end
end

% Plot the conduction velocity data on a colormap
f1 = figure('color', 'white', 'Position', [40 40 400 300]);
C = jet(256);
C(1,1:3) = [1 1 1];
J = AvgSpeed_O(:);
Y = prctile(nonzeros(J),[5 95],'all');
J = imrotate(AvgSpeed_O,90);
teto = Y(2) - Y(1);
imagesc(J - Y(1), [0, Y(2) - Y(1)]);colormap(C);
hBar1 = colorbar('eastoutside');
ylabel(hBar1, 'Conduction velocity [cm/s]', 'FontSize', 18);
hold on;
d = 200; color = 'black'; tam = 14;
set(gca, 'fontsize', 18);
ylabel('Pixels');
xlabel('Pixels');
axis off
% title('Cam 2');
% xlim([30, 100]);
% ylim([0, 50]);


%% LAT Electric Analysis

% Parameters
Data = D_SYNC.EL;
Fs = 4000;
lim1 = fix(3.12*4000); % Start sample index
lim2 = fix(3.29*4000); % End sample index

% Extract a segment of electrical data
Data_temp_E = Data(:,lim1:lim2);
% Calculate Local Activation Time (LAT) for electrical data
% Plot selected point for each electrode
for i = 1:size(Data_temp_E,1)
  if max(max(Data_temp_E(i,:))) ~= 0
    xs = Data_temp_E(i,:);
    [dvdt(i), position(i)] = max(-diff(xs));

    figure()
    plot(xs);
    hold on
    plot(position(i), xs(position(i)), 'o', 'markersize', 10, 'color', 'red');
    hold off
  end
end

%manually correct the wrong values (If needed)
position(1) = 256; %error corrregido a mão para este esperimento 
LAT_ms_E = (position) * 1000 / Fs; %changing to time in ms

% Create a matrix of example activation times
minimo = min(LAT_ms_E);
LAT_ms_E = LAT_ms_E - minimo;
tiempos_activacion = [LAT_ms_E(16) LAT_ms_E(12) LAT_ms_E(8) LAT_ms_E(4);
                      LAT_ms_E(15) LAT_ms_E(11) LAT_ms_E(7) LAT_ms_E(3);
                      LAT_ms_E(14) LAT_ms_E(10) LAT_ms_E(6) LAT_ms_E(2);
                      LAT_ms_E(13) LAT_ms_E(9)  LAT_ms_E(5) LAT_ms_E(1)];
LAT_16_E = tiempos_activacion;

% Create a meshgrid of coordinates for activation times
[X, Y] = meshgrid(1:size(tiempos_activacion, 2), 1:size(tiempos_activacion, 1));

% Create a meshgrid of coordinates for interpolation
[Xq, Yq] = meshgrid(1:0.1:size(tiempos_activacion, 2), 1:0.1:size(tiempos_activacion, 1));

% Perform bilinear interpolation
tiempos_interpolados = interp2(X, Y, tiempos_activacion, Xq, Yq, 'linear');
tiempos_interpolados_RA_S = tiempos_interpolados;

% Show the activation time map
tam = 12;
f1 = figure('color', 'white', 'Position', [40 40 200 200]);
imagesc(Xq(1,:), Yq(:,1), tiempos_interpolados,[0,100]);
hBar1 = colorbar('eastoutside');
colormap("parula");
ylabel(hBar1, 'Local Activation Time [ms]', 'FontSize', 16);
axis off
hold on;
text(1, 1, num2str(16), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');
text(1, 2, num2str(15), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');
text(1, 3, num2str(14), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');
text(1, 4, num2str(13), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');
text(2, 1, num2str(12), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');
text(2, 2, num2str(11), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');
text(2, 3, num2str(10), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');
text(2, 4, num2str(9), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');
text(3, 1, num2str(8), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');
text(3, 2, num2str(7), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');
text(3, 3, num2str(6), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');
text(3, 4, num2str(5), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');
text(4, 1, num2str(4), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');
text(4, 2, num2str(3), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');
text(4, 3, num2str(2), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');
text(4, 4, num2str(1), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');
set(gca, 'fontsize', 18)
title('MEA3');
ylabel('Pixels');
xlabel('Pixels');


%% CV Electric Analysis

% This code uses the Circle method using the Local activation times
% calculated before.
% r= circle radius

r = 3;
for i = (r + 1):(size(tiempos_interpolados, 1) - (r + 1))
    for j = (r + 1):(size(tiempos_interpolados, 2) - (r + 1))
        if max(max(squeeze(tiempos_interpolados(i, j, :)))) ~= 0
            [AvgSpeed_E_RA(i, j), StdSpeed_E_RA(i, j), Angle_E_RA(i, j)] = CV_CircleMethod(tiempos_interpolados, r, i, j, 0.3);
        end
    end
end

% Spatially filter the average speed
AvgSpeed_E2 = SpatTemp_Filtering(AvgSpeed_E_RA, 5, 0, 'CPU');

% Plot the results using jet colormap
f1 = figure('color', 'white', 'Position', [40 40 230 210]);
C = jet(256);
C(1, 1:3) = [1 1 1];
J = AvgSpeed_E2(:);
Y = prctile(nonzeros(J), [5 95], 'all');
J = imrotate(AvgSpeed_E2, 90);
imagesc(J - Y(1), [0, Y(2) - Y(1)]); colormap(C);
hBar1 = colorbar('eastoutside'); ylabel(hBar1, 'Conduction velocity [cm/s]', 'FontSize', 18);
set(gca, 'fontsize', 18);
title('Cam 2');
ylabel('Pixels'); xlabel('Pixels');
axis off




















