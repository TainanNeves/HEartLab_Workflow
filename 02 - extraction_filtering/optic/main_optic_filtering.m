%% Code for Filtering and export optical data

% Clear workspace and close existing figures
clear all; close all; clc;


%% Optical signals pre-processing
% Run for all 3 cameras and export data to future analysis
% Run the code with F9 part by part

% Load data and region of interest (ROI)
load('C:\Users\HEartLab\Documents\GitHub\HEartLab\00 - examples\Rec_15_23_24_Bin=8_Cam3.mat'); % Optical Data Binned
ROI_3 = roipoly(DATA(:,:,50)); % If you want to create a ROI / Double click in the center to finish

% If you already have the ROIs
load("ROI_E18.mat"); 

%First Plot
% Initialize color map and adjust first color to white
C = jet(256);
C(1,1:3) = [1 1 1];
% Extract a single frame from the data and rotate it 90 degrees
I = squeeze(DATA(:,:,50));
J = imrotate(I, 90);
% Create a figure and display the rotated image
f = figure('color', 'white', 'Position', [40 40 400 500]);
imagesc(J);
colormap('gray');
set(gca, 'fontsize', 14);
title('Cam');
ylabel('Pixels');
xlabel('Pixels');

%Run if you need it
Background = squeeze(DATA(:,:,20));
pick_up_a_trace(Background, DATA,1);    % Select a pixel in the image and shows the optical signal
                                        %Press space to stop

% Extract a pixel intensity profile along a line and plot it
row_pixel = 47;
col_pixel = 110;
p = double(squeeze(DATA(row_pixel,col_pixel,:)));
figure;
plot(p);

% Calculate deltaF/F (percentage) and plot the result
R = DELTAF(DATA);
p2 = double(squeeze(R(row_pixel,col_pixel,:)));
figure;
plot(p2);

% Remove DC baseline
Fsampling = 500; % Sampling frequency
Fpass = 0.5; % Passband frequency
Fcut = 1; % Cutoff frequency
n = 3; % Filter order
less = 500; % Samples eliminated at the beginning and final (1 second each)
% Apply baseline filter to the data
[DATA1] = -f_fil(R, Fpass, Fsampling, n, less);
% Apply a Gaussian spatial and temporal filtering to the data
%   SpatTemp_Filtering (3D Data, S = Size Gaussian filter matrix (must be
%   odd), T = Temporal filter matrix size, mode = 'CPU' or 'GPU') 
%   Using multiple times the result becomes Smooth, but You lose the
%   electrode positions
DATA1 = SpatTemp_Filtering(DATA1, 3, 500, 'GPU');
% DATA1 = SpatTemp_Filtering(DATA1, 3, 500, 'GPU');
% DATA1 = SpatTemp_Filtering(DATA1, 3, 500, 'GPU');
% DATA1 = SpatTemp_Filtering(DATA1, 7, 500, 'GPU');
% DATA1 = SpatTemp_Filtering(DATA1, 7, 500, 'GPU');

% Apply the mask to the filtered data
ROI = ROI_3; % put the correct ROI
for i = 1:size(DATA1,3)
    DATA3(:,:,i) = squeeze(DATA1(:,:,i)) .* ROI;
end

% Further process the data with a baseline function
DATA2 = BASE(DATA3);


%% Put in a variable
D_CAM3_filtered = DATA2; % Put the corect number for the camera

clear C col_pixel DATA DATA1 DATA2 DATA3 f Fcut Fpass Fs Fsampling;
clear i I J less n p p2 R row_pixel ROI ans Background;

%% Exporting the Optical data
identification = 'E14_F3_R4';
fileName = ['optic_data_', identification, '_filtered.mat'];
% Creating variables
D_OP.D_CAM1_filtered = D_CAM1_filtered;
D_OP.D_CAM2_filtered = D_CAM2_filtered;
D_OP.D_CAM3_filtered = D_CAM3_filtered;
D_OP.ROI.ROI_1 = ROI_1;
D_OP.ROI.ROI_2 = ROI_2;
D_OP.ROI.ROI_3 = ROI_3;
% Save variables to the .mat file
save(fileName, 'D_OP', '-v7.3');
disp(['Variables saved to ', fileName]);






