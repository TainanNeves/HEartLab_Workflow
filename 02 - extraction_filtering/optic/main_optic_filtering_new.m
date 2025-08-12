<<<<<<< HEAD
=======
% This code performs filtering of optical data, divided into two main sections. 
% The first section concerns providing the directory address where the folder containing 
% the optical recording of interest for analysis is located. 
% The second section specifically handles the filtering procedure for the data from the 
% three cameras of the recording selected in the first section. 
% The code processes the optical data filtering for the three cameras of the chosen 
% recording in a single execution within a FOR loop and terminates at the end. 
% A file is saved and exported containing the data from the three cameras of the 
% corresponding recording, the images from each camera, and the ROI for each camera.

>>>>>>> Development
%% Code for Filtering and Export Optical Data

% Clear workspace and close existing figures
clear; clc; close all;

<<<<<<< HEAD
%% Directory Selection for Recording Files

DS = uigetdir(); % DS is the selected directory
% Directory Check
if DS == 0 % Directory not selected or empty
    disp('No directory selected.'); 
else
    disp(['Selected directory: ', DS]);
    CDS = DS; % Save the Directory Path -> CDS is the path of the selected directory
    CPDS = dir(CDS); % List directories in the path -> CPDS is the content of the selected directory
    PDS = CPDS([CPDS.isdir]); % Filter directories, ignoring files -> PDS is the directories in the selected path
    PDS = PDS(~ismember({PDS.name}, {'.', '..'})); % Exclude '.' and '..' directories
    NP = length(PDS); % Count the number of directories
=======
%% Directory Selection for a Single Recording File

DS = uigetdir(); % DS is the selected directory

% Directory Check
if DS == 0 % Directory not selected or empty
    disp('No directory selected.'); 
    return; % Exit the script
else
    disp(['Selected directory: ', DS]);
    CDS = DS; % Save the Directory Path -> CDS is the path of the selected directory
    % Get the list of files in the directory
    files = dir(fullfile(CDS, '*.mat')); % Get all files in the directory with '.mat' format
    if isempty(files)
        disp('No .mat files found in the selected directory.');
        return; % Exit the script
    end
>>>>>>> Development
end

%% Optical Signals Pre-Processing

<<<<<<< HEAD
for i = 1:NP
    NameP = PDS(i).name; % Save the Directory Name
    CNameP = fullfile(CDS, NameP); % Path of the directory
    disp(['Entering the directory: ', CNameP]); % Display the directory name and path
    % Get the list of files in the directory
    files = dir(fullfile(CNameP, '*.mat')); % Get all files in the directory with '.mat' format
    %files = files(~[files.isdir]); % Filter to keep only files

for cam_idx = 1:3

    Namefiles = files(cam_idx).name; % Specific file name
    Cfiles = fullfile(CNameP, Namefiles); % File path
=======
for cam_idx = 1:3
    % Check if there are enough files for the cameras
    if cam_idx > length(files)
        disp(['Missing file for Camera ', num2str(cam_idx), '. Skipping.']);
        continue;
    end

    % Process the file for the current camera
    Namefiles = files(cam_idx).name; % Specific file name
    Cfiles = fullfile(CDS, Namefiles); % File path
>>>>>>> Development

    filename = fullfile(Cfiles); % Save the file address for the iteration

    load(filename); % Optical Data Binned

<<<<<<< HEAD
=======
    % Define ROI for each camera
>>>>>>> Development
    switch cam_idx
        case 1 %CAM 1
           ROI_CAM1 = roipoly(DATA(:,:,100));
        case 2 %CAM 2
           ROI_CAM2 = roipoly(DATA(:,:,100));
        case 3 %CAM 3
           ROI_CAM3 = roipoly(DATA(:,:,100));
    end

<<<<<<< HEAD
    %First Plot
=======
    % First Plot
>>>>>>> Development
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
<<<<<<< HEAD
    % Creating image to see electrodes

    %Run if you need it
=======

    % Run if you need it
>>>>>>> Development
    Background = squeeze(DATA(:,:,20));
    [x, y] = pick_up_a_trace(Background, DATA, 1);

    % Extract a pixel intensity profile along a line and plot it
    row_pixel = x(length(x));
    col_pixel = y(length(y));
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
    less = 500; % Samples eliminated at the beginning and end (1 second each)

    % Apply baseline filter to the data
    [DATA1] = -f_fil(R, Fpass, Fsampling, n, less);
    
    % Apply a Gaussian spatial and temporal filtering to the data
<<<<<<< HEAD
    % SpatTemp_Filtering (3D Data, S = Size Gaussian filter matrix (must be
    %   odd), T = Temporal filter matrix size, mode = 'CPU' or 'GPU') 
    %   Using multiple times the result becomes smooth, but you lose the
    %   electrode positions
    DATA1 = SpatTemp_Filtering(DATA1, 3, 500, 'GPU');
    % DATA1 = SpatTemp_Filtering(DATA1, 3, 500, 'GPU');
    % DATA1 = SpatTemp_Filtering(DATA1, 3, 500, 'GPU');
    % DATA1 = SpatTemp_Filtering(DATA1, 7, 500, 'GPU');
    % DATA1 = SpatTemp_Filtering(DATA1, 7, 500, 'GPU');
=======
    DATA1 = SpatTemp_Filtering(DATA1, 3, 500, 'GPU');
>>>>>>> Development

    % Apply the mask to the filtered data
    switch cam_idx
        case 1 %CAM 1
            for i = 1:size(DATA1,3)
                DATA3(:,:,i) = squeeze(DATA1(:,:,i)) .* ROI_CAM1;
            end
        case 2 %CAM 2
<<<<<<< HEAD
           for i = 1:size(DATA1,3)
                DATA3(:,:,i) = squeeze(DATA1(:,:,i)) .* ROI_CAM2;
            end
        case 3 %CAM 3
           for i = 1:size(DATA1,3)
=======
            for i = 1:size(DATA1,3)
                DATA3(:,:,i) = squeeze(DATA1(:,:,i)) .* ROI_CAM2;
            end
        case 3 %CAM 3
            for i = 1:size(DATA1,3)
>>>>>>> Development
                DATA3(:,:,i) = squeeze(DATA1(:,:,i)) .* ROI_CAM3;
            end
    end

    % Further process the data with a baseline function
    DATA2 = BASE(DATA3);

    % Save the processed data
    switch cam_idx
        case 1 %CAM 1
            D_CAM1_filtered = DATA2;
            D_CAM1_image(:,:,1) = J(:,:);
        case 2 %CAM 2
            D_CAM2_filtered = DATA2;
            D_CAM2_image(:,:,1) = J(:,:);
        case 3 %CAM 3
           D_CAM3_filtered = DATA2;
           D_CAM3_image(:,:,1) = J(:,:);
    end

<<<<<<< HEAD
=======
    % Clear variables for the next iteration
>>>>>>> Development
    clear C col_pixel DATA DATA1 DATA2 DATA3 f Fcut Fpass Fs Fsampling;
    clear i I J less n p p2 R row_pixel ROI ans Background;
end

<<<<<<< HEAD
%Exporting the Optical Data

identification_recording = NameP;
filenamerecording = ['optic_data_',identification_recording,'_filtered.mat'];
=======
% Exporting the Optical Data

[~, recording_name] = fileparts(DS); % Extract the name of the recording folder
filenamerecording = ['optic_data_', recording_name, '_filtered.mat'];
>>>>>>> Development
% Creating variables
D_OP.D_CAM1_filtered = D_CAM1_filtered;
D_OP.D_CAM1_rawimage = D_CAM1_image;
D_OP.D_CAM2_filtered = D_CAM2_filtered;
D_OP.D_CAM2_rawimage = D_CAM2_image;
D_OP.D_CAM3_filtered = D_CAM3_filtered;
D_OP.D_CAM3_rawimage = D_CAM3_image;
D_OP.ROI.ROI_1 = ROI_CAM1;
D_OP.ROI.ROI_2 = ROI_CAM2;
D_OP.ROI.ROI_3 = ROI_CAM3;

% Save variables to the .mat file
save(filenamerecording, 'D_OP', '-v7.3');
disp(['Variables saved to ', filenamerecording]);

<<<<<<< HEAD
clear D_CAM1_filtered  D_CAM2_filtered  D_CAM3_filtered
clear D_CAM1_image     D_CAM2_image     D_CAM3_image
clear ROI_CAM1 ROI_CAM2 ROI_CAM3 x y cam_idx

close all;

end

clear CDS CPDS DS NP PDS
clear Cfiles CNameP filename filenamerecording identification_recording
clear D_OP files Namefiles NameP
=======
% Final cleanup
clear D_CAM1_filtered  D_CAM2_filtered  D_CAM3_filtered;
clear D_CAM1_image     D_CAM2_image     D_CAM3_image;
clear ROI_CAM1 ROI_CAM2 ROI_CAM3 x y cam_idx;

close all;
>>>>>>> Development
