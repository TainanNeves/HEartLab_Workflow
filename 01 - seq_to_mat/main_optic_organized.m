% This script organizes the files obtained using the Conversion2MATLAB code; 
% The organization is achieved by creating folders/directories for each record; 
% The code is divided into two sections: 
%   one section is for running the code with an underscore (_) and 
%   the other section is for running the code with a hyphen (-). 
% These two different sections correspond to the naming format of the files, 
% so use the appropriate section according to the record's name. 
% You only need to provide the directory of the folder where the files are located, 
% which is usually the "Optical Mapping to MATLAB" folder.

clear; clc;

%% Run with Underline (_)

% Step 1: Select the directory
baseDir = uigetdir(pwd, 'Select the directory where the files are located');
if baseDir == 0
    disp('No directory selected.');
    return;
end

% Step 2: List all .mat and .png files in the directory
files = dir(fullfile(baseDir, '*.mat'));
files = [files; dir(fullfile(baseDir, '*.png'))]; % Add .png files to the list

% Initialize a struct to store files by subfolder
fileStruct = struct();

% Step 3: Identify similar names and organize the files
for k = 1:length(files)
    [~, name, ~] = fileparts(files(k).name);
    
    % Try to find the numbering in the format XX_XX_XX
    tokens = regexp(name, '(\d{2}_\d{2}_\d{2})', 'tokens');
    if ~isempty(tokens)
        key = tokens{1}{1}; % Get the found numbering
        
        % Ensure that the field name is a valid format
        validKey = matlab.lang.makeValidName(['REC_', key]);
        
        % Add the file to the list of files for the corresponding subfolder
        if isfield(fileStruct, validKey)
            fileStruct.(validKey) = [fileStruct.(validKey), {fullfile(baseDir, files(k).name)}];
        else
            fileStruct.(validKey) = {fullfile(baseDir, files(k).name)};
        end
    end
end

% Step 4: Create subfolders and move the files
keys = fieldnames(fileStruct);
for k = 1:length(keys)
    subfolder = fullfile(baseDir, keys{k});
    if ~exist(subfolder, 'dir')
        mkdir(subfolder);
    end
    
    % Move files to the corresponding subfolder
    for i = 1:length(fileStruct.(keys{k}))
        movefile(fileStruct.(keys{k}){i}, subfolder);
    end
end

disp('Files have been successfully organized.');

% Clear variables
clear baseDir files fileStruct i k key keys name subfolder tokens validKey

%% %% Run with Hyphen (-)

% Step 1: Select the directory
baseDir = uigetdir(pwd, 'Select the directory where the files are located');
if baseDir == 0
    disp('No directory selected.');
    return;
end

% Step 2: List all .mat and .png files in the directory
files = dir(fullfile(baseDir, '*.mat'));
files = [files; dir(fullfile(baseDir, '*.png'))]; % Add .png files to the list

% Initialize a struct to store files by subfolder
fileStruct = struct();

% Step 3: Identify similar names and organize the files
for k = 1:length(files)
    [~, name, ~] = fileparts(files(k).name);
    
    % Try to find the numbering in the format XX_XX_XX
    tokens = regexp(name, '(\d{2}-\d{2}-\d{2})', 'tokens');
    if ~isempty(tokens)
        key = tokens{1}{1}; % Get the found numbering
        
        % Ensure that the field name is a valid format
        validKey = matlab.lang.makeValidName(['REC_', key]);
        
        % Add the file to the list of files for the corresponding subfolder
        if isfield(fileStruct, validKey)
            fileStruct.(validKey) = [fileStruct.(validKey), {fullfile(baseDir, files(k).name)}];
        else
            fileStruct.(validKey) = {fullfile(baseDir, files(k).name)};
        end
    end
end

% Step 4: Create subfolders and move the files
keys = fieldnames(fileStruct);
for k = 1:length(keys)
    subfolder = fullfile(baseDir, keys{k});
    if ~exist(subfolder, 'dir')
        mkdir(subfolder);
    end
    
    % Move files to the corresponding subfolder
    for i = 1:length(fileStruct.(keys{k}))
        movefile(fileStruct.(keys{k}){i}, subfolder);
    end
end

disp('Files have been successfully organized.');

% Clear variables
clear baseDir files fileStruct i k key keys name subfolder tokens validKey



