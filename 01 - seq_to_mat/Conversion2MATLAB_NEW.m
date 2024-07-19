binSize = 8;
precision = 'uint16'; %for binning

String = {'Camera 1', 'Camera 2', 'Camera 3'};

% Define the base directory for saving files
baseDir = fullfile(pwd, 'gravação');

% Create the base directory if it does not exist
if ~exist(baseDir, 'dir')
    mkdir(baseDir);
end

for N = 1:3
    % Define the subfolder for the current iteration
    subDir = fullfile(baseDir, sprintf('Subpasta_%d', N));
    
    % Create the subfolder if it does not exist
    if ~exist(subDir, 'dir')
        mkdir(subDir);
    end
    
    cd 'Optical Mapping'
    cd(String{N})
    DIR = dir('*.seq');
    
    for i = 1:length(DIR)
        P = pwd;
        [header DATA] = Norpix2MATLAB(DIR(i).name, 0, 0);
        cd ..//..
        cd 'Optical Mapping to MATLAB'
        
        % Plot and save the trace
        figure;
        plot(squeeze(mean(mean(DATA(501:504, 501:504, :), 1), 2)));
        saveas(gcf, fullfile(subDir, [DIR(i).name, '_Trace_Cam', num2str(N), '.png']));
        close(gcf); % Close the figure to avoid having too many open figures
        
        % Plot and save the image
        figure;
        imagesc(DATA(:, :, 100)); colormap('gray'); axis square;
        saveas(gcf, fullfile(subDir, [DIR(i).name, '_Image_Cam', num2str(N), '.png']));
        close(gcf); % Close the figure to avoid having too many open figures
        
        % Perform binning and save the binned data
        DATA_old = DATA;
        DATA = Binning(DATA, binSize, precision);
        savefast(fullfile(subDir, [DIR(i).name(1:end-4), '_Bin=', num2str(binSize), '_Cam', num2str(N), '.mat']), 'DATA');
        
        % Optional: Additional binning (commented out in the original code)
        % DATA = Binning(DATA_old, binSize*2, 'single');
        % DATA = uint16(DATA / 4);
        % savefast(fullfile(subDir, [DIR(i).name(1:end-4), '_Bin=', num2str(binSize*2), '_Cam', num2str(N), '.mat']), 'DATA');
        
        cd(P);
    end
    
    % Go back to the 'Optical Mapping' directory
    cd ..//..
end