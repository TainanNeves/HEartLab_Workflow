% this script converts SEQ files from Norpix to MATLAB, and exports
% 1) representative plot trace from each recording
% 2) representative image      from each recording
% for multiple camera recordings the script converts all the recordings
% from all cameras and saves them into the same folder.
% The output folder is 'Raw converted to MATLAB'
% additionally the data is binned. For IMX425 based camera, the bining is
% with defalt factor of 4 reducing the images to under 256x256
% as for IMX425 data is recorded at 12bits, no need to rescaling to fit in
% 16-bit range with 4x4 binning

binSize = 8;
precision = 'uint16'; %for binning

String = {'Camera 1', 'Camera 2', 'Camera 3'}

for N=1:3
    cd 'Optical Mapping'
    cd(String{N})
    DIR = dir('*.seq');
    for i=1:length(DIR)
        P=pwd;
        [header DATA] = Norpix2MATLAB( DIR(i).name,0,0);
        cd ..//..
        cd  'Optical Mapping to MATLAB'
        figure,plot(squeeze(mean(mean(DATA(501:504,501:504,:),1),2)))
        saveas(gcf, [DIR(i).name,'_Trace_Cam',num2str(N),'.png'])
        figure,imagesc(DATA(:,:,100));colormap('gray'); axis square
        saveas(gcf, [DIR(i).name,'_Image_Cam',num2str(N),'.png'])
        
        DATA_old = DATA;
        DATA = Binning(DATA,binSize,precision);
        savefast([DIR(i).name(1:end-4),'_Bin=',num2str(binSize),'_Cam',num2str(N),'.mat'], 'DATA')
        
       %DATA = Binning(DATA_old,binSize*2,'single');
        %DATA = uint16(DATA/4);
        %savefast([DIR(i).name(1:end-4),'_Bin=',num2str(binSize*2),'_Cam',num2str(N),'.mat'], 'DATA')
        
        cd(P)
    end
    cd ..//..
end
