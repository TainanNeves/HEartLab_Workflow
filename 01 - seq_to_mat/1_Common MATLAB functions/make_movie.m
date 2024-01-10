

function make_movie(ts,outfilename,crange,C,FPS)
%% Generate a movie of 2-D time series
% use C for Vm and C1 for Ca
%C  = customcolormap([0 (0.5-0)/(1-0) 1], [0.6,0.4 0; 1 1 1;  0.4 0  0.9; ] ,64);%[1 0 0; 1 1 1; 0 0 1],32);
%C1 = customcolormap([0 (0.5-0)/(1-0) 1], [0.4,0.8 0; 1 1 1;  0.4 0  0.9; ] ,64);%[1 0 0; 1 1 1; 0 0 1],32);

% INPUT:    
%   ts          ... 2-D time series [N x M x time]
%   outfilename ... Filename of movie file (e.g. 'orig_movie.avi')
%   crange      ... Range of color axis (e.g. [0 1] for excitation variable; [-p pi] for phase
%
% OUTPUT:
%   Movie file (default is .avi but can be changed into any format)

% Show frames
h= figure(100);
i0 = zeros(size(ts(:,:,1)));
ih = imagesc(i0(:,:,1)); caxis(crange);
colormap(C); axis image off; 
set(gcf,'position',[200 200 512 512],'color',[1 1 1])
%ts = ts(:,:,1:5:end);
for frame=1:size(ts,3)
    set(ih,'cdata',ts(:,:,frame));colormap(C);
    colorbar
    title( [num2str(frame), ' ms'],'fontsize',16) 
    drawnow
    mov(frame) = getframe(gcf);
end

% Make a movie
writerObj = VideoWriter(outfilename,'MPEG-4');
writerObj.FrameRate = FPS;
writerObj.Quality = 75;
open(writerObj);
writeVideo(writerObj,mov);
close(writerObj);
close all
clear mov