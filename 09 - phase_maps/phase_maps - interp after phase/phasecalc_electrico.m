function [S, D] = phasecalc_electrico(A)
    % Check if A is a 2D array
    if length(size(A)) ~= 2
        error('Input A must be a 2D array with dimensions [nx, ny]');
    end
    
    % Initialize arrays
    [nx, ny] = size(A);
    S = zeros(nx, ny); % Signal data for all electrodes
    D = zeros(nx, ny); % Phase map for all electrodes
    
    % Loop through all electrodes
    for i = 1:nx
        disp(i);
        % Get signal and remove the DC component (mean)
        s = A(i, :) - mean(A(i, :));
        
        % Calculate the Hilbert transform and phase
        hilbert_transform = imag(hilbert(s));
        phase = -atan2(hilbert_transform, s);
        
        % Store the signal and phase
        S(i, :) = s;
        D(i, :) = phase;
    end
end



%%
% % this function need twoconsecutive beats to calculate the CL 
% % A= aaray signal with 16 electrodes one in each raw
% % frame= sample that we are going to create the map 
% % Fs = sample frequency for electrical data
% % mycmap= frequency colormap
% % MEA= MEA number to be plotted 
% % dib ...
% 
% function [S,D]=phasemap_electrico_nocorregido(A, frame, Fs,  mycmap, MEA)
% 
% 
% % phase calculation 
% if length(size(A))==2
%      [nx,ny]=size(A);D=zeros(size(A));
%      for i=1:nx
%         s=A(i,:);
%         s=s-mean(s);
%         %não corregido 
% 
%         ha=imag(hilbert(s));
%         k=-atan2(ha,s);
%         S(i,:)=s;
%         D(i,:)=k;
% 
%         f2=figure('color','white','Position', [40 40 200 350]);
%         subplot(2,1,1);
%         plot(s,'LineWidth', 1,'Color', 'black');hold on;  
%         line([frame,frame], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
%         title('Señal Original');
%         xlabel('samples');
%         ylabel('Amplitud');
%         legend('Real Part','Imaginary Part');set(gca,'fontsize', 14);
%         
%         subplot(2,1,2);
%         plot(k,'LineWidth', 1,'Color', 'black');hold on;%k
%         line([frame,frame], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
%         title('Fase de la Señal');
%         xlabel('samples');
%         ylabel('Fase (radianes)');
%         legend('-atan','phase');set(gca,'fontsize', 14);
%         linkaxes([subplot(2, 1, 1), subplot(2, 1, 2)], 'x');
%         
%      end
% end
% 
% array=D(:,frame);
% MEA_LP_PHASE(array,[-pi pi],frame,MEA,mycmap);
% 
% end
%   