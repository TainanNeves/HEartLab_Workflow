%% Tainan Function
function [S, D] = phasecalc_opticos(A)

    % Ensure A is a 3D array (x, y, t)
    assert(ndims(A) == 3, 'A must be a 3D array with dimensions [nx, ny, nt]');
    
    % Phase map for all pixels in A
    [nx, ny, nt] = size(A);
    D = zeros(nx, ny, nt); % phase map for all pixels
    S = zeros(nx, ny, nt); % signal data for all pixels
    
    % Loop through all pixels
    for i = 1:nx
        disp(i);
        for j = 1:ny
            % Get signal and remove DC component (mean)
            signal = squeeze(A(i, j, :)) - mean(squeeze(A(i, j, :)));
            % Calculate the Hilbert transform and phase
            hilbert_transform = imag(hilbert(signal));
            phase = -atan2(hilbert_transform, signal);
            
            % Store in the phase map and signal map
            D(i, j, :) = phase;
            S(i, j, :) = signal;
        end
    end

end






%% Old function
% A= array 3D como todos os sinais opticos
% A16= aaray signal with 16 optical electrodes one in each raw
% frame= sample that we are going to create the map 
% Fsampling = sample frequency for electrical data
% mycmap= frequency colormap
% ROI = A silhueta da imagem para fazer o mapa
% mycmap= colormap para frequencya

% function [S,D_16,D]=phasemap_opticos_nocorregido(A,A16,Fsampling, frame, ROI, mycmap)
% 
% % cycle lenght detection 
% %keyboard
% %16 array 
%  
%      [nx,ny]=size(A16);D=zeros(size(A16));
%      for i=1:nx
%         %keyboard
%         %não corregido
%         s=A16(i,:);
%         s=s-mean(s);ha=imag(hilbert(s));
%         k=-atan2(ha,s);
%         D_16(i,:)=k;
%         S(i,:)=s;
% 
% %         %corregido
% %         s = A16(i, :); % Extraer la fila i de la matriz A16
% %         s = s - mean(s); % Centrar la señal alrededor de cero
% %         ha=imag(hilbert(s));
% %         phase_original = -atan2(ha, s);
% %         signal_corregida = real(s) - ha;
% %         hilbert_transform_corregida = hilbert(signal_corregida);
% %         ha_corregida = imag(hilbert_transform_corregida);
% %         phase_corregida = -atan2(ha_corregida, signal_corregida);
% %         S(i,:)=s;
% %         D_16(i,:)=phase_corregida;
% %         k=phase_corregida;
% 
% 
%         f2=figure('color','white','Position', [40 40 200 350]);
%         subplot(2,1,1);
%         plot(s,'LineWidth', 1,'Color', 'black');%hold on; plot(signal_corregida,'Color','red'); 
%         hold on;
%         line([frame,frame], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
%         title('Señal Original');
%         xlabel('samples');
%         ylabel('Amplitud');
%         %legend('Real Part','Imaginary Part');set(gca,'fontsize', 14);
%         
%         subplot(2,1,2);
%         plot(k,'LineWidth', 1,'Color', 'black');%hold on; plot(phase_corregida,'Color','red'); %k
%         line([frame,frame], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
%         title('Fase de la Señal');
%         xlabel('samples');
%         ylabel('Fase (radianes)');
%        
%       end
% 
% 
% %all atrium  
% if length(size(A))==3
%               [nx,ny,nt]=size(A);
%               D=zeros(size(A));
%              for i=1:nx
%                  i;
%                  for j=1:ny
%                      s=squeeze(A(i,j,:));
%                      if(mean(s)==0)
%                      D(i,j,:)=[zeros(1, length(s))];    
%                      else
%                         %não corregido
%                         s=s-mean(s);
%                         ha=imag(hilbert(s));
%                         k=-atan2(ha,s);                   
%                         D(i,j,:)=k;
% 
% %                         %corregido
% %                         ha=imag(hilbert(s));
% %                         phase_original = -atan2(ha, s);
% %                         signal_corregida = real(s) - ha;
% %                         hilbert_transform_corregida = hilbert(signal_corregida);
% %                         ha_corregida = imag(hilbert_transform_corregida);
% %                         phase_corregida = -atan2(ha_corregida, signal_corregida);
% %                         D(i,j,:)=phase_original;
%                      end
%                  end
%              end
% end   
% 
% f1=figure('color','white','Position', [40 40 450 350]); I=squeeze(D(:,:,frame)),J=imrotate(I,90);
% ROI = imrotate(ROI,90);J2=J; J2(1,1)=+pi; J2(1,2)=-pi;J2 = mat2gray(J2); J2 = gray2ind(J2);J2 = ind2rgb(J2,colormap(mycmap));
%  for i=1:size(J2,1)
%      for j=1:size(J2,2)
%          if ROI(i,j) == 0
%             J2(i,j,1:3) = [1,1,1];
%          end
%      end
%  end
% imagesc(J2,[-pi pi]);h1 = colorbar;colorbar;
% set(gca,'fontsize', 14);
% 
%  
% 
% end
% 
