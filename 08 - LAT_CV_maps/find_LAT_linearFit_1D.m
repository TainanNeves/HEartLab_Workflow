function [T_LAT_ms] = find_LAT_linearFit_1D(y, fr, L, PCL, debug)
    % FIND_LAT_LINEARFIT_1D Find LAT in a single signal trace using linear fit.
    % 
    % INPUTS:
    %   y:      1D array representing the signal trace.
    %   fr:     Frame rate (number of images per second).
    %   L:      Length of the fit line around mid 50% point.
    %   PCL:    (Not used in this version)
    %   debug:  Boolean indicating whether to visualize debugging information.
    % 
    % OUTPUT:
    %   T_LAT_ms:   Activation Time (LAT) in milliseconds.
    % 
    % Note: This function requires that the trace always starts before the activation time.
    %

    % Smooth the signal
    temp = y;
    yy = smooth(y, 5);
    
    % Find max and min values in the smoothed signal
    [M, indexM] = max(yy);
    [m, indexm] = min(yy);
    
    % Ensure minimum value occurs before maximum value
    while indexm > indexM
        yy(1) = [];
        [M, indexM] = max(yy);
        if isempty(M)
            break;
        end
        indexM = indexM(1);
    end
    
    % Adjust indexM and indexm if necessary
    if indexm < indexM
        indexM = indexM + length(y) - length(yy);
        indexm = indexm + length(y) - length(yy);
    else
        temp = y;
        yy = smooth(y, 5);
        [M, indexM] = max(yy);
        [m, indexm] = min(yy);
        
        while indexm > indexM
            yy(end) = [];
            [M, indexM] = max(yy);
            [m, indexm] = min(yy);
        end
    end
    
    % Find the 50% point (LAT_v) in the signal
    LAT_v = ((M - m) / 2) + m;
    temp = y(indexm:indexM);
    [~, position] = min(abs(temp - LAT_v));
    T_LAT = (indexm + position);  % Location of LAT in y over the whole trace
  
    % Find the linear fit line
    if T_LAT - L <= 0
        x = 1:T_LAT + L;
    else
        x = T_LAT - L:min(T_LAT + L, length(y));
    end
    
    P = polyfit(x, y(x), 1);
    
    % Solve for T_LAT
    if P(1) < 0
        plot(y);
        disp('Correction: Slope needs to be greater than 0');
        T_LAT_ms = T_LAT;
    else
        x = (LAT_v - P(2)) / P(1);
        T_LAT = x;
    end
    
    % Check if T_LAT is within the range of indexm and indexM
    if T_LAT > indexM || T_LAT < indexm
        plot(y);
        disp('Correction: LAT outside min - max range');
        T_LAT_ms = T_LAT;
    end
    
    % Convert T_LAT to milliseconds
    T_LAT_ms = (T_LAT) * 1000 / fr;
    
    % Debug visualization
    if debug
        figure(101);
        plot(y);
        hold on;
        valid_x = max(1, min(length(y), round(x)));
        plot(valid_x, y(valid_x), 'o', 'markersize', 10, 'color', 'red');
        hold off;
    end
end



%% Original function From Jimena
% function [T_LAT_ms]= find_LAT_linearFit_1D(y,fr,L, PCL, debug)
%     % y - 1D array
%     % fr - frame rate, number of images per 1 secons
%     % +-L fit line linear lengthg
%     % finds LAT in a single signal trace with linear fit around mid 50% point as analytical solution
%     % of the fit function
%     % this function requires that trace always starts before the activation
%     % time
%     % debug is for visualization purposes
%     %
%  
% % y=squeeze(DATA3(80,80,:));;
% % fr=500;
% % L=7;
% % debug=1;
% %keyboard
%     temp = y;
%     yy = smooth(y,5);
%     
%     M=max(yy); indexM = find(yy==M);indexM = indexM(1);% find the location of the maximum valiu 
%     m=min(yy); indexm = find(yy==m);indexm = indexm(1);% find the location of the minimum valiu 
%     
%     %this part only confirm that the minimum value is before the maximum
%     %value
%     while indexm>indexM
%         yy(1) = [];
%         M=max(yy); 
%         if isempty(M)
%             break
%         end          
%         indexM = find(yy==M);
%         indexM = indexM(1);
%     end
%     if indexm < indexM
%         indexM = indexM + length(y)-length(yy);
%         indexm = indexm + length(y)-length(yy);
%     else
%         temp = y;
%         yy = smooth(y,5);
%         M=max(yy); indexM = find(yy==M);indexM = indexM(1);
%         m=min(yy); indexm = find(yy==m);indexm = indexm(1);
%         
%         while indexm>indexM
%             yy(end)=[];
%             M=max(yy); indexM = find(yy==M);indexM = indexM(1);
%             m=min(yy); indexm = find(yy==m);indexm = indexm(1);
%         end
%     end
% % keyboard;   
% %agarrar solo la mitad para encontrar el 50% 
%     LAT_v=((M-m)/2)+m;%encuentro el valor del 50 porciento en el tiempo
%     temp=y(indexm:indexM);
%     l=length(temp);
%     [minimo,position]=min(abs(temp - LAT_v));
%     T_LAT = (indexm+position);%-1;% locacion de lat en y todo el tramo
%   
%  
% %find fit line linear
%     if T_LAT-L<=0
%        x = [1:T_LAT+L];
%     else
%        x = [T_LAT-L:min(T_LAT+L, length(y))];
%     end
%     P = polyfit( x,y(x),1);
%     %solve for T_LAT
%     %yfit = P(1)*x+P(2);
%      if P(1)<0 
%         plot(y)
%         display('correction slope need to be greater than 0')
%         T_LAL = T_LAT    
%      else
%         x = (LAT_v - P(2)) / P(1);
%         T_LAT = x;
%      end
%      
%      % check if point is in betwerern min and max, if not take the old
%      % value
%      if T_LAT>indexM | T_LAT<indexm
%         %T_LAT = (indexm+position) ;
%         plot(y)
%         display('correction LAT outside min - max range')
%         T_LAL = T_LAT    
%      end
%      
% %convertir a milisegundos
% T_LAT_ms=(T_LAT)*1000/fr;
%     %T_LAT_ms=(T_LAT-indexm)*1000/fr;
%     
%     if debug
%         figure(101)
%         plot(squeeze(y))
%         hold on
%         plot(round(x), y(round(x)),'o','markersize',10,'color','red')
%         hold off
%     end