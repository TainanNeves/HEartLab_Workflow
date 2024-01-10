function [T_LAT_ms]= find_LAT_linearFit_1D(y,fr,L, PCL, debug)
    % y - 1D array
    % fr - frame rate, number of images per 1 secons
    % +-L fit line linear lengthg
    % finds LAT in a single signal trace with linear fit around mid 50% point as analytical solution
    % of the fit function
    % this function requires that trace always starts before the activation
    % time
    % debug is for visualization purposes
    %
    temp = y;
    yy = smooth(y,5);
    
    M=max(yy); indexM = find(yy==M);indexM = indexM(1);
    m=min(yy); indexm = find(yy==m);indexm = indexm(1);
    
    while indexm>indexM
        yy(1) = [];
        M=max(yy); 
        if isempty(M)
            break
        end          
        indexM = find(yy==M);
        indexM = indexM(1);
    end
    if indexm < indexM
        indexM = indexM + length(y)-length(yy);
        indexm = indexm + length(y)-length(yy);
    else
        temp = y;
        yy = smooth(y,5);
        M=max(yy); indexM = find(yy==M);indexM = indexM(1);
        m=min(yy); indexm = find(yy==m);indexm = indexm(1);
        
        while indexm>indexM
            yy(end)=[];
            M=max(yy); indexM = find(yy==M);indexM = indexM(1);
            m=min(yy); indexm = find(yy==m);indexm = indexm(1);
        end
    end
    
%agarrar solo la mitad para encontrar el 50% 
    LAT_v=((M-m)/2)+m;%encuentro el valor del 50 porciento 
    temp=y(indexm:indexM);
    l=length(temp);
    [minimo,position]=min(abs(temp - LAT_v));
    T_LAT = (indexm+position);%-1;% locacion de lat en y todo el tramo
  
 
%find fit line linear
    if T_LAT-L<=0
       x = [1:T_LAT+L];
    else
       x = [T_LAT-L:min(T_LAT+L, length(y))];
    end
    P = polyfit( x,y(x),1);
    %solve for T_LAT
    %yfit = P(1)*x+P(2);
     if P(1)<0 
        plot(y)
        display('correction slope need to be greater than 0')
        T_LAL = T_LAT    
     else
        x = (LAT_v - P(2)) / P(1);
        T_LAT = x;
     end
     
     % check if point is in betwerern min and max, if not take the old
     % value
     if T_LAT>indexM | T_LAT<indexm
        %T_LAT = (indexm+position) ;
        plot(y)
        display('correction LAT outside min - max range')
        T_LAL = T_LAT    
     end
     
%convertir a milisegundos
    T_LAT_ms=T_LAT*1000/fr;
    
    if debug
        figure(101)
        plot(squeeze(y))
        hold on
        plot(round(x), y(round(x)),'o','markersize',10,'color','red')
        hold off
    end