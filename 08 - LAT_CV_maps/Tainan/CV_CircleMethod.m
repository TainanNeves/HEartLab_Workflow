function [AvgSpeed,  StdSpeed, Angle]= CV_CircleMethod(LAT, r, row, col, SpaceScale)
% Get dimensions of the input matrix
[maxR, maxC] = size(LAT);

% define the circle
    theta = linspace(0, 2*pi,180); L = length(theta);
    col2 = col + r*cos(theta);
    row2 = row + r*sin(theta);
    
% evaluate LAT along the circle
    LAT_Circle = zeros(1, L); % Pre-allocate for speed
    for k=1:L
        % 1. Determine local floor coordinates
        C = floor(col2(k));
        R = floor(row2(k));
        
        % 2. BOUNDARY CHECK: Ensure R and C are within [1, max-1]
        % We use max-1 because interpolation looks at R+1 and C+1
        R = max(1, min(R, maxR - 1));
        C = max(1, min(C, maxC - 1));
        
        temp1 = (row2(k) - floor(row2(k)));
        temp2 = (col2(k) - floor(col2(k)));

        % Now this line is safe from "Index out of bounds"
        LAT_Circle(k) = ( (1-temp1)*LAT(R,C) + temp1*LAT(R+1,C) + (1-temp2)*LAT(R,C) +  temp2*LAT(R,C+1) ) / 2;
    end

% evaluate CV along the circle
    Distance = 2*r*SpaceScale; % [mm]
    CV_f = Distance./( LAT_Circle(1:L/2) - LAT_Circle(L/2+1:L) ); %[mm/ms]
    CV_f = CV_f*100; % converting to [cm/s]
    
% smoothing data to detect propagation direction (index value)
    sigma = 5;
    % Handle potential NaNs or Infs in CV_f before smoothing
    CV_f(isinf(CV_f)) = 0; 
    CV_f(isnan(CV_f)) = 0;
    
    CV_f_smooth = smoothdata( abs([CV_f(end-sigma+1:end), CV_f, CV_f(1:sigma)]),'gaussian',sigma);
    [~,loc] = min(CV_f_smooth(sigma+1:end-sigma));

% rotate un-smoothed CV 
    CV_oriented = circshift(abs(CV_f),-loc-L/4);

% average around the index value 
    w=floor(15/(360/L));
    % Ensure indices don't go out of bounds during mean calculation
    idx_range = (L/4-w):(L/4+w);
    idx_range(idx_range < 1) = 1;
    idx_range(idx_range > length(CV_oriented)) = length(CV_oriented);
    
    AvgSpeed = abs(mean(CV_oriented(idx_range)./cosd( [-w:+w]*360/L)  ));
    StdSpeed = std(CV_oriented(idx_range)./cosd( [-w:+w]*360/L)  );
    
% find angle
    if CV_f(loc)<0
        Angle =  loc*(360/(L)) + 180 ; 
    else
        Angle = loc*(360/(L)); 
    end
end