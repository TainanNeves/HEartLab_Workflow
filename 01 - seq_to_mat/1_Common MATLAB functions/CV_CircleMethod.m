function [AvgSpeed,  StdSpeed, Angle]= CV_CircleMethod(LAT, r, row, col, SpaceScale)
% Calculate conduction velocity via the Circle Method
%
% [AvgSpeed, Angle]= CV_CircleMethod( LAT, r, row, col, SpaceScale)
%
% Inputs
%          LAT     2-D matrix of Local Activation Time in [ms]
%            r     radius of the circle used to calulate time difference in [px]
%          row     row at which CV is being calculated
%          col     column at which CV is being calculated
%   SpaceScale     spatial resolution in [mm/px]
%
% Outputs
%     AvgSpeed     speed of the wave front
%        Angle     angular orientation of the wave front propagation 

% define the circle
	theta = linspace(0, 2*pi,180); L = length(theta);
    col2 = col + r*cos(theta);
    row2 = row + r*sin(theta);
    
% evaluate LAT along the circle
    for k=1:L
        temp1 = (row2(k) - floor(row2(k)));
        temp2 = (col2(k) - floor(col2(k)));
        C = floor(col2(k));
        R = floor(row2(k));
        LAT_Circle(k) = ( (1-temp1)*LAT(R,C) + temp1*LAT(R+1,C) + (1-temp2)*LAT(R,C) +  temp2*LAT(R,C+1) ) / 2;
    end
    
% evaluate CV along he circle
    Distance = 2*r*SpaceScale; % [mm]
    CV_f = Distance./( LAT_Circle(1:L/2) - LAT_Circle(L/2+1:L) ); %[mm/ms]
    CV_f = CV_f*100; % converting to [cm/s]

% smoothing data to detect propagation direction (index value)
    sigma = 5;
    CV_f_smooth = smoothdata( abs([CV_f(end-sigma+1:end), CV_f, CV_f(1:sigma)]),'gaussian',sigma);
    [~,loc] = min(CV_f_smooth(sigma+1:end-sigma));

% rotate un-smoothed CV to reorient to make the propagation direction at the center of array
    CV_oriented = circshift(abs(CV_f),-loc-L/4);

% average around the index value in the range of -15 to 15 degrees
    w=floor(15/(360/L));
    AvgSpeed = abs(mean(CV_oriented((L/4-w):(L/4+w))./cosd( [-w:+w]*360/L)  ));
    StdSpeed = (std(CV_oriented((L/4-w):(L/4+w))./cosd( [-w:+w]*360/L)  ));
    
% find angle from smoothed data
if CV_f(loc)<0
    Angle =  loc*(360/(L)) + 180 ; % [degrees]
else
    Angle = loc*(360/(L)); % [degrees]
end

end
    