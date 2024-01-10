function h=ndgaussian(varargin,sigma)
% the last argument is specified is sigma
% use sigma=1 for 3x3, sigma 1.5 for 5x5 spat filtering
% Creates gaussian 1D or 2D map, 3D map needs work
% varargin can be 1 or 2 numbers specifying dimensions
% Examples: general expression is ndgaussian(m,n) outputing 1D or 2D
% gaussian map of a final dimesion m,n
m=0; % This is for ndgrid creation so that in the case of gaussian dim <3 proper mesh is constructed
n=0;
p=0;
sigma_m=0.5; % so that in gaussian curve calculation factors in the case dim<3 enter as multiplication by 1
sigma_n=0.5;
sigma_p=0.5;

dim=length(varargin); % 

switch dim
    
    case 1
        m=varargin;
        [x1] = ndgrid( -(m-1)/2:(m-1)/2);
        sigma_m= sigma;
        h=exp(-(x1.^2)/2/sigma_m^2);
    
    case 2
        m=varargin(1);
        n=varargin(2);
        [x1,x2] = ndgrid( -(m-1)/2:(m-1)/2, -(n-1)/2:(n-1)/2);
        sigma_m= sigma;
        sigma_n= sigma;
        h=exp(-(x1.^2)/2/sigma_m^2).*exp(-(x2.^2)/2/sigma_n^2);
    
    case 3
        m=varargin{1};
        n=varargin{2};
        p=varargin{3};
        [x1,x2,x3] = ndgrid( (-(m-1)/2):((m-1)/2), -(n-1)/2:(n-1)/2, -(p-1)/2:(p-1)/2);
        sigma_m=m/(2*sqrt(2*log(2)));
        sigma_n=n/(2*sqrt(2*log(2)));
        sigma_p=p/(2*sqrt(2*log(2)));
        h=exp(-(x1.^2)/2/sigma_m^2).*exp(-(x2.^2)/2/sigma_n^2);
end

h=h/sum(sum(sum(h))); % so that h is normalized