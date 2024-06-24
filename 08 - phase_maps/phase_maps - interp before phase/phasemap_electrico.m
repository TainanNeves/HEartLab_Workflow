function [D]=phasemap_electrico(A)

  if length(size(A))==2
              [nx,ny]=size(A);
              D=zeros(size(A));
             for i=1:nx
                 disp(i);
                 
                    s=A(i,:);s=s-mean(s);
                    ha=imag(hilbert(s));
                    D(i,:)=-atan2(ha,s);
             
             end
  end      

% Hilbert Transformation:
% Hilbert transform is a mathematical operation that takes 
% a real-valued signal and produces a complex-valued signal.
% The code select just the imaginary part to work

% In summary:
% The code utilizes the Hilbert transform to 
% obtain the analytic representation of each row in the 
% input array and then applies atan2 to extract the 
% phase information, which is stored in the output array.
