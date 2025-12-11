function [D]=phasemap_opticos(A,dib)
  if length(size(A))==3
              [nx,ny,nt]=size(A);
              D=zeros(size(A));
             for i=1:nx
                 disp(i);
                 for j=1:ny
                    s=squeeze(A(i,j,:));
                    s=s-mean(s);
                    ha=imag(hilbert(s));
                    D(i,j,:)=-atan2(ha,s);
                 end
             end
  end    
end
  