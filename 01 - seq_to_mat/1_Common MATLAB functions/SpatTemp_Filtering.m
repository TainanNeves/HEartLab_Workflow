function DATA = SpatTemp_Filtering(DATA,S,T,mode)
% syntax SpatTemp_Filtering(DATA,S,T,'mode')
% S,T specify the size of the filter kernel matices. 
% For example S = 5 creates 5x5 matrix that for each pixel in the image 
% filters (convolutes the filter kernel) of the data in the radius of 2 pixels (5-1)/2 = 2. Use always odd numbers. 
% T means the total width -> the range is +-T/2  
%Examples:
% SpatTemp_Filtering(DATA,3,5,'CPU')
% SpatTemp_Filtering(DATA,3,11,'CPU')
% Try first to filter DATA with increasing the T size, but if T >15 doesn't
% smooth it enoug then increase the S to 5. 
% The upper limits should be S = 9, T =21, 
%
% mode,  'GPU' or 'CPU' processing 

sigma=1;

DATA = single(DATA);
if T>1
            if T==3, sigma = 1; end
            if T==5, sigma = 1.5; end
            if T==7, sigma = 1.5; end
            if T==9, sigma = 2; end
    
        if mode == 'GPU' 
            
        H = gpuArray( ndgaussian(T,sigma));
        for i=1:size(DATA,1)
            X = gpuArray( single(squeeze(DATA(i,:,:)) ) );
            X2 = X;
        for j=1:size(X,1)
            X2(j,:) = conv( X(j,:),H,'same');
        end
        DATA(i,:,:) = gather(X2);
        end

        else
        
        H = single(ndgaussian(T,sigma));
        for i=1:size(DATA,1)
            for j=1:size(DATA,2)
                X = single(squeeze(DATA(i,j,:)));
                DATA(i,j,:) = conv( X,H,'same');
            end
        end
    end

end

if S>1
            if S==3, sigma = 1; end
            if S==5, sigma = 1.5; end
            if S==7, sigma = 1.5; end
            if S==9, sigma = 2; end
    H = ndgaussian([S,S],sigma);
    for i=1:size(DATA,3)
        DATA(:,:,i) = conv2(DATA(:,:,i),H,'same');
    end
end


    
