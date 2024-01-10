function DATA = Binning(DATA,binSize,precision)
% with bining, tha value taken at the end should fit in 16-bit range
% data is assumed take with EM-HB1800s camera in 12 bits so with 4x4
% binning there is room for 16-bit range but in larger bining there is not

ScaleFactor = binSize^2/16;
if binSize<=4
     temp = zeros(size(DATA,1)/binSize,size(DATA,2)/binSize,size(DATA,3),'uint16');
else
     temp = zeros(size(DATA,1)/binSize,size(DATA,2)/binSize,size(DATA,3),'single');
end

    for i=1:size(DATA,3)
        A = DATA(:,:,i);
        C = sum(reshape(A,binSize,[]));
        C = reshape(C,size(A,1) / binSize,[])';
        C = sum(reshape(C,binSize,[]));
        C = reshape(C,size(A,2) / binSize,[])';
        temp(:,:,i) = C;
    end
    DATA = uint16(temp/ScaleFactor);