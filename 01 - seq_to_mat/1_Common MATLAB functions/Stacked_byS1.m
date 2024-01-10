function Stacked = Stacked_byS1(DATA,S1,pacing)
% syntax, Stacked_byS1(DATA,S1,pacing)
% DATA, time vector or series of images for ensable averaging  (staking)m 
% S1 - 1D vector of time (sample) point when each next period (beat) begins
% PCL - pacing cycle length in number of samples

if ndims(DATA) == 2;
    Stacked = zeros(1,pacing,'single');
    for i = 1:length(S1)
    Stacked = Stacked + single(DATA( S1(i):S1(i)+pacing-1))';
    end
end


if ndims(DATA) == 3
    Stacked = zeros(size(DATA,1),size(DATA,2),pacing,'single');
    for i = 1:length(S1)
        Stacked(:,:,1:pacing) = Stacked(:,:,1:pacing) + single(DATA(:,:,S1(i):S1(i)+pacing-1));
    end
end

Stacked = Stacked/(length(S1)-1);