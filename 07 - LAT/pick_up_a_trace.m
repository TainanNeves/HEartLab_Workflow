function [x,y,Trace1, Trace2] = pick_up_a_trace(I,DATA,N)

button =0;
   
while button ~=32
    
   
    figure(99),imagesc(I),colormap('gray'), axis square
    [y,x button]=ginput(1);
    x=round(x);y=round(y);

    if button == 3 
        I(x,y) = max(max(I));
        imagesc(I);
    end
    if x<=0 || y<=0, break, end
    X = squeeze(mean(mean(DATA(x-N:x+N, y-N:y+N, :),1),2));
    %X2 = squeeze(mean(mean(DATA2(x-N:x+N, y-N:y+N, :),1),2));
    Trace1 = X;
    %Trace2 = X2;
    figure(103),plot(X); title(strcat('row: ', num2str(x), ' col: ', num2str(y)))
    %figure(104),plot(X2); title(strcat('row: ', num2str(x), ' col: ', num2str(y)))
    display((max(X)-min(X))*100);
end

