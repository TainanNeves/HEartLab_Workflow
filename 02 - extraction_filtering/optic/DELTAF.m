function [R]=DELTAF(DATA)
[a b c]=size(DATA);
for i=1:a
    for j=1:b
        p=double(squeeze(DATA(i,j,:)));
        m=mean(p);
        DF=((p - m)./m)*100;
        R(i,j,:)=DF;
    end
end