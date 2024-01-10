function [R]=BASE(DATA)
[a b c]=size(DATA);
for i=1:a
    for j=1:b
        x=squeeze(DATA(i,j,:));
        R(i,j,:)=x-(min(x));
    end 
end