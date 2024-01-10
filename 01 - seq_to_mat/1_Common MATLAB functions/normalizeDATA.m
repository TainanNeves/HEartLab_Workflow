function DATA = normalizeDATA(DATA)

tf = isa(DATA,'single');
if ~tf
    tf = isa(DATA,'double');
    if ~tf
        display('input variable has to be in single or double precision')
        return
    end
end

R = size(DATA,1);
C = size(DATA,2);
T = size(DATA,3);

DATA = (reshape(DATA, [R*C, T] ))';

DATA2 = DATA;
for i = 1:10
    clear DATA3
    [Max(i,:), Loc] = max ( DATA2 );
    for j=1:size(DATA,2)
        X = DATA2(:,j); X(Loc(j)) = [];
        DATA3(:,j) = X;
    end
    DATA2 = DATA3;
end

DATA2 = DATA;
for i = 1:10
    clear DATA3
    [Min(i,:), Loc] = min ( DATA2 );
    for j=1:size(DATA,2)
        X = DATA2(:,j); X(Loc(j)) = [];
        DATA3(:,j) = X;
    end
    DATA2 = DATA3;
end

Max = mean(Max,1);
Min = mean(Min,1);

Min = min ( DATA );
DATA = (DATA - Min)./ (Max - Min);
DATA = (reshape(DATA', [R, C, T] ));

%taking care of NaNs
DATA(isnan(DATA)) = 0;