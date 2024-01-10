close all;
clear all;

D = load('aas.mat');

XLimInf = floor(D.Timestamps(end)/2);
XLimSup = XLimInf + 1;

[~,idxInf]= min(abs(D.Timestamps - XLimInf));
[~,idxSup]= min(abs(D.Timestamps - XLimSup));

a = min(D.Data(1, idxInf:idxSup));
b = max(D.Data(1, idxInf:idxSup));

YLimInf = a - a*0.1;
YLimSup = b + b*0.1;

clear