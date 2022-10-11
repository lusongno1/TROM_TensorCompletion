clc
clear
close all
%%
%load('../../data274D5333T100');
load('../../data574D9555');
clearvars -except Phi
Y_true = Phi;
%%
sz = size(Y_true);
missingRate = 0.5;
make_missing;
X0 = Y_true;
%W = [1 0;0 1];
X0 = X0.*W;
X0 = tensor(X0);
[Tind,Tval] = find(X0);% non zero positions and values
[Tind0,Tval0] = find(X0==0);% non zero positions and values
nnz(X0)./prod(sz)
if(size(Tind0)~=0)
    X0init = mean(Tval);
    X0(Tind0) = X0init;
end
%%
R = 100;
for i=1:10
    Y = cp_als(X0,R);
    cal_acc(double(Y),double(Y_true))
    cal_acc_avail(double(Y),double(Y_true),double(W))
    X0 = tensor(Y);
    X0(Tind) = Tval;
end