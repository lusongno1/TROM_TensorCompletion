clc
clear
close all
global Y_true Tind Tind0
%%
load('../../data274D5333T100');
%load('../../data574D9555');
%load('../../data274D8444T100');
%Phi = double(Phi_tucker);
clearvars -except Phi
%Phi = [1 2; 3 4];
Y_true = Phi;%tensor(Phi);
%%
sz = size(Y_true);
missingRate = 0.85;
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
tic
Y = cp_tcals(X0,R,Tind,Tval);
toc
cal_acc(double(Y),double(Y_true))
%cal_acc(double(Y),double(X0))

