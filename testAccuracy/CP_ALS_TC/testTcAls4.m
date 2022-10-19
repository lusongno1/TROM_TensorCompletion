%% for random missing
clc
clear
close all
global Y_true Tind Tind0
%%
load('../../data274D5333T100');
clearvars -except Phi
Y_true = Phi;
missingRate = 0.85;
X0 = Phi;
II = size(X0);
N  = prod(II);
missing_rate = missingRate;
idd = (randperm(N) > N*missing_rate);
Q   = reshape(idd,II);
T   = zeros(II);
T(Q)= X0(Q);
X0 = tensor(T);
%%
sz = size(Y_true);
[Tind,Tval] = find(X0);% non zero positions and values
[Tind0,Tval0] = find(X0==0);% zero positions and values
nnz(X0)./prod(sz)
%%
%X0 = sptensor(Tind,Tval,sz);%dont't use sparse tensor
%%
R = 100;
tic
%Y = cp_tcals2(X0,R,Tind,Tval,0.01,'tol',1e-12,'maxiters',10);
Y = cp_tcals2(X0,R,Tind,'maxiters',10);
%Y = cp_tcals(X0,R,Tind,Tval,'maxiters',10);
toc
cal_acc(double(Y),double(Y_true))