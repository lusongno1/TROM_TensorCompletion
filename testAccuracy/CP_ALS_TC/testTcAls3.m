%% for random missing
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

missingRate = 0.85;
X0 = Phi;
%%% make missing entry
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
if(size(Tind0)~=0)
    X0init = mean(Tval);
    X0(Tind0) = X0init;
end
%%
%tic
%X0 = sptensor(Tind,Tval,sz);
%toc
R = 100;
tic
%Y = cp_als(X0,R,'maxiters',30);
Y = cp_tcals(X0,R,Tind,Tval,0.01,'tol',1e-12,'maxiters',50);
toc
cal_acc(double(Y),double(Y_true))
%cal_acc(double(Y),double(X0))