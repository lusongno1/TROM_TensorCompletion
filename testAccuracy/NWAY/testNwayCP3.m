%% Last modified on 10/11/2022
%  lusong@lsec.cc.ac.cn
clc
clear
close all
load('../../data274D5333T100');
clearvars -except Phi
Y = Phi;
%%
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
X0 = T;%tensor(T);
%%
X0(X0==0)=NaN;
%%
%X = Phi;
RCP = 30;
Options(1) = 1e-3;
%Options(4) = 2;
%const = [0 0 0];
tic
model1 = parafac(X0,RCP,Options);
s1 = toc
Xkrec = ktensor(model1);
Xm = double(Xkrec);
cal_acc(double(Y),Xm)