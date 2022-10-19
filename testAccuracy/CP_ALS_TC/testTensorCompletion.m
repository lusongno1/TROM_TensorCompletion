clc
clear
close all
%%
global Y_true ERR W
Y_true = [];
ERR = [];

%%
data_path = '../../';
load([data_path 'data274D4222T100.mat']);
Phi = double(Phi_tucker);
sz  = size(Phi);
N = prod(sz);
Y_true = Phi;
clear Phi;
%%
missingRate = 0.0;
make_missing;
%%
Q = logical(double(W));
X0 = Y_true.*W;
clearvars -except X0 Y_true;
[Tind,Tval] = find(X0);% non zero positions and values
%%
R = 100;
tic
Y = cp_tcals(X0,R,Tind,Tval);
toc
cal_acc(double(Y),double(Y_true))