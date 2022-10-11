clc
clear
close all
%%
load('../../data274D5333T100.mat');
clearvars -except Phi
Y_true = Phi;
%%
missingRate = 0.85;
sz = size(Phi);
make_missing;
X = Y_true.*W;
Rcp = 10;
options.factr = 1e7;
options.maxIts = 300;
tic
[Y,~,output] = cp_wopt(X, W, Rcp,'opt_options',options);
%[Y,~,output] = cp_wopt(X, W, Rcp);
s1 = toc
err = cal_acc(Y_true,double(Y))
%%
tic
%[Y,~,output] = cp_wopt(X, W, Rcp,'opt_options',options);
Yt = tensor(Y_true);
options2.maxIts = 300;
[Y,~,output] = cp_opt(Yt, Rcp,'opt_options',options2);
s2 = toc
err = cal_acc(Y_true,double(Y))