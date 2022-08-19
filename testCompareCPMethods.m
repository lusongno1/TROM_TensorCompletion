% Last modified on 8/16/2022
clc
clear
close all
load ./data574.mat
addpath(genpath('./wopt')) 
Phits = tensor(Phi);
sz  = size(Phi);
Y_true = Phi;
%Acc_TIME;
maxItersSet = [2]
Acc_TIME = ones(4,length(maxItersSet));
ii = 0;
rng('default')
rng(0);
for maxIters = maxItersSet
    ii = ii+1;
    %% set parameters
    Rcp = 10;
    %maxIters = 5;
    missingRate = 0.75;
    %% creat missing
    creat_missing;
    Wd = double(W);
    X = Y_true.*W;
    %% cp_wopt
    options.factr = 1e7;
    options.maxIts = maxIters;
    tic
    [Y,~,output] = cp_wopt(X, W, Rcp, 'opt_options',options);
    cpu_time_wopt = toc;
    Yd = double(Y);
    acc_wopt = cal_acc(Y_true,Yd);
    %% nway toolbox
    Options(6) = maxIters;
    %Options(2) = 10;
    XX = double(X);
    XX(Wd==0)=NaN;
    RCP = Rcp;
    tic
    model1 = parafac(XX,RCP,Options);
    cpu_time_nway = toc;
    Xkrec = ktensor(model1);
    Xm = double(Xkrec);
    acc_nway = cal_acc(Y_true,Xm);
    %% BCPF
    
    %%
    acc_wopt
    acc_nway
    cpu_time_wopt
    cpu_time_nway
    Acc_TIME(1:4,ii) = [acc_wopt acc_nway cpu_time_wopt cpu_time_nway]
end
%% SPC



%%





