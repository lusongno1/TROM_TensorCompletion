%% test cp_als %cp_wopt cp_als and spc time
clc
clear
close all
%%
load('../../data274D5333T100.mat');
clearvars -except Phi
R = 100;
Xt = tensor(Phi);
Xk = cp_als(Xt,100);
Xk([1,1,1,1,1,1])
%% cp_als
tic
Y = cp_als(Xk,R);
toc
cal_acc(double(Y),Phi)