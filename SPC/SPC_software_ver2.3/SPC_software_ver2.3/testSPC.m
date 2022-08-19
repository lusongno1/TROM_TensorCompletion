%%
% Last modifed on 08/19/2022
% lusong@lsec.cc.ac.cn
clc
clear
close all
%%
global Y_true
load ../../../data574.mat
Phits = tensor(Phi);
Rcp = 12;
sz  = size(Phi);
Y_true = Phi;
X = Phi;
%%
missingRate = 0.85;
creat_missing;
Wd = double(W);
X = Y_true.*W;
%checkRate = nnz(X)/prod(sz);
%%
addpath Function_SPC
addpath plotting_function

%%% make synthetic data
N = 50;

X0 = Phi;

%%% make missing entry
II = size(X0);
N  = prod(II);
missing_rate = missingRate;

idd = (randperm(N) > N*missing_rate);
Q   = reshape(idd,II);
T   = zeros(II);
T(Q)= X0(Q);

%% hyperparameters and run SPC-TV

TVQV    = 'tv';        % 'tv' or 'qv' ;
rho     = [0.01 0.01 0.01 0.01 0.01 0.01]; % smoothness (0.1 - 1.0) for 'qv' and (0.01 - 0.5) for 'tv' is recommended.
K       = 10;          % Number of components which are updated in one iteration.
SNR     = 50;          % error bound
nu      = 0.01;        % threshold for R <-- R + 1.
maxiter = 300;       % maximum number of iteration
tol     = 1e-7;        % tolerance
out_im  = 0;           % you can monitor the process of 'image' completion if out == 1.

[Xtv Z G U histo histo_R] = SPC(T,Q,TVQV,rho,K,SNR,nu,maxiter,tol,out_im);


%% hyperparameters and run SPC-QV

% TVQV    = 'qv';        % 'tv' or 'qv' ;
% rho     = [1.0 1.0 1.0 1.0 1.0 1.0]; % smoothness (0.1 - 1.0) for 'qv' and (0.01 - 0.5) for 'tv' is recommended.
% K       = Rcp;          % Number of components which are updated in one iteration.
% SNR     = 50;          % error bound
% nu      = 0.01;        % threshold for R <-- R + 1.
% maxiter = 10;       % maximum number of iteration
% tol     = 1e-7;        % tolerance
% out_im  = 0;           % you can monitor the process of 'image' completion if out == 1.
% 
% [Xqv Z G U histo histo_R] = SPC(T,Q,TVQV,rho,K,SNR,nu,maxiter,tol,out_im);

err = cal_acc(Y_true,Xtv)



