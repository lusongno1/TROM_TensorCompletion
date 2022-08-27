%%
% Last modifed on 08/19/2022
% lusong@lsec.cc.ac.cn
clc
clear
close all
%%
global Y_true ERR
Y_true = [];
ERR = [];
load ../../data574.mat
Phits = tensor(Phi);
sz  = size(Phi);
Y_true = Phi;
%%
missingRate = 0.85;
%creat_missing;
%Wd = double(W);
%X = Y_true.*W;
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

checkRate = nnz(T)/N;

%% hyperparameters and run SPC-TV

TVQV    = 'qv';        % 'tv' or 'qv' ;
rho     = [0.01 0.01 0.01 0.01 0.01 0.01]; % smoothness (0.1 - 1.0) for 'qv' and (0.01 - 0.5) for 'tv' is recommended.
K       = 5;          % Number of components which are updated in one iteration.
SNR     = 50;          % error bound
nu      = 0.01;        % threshold for R <-- R + 1.
maxiter = 50;       % maximum number of iteration
tol     = 1e-7;        % tolerance

out_im  = 0;           % you can monitor the process of 'image' completion if out == 1.
tic
[Xtv Z G U histo histo_R] = SPC(T,Q,TVQV,rho,K,SNR,nu,maxiter,tol,out_im);
%iter:  objective :: epsilon :: conv. speed :: number of components 

%% Err
err = cal_acc(Y_true,Xtv)
%err = cal_acc_avail_std(Y_true,Xtv,Q)
%checkRate = nnz(Q)/N;
ERR(:,4) = sqrt(histo./sum((Q.*Y_true).^2,'all'));
close all;
%%
Accs = ERR(:,1);
AccsAvail = ERR(:,4);
Ranks = ERR(:,2);
Time = ERR(:,3);
iters = 1:size(ERR,1);
%%
figure;
yyaxis left 
%plot(iters,Accs,'Marker','.');
semilogy(iters,Accs);
hold on;
semilogy(iters,AccsAvail);
ylabel('Accuracy')
yyaxis right 
plot(iters,Ranks);
ylabel('Ranks')
legend('Accuracy $\varepsilon$','Accuracy $\tilde \varepsilon$','Ranks','Interpreter','LaTex');
xlabel('Iterations')
title(['SPC method: Missing = ' num2str(missingRate)])

%%
figure;
yyaxis left 
plot(iters,Time,'Marker','.');
ylabel('Time Costs(s)')
yyaxis right 
plot(iters,Ranks);
ylabel('Ranks')
legend('Time Per Iter','Ranks','Interpreter','LaTex');
xlabel('Iterations')
title(['SPC method: Missing = ' num2str(missingRate)])


%%
figure
[Runq,Inds] = unique(Ranks);
%ed = Inds(end);
Inds = Inds-1;
Inds(1) = [];
Inds(end+1) = length(Ranks);
%Inds = Inds(180:end);
accs = ERR(Inds,1);
accsAvail = ERR(Inds,4);
ranks = ERR(Inds,2);
time = ERR(Inds,3);
yyaxis left 
%plot(ranks,accs,'Marker','.');
semilogy(ranks,accs,'Marker','.');
hold on;
semilogy(ranks,accsAvail,'Marker','.');
ylabel('Accuracy')
yyaxis right 
plot(ranks,time,'Marker','.');
ylabel('Time Costs(s)')
legend('Accuracy $\varepsilon$','Accuracy  $\tilde \varepsilon$','Time Cost Per Iter(s)','Interpreter','LaTex');
xlabel('Ranks')
title(['SPC method: Missing = ' num2str(missingRate)])





