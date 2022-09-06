clc
clear
close all
%%
load ../../data.mat
T = Phi;
%T = ones(3,3,3,3,3,3);
Y_true = T;
%% creat missing
missingRate = 0.8;
II = size(Y_true);
N  = prod(II);
idd = (randperm(N) > N*missingRate);
W   = reshape(idd,II);
X0   = zeros(II);
X0(W)= Y_true(W);
checkRate = nnz(X0)/N;
%%
T = X0;
Q = W;
%[Xest, histo, histoR, G, U, S, D] = MDT_Tucker_incR(T,Q,tau,param)
[Xest, histo, histoR, G, U, S, D] = MDT_Tucker_incR(T,Q,ones(1,6)*1)
