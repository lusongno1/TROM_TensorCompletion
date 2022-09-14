clc
clear
close all
%%
global Y_true Err
Y_true = [];
Err = [];
load ../../data.mat
T = Phi;
%T = ones(3,3,3,3,3,3);
Y_true = T;
%% creat missing
sz = size(T);
N = prod(sz);
missingRate = 0.8;
make_missing;
% II = size(Y_true);
% N  = prod(II);
% idd = (randperm(N) > N*missingRate);
% W   = reshape(idd,II);
% X0   = zeros(II);
X0 = Y_true.*W;
checkRate = nnz(X0)/N;
%%
T = double(X0);
Q = double(W);
%[Xest, histo, histoR, G, U, S, D] = MDT_Tucker_incR(T,Q,tau,param)
[Xest, histo, histoR, G, U, S, D] = MDT_Tucker_incR(T,Q,[1 3 3 3 3 3]);
%Err = sqrt(sum((Xest-Y_true).^2,'all')/sum(Xest.^2,'all'))
Accs = Err(:,1);
AccsAvail = Err(:,2);
iters = 1:size(Accs,1);
semilogy(iters,Accs,'Marker','.');
hold on;
semilogy(iters,AccsAvail,'Marker','.');
ylabel('Accuracy')
legend('Accuracy $\varepsilon$','Accuracy $\tilde \varepsilon$','Interpreter','LaTex','Location','best');
xlabel('Iterations')
title(['MDT method: Missing = ' num2str(missingRate)])
