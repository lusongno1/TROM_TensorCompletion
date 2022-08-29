% last modified on 08/22/2022
% lusong@lsec.cc.ac.cn
clc
clear
close all
load ../../../../data574.mat
Y_true = Phi;
%%
sz  = size(Phi);
missingRate = 0.75;
creat_missing;
Wd = double(W);
X = Y_true.*W;
checkRate = nnz(X)/prod(sz);
%%

[Omega0,X_Omega0] = find(X);%size(X_Omega0,1)./prod(sz)
sizeOmega0 = size(Omega0,1);
%% test set
[Gamma,X_Gamma] = find(tensor(Y_true));
%checkRate = nnz(X_Gamma)/prod(sz);
%% training set and control set
Omega = Omega0;
X_Omega = X_Omega0;
sizeOmega_C = round(sizeOmega0*0.0001);
Omega_C_ind = randperm( sizeOmega0, sizeOmega_C );
Omega_C = Omega( Omega_C_ind, : );
Omega( Omega_C_ind, : ) = [];
X_Omega_C = X_Omega0(Omega_C_ind,:);
X_Omega( Omega_C_ind, : ) = [];
%% completion
d = ndims(X);
n = size(X);
r = [1, 1*ones(1,d-1), 1];
X0 = TTeMPS_rand( r, n );
X0 = orthogonalize( X0, X0.order );
maxrank = 10;%7
opts_cg = struct('maxiter', 3,'maxiter_final',3, 'tol', 1e-6, ...
    'reltol', 1e-6, 'gradtol', 0, 'maxrank', maxrank,'epsilon',1e-8);
[X,cost,test,stats] = completion_rankincrease( 'GeomCG', X_Omega, ...
    Omega, X_Omega_C, Omega_C, X_Gamma, Gamma, X0, opts_cg );
%checkRate = 1951600/prod(sz);
%% plot
stats.rankidx = cumsum(stats.rankidx)
subplot(1,2,1)
semilogy( 1:length(cost), cost,'Markersize',8);
hold on
line = [1e-6,1e0];
for i=1:length(stats.rankidx)
    semilogy( [stats.rankidx(i), stats.rankidx(i)], line, '--','color',[0.7,0.7,0.7]);
end
title('Reduction of cost function')
xlabel('Number of individual RTTC iterations performed')
ylabel('Cost function')
legend('Cost function','rank increase in one mode')
set(gca,'fontsize',16)

subplot(1,2,2)
semilogy( 1:length(test), test,'Markersize',8);
title('Reduction of rel. error on test set')
xlabel('Number of full RTTC runs')
ylabel('Rel. error after one RTTC run for a certain rank')
set(gca,'fontsize',16)

%%
Yd  = full(X);
err = cal_acc(Y_true,Yd)
err2 = cal_acc_avail_std(Y_true,Yd,Wd)

