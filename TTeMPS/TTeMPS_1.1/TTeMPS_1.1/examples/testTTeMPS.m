% last modified on 08/08/2022
% lusong@lsec.cc.ac.cn
clc
clear
close all
load ../../../../data_missing
%X = X(1:1000,:,:,:,:,:);
[Omega0,X_Omega0] = find(X);
sizeOmega0 = size(Omega0,1);
%% test set
Gamma = Omega0;
X_Gamma = X_Omega0;
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
maxrank = 3;%7
opts_cg = struct('maxiter', 10,'maxiter_final',10, 'tol', 1e-6, ...
    'reltol', 1e-6, 'gradtol', 0, 'maxrank', maxrank,'epsilon',1e-8);
[X,cost,test,stats] = completion_rankincrease( 'GeomCG', X_Omega, ...
    Omega, X_Omega_C, Omega_C, X_Gamma, Gamma, X0, opts_cg );
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
