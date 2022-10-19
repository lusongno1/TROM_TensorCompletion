%% for random missing
clc
clear
close all
global Y_true Tind Tind0 W;
%%
load('../../data274D5333T100');
%clearvars -except Phi
Y_true = Phi;
missingRate = 0.85;
X0 = Phi;
II = size(X0);
N  = prod(II);
missing_rate = missingRate;
idd = (randperm(N) > N*missing_rate);
Q   = reshape(idd,II);
T   = zeros(II);
T(Q)= X0(Q);
X0 = tensor(T);
W = Q;%double(Q);
%%
sz = size(Y_true);
[Tind,Tval] = find(X0);% non zero positions and values
[Tind0,Tval0] = find(X0==0);% zero positions and values
nnz(X0)./prod(sz)
%%
R = 100;

%%
if 0
tic
%Y = cp_tcals(X0,R,Tind,Tval,'tilde_epsilon',0.001,'maxiters',inf,'tol',0,'monitor',true);
Y = cp_tcals(X0,R,Tind,Tval,'tilde_epsilon',0.001,'maxiters',26,'tol',0);

toc
cal_acc_avail_std(double(Y),double(Y_true),W)
cal_acc(double(Y),double(Y_true))
end
%%
if 1
TVQV    = 'qv';        % 'tv' or 'qv' ;
rho     = zeros(1,6);
rh = 0.2;
rho(1) = 0.0; %not quit sensitive, set to zero
rho(2) = rh; %not sensitive
rho(3) = rh; %sensitive to this parameter. set to zero
rho(4) = rh; %not quit sensitive to this parameter. set to zero
rho(5) = rh; %sensitive, set to zero
rho(6) = rh; %not sensive, set to zere
K       = 10;          % Number of components which are updated in one iteration.
tilde_epsilon = 0.000001; 
SNR = -log10((tilde_epsilon^2))*10;
nu      = 0.001;%0.2;        % threshold for R <-- R + 1.
maxR = R;
maxiter = inf;       % maximum number of iteration
tol     = 0;%1e-15;        % tolerance
out_im  = 0;           % you can monitor the process of 'image' completion if out == 1.
sall = tic;
[Z G U histo histo_R] = SPC_Sparse(T,Q,TVQV,rho,K,SNR,nu,maxiter,tol,out_im,maxR,1);
ssall = toc(sall)
cal_acc_avail_std(double(Z),double(Y_true),W)
cal_acc(double(Z),double(Y_true))
end