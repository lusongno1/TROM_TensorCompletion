%%
% Last modifed on 09/07/2022
% lusong@lsec.cc.ac.cn
clc
clear
close all
rng(0);
%%
global Y_true ERR
Y_true = [];
ERR = [];
%%
load ../../data574.mat
%load ../../data574D9555.mat
%load ../../data3072D9555.mat
Phits = tensor(Phi);
sz  = size(Phi);
N = prod(sz);
Y_true = Phi;
%%
missingRate = 0.85;
%nnzW = nnz(W);
%checkRate = nnzW/N;
make_missing;
%make_missing_diag;
%make_missing_alldim;

%%
Q = logical(double(W));
T = double(Y_true.*W);
%%
addpath Function_SPC
addpath plotting_function
%% hyperparameters and run SPC-TV

TVQV    = 'qv';        % 'tv' or 'qv' ;
%rho     = [0.1,0.01,0.01,0.01,0.01,0.01]; % smoothness (0.1 - 1.0) for 'qv' and (0.01 - 0.5) for 'tv' is recommended.
%rho     = [0.001,0.001,0.001,0.001,0.001,0.001];
rho     = zeros(1,6);
rh = 0.5;
rho(1) = 0.0; %not quit sensitive, set to zero
rho(2) = rh; %not sensitive
rho(3) = rh; %sensitive to this parameter. set to zero
rho(4) = rh; %not quit sensitive to this parameter. set to zero
rho(5) = rh; %sensitive, set to zero
rho(6) = rh; %not sensive, set to zere
K       = 10;          % Number of components which are updated in one iteration.
%SNR     = 50;          % error bound
tilde_epsilon = 0.0001; 
SNR = -log10((tilde_epsilon^2))*10;
nu      = 0.2;%0.2;        % threshold for R <-- R + 1.

maxR = 10;
maxiter = inf;       % maximum number of iteration
tol     = 0;%1e-15;        % tolerance
out_im  = 0;           % you can monitor the process of 'image' completion if out == 1.
[Xtv Z G U histo histo_R] = SPC(T,Q,TVQV,rho,K,SNR,nu,maxiter,tol,out_im,maxR);
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
subplot(2,2,1)
yyaxis left 
%plot(iters,Accs,'Marker','.');
semilogy(iters,Accs);
hold on;
semilogy(iters,AccsAvail);
ylabel('Accuracy')
yyaxis right 
plot(iters,Ranks);
ylabel('Ranks')
legend('Accuracy $\varepsilon$','Accuracy $\tilde \varepsilon$','Ranks','Interpreter','LaTex','Location','best');
xlabel('Iterations')
title(['SPC method: Missing = ' num2str(missingRate)])

%%
%figure;
subplot(2,2,2)
yyaxis left 
plot(iters,Time,'Marker','.');
ylabel('Time Costs(s)')
yyaxis right 
plot(iters,Ranks);
ylabel('Ranks')
legend('Time Per Iter','Ranks','Interpreter','LaTex','Location','best');
xlabel('Iterations')
title(['SPC method: Missing = ' num2str(missingRate)])


%%
%figure
subplot(2,2,3)
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
semilogy(ranks,accs);
hold on;
semilogy(ranks,accsAvail);
ylabel('Accuracy')
yyaxis right 
plot(ranks,time,'Marker','.');
ylabel('Time Costs(s)')
legend('Accuracy $\varepsilon$','Accuracy  $\tilde \varepsilon$','Time Cost Per Iter(s)',...
    'Interpreter','LaTex','Location','best');
xlabel('Ranks')
title(['SPC method: Missing = ' num2str(missingRate)])

%%
subplot(2,2,4)
%close;
title('Parameters used and result show!!!');
% tensorsize;
% TVQV    = 'qv';        % 'tv' or 'qv' ;
% rho     = [0.01 0.01 0.01 0.01 0.01 0.01]; % smoothness (0.1 - 1.0) for 'qv' and (0.01 - 0.5) for 'tv' is recommended.
% K       = 5;          % Number of components which are updated in one iteration.
% SNR     = 50;          % error bound
% nu      = 0.01;        % threshold for R <-- R + 1.
% maxiter = 20;       % maximum number of iteration
% tol     = 1e-7;        % tolerance
text(0.1,0.9,['TVQV = ' TVQV],'Interpreter','LaTex');
text(0.1,0.85,['$\rho =$ [' num2str(rho) ']'],'Interpreter','LaTex');
text(0.1,0.8,['K = ' num2str(K)],'Interpreter','LaTex');
text(0.1,0.75,['SNR = ' num2str(SNR)],'Interpreter','LaTex');
text(0.1,0.7,['$\nu =$ ' num2str(nu)],'Interpreter','LaTex');
text(0.1,0.65,['maxiter = ' num2str(maxiter)],'Interpreter','LaTex');
text(0.1,0.6,['total iters = ' num2str(size(ERR,1))],'Interpreter','LaTex');
text(0.1,0.55,['tol = ' num2str(tol)],'Interpreter','LaTex');

text(0.1,0.5,['input data size = [' num2str(sz) ']'],'Interpreter','LaTex');
text(0.1,0.45,['missing rate = ' num2str(missingRate)],'Interpreter','LaTex');
text(0.1,0.4,['min $\varepsilon$ = ' num2str(min(ERR(:,1)))],'Interpreter','LaTex');
text(0.1,0.35,['min $\tilde \varepsilon$ = ' num2str(min(ERR(:,4)))],'Interpreter','LaTex');
text(0.1,0.3,['max Rank = ' num2str(max(ERR(:,2)))],'Interpreter','LaTex');
%sss = toc;
text(0.1,0.25,['Total time cost = ' num2str(sum(ERR(:,3))) 's'],'Interpreter','LaTex');

%%
saveas(gcf,'result.png')
save('vars.mat')




