%% Last modified on 10/11/2022
%  lusong@lsec.cc.ac.cn
clc
clear
close all
load('../../data274D5333T100');
clearvars -except Phi
Y = Phi;
%%
X = Phi;
RCP = 10;
Options(1) = 1e-2;
%Options(4) = 2;
%const = [0 0 0];
tic
model1 = parafac2(X,RCP,Options);
s1 = toc
Xkrec = ktensor(model1);
Xm = double(Xkrec);
cal_acc(double(Y),Xm)
% %% for missing
% missingRate = 0.85;
% sz = size(Y);
% make_missing;
% X = double(Y.*W);
% sum(X)./prod(sz)
% X(X==0)=NaN;
% sum(isnan(X),'all')./prod(sz)
% tic
% model2 = parafac(X,RCP,Options);
% s2 = toc
% Xkrec = ktensor(model2);
% Xm = double(Xkrec);
% cal_acc(double(Y),Xm)
% %%
% II = size(X);
% N  = prod(II);
% idd = (randperm(N) > N*missingRate);
% Q   = reshape(idd,II);
% T   = zeros(II);
% T(Q)= X(Q);
% checkRate = nnz(T)/N;
% T(T==0)=NaN;
% sum(isnan(T),'all')./prod(sz)
% tic
% model3 = parafac(T,RCP,Options);
% s2 = toc
% Xkrec = ktensor(model3);
% Xm = double(Xkrec);
% cal_acc(double(Y),Xm)
