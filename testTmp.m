% Last modifed on 08/10/202
% lusong@lsec.cc.ac.cn
clc
clear
close all
%%
load ./data574.mat
Phits = tensor(Phi);
Rcp = 12;
Phik = cp_als(Phits,Rcp);
ndim = ndims(Phik);
ncpt = ncomponents(Phik);
sz  = size(Phi);
Y_true = Phi;
X = Phi;
%%
missingRate = 0.75;
creat_missing;
Wd = double(W);
X = Y_true.*W;
%checkRate = nnz(X)/prod(sz);
%%
Rcp = 10;
options.factr = 1e7;
options.maxIts = 10;
tic
[Y,~,output] = cp_wopt(X, W, Rcp, 'opt_options',options);
cputime = toc;


