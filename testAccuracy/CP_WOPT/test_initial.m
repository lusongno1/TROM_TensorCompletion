% Last modifed on 08/10/202
% lusong@lsec.cc.ac.cn
clc
clear
close all
global Y_true Yd Wd Err rng_s
rng_seed = 5;
rng(rng_seed);
Err = [];
Y_true = [];
Yd = [];
Wd = [];
%%
load ../../data574.mat
Phits = tensor(Phi);
Rcp = 12;
Phik = cp_als(Phits,Rcp);
ndim = ndims(Phik);
ncpt = ncomponents(Phik);
sz  = size(Phi);
Y_true = Phi;
X = Phi;
%%
missingRate = 0.5;
creat_missing;
Wd = double(W);
X = Y_true.*W;
%checkRate = nnz(X)/prod(sz);
%%
Rcp = 12;
options.factr = 1e7;
options.maxIts = 30;
tic
[Y,~,output] = cp_wopt(X, W, Rcp, 'opt_options',options);
cputime = toc;
%% calculate accuracy
Yd = double(Y);
errW = cal_acc_avail(Y_true,Yd,Wd)
err = cal_acc(Y_true,Yd)
%%
close all
Err2 = Err;
Err2(:,1) = sqrt(2*Err(:,1))./sqrt(sum((Wd.*Y_true).^2,'all'));
plot(Err2(:,1),'Marker','.');
hold on;
plot(Err2(:,2),'Marker','.');
xlabel('Iteraion');
ylabel('Accuracy');
title(['Missing:' num2str(missingRate) ',  ' 'randn seed:' num2str(rng_seed)]);
[~,ind] = min(Err2(:,2));
text_point(ind,Err2(ind,2))
legend('$\hat \varepsilon$','$\varepsilon$','minimum','Interpreter','LaTex')

