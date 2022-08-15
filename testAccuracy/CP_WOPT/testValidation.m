% Last modifed on 08/13/202
% lusong@lsec.cc.ac.cn
clc
clear
close all
global Y_true Yd Wd Err rng_s  Wdval Wdmiss
rng('default')
rng_seed = 0;
rng(rng_seed);
Err = [];
Y_true = [];
Yd = [];
Wd = [];
Wdval = [];
Wdmiss = [];
%%
load ../../data574.mat
Phits = tensor(Phi);
Rcp = 12;
Phik = cp_als(Phits,Rcp);
ndim = ndims(Phik);
ncpt = ncomponents(Phik);
sz  = size(Phi);
Y_true = double(Phi);
X = Phi;
%%
missingRate = 0.6;
creat_missing;
%% divide available into training set and validation set
valid_rate = 0.1;
%rng(2);
split_avail_data;
%% left part
% Wmiss = W;
% Wmiss = (Wmiss==0);
%rng(3);
get_missing_inds;
%%
checkRateVal = nnz(Wval)/prod(sz);
checkRateTr = nnz(Wtr)/prod(sz);
checkRateAvail = nnz(W)/prod(sz);
checkRateMiss = nnz(Wmiss)/prod(sz);
%%
W = Wtr;
Wd = double(W);
Wdval = double(Wval);
Wdmiss = double(Wmiss);
X = Y_true.*W;
%checkRate = nnz(X)/prod(sz);
%%
Rcp = 12;
options.factr = 1e7;
options.maxIts = 10;
tic
[Y,~,output] = cp_wopt(X, W, Rcp, 'opt_options',options);
cputime = toc;
%%
close all
Err2 = Err;
Err2(:,1) = sqrt(2*Err(:,1))./sqrt(sum((Wd.*Y_true).^2,'all'));
Err2(:,3) = sqrt(2*Err(:,3))./sqrt(sum((Wdval.*Y_true).^2,'all'));
Err2(:,4) = sqrt(2*Err(:,4))./sqrt(sum((Wdmiss.*Y_true).^2,'all'));
plot(Err2(:,1),'Marker','.');
hold on;
plot(Err2(:,2),'Marker','.');
plot(Err2(:,3),'Marker','.');
plot(Err2(:,4),'Marker','.');
xlabel('Iteraion');
ylabel('Accuracy');
title(['Missing:' num2str(missingRate) ',  '  'Validation of avail:' num2str(valid_rate)]);
[~,ind] = min(Err2(:,2));
text_point(ind,Err2(ind,2))
legend('$\hat \varepsilon_{tr}$','$\varepsilon$','$\hat \varepsilon_{val}$','$\tilde \varepsilon$','minimum','Interpreter','LaTex')

