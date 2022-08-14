%% last version modified on 06/08/2022
%  lusong@lsec.cc.ac.cn
clc
clear
close all
%% pre-process of data
load ../../data574.mat
Phits = tensor(Phi);
Rcp = 12;
Phik = cp_als(Phits,Rcp);%low rank cp approximation
ndim = ndims(Phik);
ncpt = ncomponents(Phik);
sz  = size(Phi);
%% creat indicator tensor W by create_problem
sizeParas = sz(2:end-1);
missingRate = 0.85;
info = create_problem('Lambda_Generator', @ones,'Size', sizeParas,'Num_Factors', 1, ...
    'M', missingRate, 'Noise', 0.00);
Ws = info.Pattern;
U = tenones(sz(1),1);
V = tenones(sz(end),1);
W  = ttt(U,Ws);
W = ttt(W,V);
W = squeeze(W);
nnzW = nnz(W);
nprod = prod(sz);
checkRate = nnzW/nprod;
%% make incomplete tensor with W
%M = Phits;
Y_true = Phik;
X = Y_true.*W;
%checkRate = nnz(X)/prod(sz);
%% choose initial guess Y0
% pertubation of true soluiton
Y0 = create_guess('Soln',Y_true, 'Factor_Generator', 'pertubation','Pertubation',0.1);
Y0{1} = Y0{1}.*repmat(Y_true.lambda.',size(Y0{1},1),1); %fix a bug of tensor toolbox
% tmp = tt_fac_to_vec(Y0)
% Y0ts = ktensor(tt_cp_vec_to_fac(tmp,Y_true))
% score(Y_true,Y0ts)

% generate initial value by HOSVD
% Y0 = create_guess('Data', X, 'Num_Factors', Rcp, ...
%     'Factor_Generator', 'nvecs');%should set Rcp<=3

% generate initial value randomly
% Y0 = create_guess('Data', X, 'Num_Factors', Rcp, ...
%     'Factor_Generator', 'rand');
%% recover Phi by cp_wopt
%tol = 1e7^2*nnzW*0.5;
tol = 1e9; % set this value to st   op early tol*epsm
options.factr = tol;
options.maxIts = 10;
[Y,~,output] = cp_wopt(X, W, Rcp, 'init', Y0,'opt_options',options);
%[Y,~,output] = cp_wopt(X, W, Rcp, 'init', Y0);
%% score
exitmsg = output.ExitMsg
scr1 = score(Y,Y_true)
scr2 = norm(full(Y_true)-full(Y))/sqrt(nprod)
scr3 = norm(full(Y_true.*W)-full(Y.*W))/sqrt(nnzW)
viz(Y)
viz(Y_true)


