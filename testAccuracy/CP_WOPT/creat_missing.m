sizeParas = sz(2:end-1);
%% not a good choice, since it is not real random
% info = create_problem('Lambda_Generator', @ones,'Size', sizeParas,'Num_Factors', 1, ...
%     'M', missingRate, 'Noise', 0.00);
% Ws = info.Pattern;
%% real random? not important
avail_inds_flat = randperm(prod(sizeParas),round(prod(sizeParas)*(1-missingRate)));
Ws = tenzeros(sizeParas);
Ws = double(Ws);
Ws(avail_inds_flat) = 1;
Ws = tensor(Ws);
%%
%1-nnz(Ws<1e-1)/prod(sizeParas)
U = tenones(sz(1),1);
V = tenones(sz(end),1);
W  = ttt(U,Ws);
W = ttt(W,V);
W = squeeze(W);
nnzW = nnz(W);
nprod = prod(sz);
checkRate = nnzW/nprod;