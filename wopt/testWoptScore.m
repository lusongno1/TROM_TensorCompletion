clc
clear
close all
R = 2;
info = create_problem('Size', [3 5], 'Num_Factors', R, ...
    'M', 0.5, 'Noise', 0.00);
X = info.Data;
P = info.Pattern;
M_true= info.Soln;
M_init = create_guess('Data', X, 'Num_Factors', R, ...
    'Factor_Generator', 'nvecs');
[M,~,output] = cp_wopt(X, P, R, 'init', M_init);
exitmsg = output.ExitMsg
scr = score(M,M_true)

M = normalize(M);
M_true = normalize(M_true);
Err = [];
for k=1:R
scr2 = 1;
%score = penalty * (a1'*b1) * (a2'*b2) * ... * (aR'*bR),
%penalty = 1 - abs(lambda_a - lambda_b) / max(lamdba_a, lambda_b).
for i=1:2
    ai = M.U{i}(:,k);
    bi = M_true.U{i}(:,k);
    scr2 = scr2*ai.'*bi;
end
lambda_a = M.lambda(k);
lambda_b = M_true.lambda(k);
penalty = 1 - abs(lambda_a - lambda_b) / max(lambda_a, lambda_b);
scr2 = scr2*penalty
Err = [Err scr2];
end
scr22 = norm(Err,2)