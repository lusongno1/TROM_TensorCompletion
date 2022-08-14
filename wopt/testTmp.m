clc
clear
close all
R = 2;
info = create_problem('Lambda_Generator', @ones,'Size', [2 3], 'Num_Factors', R, ...
    'M', 0.5, 'Noise', 0.00);
% info = create_problem('Size', [2 3], 'Num_Factors', R, ...
%     'M', 0.5, 'Noise', 0.00);
X = info.Data;
P = info.Pattern;
M_true= info.Soln;
% M_init = create_guess('Data', X, 'Num_Factors', R, ...
%     'Factor_Generator', 'nvecs');
M_init = create_guess('Soln',M_true, 'Factor_Generator', 'pertubation','Pertubation',0.00);
tmp = tt_fac_to_vec(M_init)
A = ktensor(tt_cp_vec_to_fac(tmp,M_true))
full(A)
full(M_true)

[M,~,output] = cp_wopt(X, P, R, 'init', M_init);
exitmsg = output.ExitMsg
scr = score(M,M_true)
full(M)
full(M_true)
tmp = P .* (X - M) 