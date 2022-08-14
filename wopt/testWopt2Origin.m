clc
clear
close all
R = 2;
info = create_problem('Size', [150 100 50], 'Num_Factors', R, ...
    'M', 0.95, 'Sparse_M', true, 'Noise', 0.00);
X = info.Data;
P = info.Pattern;
M_true= info.Soln;
M_init = create_guess('Data', X, 'Num_Factors', R, ...
    'Factor_Generator', 'nvecs');
[M,~,output] = cp_wopt(X, P, R, 'init', M_init);
exitmsg = output.ExitMsg
scr = score(M,M_true)
