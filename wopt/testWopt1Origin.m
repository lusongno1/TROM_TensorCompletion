clc
clear
close all
R = 2;
info = create_problem('Size', [15 10 5], 'Num_Factors', R, ...
    'M', 0.25, 'Noise', 0.10);
X = info.Data;
P = info.Pattern;
M_true= info.Soln;
M_init = create_guess('Data', X, 'Num_Factors', R, ...
    'Factor_Generator', 'nvecs');
tt_fac_to_vec(M_init,M_true)
[M,~,output] = cp_wopt(X, P, R, 'init', M_init);
exitmsg = output.ExitMsg
scr = score(M,M_true)