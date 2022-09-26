
%ttbpath = '../tensor_toolbox-v3.2.1/';
%tttpath = '../TT-Toolbox-2.3/';

%addpath(ttbpath);
%savepath = pwd;
%cd(tttpath); setup; cd(savepath);

rngseed = 0;
rng(rngseed);

Phi = rand(3, 4, 5);
R = 2;

[onlineCP, UCP, PhiCPk] = tromcpoffline(Phi, R, 100, 1e-6);

ea{1} = ones(4, 1);

Uc = tromcponline(onlineCP, ea);

