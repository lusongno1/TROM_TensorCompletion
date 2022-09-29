%==========================================================================
% Snapshot generation for the
% heat equation model with 3 holes and 4 parameters: 
% Biot number on left boundary Omega_out
% and temperatures on the holes' boundaries
%
% Alexander Mamonov, University of Houston, 2022
%==========================================================================

%clear all; 
%clear model;
clc
clear
close all

fignum = 100;

%ttbpath = '../tensor_toolbox-v3.2.1/';
%tttpath = '../TT-Toolbox-2.3/';

%addpath(ttbpath);
%savepath = pwd;
%cd(tttpath); setup; cd(savepath);

mstr = 'gensnap3hole4par'; % method name string

%==========================================================================
% flags

isrerunsnap = true; % rerun snapshot generation
%isrerunsnap = false;

iscartsample = true; % Cartesian sampling
%iscartsample = false;

%==========================================================================
% parameters

% parameters: Biot numbers and ranges
Biomin = 0.01;
Biomax = 0.5;
Bii = 0.5; % fixed

% parameters: external temperatures and ranges
uo = 1; % fixed
uimin = 0;
uimax = 0.9;

if iscartsample
    D = 4; % dimension of parameter space 
    ni = zeros(D, 1); % sampling size in each parameter direction     
    ni(1) = 8; 
    ni(2) = 4;
    ni(3) = 4; 
    ni(4) = 4;
else
    error('Not implemented');
end

K = prod(ni); % total number of parameter samples

% time interval
tmax = 20; 
N = 100; % number of time steps
dt = tmax / N;

%==========================================================================
% parameter sampling

alphahat = zeros(K, D);
aind = zeros(K, D);

% parameter samples in each direction
alspace{1} = linspace(Biomin, Biomax, ni(1));
alspace{2} = linspace(uimin, uimax, ni(2));
alspace{3} = linspace(uimin, uimax, ni(3));
alspace{4} = linspace(uimin, uimax, ni(4));

if iscartsample
    % Cartesian parameter sampling
    k = 0;
    
    for i1 = 1:ni(1), for i2 = 1:ni(2)
    for i3 = 1:ni(3), for i4 = 1:ni(4)
        k = k + 1;

        alphahat(k, :) = [alspace{1}(i1), alspace{2}(i2), ...
                          alspace{3}(i3), alspace{4}(i4)];  
        aind(k, :) = [i1, i2, i3, i4];
    end; end; end; end
end

%==========================================================================
% Geometry and mesh
% generate the rectangular domain with three square holes

model = createpde;

% outer rectangle
R0 = [3, 4, 0, 10, 10, 0, 4, 4, 0, 0]'; 
% three holes
R1 = [3, 4, 1, 3, 3, 1, 3, 3, 1, 1]';
R2 = [3, 4, 4, 6, 6, 4, 3, 3, 1, 1]';
R3 = [3, 4, 7, 9, 9, 7, 3, 3, 1, 1]';

gm = [R0, R1, R2, R3];
sf = 'R0 - R1 - R2 - R3';
ns = char('R0', 'R1', 'R2', 'R3')';
dl = decsg(gm, sf, ns);

geometryFromEdges(model, dl);
generateMesh(model, 'Hmax', 1);

% set up the PDE
% Matlab PDE: m(∂2u/∂t2)+d(∂u/∂t)−∇·(c∇u)+au=f
% Heat equation (∂u/∂t)−∇·(∇u) = 0
specifyCoefficients(model, 'm', 0, 'd', 1, 'c', 1, 'a', 0, 'f', 0);

% zero initial condition
setInitialConditions(model, 0);

% times at which to evaluate the solution
tlist = linspace(0, tmax, N+1);

%==========================================================================
% generate the snapshots corresponding to all parameter samples

if isrerunsnap

    % assemble FEM matrices to extract the number 
    % M of spatial degrees of freedom
    FEM = fem3hole(model, uo, uimin, uimin, uimin, Biomin, Bii);
    
    M = size(FEM.K, 1);

    % multi-array of solution snapshots
    Phi = zeros([M, ni', N]);
   
    fprintf('Generating parameter sample snapshots...\n'); tic;

    if iscartsample
        for k = 1:K
            if mod(k, 25) == 0
                fprintf('Snapshot %d out of %d\n', k, K);
            end

            Phi(:, aind(k, 1), aind(k, 2), aind(k, 3), aind(k, 4), :) = ...
                solve3hole(model, uo, ...
                alphahat(k, 2), alphahat(k, 3), alphahat(k, 4), ...
                alphahat(k, 1), Bii, tlist);
        end
    end 
    fprintf('done in %f sec\n', toc);
end
Phi_t = tensor(Phi);
Phi_tucker = hosvd(Phi_t,1e-15);
relerr = norm(Phi_t-full(Phi_tucker))/norm(Phi_t);
sz = size(Phi)
clear Phi Phi_t
data_path = '../../';
save([data_path 'data' num2str(sz(1)) 'D' num2str(sz(2)) num2str(sz(3)) num2str(sz(4))...
    num2str(sz(5)) 'T' num2str(sz(6)) '.mat'])
