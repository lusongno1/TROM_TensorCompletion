% Interpolatory tensor ROM testing
%==========================================================================
%
% Heat equation model with 3 holes and 2 parameters: 
% Biot number on left boundary Omega_out
% and the temperature on the holes' boundaries
%
% Alexander Mamonov, University of Houston, 2022
%==========================================================================

%clear all; 
%clear model;

clc
clear
close all
fignum = 100;

% ttbpath = '../tensor_toolbox-v3.2.1/';
% tttpath = '../TT-Toolbox-2.3/';
% 
% addpath(ttbpath);
% savepath = pwd;
% cd(tttpath); setup; cd(savepath);

mstr = 'itrom3hole2par'; % method name string

%==========================================================================
% flags

isrerunsnap = true; % rerun snapshot generation
%isrerunsnap = false;

isrerunoffline = true; % rerun the offline part
%isrerunoffline = false;

isrerunonline = true; % rerun the online part
%isrerunonline = false;

% flags to run each decomposition
%iscp    = true; 
iscp    = false; 

ishosvd = true; 
%ishosvd = false;

%istt    = true; 
istt    = false;

ispod   = true;

iscartsample = true; % Cartesian sample
%iscartsample = false;

%==========================================================================
% parameters

% parameters: Biot numbers and ranges
Biomin = 0.01;
Biomax = 0.5;
Bii = 0.5; % fixed

% parameters: external temperature and ranges
uo = 1; % fixed
uimin = 0;
uimax = 0.9;

if iscartsample
    D = 2; % dimension of parameter space 
    ni = zeros(D, 1); % sampling size in each parameter direction     
    ni(1) = 5; 
    ni(2) = 5;
else
    error('Not implemented');
end

K = prod(ni); % total number of parameter samples

% time interval
tmax = 20; 
N = 100; % number of time steps
dt = tmax / N;

% POD rank
nPOD = 20; 

% tensor decomposition accuracy
epsa = 1e-6;

p = 2; % number of polynomial interpolation points

% CP parameters
epsCP  = epsa; % a priori accuracy
nCP = nPOD; % number of reduced basis vectors for CP-TROM
RCP = 100; %250; % target canonical rank
maxitCP = 200; % maximum number of iterations for CP decomposition 

% Tucker (HOSVD) parameters
epsHOSVD  = epsa; % a priori decomposition accuracy
nHOSVD = nPOD; % number of reduced basis vectors for HOSVD-TROM
%ranksHOSVD = [4, 4, 4]; % a priori ranks 

% Tensor train parameters
epsTT  = epsa; % a priori decomposition accuracy
nTT = nPOD; % number of reduced basis vectors for TT-TROM
%ranksTT = [4, 4, 4]; % a priori ranks 

rngseed = 0; % seed for random parameter sampling (if chosen)

%==========================================================================
% parameter sampling

alphahat = zeros(K, D);
aind = zeros(K, D);

% parameter samples in each direction
alspace{1} = linspace(Biomin, Biomax, ni(1));
alspace{2} = linspace(uimin, uimax, ni(2));

if iscartsample
    % Cartesian parameter sampling
    k = 0;
    
    for i1 = 1:ni(1),  for i2 = 1:ni(2) 
        k = k + 1;
        alphahat(k, :) = [alspace{1}(i1), alspace{2}(i2)];  
        aind(k, :) = [i1, i2];
    end; end
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
generateMesh(model, 'Hmax', 0.2);

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
   
    fprintf('Generating parameter sample snapshots...'); tic;

    if iscartsample
        for k = 1:K
            Phi(:, aind(k, 1), aind(k, 2), :) = ...
                solve3hole(model, uo, ...
                alphahat(k, 2), alphahat(k, 2), alphahat(k, 2), ...
                alphahat(k, 1), Bii, tlist);
        end
    end
    
    fprintf('done in %f sec\n', toc);
end

%==========================================================================
% Offline stage: compute tensor ROMs

if isrerunoffline

fprintf('Offline stage:\n'); tic;
    
Phim = reshape(Phi, M, []); % big matrix of all snapshots
Phit = tensor(Phi); % tensor object from multi-array

% POD: SVD of the big snapshot matrix
[UPhi, SPhi, ~] = svd(Phim, 'econ'); 

% CP decomposition
if iscp
    [onlineCP, UCP, PhiCPk] = tromcpoffline(Phi, RCP, maxitCP);
    
    % project FEM matrices onto UCP-basis for use at the online stage
    % only mass and stiffness matrix projections can be pre-computed, 
    % since Q and G depend on parameters that enter via
    % the boundary conditions
    FEMCP = [];
    FEMCP.M = UCP' * FEM.M * UCP;
    FEMCP.K = UCP' * FEM.K * UCP;
end

% Tucker decomposition
if ishosvd
    tol = 1e-6;
    [onlineHOSVD, U, PhiHOSVD] = tromhosvdoffline(Phi, tol);
    [Uc,Sc] = tromhosvdonline(onlineHOSVD, ea);
    projectionAll;
end

% Tensor train
if istt
    tol = 1e-6;
    [onlineTT, U, PhiTT] = tromttoffline(Phi, tol);
    [Uc,Sc] = tromttonline(onlineTT, ea);
    projectionAll;
end

fprintf('Offline stage completed in %f sec\n', toc);
    
end

%==========================================================================
% online stage

if isrerunonline

% random test parameters from the intervals 
alpha = [(Biomax - Biomin)*rand(1) + Biomin, ...
         (uimax - uimin)*rand(1) + uimin];

fprintf('Test parameter values: '); 
fprintf(' %g', alpha); fprintf('\n');

% assemble FEM matrices for test parameters alpha
FEMa = fem3hole(model, uo, alpha(2), alpha(2), alpha(2), alpha(1), Bii);

% true snapshots for test parameters alpha
Phia = tstepfemcn0qr(FEMa, tlist);

% interpolation coefficients
ea = cell(D, 1);
for i = 1:D
    ea{i} = lagrangeinterp(alpha(i), alspace{i}, p);
end

% local reduced basis coordinate computation

% CP decomposition
if iscp
    [UcCP] = tromcponline(onlineCP, ea);
    
    % only the first nCP columns of UcCP are used as 
    % the coordinates of the local basis
    UcCPn = UcCP(:, 1:nCP);
    
    % online ROM FEM matrices
    FEMCPc = [];
    
    % pre-computed in the offline stage, parameter-independent
    FEMCPc.M = UcCPn' * FEMCP.M * UcCPn;
    FEMCPc.K = UcCPn' * FEMCP.K * UcCPn;

    % parameter-dependent, need to project on UCP
    FEMCPc.Q = UcCPn' * (UCP' * FEMa.Q * UCP) * UcCPn;
    FEMCPc.G = UcCPn' * (UCP' * FEMa.G);
    
    % ROM snapshot computation
    PhiCPc = UCP * (UcCPn * tstepfemcn0qr(FEMCPc, tlist));
    
    errPhiCP = norm(PhiCPc - Phia) / norm(Phia);
    fprintf('CP ROM snapshot relative error: %e\n', errPhiCP);
end

% Tucker decomposition
% if ishosvd
% 
% end

% Tensor train
% if istt
% 
% end

end

return;
