clc
clear
close all
%% load Tensor Decomposition Data
load ../../TC_RES274D5333T100.mat
load ../../data274D5333T100.mat
Phi = Z;
RCP = size(U{1},2);
PhiCPk = ktensor(G.',U);
isTC = 1;
%%
iscp = true;
M = size(U{1},1);
p = 2; % number of polynomial interpolation points
nCP = 20; % number of reduced basis vectors for CP-TROM
% Biomax = 0.5;
% Biomin = 0.01;
% uimax = 0.9;
% uimin = 0;
% uo = 1; % fixed
% Bii = 0.5; % fixed
% tlist = linspace(0, tmax, N+1);
%% offline stage
Phim = reshape(Phi, M, []); % big matrix of all snapshots
Phit = tensor(Phi); % tensor object from multi-array

% POD: SVD of the big snapshot matrix
%[UPhi, SPhi, ~] = svd(Phim, 'econ'); 

% CP decomposition
if iscp
if isTC    
    [onlineCP, UCP] = tromcpofflineTC(Phi, PhiCPk, RCP);
else
    maxitCP = 200;  
    [onlineCP, UCP, PhiCPk] = tromcpoffline(Phi, RCP, maxitCP);  
end
    
    % project FEM matrices onto UCP-basis for use at the online stage
    % only mass and stiffness matrix projections can be pre-computed, 
    % since Q and G depend on parameters that enter via
    % the boundary conditions
    FEMCP = [];
    FEMCP.M = UCP' * FEM.M * UCP;
    FEMCP.K = UCP' * FEM.K * UCP;
end
%% oneline stage
% random test parameters from the intervals 
alpha = [(Biomax - Biomin)*rand(1) + Biomin, ...
         (uimax - uimin)*rand(1) + uimin...
         (Biomax - Biomin)*rand(1) + Biomin, ...
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
