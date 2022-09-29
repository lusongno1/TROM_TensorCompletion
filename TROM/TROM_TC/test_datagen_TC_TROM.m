clc
clear
close all
%%
data_path = '../../';
addpath(genpath(data_path));
addpath(genpath([data_path 'testAccuracy']));
mesh_h = 1;
DD = [5 3 3 3];
TT = 100;
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
    ni(1) = DD(1); 
    ni(2) = DD(2);
    ni(3) = DD(3); 
    ni(4) = DD(4);
else
    error('Not implemented');
end

K = prod(ni); % total number of parameter samples

% time interval
tmax = 20; 
N = TT; % number of time steps
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
generateMesh(model, 'Hmax', mesh_h);

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
%%
%%
% Last modifed on 09/07/2022
% lusong@lsec.cc.ac.cn
%clc
%clear
%close all
rng(0);
%%
global Y_true ERR
Y_true = [];
ERR = [];

%%
%data_path = '../../';
load([data_path 'data' path_str '.mat']);
%load ../../data274.mat
%load ../../data274D8444T100.mat
Phi = double(Phi_tucker);
%load ../../data3072D9555.mat
%Phits = tensor(Phi);
sz  = size(Phi);
N = prod(sz);
Y_true = Phi;
clear Phi;
%%
missingRate = 0.85;
%nnzW = nnz(W);
%checkRate = nnzW/N;
make_missing;
%make_missing_diag;
%make_missing_alldim;

%%
Q = logical(double(W));
T = double(Y_true.*W);
clear W;
%%
addpath Function_SPC
addpath plotting_function
%% hyperparameters and run SPC-TV

TVQV    = 'qv';        % 'tv' or 'qv' ;
%rho     = [0.1,0.01,0.01,0.01,0.01,0.01]; % smoothness (0.1 - 1.0) for 'qv' and (0.01 - 0.5) for 'tv' is recommended.
%rho     = [0.001,0.001,0.001,0.001,0.001,0.001];
rho     = zeros(1,6);
rh = 0.2;
rho(1) = 0.0; %not quit sensitive, set to zero
rho(2) = rh; %not sensitive
rho(3) = rh; %sensitive to this parameter. set to zero
rho(4) = rh; %not quit sensitive to this parameter. set to zero
rho(5) = rh; %sensitive, set to zero
rho(6) = rh; %not sensive, set to zere
K       = 10;          % Number of components which are updated in one iteration.
%SNR     = 50;          % error bound
tilde_epsilon = 0.0001; 
SNR = -log10((tilde_epsilon^2))*10;
nu      = 0.2;%0.2;        % threshold for R <-- R + 1.

maxR = 100;
maxiter = inf;       % maximum number of iteration
tol     = 0;%1e-15;        % tolerance
out_im  = 0;           % you can monitor the process of 'image' completion if out == 1.
[Xtv Z G U histo histo_R] = SPC(T,Q,TVQV,rho,K,SNR,nu,maxiter,tol,out_im,maxR);
%iter:  objective :: epsilon :: conv. speed :: number of components 

%% Err
err = cal_acc(Y_true,Xtv)
%err = cal_acc_avail_std(Y_true,Xtv,Q)
%checkRate = nnz(Q)/N;
ERR(:,4) = sqrt(histo./sum((Q.*Y_true).^2,'all'));
close all;
%%
Accs = ERR(:,1);
AccsAvail = ERR(:,4);
Ranks = ERR(:,2);
Time = ERR(:,3);
iters = 1:size(ERR,1);
%%
figure;
subplot(2,2,1)
yyaxis left 
%plot(iters,Accs,'Marker','.');
semilogy(iters,Accs);
hold on;
semilogy(iters,AccsAvail);
ylabel('Accuracy')
yyaxis right 
plot(iters,Ranks);
ylabel('Ranks')
legend('Accuracy $\varepsilon$','Accuracy $\tilde \varepsilon$','Ranks','Interpreter','LaTex','Location','best');
xlabel('Iterations')
title(['SPC method: Missing = ' num2str(missingRate)])

%%
%figure;
subplot(2,2,2)
yyaxis left 
plot(iters,Time,'Marker','.');
ylabel('Time Costs(s)')
yyaxis right 
plot(iters,Ranks);
ylabel('Ranks')
legend('Time Per Iter','Ranks','Interpreter','LaTex','Location','best');
xlabel('Iterations')
title(['SPC method: Missing = ' num2str(missingRate)])


%%
%figure
subplot(2,2,3)
[Runq,Inds] = unique(Ranks);
%ed = Inds(end);
Inds = Inds-1;
Inds(1) = [];
Inds(end+1) = length(Ranks);
%Inds = Inds(180:end);
accs = ERR(Inds,1);
accsAvail = ERR(Inds,4);
ranks = ERR(Inds,2);
time = ERR(Inds,3);
yyaxis left 
%plot(ranks,accs,'Marker','.');
semilogy(ranks,accs);
hold on;
semilogy(ranks,accsAvail);
ylabel('Accuracy')
yyaxis right 
plot(ranks,time,'Marker','.');
ylabel('Time Costs(s)')
legend('Accuracy $\varepsilon$','Accuracy  $\tilde \varepsilon$','Time Cost Per Iter(s)',...
    'Interpreter','LaTex','Location','best');
xlabel('Ranks')
title(['SPC method: Missing = ' num2str(missingRate)])

%%
subplot(2,2,4)
%close;
title('Parameters used and result show!!!');
% tensorsize;
% TVQV    = 'qv';        % 'tv' or 'qv' ;
% rho     = [0.01 0.01 0.01 0.01 0.01 0.01]; % smoothness (0.1 - 1.0) for 'qv' and (0.01 - 0.5) for 'tv' is recommended.
% K       = 5;          % Number of components which are updated in one iteration.
% SNR     = 50;          % error bound
% nu      = 0.01;        % threshold for R <-- R + 1.
% maxiter = 20;       % maximum number of iteration
% tol     = 1e-7;        % tolerance
text(0.1,0.9,['TVQV = ' TVQV],'Interpreter','LaTex');
text(0.1,0.85,['$\rho =$ [' num2str(rho) ']'],'Interpreter','LaTex');
text(0.1,0.8,['K = ' num2str(K)],'Interpreter','LaTex');
text(0.1,0.75,['SNR = ' num2str(SNR)],'Interpreter','LaTex');
text(0.1,0.7,['$\nu =$ ' num2str(nu)],'Interpreter','LaTex');
text(0.1,0.65,['maxiter = ' num2str(maxiter)],'Interpreter','LaTex');
text(0.1,0.6,['total iters = ' num2str(size(ERR,1))],'Interpreter','LaTex');
text(0.1,0.55,['tol = ' num2str(tol)],'Interpreter','LaTex');

text(0.1,0.5,['input data size = [' num2str(sz) ']'],'Interpreter','LaTex');
text(0.1,0.45,['missing rate = ' num2str(missingRate)],'Interpreter','LaTex');
text(0.1,0.4,['min $\varepsilon$ = ' num2str(min(ERR(:,1)))],'Interpreter','LaTex');
text(0.1,0.35,['min $\tilde \varepsilon$ = ' num2str(min(ERR(:,4)))],'Interpreter','LaTex');
text(0.1,0.3,['max Rank = ' num2str(max(ERR(:,2)))],'Interpreter','LaTex');
%sss = toc;
text(0.1,0.25,['Total time cost = ' num2str(sum(ERR(:,3))) 's'],'Interpreter','LaTex');

%%
saveas(gcf,'result.png')
%save('vars.mat')
%data_path = '../../';
save([data_path 'TC_RES' num2str(sz(1)) 'D' num2str(sz(2)) num2str(sz(3)) num2str(sz(4))...
    num2str(sz(5)) 'T' num2str(sz(6)) '.mat'],'G','U')



%%
%clc
%clear
%close all
%% load Tensor Decomposition Data
%data_path = '../../';
load( [data_path 'TC_RES' path_str '.mat']);
load( [data_path 'data' path_str '.mat']);
Phi = double(Phi_tucker);
%Phi = Z;
RCP = size(U{1},2);
PhiCPk = ktensor(G.',U);
isTC = 0; %change this one to be SPC or ALS
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
         (uimax - uimin)*rand(1) + uimin, ...
         (uimax - uimin)*rand(1) + uimin];
alpha = [0.25 0.45 0.45 0.45];

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






