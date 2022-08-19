clc
clear
close all;
randn('state',1); rand('state',1); %#ok<RAND>
loadpath;
%%
load ../data574.mat

%% Generate a low-rank tensor
DIM = size(Phi);     % Dimensions of data
R = 8;                % True CP rank
DataType = 1;         % 1: Random factors   2: The deterministic factors (sin, cos, square)
% Generate tensor by factor matrices
X = Phi;

%% Random missing values
ObsRatio = 0.2;               % Observation rate: [0 ~ 1]
Omega = randperm(prod(DIM)); 
Omega = Omega(1:round(ObsRatio*prod(DIM)));
O = zeros(DIM); 
O(Omega) = 1;


%% Generate observation tensor Y
Y = O.*X;

%% Run BayesCP
fprintf('------Bayesian CP Factorization---------- \n');
ts = tic;
if ObsRatio~=1 
    % Bayes CP algorithm for incomplete tensor and tensor completion    
    [model] = BCPF_TC(Y, 'obs', O, 'init', 'ml', 'maxRank', max([DIM 2*R]), 'dimRed', 1, 'tol', 1e-6, 'maxiters', 100, 'verbose', 2);
else
    % Bayes CP algorithm for fully observed tensor 
    [model] = BCPF(Y, 'init', 'ml', 'maxRank', max([DIM 2*R]), 'dimRed', 1, 'tol', 1e-6, 'maxiters', 200, 'verbose', 2);
end
t_total = toc(ts);

% Performance evaluation
X_hat = double(model.X);
err = X_hat(:) - X(:);
rmse = sqrt(mean(err.^2));
rrse = sqrt(sum(err.^2)/sum(X(:).^2));

% Report results
fprintf('\n------------Bayesian CP Factorization-----------------------------------------------------------------------------------\n')
fprintf('Observation ratio = %g, SNR = %g, True Rank=%d\n', ObsRatio, SNR, R);
fprintf('RRSE = %g, RMSE = %g, Estimated Rank = %d, \nEstimated SNR = %g, Time = %g\n', ...
    rrse, rmse, model.TrueRank, model.SNR, t_total);
fprintf('--------------------------------------------------------------------------------------------------------------------------\n')

%% Visualization of data and results
plotYXS(Y, X_hat);
factorCorr = plotFactor(Z,model.X.U);



