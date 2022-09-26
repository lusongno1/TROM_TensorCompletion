function [onlineCP, UCP, PhiCPk] = tromcpoffline(Phi, R, maxit, tol)
%TROMCPOFFLINE Offline stage of CP-TROM computation
% Input:  
%  Phi = snapshot multi-array
%  R = target canonical rank
%  maxit = maximum number of cp_als iterations
%  tol = cp_als tolerance
% Output:
%  PhiCPk = CP approximation to Phi in Kruskal tensor format
%  onlineCP.lambda = Kruskal tensor weights of PhiCPk
%  onlineCP.sigma  = cell array of \sigma_i, each having R columns
%  onlineCP.RU = upper triangular factor in thin QR factorization of \hat{U}
%  onlineCP.RV = upper triangular factor in thin QR factorization of \hat{V}
%  UCP = MxR orthogonal matrix with columns forming and orthonormal basis 
%        for the universal space
%
% Alexander Mamonov, University of Houston, 2022
%==========================================================================

onlineCP = [];

sizePhi = size(Phi);
onlineCP.R  = R;
onlineCP.D  = length(sizePhi) - 2;
onlineCP.M  = sizePhi(1);
onlineCP.N  = sizePhi(onlineCP.D + 2);
onlineCP.ni = sizePhi(2:(onlineCP.D + 1));

Phit = tensor(Phi); % tensor object from Phi multi-array

% default CP algorithm parameter values
if nargin < 3, maxit = 100; end
if nargin < 4, tol   = 1e-4; end

[PhiCPk, ~, cpout] = cp_als(Phit, R, 'maxiters', maxit, 'tol', tol);

errcp = sqrt( norm(Phit)^2 + norm(PhiCPk)^2 ...
            - 2 * innerprod(Phit, PhiCPk) ) / norm(Phit);

fprintf('CP-TROM offline stage:\n');
fprintf('cp_als converged in %d iterations\n', cpout.iters);
fprintf('CP decomposition relative error %e\n', errcp);

Uhat = PhiCPk.U{1};
Vhat = PhiCPk.U{end};
onlineCP.lambda = PhiCPk.lambda;
onlineCP.sigma  = PhiCPk.U(2:end-1);

% thin QR factorizations
[UCP, onlineCP.RU] = qr(Uhat, 0);
[~,   onlineCP.RV] = qr(Vhat, 0);

end
