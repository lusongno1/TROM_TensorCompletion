function [onlineTT, U, PhiTT] = tromttoffline(Phi, tol)
%TROMTTOFFLINE Offline stage of TT-TROM computation
% Input:  
%  Phi = snapshot multi-array
%  tol = TT tolerance
% Output:
%  PhiTT = Low-lank approximation to Phi in TT tensor format
%  onlineTT.S  = cell array of tensors S_i, 
%  onlineTT.Wc = scaling matrix Wc
%  onlineTT.r = TT ranks.
%  onlineTT.D - dimension of the parameter space
%  U = M x\tilde{r_1} orthogonal matrix with columns forming and orthonormal basis 
%        for the universal space
%
% Alexander Mamonov, University of Houston, 2022
%==========================================================================

onlineTT = [];

sizePhi = size(Phi);
onlineTT.D  = length(sizePhi) - 2;

% default TT algorithm parameter values
if nargin < 2, tol   = 1e-5; end

fprintf('TT-TROM offline stage:\n');
    PhiTT = tt_tensor(Phi, tol);
    
    errTT = norm(Phi(:) - full(PhiTT)) / norm(Phi(:));
    
    fprintf('TT decomposition with ranks:'); 
    fprintf('[%d]', PhiTT.r); 
    fprintf('and accuracy = %g\n', errTT);
    
LASTdim=size(size(PhiTT),2);
Vtmp = core(PhiTT, LASTdim); 
onlineTT.Wc = sqrt(sum(Vtmp.^2, 2));
onlineTT.r = PhiTT.r;

U = squeeze(core(PhiTT, 1));
for i = 2:LASTdim-1
    onlineTT.S{i-1} = tensor(core(PhiTT, i));
end

end