function [onlineHOSVD, U, PhiHOSVD] = tromhosvdoffline(Phi, tol)
%TROMHOSVDOFFLINE Offline stage of HOSVD-TROM computation
% Input:  
%  Phi = snapshot multi-array
%  tol = hosvd tolerance
% Output:
%  PhiHOSVD = Low-lank approximation to Phi in Tuker tensor format
%  onlineHOSVD.S  = cell array of matrices S_i, 
%  onlineHOSVD.C = core C tensor of the Tuker format
%  onlineHOSVD.ranks = Tuker ranks.
%  U = Mx\tile{M} orthogonal matrix with columns forming and orthonormal basis 
%        for the universal space
%
% Alexander Mamonov, University of Houston, 2022
%==========================================================================

onlineHOSVD = [];

sizePhi = size(Phi);
onlineHOSVD.D  = length(sizePhi) - 2;

Phit = tensor(Phi); % tensor object from Phi multi-array

% default hosvd algorithm parameter values
if nargin < 2, tol   = 1e-5; end

fprintf('HOCVD-TROM offline stage:\n');
PhiHOSVD = hosvd(Phit, tol); 
   
U = PhiHOSVD.U{1};
onlineHOSVD.S  = PhiHOSVD.U(2:end-1);
onlineHOSVD.ranks = size(PhiHOSVD.core);
onlineHOSVD.C=PhiHOSVD.core;
end
