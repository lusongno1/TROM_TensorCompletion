function [Uc,Sc] = tromhosvdonline(onlineHOSVD, ea)
%TROMHOSVDONLINE Online stage of HOSBD-TROM computation
% Input:  
%  onlineHOSVD = struct computed by TROMHOSVDOFFLINE
%  ea = cell array of interpolation vectors e^i(\alpha), i=1,...,D
% Output:
%  Uc = \tilde M x \tilde M  (orthogonal) matrix of coordinates of local reduced basis
%
% Alexander Mamonov, University of Houston, 2022
%==========================================================================

sea = cell(onlineHOSVD.D, 1);

for i = 1:onlineHOSVD.D
	sea{i} = onlineHOSVD.S{i}'*ea{i};
end

% core matrix 
Ca = ttv(onlineHOSVD.C, sea, 2:onlineHOSVD.D+1);

% coordinates of local basis in the universal one
[Uc, Sc, ~] = svd(double(Ca), 'econ');

end