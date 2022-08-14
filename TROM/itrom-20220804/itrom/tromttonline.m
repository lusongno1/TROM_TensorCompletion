function [Uc,Sc] = tromttonline(onlineTT, ea)
%TROMTTONLINE Online stage of HOSBD-TROM computation
% Input:  
%  onlineTT = struct computed by TROMTTOFFLINE
%  ea = cell array of interpolation vectors e^i(\alpha), i=1,...,D
% Output:
%  Uc = \tilde M x \tilde M  (orthogonal) matrix of coordinates of local reduced basis
%
% Alexander Mamonov, University of Houston, 2022
%==========================================================================
Ctmp = cell(onlineTT.D, 1);

for i = 1:onlineTT.D
    Ctmp{i} = double(ttv(onlineTT.S{i}, ea{i}, 2));
end

Ca = Ctmp{onlineTT.D} * diag(onlineTT.Wc); 
for i = onlineTT.D-1:-1:1
    Ca = Ctmp{i} * Ca;
end

% coordinates of local basis in the universal one
[Uc, Sc, ~] = svd(double(Ca), 'econ');

end