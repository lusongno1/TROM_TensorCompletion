function [Uc] = tromcponline(onlineCP, ea)
%TROMCPONLINE Online stage of CP-TROM computation
% Input:  
%  onlineCP = struct computed by TROMCPOFFLINE
%  ea = cell array of interpolation vectors e^i(\alpha), i=1,...,D
% Output:
%  Uc = RxR (orthogonal) matrix of coordinates of local reduced basis
%
% Alexander Mamonov, University of Houston, 2022
%==========================================================================

% inner products of \sigma_i and e^i(\alpha), for all R at once
sea = zeros(onlineCP.D, onlineCP.R);

for i = 1:onlineCP.D
	sea(i, :) = ea{i}' * onlineCP.sigma{i};
end

% row vector, diagonal of S(\alpha)
sa = onlineCP.lambda' .* prod(sea, 1);

% sanity check
if any(size(sa) ~= [1, onlineCP.R]) 
    error('sa must be a row vector'); 
end

% sa must be row vector to have an effect 
% of multiplying by diag(sa) on the right
Ca = (onlineCP.RU .* sa) * onlineCP.RV'; % core matrix 

[Uc, ~, ~] = svd(Ca);

end
