function [eia, indknn] = lagrangeinterp(alphai, alphahati, k)
%LAGRANGEINTERP scalar polynomial interpolation of degree k-1
% computes interpolation coefficients at a point alphai via
% k-NN of alphai among the samples in alphahati
% 
% Input:
%  alphai = scalar value of i-th parameter to be interpolated
%  alphahati = samples of i-th parameter, vector of length n_i
%  k = number of nearest neighbors to use in interpolation
% 
% Output:
%  eia = interpolation vector of length n_i
%  indknn = indices of alphahati containing k-NN of alphai,
%           also indices of nonzero entries of eia
%
% Alexander Mamonov, University of Houston, 2022
%==========================================================================

ni = length(alphahati);

% find k-NN of parval
[~, ind] = sort(abs(alphahati - alphai));
indknn = sort(ind(1:k));

% interpolation points, reshape into a row vector
alphaint = reshape(alphahati(indknn), 1, k);

% interpolation coefficients
eia = zeros(ni, 1);

for j = 1:k
    ej = zeros(1, k);
    ej(j) = 1;
    % j-th Lagrange polynomial of degree k-1 evaluated at parval
    eia(indknn(j)) = lagrangepoly(alphaint, ej, alphai);
end

end
