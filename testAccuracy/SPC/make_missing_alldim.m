II = size(Y_true);
N  = prod(II);
idd = (randperm(N) > N*missingRate);
W   = reshape(idd,II);