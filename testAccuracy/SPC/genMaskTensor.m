function W = genMaskTensor(para_size,missing_rate)
N = prod(para_size);
idd = (randperm(N) > N*missing_rate);
W   = reshape(idd,para_size);
end

