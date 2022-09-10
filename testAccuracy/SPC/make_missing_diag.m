%make sure massingRate>0.6
sz1 = [9,5,5,5];
W1 = zeros(5,5,5,5);
for i=1:5
    W1(i,i,i,i) = 1;
end

sz2 = [4,5,5,5];
Navail = prod(sz1)*(1-missingRate);
missingRateTemp = 1-(Navail - 5)/prod(sz2);
W2 = genMaskTensor(sz2,missingRateTemp);
Ws = cat(1,W1,W2);
%checkRate = nnz(Ws)/prod(sz1);

U = tenones(sz(1),1);
V = tenones(sz(end),1);
W  = ttt(U,tensor(Ws));
W = ttt(W,V);
W = squeeze(W);