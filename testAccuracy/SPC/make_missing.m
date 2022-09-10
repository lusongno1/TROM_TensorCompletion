sizeParas = sz(2:end-1);
Ws = genMaskTensor(sizeParas,missingRate);
U = tenones(sz(1),1);
V = tenones(sz(end),1);
W  = ttt(U,tensor(Ws));
W = ttt(W,V);
W = squeeze(W);