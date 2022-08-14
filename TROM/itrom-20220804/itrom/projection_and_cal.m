% only the first nCP columns of UcCP are used as
% the coordinates of the local basis
UcCPn = UcCP(:, 1:nCP);

% online ROM FEM matrices
FEMCPc = [];

% pre-computed in the offline stage, parameter-independent
FEMCPc.M = UcCPn' * FEMCP.M * UcCPn;
FEMCPc.K = UcCPn' * FEMCP.K * UcCPn;

% parameter-dependent, need to project on UCP
FEMCPc.Q = UcCPn' * (UCP' * FEMa.Q * UCP) * UcCPn;
FEMCPc.G = UcCPn' * (UCP' * FEMa.G);

% ROM snapshot computation
PhiCPc = UCP * (UcCPn * tstepfemcn0qr(FEMCPc, tlist));

errPhiCP = norm(PhiCPc - Phia) / norm(Phia);
fprintf('CP ROM snapshot relative error: %e\n', errPhiCP);