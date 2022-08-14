function [Phi, FEM] = solve3hole(model, uo, ui1, ui2, ui3, Bio, Bii, tlist)
%MSOLVE3HOLE Generate the snapshots 
% for the rectangular 3-hole rectangular model
%
% Alexander Mamonov, University of Houston, 2021
%==========================================================================

% assemble FEM matrices
FEM = fem3hole(model, uo, ui1, ui2, ui3, Bio, Bii);

% Solve using custom Crank – Nicolson scheme
Phi = tstepfemcn0qr(FEM, tlist);

end
