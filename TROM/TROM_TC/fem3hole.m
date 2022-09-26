function [FEM] = fem3hole(model, uo, ui1, ui2, ui3, Bio, Bii)
%FEM3HOLE Assemble FEM matrices for given parameters
% for the rectangular 3-hole rectangular model
%
% Alexander Mamonov, University of Houston, 2021
%==========================================================================

% convection bdry cond: n·∇u + Bi*u = Bi*u_ext
% Matlab generalized Neumann, i.e., Robin: n·(c×∇u) + qu = g

% left boundary: Omega_out
applyBoundaryCondition(model, 'neumann', 'Edge', 4, ...
                       'q', Bio, 'g', Bio*uo);

% internal boundaries: Omega_in
% (normalized) ui1, ui2, ui3-temperature air flow
applyBoundaryCondition(model, 'neumann', 'Edge', [6,7,9,13], ...
                       'q', Bii, 'g', Bii*ui1); 
                   
applyBoundaryCondition(model, 'neumann', 'Edge', [8,10,12,14], ...
                       'q', Bii, 'g', Bii*ui2);
                   
applyBoundaryCondition(model, 'neumann', 'Edge', [3,11,15,16], ...
                       'q', Bii, 'g', Bii*ui3);
                   
% right, top and bottom: insulated, zero Neumann
applyBoundaryCondition(model, 'neumann', 'Edge', [1, 2, 5], ...
                       'q', 0, 'g', 0);

% zero initial condition
setInitialConditions(model, 0);

% Assembling FEM matrices
FEM = assembleFEMatrices(model); 
