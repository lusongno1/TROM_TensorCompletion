% project FEM matrices onto UCP-basis for use at the online stage
% only mass and stiffness matrix projections can be pre-computed,
% since Q and G depend on parameters that enter via
% the boundary conditions
FEMCP = [];
FEMCP.M = UCP' * FEM.M * UCP;
FEMCP.K = UCP' * FEM.K * UCP;