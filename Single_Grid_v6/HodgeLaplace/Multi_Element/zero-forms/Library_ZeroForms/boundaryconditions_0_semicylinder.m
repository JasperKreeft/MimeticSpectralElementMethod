% RAAR POTENTIAL vs STOKES

function [PHIbc,Ubc,boundary_points,interior_points] = boundaryconditions_0_semicylinder()

% int a^0 wedge star b^1     or     int a^{n-2} wedge star b^{n-1}

global N numElements
global globalnr_0
global nr_0


boundary_points = globalnr_0(1:N+1:(N+1)^2,1);

interior_points = (1:nr_0)';
interior_points(boundary_points) = [];

PHIbc = ones(N+1,1)*3/8;
Ubc = zeros(nr_0,1);
ind = globalnr_0(N+1:N+1:(N+1)^2,6);
Ubc(ind) = 3/8;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%