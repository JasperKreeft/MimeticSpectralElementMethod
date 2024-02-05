function [Ubc] = boundaryconditions_0_stokes_triangle(Mesh)

% int a^0 wedge star b^1     or     int a^{n-2} wedge star b^{n-1}

global N
global nr_0

Ubc = zeros(nr_0,1);

Ubc(1:N+1) = -1/2;
Ubc((N+1)^2+(1:N)) = -1/2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%